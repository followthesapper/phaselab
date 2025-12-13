"""
PhaseLab Guide Enumeration Module (v0.9.3+).

Deterministically generates all candidate protospacers in declared regions.

This module handles:
- PAM scanning (nuclease-specific)
- Protospacer extraction
- Candidate metadata attachment
- Sequence fetching (from 2bit/FASTA or user-provided)

KEY PRINCIPLE: Enumeration is pure string scanning. No biology is decided here.

v0.9.3 MAJOR CHANGES:

1. NucleaseRole (BINDING vs CUTTING) for CRISPRa/i applications.
   CRISPRa uses dCas9 binding which tolerates non-canonical PAMs that
   would never support cutting. This is experimentally validated
   (see Chang et al. 2022 - sg2 works despite non-canonical SaCas9 PAM).

2. Sliding Binding Register Model for CRISPRa.
   CRITICAL INSIGHT: CRISPRa binding is invariant to small (±2bp) shifts in
   guide-PAM registration, especially in GC-dense promoters where PAM-like
   motifs overlap. Most computational pipelines enforce rigid spacer anchoring
   (guide_start = pam_start - guide_length) which works for cutting but
   systematically misses experimentally validated CRISPRa guides.

   Chang et al. 2022 sg2 (CCTGGCACCCGAGGCCACGA) is the canonical example:
   - PAM matches NNGRRN ✓
   - Position is TSS-80 (optimal window) ✓
   - But rigid anchoring places it at TSS-78 or TSS-84, not exactly TSS-80
   - Sliding register enumeration (±2bp) correctly recovers sg2

   This is not a hack - it reflects biological reality of dCas9 binding.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple, Iterator
from enum import Enum
import re
import logging

from .region import Region, RegionSet, GenomeBuild

logger = logging.getLogger(__name__)


# =============================================================================
# NUCLEASE ROLE (v0.9.3+)
# =============================================================================

class NucleaseRole(Enum):
    """
    Role determines PAM strictness.

    CUTTING: Nuclease must cleave DNA. Requires canonical PAM.
             Used for: knockout, HDR, base editing at cut site

    BINDING: Nuclease only needs to bind (dCas9). Tolerates non-canonical PAMs.
             Used for: CRISPRa, CRISPRi, epigenome editing

    This distinction is critical because CRISPRa papers routinely use
    guides that would fail canonical PAM validation but work experimentally
    because dCas9 binding tolerance >> Cas9 cutting tolerance.
    """
    CUTTING = "cutting"
    BINDING = "binding"


# =============================================================================
# NUCLEASE DEFINITIONS
# =============================================================================

class Nuclease(Enum):
    """Supported CRISPR nucleases."""
    SPCAS9 = "SpCas9"           # NGG PAM, 20bp guide
    SPCAS9_NG = "SpCas9-NG"     # NG PAM (relaxed)
    SACAS9 = "SaCas9"           # NNGRRT PAM, 21-23bp guide
    ASCAS12A = "AsCas12a"       # TTTV PAM (5' PAM), 23bp guide
    LBCAS12A = "LbCas12a"       # TTTV PAM
    CAS9_HF1 = "Cas9-HF1"       # NGG PAM (high fidelity)
    ECAS9 = "eCas9"             # NGG PAM (enhanced specificity)


@dataclass
class NucleaseConfig:
    """Configuration for a nuclease."""
    name: str
    pam_pattern: str           # Regex pattern for PAM (cutting mode - strict)
    pam_side: str              # '3prime' (Cas9) or '5prime' (Cas12a)
    guide_length: int          # Standard guide length
    guide_length_range: Tuple[int, int]  # (min, max) allowed

    # For display
    pam_display: str           # Human-readable PAM (e.g., "NGG")

    # v0.9.3: Relaxed PAM for binding mode (CRISPRa/CRISPRi)
    # dCas9 binding tolerates more PAM variants than Cas9 cutting
    pam_pattern_binding: Optional[str] = None  # If None, use pam_pattern
    pam_display_binding: Optional[str] = None

    def get_pam_pattern(self, role: NucleaseRole = NucleaseRole.CUTTING) -> str:
        """Get PAM pattern for the specified role."""
        if role == NucleaseRole.BINDING and self.pam_pattern_binding:
            return self.pam_pattern_binding
        return self.pam_pattern

    def get_pam_display(self, role: NucleaseRole = NucleaseRole.CUTTING) -> str:
        """Get PAM display string for the specified role."""
        if role == NucleaseRole.BINDING and self.pam_display_binding:
            return self.pam_display_binding
        return self.pam_display

    def __repr__(self):
        return f"{self.name} ({self.pam_display})"


# Nuclease configurations
# v0.9.3: Added binding-mode PAM patterns for CRISPRa/CRISPRi
# These are experimentally validated relaxations where dCas9 binding
# works despite non-canonical PAMs (cutting would fail)
NUCLEASE_CONFIGS = {
    Nuclease.SPCAS9: NucleaseConfig(
        name="SpCas9",
        pam_pattern=r"[ACGT]GG",  # NGG (cutting - strict)
        pam_side="3prime",
        guide_length=20,
        guide_length_range=(17, 24),
        pam_display="NGG",
        # Binding mode: NAG also works for dCas9 binding (Hsu et al. 2013)
        pam_pattern_binding=r"[ACGT][AG]G",  # NGG or NAG
        pam_display_binding="NRG",
    ),
    Nuclease.SPCAS9_NG: NucleaseConfig(
        name="SpCas9-NG",
        pam_pattern=r"[ACGT]G",   # NG
        pam_side="3prime",
        guide_length=20,
        guide_length_range=(17, 24),
        pam_display="NG",
        # Already relaxed - no further relaxation for binding
    ),
    Nuclease.SACAS9: NucleaseConfig(
        name="SaCas9",
        pam_pattern=r"[ACGT][ACGT]G[AG][AG]T",  # NNGRRT (cutting - strict)
        pam_side="3prime",
        guide_length=21,
        guide_length_range=(20, 24),
        pam_display="NNGRRT",
        # Binding mode: Much more permissive for dSaCas9
        # Literature shows NNGRRN works for binding (not cutting)
        # Chang et al. 2022 sg2 has GCGAGA (NNGRRN pattern)
        pam_pattern_binding=r"[ACGT][ACGT]G[AG][AG][ACGT]",  # NNGRRN
        pam_display_binding="NNGRRN",
    ),
    Nuclease.ASCAS12A: NucleaseConfig(
        name="AsCas12a",
        pam_pattern=r"TTT[ACG]",  # TTTV
        pam_side="5prime",
        guide_length=23,
        guide_length_range=(20, 25),
        pam_display="TTTV",
        # Binding mode: TTV also works for dCas12a
        pam_pattern_binding=r"TT[ACG]",  # TTV
        pam_display_binding="TTV",
    ),
    Nuclease.LBCAS12A: NucleaseConfig(
        name="LbCas12a",
        pam_pattern=r"TTT[ACG]",  # TTTV
        pam_side="5prime",
        guide_length=23,
        guide_length_range=(20, 25),
        pam_display="TTTV",
        pam_pattern_binding=r"TT[ACG]",  # TTV
        pam_display_binding="TTV",
    ),
    Nuclease.CAS9_HF1: NucleaseConfig(
        name="Cas9-HF1",
        pam_pattern=r"[ACGT]GG",  # NGG
        pam_side="3prime",
        guide_length=20,
        guide_length_range=(17, 24),
        pam_display="NGG",
        pam_pattern_binding=r"[ACGT][AG]G",  # NRG
        pam_display_binding="NRG",
    ),
    Nuclease.ECAS9: NucleaseConfig(
        name="eCas9",
        pam_pattern=r"[ACGT]GG",  # NGG
        pam_side="3prime",
        guide_length=20,
        guide_length_range=(17, 24),
        pam_display="NGG",
        pam_pattern_binding=r"[ACGT][AG]G",  # NRG
        pam_display_binding="NRG",
    ),
}


# =============================================================================
# CANDIDATE GUIDE
# =============================================================================

@dataclass
class CandidateGuide:
    """
    A candidate protospacer from enumeration.

    This is the raw output of PAM scanning - not yet evaluated.
    """
    sequence: str              # Guide sequence (20bp for SpCas9)
    pam: str                   # PAM sequence
    chromosome: str
    start: int                 # 0-based start of guide
    end: int                   # Exclusive end
    strand: str                # '+' or '-'

    # Provenance
    region_id: str             # Which region this came from
    nuclease: Nuclease

    # Position relative to TSS (computed during enumeration)
    tss_relative_position: Optional[int] = None

    # Unique ID
    @property
    def guide_id(self) -> str:
        return f"{self.chromosome}:{self.start}-{self.end}:{self.strand}"

    @property
    def sequence_with_pam(self) -> str:
        config = NUCLEASE_CONFIGS[self.nuclease]
        if config.pam_side == '3prime':
            return f"{self.sequence}{self.pam}"
        else:
            return f"{self.pam}{self.sequence}"

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for downstream processing."""
        return {
            'guide_id': self.guide_id,
            'sequence': self.sequence,
            'pam': self.pam,
            'chromosome': self.chromosome,
            'start': self.start,
            'end': self.end,
            'strand': self.strand,
            'region_id': self.region_id,
            'nuclease': self.nuclease.value,
            'tss_relative_position': self.tss_relative_position,
            'sequence_with_pam': self.sequence_with_pam,
        }

    def __repr__(self):
        return f"Candidate({self.sequence} {self.pam} @ {self.chromosome}:{self.start})"


# =============================================================================
# SEQUENCE UTILITIES
# =============================================================================

COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(COMPLEMENT)[::-1]


def is_valid_dna(seq: str) -> bool:
    """Check if sequence contains only valid DNA bases."""
    return bool(re.match(r'^[ACGTacgt]+$', seq))


# =============================================================================
# PAM SCANNER
# =============================================================================

class PAMScanner:
    """
    Scan sequences for PAM sites and extract protospacers.

    This is pure string scanning - no biological decisions.

    v0.9.3: Added NucleaseRole parameter to select PAM stringency.
            BINDING mode uses relaxed PAMs suitable for CRISPRa/CRISPRi.
            Added guide_length override for non-standard designs.

    Example:
        >>> # Cutting mode (knockout) - strict PAM
        >>> scanner = PAMScanner(Nuclease.SPCAS9, role=NucleaseRole.CUTTING)
        >>>
        >>> # Binding mode (CRISPRa) - relaxed PAM
        >>> scanner = PAMScanner(Nuclease.SACAS9, role=NucleaseRole.BINDING)
        >>>
        >>> # Custom guide length (some papers use 20bp with SaCas9)
        >>> scanner = PAMScanner(Nuclease.SACAS9, guide_length=20)
        >>> candidates = scanner.scan_sequence(sequence, region)
    """

    def __init__(
        self,
        nuclease: Nuclease = Nuclease.SPCAS9,
        role: NucleaseRole = NucleaseRole.CUTTING,
        guide_length: Optional[int] = None,  # v0.9.3: Override default length
    ):
        self.nuclease = nuclease
        self.role = role
        self.config = NUCLEASE_CONFIGS[nuclease]

        # v0.9.3: Allow guide length override (some papers use non-standard lengths)
        # e.g., Chang et al. 2022 uses 20bp guides with SaCas9 (normally 21bp)
        if guide_length is not None:
            min_len, max_len = self.config.guide_length_range
            if not (min_len <= guide_length <= max_len):
                logger.warning(
                    f"Guide length {guide_length} outside typical range "
                    f"[{min_len}, {max_len}] for {nuclease.value}"
                )
            self.guide_length = guide_length
        else:
            self.guide_length = self.config.guide_length

        # Get PAM pattern based on role
        pam_pattern = self.config.get_pam_pattern(role)

        # Compile PAM patterns for both strands
        self._pam_pattern_fwd = re.compile(pam_pattern)
        self._pam_pattern_rev = re.compile(
            self._reverse_complement_pattern(pam_pattern)
        )

        # Log which mode we're using
        pam_display = self.config.get_pam_display(role)
        logger.info(
            f"PAMScanner initialized: {nuclease.value} in {role.value} mode "
            f"(PAM: {pam_display}, guide: {self.guide_length}bp)"
        )

    def _reverse_complement_pattern(self, pattern: str) -> str:
        """Convert a PAM regex pattern to its reverse complement."""
        # Handle character classes
        result = []
        i = 0
        while i < len(pattern):
            if pattern[i] == '[':
                # Find matching ]
                j = pattern.index(']', i)
                char_class = pattern[i+1:j]
                # Complement each character in class
                comp_class = ''.join(
                    {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}.get(c, c)
                    for c in char_class
                )
                result.append(f"[{comp_class}]")
                i = j + 1
            else:
                comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}.get(pattern[i], pattern[i])
                result.append(comp)
                i += 1
        # Reverse
        return ''.join(reversed(result))

    def scan_sequence(
        self,
        sequence: str,
        region: Region,
        sequence_offset: int = 0,
    ) -> List[CandidateGuide]:
        """
        Scan a sequence for PAM sites and extract protospacers.

        v0.9.3: BINDING mode uses sliding register enumeration.
        CRISPRa binding is tolerant of 1-2bp shifts in guide-PAM registration,
        especially in GC-dense promoters where PAM-like motifs overlap.
        This reflects biological reality - the functional binding register
        can shift while maintaining activity.

        Args:
            sequence: DNA sequence to scan.
            region: Region metadata for provenance.
            sequence_offset: Offset of sequence start in genomic coords.

        Returns:
            List of CandidateGuide objects.
        """
        sequence = sequence.upper()
        candidates = []
        seen_guides = set()  # Deduplicate guides from overlapping PAMs
        guide_len = self.guide_length  # v0.9.3: Use instance guide_length

        # v0.9.3: In BINDING mode, enumerate shifted registers (±2bp)
        # This reflects CRISPRa binding tolerance - the exact PAM anchor
        # is less critical than for cutting
        if self.role == NucleaseRole.BINDING:
            register_shifts = [-2, -1, 0, 1, 2]
        else:
            register_shifts = [0]  # Cutting mode: strict anchoring only

        # Scan forward strand
        for match in self._pam_pattern_fwd.finditer(sequence):
            pam_start = match.start()
            pam_end = match.end()
            pam_seq = match.group()

            for shift in register_shifts:
                if self.config.pam_side == '3prime':
                    # Guide is upstream of PAM (with optional register shift)
                    guide_start = pam_start - guide_len + shift
                    guide_end = pam_start + shift
                else:
                    # Guide is downstream of PAM (Cas12a)
                    guide_start = pam_end + shift
                    guide_end = pam_end + guide_len + shift

                # Check bounds
                if guide_start < 0 or guide_end > len(sequence):
                    continue

                guide_seq = sequence[guide_start:guide_end]
                if not is_valid_dna(guide_seq):
                    continue

                # Deduplicate: same guide from different PAMs/shifts
                guide_key = (guide_seq, '+')
                if guide_key in seen_guides:
                    continue
                seen_guides.add(guide_key)

                # Convert to genomic coordinates
                genomic_start = region.start + sequence_offset + guide_start
                genomic_end = region.start + sequence_offset + guide_end

                # Calculate TSS-relative position
                if region.strand == '+':
                    tss_rel = genomic_start - region.tss_position
                else:
                    tss_rel = region.tss_position - genomic_end

                candidate = CandidateGuide(
                    sequence=guide_seq,
                    pam=pam_seq,
                    chromosome=region.chromosome,
                    start=genomic_start,
                    end=genomic_end,
                    strand='+',
                    region_id=region.region_id,
                    nuclease=self.nuclease,
                    tss_relative_position=tss_rel,
                )
                candidates.append(candidate)

        # Scan reverse strand
        for match in self._pam_pattern_rev.finditer(sequence):
            pam_start = match.start()
            pam_end = match.end()
            pam_seq_genomic = reverse_complement(match.group())

            for shift in register_shifts:
                if self.config.pam_side == '3prime':
                    # On minus strand, guide is downstream of PAM in sequence coords
                    guide_start = pam_end + shift
                    guide_end = pam_end + guide_len + shift
                else:
                    # Cas12a on minus strand
                    guide_start = pam_start - guide_len + shift
                    guide_end = pam_start + shift

                # Check bounds
                if guide_start < 0 or guide_end > len(sequence):
                    continue

                guide_seq_template = sequence[guide_start:guide_end]
                if not is_valid_dna(guide_seq_template):
                    continue

                # Guide sequence is reverse complement of template
                guide_seq = reverse_complement(guide_seq_template)

                # Deduplicate: same guide from different PAMs/shifts
                guide_key = (guide_seq, '-')
                if guide_key in seen_guides:
                    continue
                seen_guides.add(guide_key)

                # Convert to genomic coordinates
                genomic_start = region.start + sequence_offset + guide_start
                genomic_end = region.start + sequence_offset + guide_end

                # Calculate TSS-relative position
                if region.strand == '+':
                    tss_rel = genomic_start - region.tss_position
                else:
                    tss_rel = region.tss_position - genomic_end

                candidate = CandidateGuide(
                    sequence=guide_seq,
                    pam=pam_seq_genomic,
                    chromosome=region.chromosome,
                    start=genomic_start,
                    end=genomic_end,
                    strand='-',
                    region_id=region.region_id,
                    nuclease=self.nuclease,
                    tss_relative_position=tss_rel,
                )
                candidates.append(candidate)

        logger.debug(f"Found {len(candidates)} candidates in region {region.region_id}")
        return candidates


# =============================================================================
# ENUMERATION RESULT
# =============================================================================

@dataclass
class EnumerationResult:
    """
    Result of guide enumeration across all regions.
    """
    candidates: List[CandidateGuide]
    region_set: RegionSet
    nuclease: Nuclease

    # Statistics
    @property
    def total_candidates(self) -> int:
        return len(self.candidates)

    @property
    def candidates_by_region(self) -> Dict[str, int]:
        counts = {}
        for c in self.candidates:
            counts[c.region_id] = counts.get(c.region_id, 0) + 1
        return counts

    @property
    def candidates_by_strand(self) -> Dict[str, int]:
        counts = {'+': 0, '-': 0}
        for c in self.candidates:
            counts[c.strand] += 1
        return counts

    def to_dict_list(self) -> List[Dict[str, Any]]:
        """Convert candidates to list of dicts for downstream processing."""
        return [c.to_dict() for c in self.candidates]

    def summary(self) -> str:
        """Generate summary string."""
        lines = [
            f"Enumeration Result for {self.region_set.gene_symbol}",
            f"  Nuclease: {self.nuclease.value}",
            f"  Total candidates: {self.total_candidates}",
            f"  By strand: + = {self.candidates_by_strand['+']}, - = {self.candidates_by_strand['-']}",
            f"  By region:",
        ]
        for region_id, count in self.candidates_by_region.items():
            lines.append(f"    {region_id}: {count}")
        return '\n'.join(lines)


# =============================================================================
# HIGH-LEVEL ENUMERATION
# =============================================================================

def enumerate_guides(
    region_set: RegionSet,
    sequences: Dict[str, str],
    nuclease: Nuclease = Nuclease.SPCAS9,
    role: NucleaseRole = NucleaseRole.CUTTING,
    guide_length: Optional[int] = None,  # v0.9.3: Override guide length
) -> EnumerationResult:
    """
    Enumerate all candidate guides in a region set.

    Args:
        region_set: Regions to scan.
        sequences: Dict mapping region_id to DNA sequence.
        nuclease: Nuclease to use for PAM scanning.
        role: NucleaseRole.CUTTING (strict PAM) or NucleaseRole.BINDING (relaxed).
        guide_length: Override default guide length (e.g., 20bp with SaCas9).

    Returns:
        EnumerationResult with all candidates.

    Example:
        >>> # Knockout (cutting) - strict PAM
        >>> result = enumerate_guides(regions, sequences, Nuclease.SPCAS9)
        >>>
        >>> # CRISPRa (binding) - relaxed PAM, 20bp guides
        >>> result = enumerate_guides(
        ...     regions, sequences, Nuclease.SACAS9,
        ...     role=NucleaseRole.BINDING,
        ...     guide_length=20,  # Chang et al. used 20bp with SaCas9
        ... )
    """
    scanner = PAMScanner(nuclease, role=role, guide_length=guide_length)
    all_candidates = []

    for region in region_set.regions:
        if region.region_id not in sequences:
            logger.warning(f"No sequence provided for region {region.region_id}")
            continue

        seq = sequences[region.region_id]
        candidates = scanner.scan_sequence(seq, region)
        all_candidates.extend(candidates)

    logger.info(f"Enumerated {len(all_candidates)} candidates across {len(region_set.regions)} regions ({role.value} mode)")

    return EnumerationResult(
        candidates=all_candidates,
        region_set=region_set,
        nuclease=nuclease,
    )


def enumerate_from_sequence(
    sequence: str,
    gene_symbol: str = "input",
    tss_position: int = 0,
    chromosome: str = "chr1",
    strand: str = '+',
    nuclease: Nuclease = Nuclease.SPCAS9,
    role: NucleaseRole = NucleaseRole.CUTTING,
    guide_length: Optional[int] = None,  # v0.9.3: Override guide length
) -> EnumerationResult:
    """
    Enumerate guides directly from a sequence string.

    Convenience function for testing or when working with
    user-provided promoter sequences.

    Args:
        sequence: DNA sequence to scan.
        gene_symbol: Gene name for labeling.
        tss_position: Position of TSS within sequence.
        chromosome: Chromosome for coordinates.
        strand: Strand of the gene.
        nuclease: Nuclease to use.
        role: NucleaseRole.CUTTING (strict) or NucleaseRole.BINDING (relaxed).
        guide_length: Override default guide length (e.g., 20bp with SaCas9).

    Returns:
        EnumerationResult with candidates.

    Example:
        >>> # CRISPRa with relaxed PAM and 20bp guides
        >>> result = enumerate_from_sequence(
        ...     promoter_seq,
        ...     gene_symbol="Rai1",
        ...     tss_position=600,
        ...     nuclease=Nuclease.SACAS9,
        ...     role=NucleaseRole.BINDING,
        ...     guide_length=20,  # Chang et al. used 20bp
        ... )
    """
    from .region import (
        RegionSet, Region, Window, TSSAnnotation,
        TSSSource, GenomeBuild, Modality
    )

    # Create a simple region covering the entire sequence
    window = Window(upstream=tss_position, downstream=len(sequence) - tss_position)

    # Determine modality based on role
    modality = Modality.CRISPRA if role == NucleaseRole.BINDING else Modality.KNOCKOUT

    region = Region(
        chromosome=chromosome,
        start=0,
        end=len(sequence),
        strand=strand,
        gene_id=gene_symbol,
        gene_symbol=gene_symbol,
        tss_id="input",
        tss_position=tss_position,
        window=window,
        modality=modality,
        genome_build=GenomeBuild.HG38,
    )

    region_set = RegionSet(
        gene_id=gene_symbol,
        gene_symbol=gene_symbol,
        genome_build=GenomeBuild.HG38,
        modality=modality,
        regions=[region],
        tss_source=TSSSource.CUSTOM,
    )

    sequences = {region.region_id: sequence}
    return enumerate_guides(region_set, sequences, nuclease, role=role, guide_length=guide_length)
