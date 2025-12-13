"""
PhaseLab Region Declaration Module (v0.9.2+).

Handles promoter/region definition for guide enumeration.

KEY PRINCIPLE: PhaseLab evaluates guides *conditional on* a declared
promoter hypothesis. It does NOT decide which TSS is "correct" - that
is a biological declaration made by the user.

This module provides:
- TSS coordinate loading from GENCODE/RefSeq
- Multi-TSS hypothesis support
- Window-based region definition
- Region validation and metadata
"""

from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple, Union
from enum import Enum
import logging

logger = logging.getLogger(__name__)


class GenomeBuild(Enum):
    """Supported genome builds."""
    HG38 = "hg38"
    HG19 = "hg19"
    MM10 = "mm10"
    MM39 = "mm39"


class TSSSource(Enum):
    """TSS annotation sources."""
    GENCODE = "gencode"
    REFSEQ = "refseq"
    ENSEMBL = "ensembl"
    CUSTOM = "custom"  # User-provided coordinates


class Modality(Enum):
    """CRISPR modality - affects window positioning."""
    CRISPRA = "CRISPRa"      # Activation: typically -400 to -50 from TSS
    CRISPRI = "CRISPRi"      # Interference: typically -50 to +300 from TSS
    KNOCKOUT = "knockout"    # Cutting: typically exon-focused
    PRIME = "prime"          # Prime editing: position-specific
    BASE = "base"            # Base editing: position-specific


@dataclass
class TSSAnnotation:
    """
    Transcription Start Site annotation.

    Represents a single TSS from an annotation source.
    """
    gene_id: str
    gene_symbol: str
    transcript_id: str
    chromosome: str
    position: int           # 0-based genomic coordinate
    strand: str             # '+' or '-'
    source: TSSSource
    source_version: Optional[str] = None

    # Confidence/evidence (if available)
    support_level: Optional[int] = None  # GENCODE support level (1-5)
    is_canonical: bool = False           # MANE/canonical transcript

    def __repr__(self):
        return f"TSS({self.gene_symbol}:{self.transcript_id} @ {self.chromosome}:{self.position}{self.strand})"


@dataclass
class Window:
    """
    Relative window around a TSS.

    Coordinates are relative to TSS (negative = upstream).
    For minus-strand genes, upstream/downstream are automatically flipped.
    """
    upstream: int = 400     # bp upstream of TSS (negative direction)
    downstream: int = 50    # bp downstream of TSS (positive direction)
    name: Optional[str] = None

    def to_genomic(self, tss_position: int, strand: str) -> Tuple[int, int]:
        """
        Convert relative window to absolute genomic coordinates.

        Args:
            tss_position: Absolute TSS position.
            strand: '+' or '-'.

        Returns:
            (start, end) as 0-based half-open interval.
        """
        if strand == '+':
            start = tss_position - self.upstream
            end = tss_position + self.downstream
        else:
            # Minus strand: upstream is in positive direction
            start = tss_position - self.downstream
            end = tss_position + self.upstream
        return (start, end)


# Default windows for each modality
DEFAULT_WINDOWS = {
    Modality.CRISPRA: Window(upstream=400, downstream=50, name="CRISPRa_standard"),
    Modality.CRISPRI: Window(upstream=50, downstream=300, name="CRISPRi_standard"),
    Modality.KNOCKOUT: Window(upstream=0, downstream=500, name="exon1"),
}


@dataclass
class Region:
    """
    A genomic region for guide enumeration.

    Represents a specific interval to scan for candidate guides,
    with full provenance tracking.
    """
    chromosome: str
    start: int              # 0-based
    end: int                # exclusive
    strand: str

    # Provenance
    gene_id: str
    gene_symbol: str
    tss_id: str             # Which TSS this region is derived from
    tss_position: int
    window: Window
    modality: Modality

    # Metadata
    genome_build: GenomeBuild = GenomeBuild.HG38

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def region_id(self) -> str:
        """Unique identifier for this region."""
        return f"{self.gene_symbol}_{self.tss_id}_{self.window.name or 'custom'}"

    def __repr__(self):
        return f"Region({self.chromosome}:{self.start}-{self.end} [{self.gene_symbol}])"


@dataclass
class RegionSet:
    """
    Collection of regions for a single target.

    Supports multiple TSS hypotheses and multiple windows per TSS.
    """
    gene_id: str
    gene_symbol: str
    genome_build: GenomeBuild
    modality: Modality

    regions: List[Region] = field(default_factory=list)
    tss_annotations: List[TSSAnnotation] = field(default_factory=list)

    # Metadata for reproducibility
    tss_source: TSSSource = TSSSource.GENCODE
    tss_source_version: Optional[str] = None

    def add_region(self, region: Region):
        """Add a region to the set."""
        self.regions.append(region)

    @property
    def total_length(self) -> int:
        """Total bp covered (may overlap)."""
        return sum(r.length for r in self.regions)

    def to_dict(self) -> Dict[str, Any]:
        """Serialize for manifest."""
        return {
            'gene_id': self.gene_id,
            'gene_symbol': self.gene_symbol,
            'genome_build': self.genome_build.value,
            'modality': self.modality.value,
            'tss_source': self.tss_source.value,
            'tss_source_version': self.tss_source_version,
            'num_regions': len(self.regions),
            'total_length': self.total_length,
            'regions': [
                {
                    'region_id': r.region_id,
                    'chromosome': r.chromosome,
                    'start': r.start,
                    'end': r.end,
                    'strand': r.strand,
                    'tss_id': r.tss_id,
                    'tss_position': r.tss_position,
                    'window': {
                        'upstream': r.window.upstream,
                        'downstream': r.window.downstream,
                    },
                }
                for r in self.regions
            ],
        }


# =============================================================================
# REGION BUILDER
# =============================================================================

class RegionBuilder:
    """
    Build regions from TSS annotations and window specifications.

    This is the main interface for region declaration.

    Example:
        >>> builder = RegionBuilder(GenomeBuild.HG38, Modality.CRISPRA)
        >>> builder.add_tss(TSSAnnotation(...))
        >>> builder.add_window(Window(upstream=400, downstream=50))
        >>> region_set = builder.build("RAI1")
    """

    def __init__(
        self,
        genome_build: GenomeBuild = GenomeBuild.HG38,
        modality: Modality = Modality.CRISPRA,
    ):
        self.genome_build = genome_build
        self.modality = modality
        self.tss_list: List[TSSAnnotation] = []
        self.windows: List[Window] = []
        self.tss_source = TSSSource.CUSTOM
        self.tss_source_version = None

    def add_tss(self, tss: TSSAnnotation) -> 'RegionBuilder':
        """Add a TSS annotation."""
        self.tss_list.append(tss)
        if tss.source != TSSSource.CUSTOM:
            self.tss_source = tss.source
            self.tss_source_version = tss.source_version
        return self

    def add_tss_manual(
        self,
        gene_symbol: str,
        chromosome: str,
        position: int,
        strand: str,
        transcript_id: str = "manual",
        gene_id: Optional[str] = None,
    ) -> 'RegionBuilder':
        """Add a TSS with manual coordinates."""
        tss = TSSAnnotation(
            gene_id=gene_id or gene_symbol,
            gene_symbol=gene_symbol,
            transcript_id=transcript_id,
            chromosome=chromosome,
            position=position,
            strand=strand,
            source=TSSSource.CUSTOM,
        )
        return self.add_tss(tss)

    def add_window(self, window: Window) -> 'RegionBuilder':
        """Add a window specification."""
        self.windows.append(window)
        return self

    def use_default_window(self) -> 'RegionBuilder':
        """Use the default window for the current modality."""
        if self.modality in DEFAULT_WINDOWS:
            self.windows.append(DEFAULT_WINDOWS[self.modality])
        else:
            # Fallback
            self.windows.append(Window(upstream=400, downstream=50, name="default"))
        return self

    def build(self, gene_symbol: Optional[str] = None) -> RegionSet:
        """
        Build the region set.

        Generates one region per (TSS, window) combination.

        Args:
            gene_symbol: Override gene symbol (uses first TSS if not provided).

        Returns:
            RegionSet with all declared regions.
        """
        if not self.tss_list:
            raise ValueError("No TSS annotations provided. Use add_tss() first.")

        if not self.windows:
            logger.warning("No windows specified, using default for modality")
            self.use_default_window()

        # Determine gene info
        first_tss = self.tss_list[0]
        gene_symbol = gene_symbol or first_tss.gene_symbol
        gene_id = first_tss.gene_id

        region_set = RegionSet(
            gene_id=gene_id,
            gene_symbol=gene_symbol,
            genome_build=self.genome_build,
            modality=self.modality,
            tss_annotations=self.tss_list.copy(),
            tss_source=self.tss_source,
            tss_source_version=self.tss_source_version,
        )

        # Generate regions for each (TSS, window) combination
        for tss in self.tss_list:
            for window in self.windows:
                start, end = window.to_genomic(tss.position, tss.strand)

                region = Region(
                    chromosome=tss.chromosome,
                    start=start,
                    end=end,
                    strand=tss.strand,
                    gene_id=tss.gene_id,
                    gene_symbol=tss.gene_symbol,
                    tss_id=tss.transcript_id,
                    tss_position=tss.position,
                    window=window,
                    modality=self.modality,
                    genome_build=self.genome_build,
                )
                region_set.add_region(region)

        logger.info(f"Built {len(region_set.regions)} regions for {gene_symbol}")
        return region_set


# =============================================================================
# TSS LOOKUP (from annotations)
# =============================================================================

# Common TSS data (built-in for key genes)
# This is a fallback when annotation files aren't available
BUILTIN_TSS = {
    ('RAI1', 'hg38'): [
        TSSAnnotation(
            gene_id='ENSG00000108557',
            gene_symbol='RAI1',
            transcript_id='ENST00000225568.10',
            chromosome='chr17',
            position=17584786,  # GENCODE v44
            strand='+',
            source=TSSSource.GENCODE,
            source_version='v44',
            is_canonical=True,
        ),
    ],
    ('RAI1', 'hg19'): [
        TSSAnnotation(
            gene_id='ENSG00000108557',
            gene_symbol='RAI1',
            transcript_id='ENST00000225568',
            chromosome='chr17',
            position=17672279,
            strand='+',
            source=TSSSource.GENCODE,
            source_version='v19',
            is_canonical=True,
        ),
    ],
    ('Rai1', 'mm10'): [
        TSSAnnotation(
            gene_id='ENSMUSG00000062475',
            gene_symbol='Rai1',
            transcript_id='ENSMUST00000078364',
            chromosome='chr11',
            position=60104470,  # From Chang et al. paper context
            strand='+',
            source=TSSSource.GENCODE,
            source_version='vM25',
            is_canonical=True,
        ),
    ],
}


def get_tss_for_gene(
    gene_symbol: str,
    genome_build: GenomeBuild,
    source: TSSSource = TSSSource.GENCODE,
) -> List[TSSAnnotation]:
    """
    Look up TSS annotations for a gene.

    Currently uses built-in data. In production, would query
    GENCODE/RefSeq annotation files.

    Args:
        gene_symbol: Gene name (e.g., "RAI1").
        genome_build: Genome build.
        source: Annotation source to use.

    Returns:
        List of TSSAnnotation objects.
    """
    key = (gene_symbol, genome_build.value)

    if key in BUILTIN_TSS:
        return BUILTIN_TSS[key]

    # Also try case variations
    for (g, b), tss_list in BUILTIN_TSS.items():
        if g.upper() == gene_symbol.upper() and b == genome_build.value:
            return tss_list

    logger.warning(f"No built-in TSS data for {gene_symbol} ({genome_build.value})")
    return []


def build_regions_for_gene(
    gene_symbol: str,
    genome_build: GenomeBuild = GenomeBuild.HG38,
    modality: Modality = Modality.CRISPRA,
    windows: Optional[List[Window]] = None,
    use_canonical_only: bool = True,
) -> RegionSet:
    """
    Convenience function to build regions for a gene using built-in TSS data.

    Args:
        gene_symbol: Gene name.
        genome_build: Genome build.
        modality: CRISPR modality.
        windows: Custom windows (uses default if None).
        use_canonical_only: Only use canonical/MANE transcripts.

    Returns:
        RegionSet for the gene.

    Example:
        >>> regions = build_regions_for_gene("RAI1", GenomeBuild.HG38, Modality.CRISPRA)
        >>> print(regions.regions[0])
    """
    tss_list = get_tss_for_gene(gene_symbol, genome_build)

    if not tss_list:
        raise ValueError(
            f"No TSS data found for {gene_symbol} ({genome_build.value}). "
            f"Use RegionBuilder.add_tss_manual() to provide coordinates."
        )

    if use_canonical_only:
        canonical = [t for t in tss_list if t.is_canonical]
        if canonical:
            tss_list = canonical

    builder = RegionBuilder(genome_build, modality)

    for tss in tss_list:
        builder.add_tss(tss)

    if windows:
        for w in windows:
            builder.add_window(w)
    else:
        builder.use_default_window()

    return builder.build(gene_symbol)
