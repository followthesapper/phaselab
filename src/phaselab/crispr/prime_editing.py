"""
PhaseLab Prime Editing: pegRNA design for precise genome editing.

Implements:
- pegRNA (prime editing guide RNA) design
- PBS (Primer Binding Site) optimization
- RT (Reverse Transcriptase) template design
- Secondary structure risk assessment
- Off-target prime editing site screening
- IR coherence validation
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass

from .pam_scan import find_pam_sites, PAMHit
from .scoring import (
    gc_content,
    delta_g_santalucia,
    mit_specificity_score,
    cfd_score,
    max_homopolymer_run,
    sequence_complexity,
)
from ..core.coherence import coherence_score, go_no_go
from ..core.hamiltonians import build_grna_hamiltonian


# PBS length optimization parameters
PBS_OPTIMAL_LENGTH = (13, 15)  # bp
PBS_MIN_LENGTH = 8
PBS_MAX_LENGTH = 17

# RT template length parameters
RT_OPTIMAL_LENGTH = (10, 16)  # bp
RT_MIN_LENGTH = 7
RT_MAX_LENGTH = 34

# Nick-to-edit distance (for nicking guides)
NICK_TO_EDIT_OPTIMAL = (40, 100)  # bp

# Secondary structure energy thresholds (kcal/mol)
SECONDARY_STRUCTURE_WARNING = -5.0
SECONDARY_STRUCTURE_FAIL = -10.0


@dataclass
class PrimeEditConfig:
    """Configuration for prime editing guide design."""

    # PAM settings
    pam: str = "NGG"
    guide_length: int = 20

    # PBS parameters
    pbs_length_min: int = PBS_MIN_LENGTH
    pbs_length_max: int = PBS_MAX_LENGTH
    pbs_optimal_length: Tuple[int, int] = PBS_OPTIMAL_LENGTH
    pbs_gc_min: float = 0.35
    pbs_gc_max: float = 0.65

    # RT template parameters
    rt_length_min: int = RT_MIN_LENGTH
    rt_length_max: int = RT_MAX_LENGTH
    rt_optimal_length: Tuple[int, int] = RT_OPTIMAL_LENGTH

    # Edit type
    edit_type: str = "substitution"  # "substitution", "insertion", "deletion"

    # Filtering
    min_gc: float = 0.35
    max_gc: float = 0.75
    max_homopolymer: int = 4

    # Secondary structure
    check_secondary_structure: bool = True
    max_secondary_structure_dg: float = -5.0

    # Scoring weights
    weight_mit: float = 1.0
    weight_cfd: float = 1.0
    weight_pbs_score: float = 1.2
    weight_rt_score: float = 1.2
    weight_nick_distance: float = 0.8
    weight_coherence: float = 1.0

    # Coherence
    compute_coherence: bool = True
    coherence_shots: int = 2000

    # Output
    top_n: int = 10


@dataclass
class PegRNACandidate:
    """A candidate pegRNA design."""
    spacer: str              # 20bp spacer sequence
    pam: str                 # PAM sequence
    pbs: str                 # Primer binding site
    rt_template: str         # RT template (includes edit)
    edit_position: int       # Position of edit in genome
    edit_type: str           # substitution, insertion, deletion
    edit_from: str           # Original sequence
    edit_to: str             # New sequence

    @property
    def full_pegrna(self) -> str:
        """Full pegRNA sequence (spacer + scaffold + RT + PBS)."""
        # Note: This is a simplified representation
        # Real pegRNA includes scaffold sequence between spacer and extension
        return f"{self.spacer}[scaffold]{self.rt_template}{self.pbs}"


def pbs_score(pbs_seq: str) -> float:
    """
    Score a PBS sequence.

    Args:
        pbs_seq: PBS sequence.

    Returns:
        Score (0 to ~1.5).
    """
    gc = gc_content(pbs_seq)
    length = len(pbs_seq)

    # GC scoring
    gc_score = 1.0
    if gc < 0.35 or gc > 0.65:
        gc_score = 0.7
    if gc < 0.25 or gc > 0.75:
        gc_score = 0.4

    # Length scoring
    length_score = 1.0
    if PBS_OPTIMAL_LENGTH[0] <= length <= PBS_OPTIMAL_LENGTH[1]:
        length_score = 1.2

    # Secondary structure
    ss_dg = estimate_hairpin_dg(pbs_seq)
    ss_score = 1.0 if ss_dg > SECONDARY_STRUCTURE_WARNING else 0.6

    return gc_score * length_score * ss_score


def rt_template_score(rt_seq: str) -> float:
    """
    Score an RT template sequence.

    Args:
        rt_seq: RT template sequence.

    Returns:
        Score (0 to ~1.5).
    """
    gc = gc_content(rt_seq)
    length = len(rt_seq)

    # GC scoring
    gc_score = 1.0
    if gc < 0.35 or gc > 0.65:
        gc_score = 0.7

    # Length scoring
    length_score = 1.0
    if RT_OPTIMAL_LENGTH[0] <= length <= RT_OPTIMAL_LENGTH[1]:
        length_score = 1.2

    # Secondary structure
    ss_dg = estimate_hairpin_dg(rt_seq)
    ss_score = 1.0 if ss_dg > SECONDARY_STRUCTURE_WARNING else 0.6

    return gc_score * length_score * ss_score


def design_pbs(
    target_sequence: str,
    nick_position: int,
    pbs_lengths: List[int] = None,
) -> List[Dict[str, Any]]:
    """
    Design Primer Binding Sites for pegRNA.

    PBS binds to the nicked DNA strand to initiate RT-mediated synthesis.

    Args:
        target_sequence: DNA sequence around target site.
        nick_position: Position where Cas9 nicks (3bp upstream of PAM).
        pbs_lengths: List of PBS lengths to try.

    Returns:
        List of PBS designs with scores.
    """
    if pbs_lengths is None:
        pbs_lengths = list(range(PBS_MIN_LENGTH, PBS_MAX_LENGTH + 1))

    pbs_designs = []

    for length in pbs_lengths:
        # PBS is complementary to sequence upstream of nick
        start = nick_position - length
        if start < 0:
            continue

        pbs_target = target_sequence[start:nick_position]
        pbs_seq = reverse_complement(pbs_target)

        # Score PBS
        gc = gc_content(pbs_seq)
        tm = estimate_melting_temp(pbs_seq)

        # Penalize extreme GC
        gc_score = 1.0
        if gc < 0.35 or gc > 0.65:
            gc_score = 0.7
        if gc < 0.25 or gc > 0.75:
            gc_score = 0.4

        # Optimal length bonus
        length_score = 1.0
        if PBS_OPTIMAL_LENGTH[0] <= length <= PBS_OPTIMAL_LENGTH[1]:
            length_score = 1.2

        # Secondary structure penalty
        ss_dg = estimate_hairpin_dg(pbs_seq)
        ss_score = 1.0 if ss_dg > SECONDARY_STRUCTURE_WARNING else 0.6

        overall_score = gc_score * length_score * ss_score

        pbs_designs.append({
            'sequence': pbs_seq,
            'length': length,
            'gc': gc,
            'tm': tm,
            'secondary_structure_dg': ss_dg,
            'score': overall_score,
        })

    return sorted(pbs_designs, key=lambda x: x['score'], reverse=True)


def design_rt_template(
    target_sequence: str,
    nick_position: int,
    edit_position: int,
    edit_from: str,
    edit_to: str,
    rt_lengths: List[int] = None,
) -> List[Dict[str, Any]]:
    """
    Design RT templates for pegRNA.

    RT template encodes the desired edit and homology arms.

    Args:
        target_sequence: DNA sequence around target site.
        nick_position: Position where Cas9 nicks.
        edit_position: Position of the desired edit.
        edit_from: Original sequence at edit site.
        edit_to: Desired new sequence.
        rt_lengths: List of RT template lengths to try.

    Returns:
        List of RT template designs with scores.
    """
    if rt_lengths is None:
        rt_lengths = list(range(RT_MIN_LENGTH, RT_MAX_LENGTH + 1))

    rt_designs = []
    edit_offset = edit_position - nick_position

    for length in rt_lengths:
        # Calculate template boundaries
        # RT synthesizes from nick toward edit and beyond
        rt_start = nick_position
        rt_end = nick_position + length

        if rt_end > len(target_sequence):
            continue
        if edit_position > rt_end:
            continue  # Edit not covered

        # Build RT template with edit incorporated
        template_seq = target_sequence[rt_start:rt_end]

        # Insert the edit
        edit_rel_pos = edit_position - rt_start
        if 0 <= edit_rel_pos < len(template_seq):
            # For substitution
            if len(edit_from) == len(edit_to):
                template_with_edit = (
                    template_seq[:edit_rel_pos] +
                    edit_to +
                    template_seq[edit_rel_pos + len(edit_from):]
                )
            # For insertion
            elif len(edit_from) < len(edit_to):
                template_with_edit = (
                    template_seq[:edit_rel_pos] +
                    edit_to +
                    template_seq[edit_rel_pos + len(edit_from):]
                )
            # For deletion
            else:
                template_with_edit = (
                    template_seq[:edit_rel_pos] +
                    edit_to +
                    template_seq[edit_rel_pos + len(edit_from):]
                )
        else:
            template_with_edit = template_seq

        # RT template is reverse complement (for RNA)
        rt_rna = reverse_complement(template_with_edit)

        # Score RT template
        gc = gc_content(rt_rna)
        homology_3prime = min(10, length - edit_offset)  # 3' homology after edit

        gc_score = 1.0
        if gc < 0.35 or gc > 0.65:
            gc_score = 0.7

        length_score = 1.0
        if RT_OPTIMAL_LENGTH[0] <= length <= RT_OPTIMAL_LENGTH[1]:
            length_score = 1.2

        # More 3' homology is better
        homology_score = min(1.0, homology_3prime / 10)

        ss_dg = estimate_hairpin_dg(rt_rna)
        ss_score = 1.0 if ss_dg > SECONDARY_STRUCTURE_WARNING else 0.6

        overall_score = gc_score * length_score * homology_score * ss_score

        rt_designs.append({
            'sequence': rt_rna,
            'length': length,
            'gc': gc,
            'homology_3prime': homology_3prime,
            'secondary_structure_dg': ss_dg,
            'score': overall_score,
        })

    return sorted(rt_designs, key=lambda x: x['score'], reverse=True)


def design_prime_edit(
    sequence: str,
    edit_position: int,
    edit_from: str,
    edit_to: str,
    config: Optional[PrimeEditConfig] = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Design pegRNAs for prime editing.

    This designs complete prime editing guide RNAs with optimized
    PBS and RT template components.

    Args:
        sequence: DNA sequence containing the edit site.
        edit_position: Position of the edit in sequence (0-based).
        edit_from: Original sequence to replace.
        edit_to: New sequence (the edit).
        config: PrimeEditConfig with parameters.
        verbose: Print progress messages.

    Returns:
        DataFrame with ranked pegRNA candidates.

    Example:
        >>> from phaselab.crispr import design_prime_edit
        >>> pegrnas = design_prime_edit(
        ...     sequence=gene_region,
        ...     edit_position=150,
        ...     edit_from="A",
        ...     edit_to="G",  # A-to-G substitution
        ... )
        >>> print(pegrnas[['spacer', 'pbs_length', 'rt_length', 'score']].head())
    """
    if config is None:
        config = PrimeEditConfig()

    sequence = sequence.upper()
    edit_from = edit_from.upper()
    edit_to = edit_to.upper()

    if verbose:
        print(f"Designing prime editing guides for {edit_from}->{edit_to} at position {edit_position}")

    # Find PAM sites near the edit
    # PE requires the nick to be within ~50bp of the edit
    search_start = max(0, edit_position - 100)
    search_end = min(len(sequence), edit_position + 100)
    search_seq = sequence[search_start:search_end]

    all_hits = find_pam_sites(
        search_seq,
        pam=config.pam,
        guide_length=config.guide_length,
        both_strands=True,
    )

    if verbose:
        print(f"Found {len(all_hits)} PAM sites near edit")

    if not all_hits:
        return _empty_prime_edit_df()

    candidates = []

    for hit in all_hits:
        # Adjust positions back to original sequence coordinates
        abs_guide_start = search_start + hit.guide_start
        abs_guide_end = search_start + hit.guide_end

        # Nick position (3bp upstream of PAM for NGG)
        if hit.strand == '+':
            nick_pos = abs_guide_end - 3
        else:
            nick_pos = abs_guide_start + 3

        # Distance from nick to edit
        nick_to_edit = edit_position - nick_pos
        if hit.strand == '-':
            nick_to_edit = -nick_to_edit

        # PE3 requires nick downstream of edit on the edited strand
        # PE2 can work with various configurations
        if nick_to_edit < -5 or nick_to_edit > 100:
            continue  # Too far from edit

        guide_seq = hit.guide

        # Quality filters
        gc = gc_content(guide_seq)
        if gc < config.min_gc or gc > config.max_gc:
            continue

        homo = max_homopolymer_run(guide_seq)
        if homo > config.max_homopolymer:
            continue

        # Design PBS options
        pbs_designs = design_pbs(
            sequence,
            nick_pos,
            list(range(config.pbs_length_min, config.pbs_length_max + 1)),
        )
        if not pbs_designs:
            continue

        best_pbs = pbs_designs[0]

        # Design RT template options
        rt_designs = design_rt_template(
            sequence,
            nick_pos,
            edit_position,
            edit_from,
            edit_to,
            list(range(config.rt_length_min, config.rt_length_max + 1)),
        )
        if not rt_designs:
            continue

        best_rt = rt_designs[0]

        # Standard scores
        delta_g = delta_g_santalucia(guide_seq)
        mit = mit_specificity_score(guide_seq)
        cfd = cfd_score(guide_seq)

        # IR coherence
        R_bar = None
        go_status = None
        if config.compute_coherence:
            R_bar = _compute_guide_coherence(guide_seq)
            go_status = go_no_go(R_bar)

        # Combined pegRNA secondary structure check
        full_extension = best_rt['sequence'] + best_pbs['sequence']
        extension_ss_dg = estimate_hairpin_dg(full_extension)

        # Combined score
        combined = _compute_prime_edit_score(
            pbs_score=best_pbs['score'],
            rt_score=best_rt['score'],
            nick_to_edit=nick_to_edit,
            mit=mit,
            cfd=cfd,
            R_bar=R_bar,
            config=config,
        )

        candidates.append({
            'spacer': guide_seq,
            'pam': hit.pam,
            'strand': hit.strand,
            'nick_position': nick_pos,
            'nick_to_edit_distance': nick_to_edit,
            'pbs_sequence': best_pbs['sequence'],
            'pbs_length': best_pbs['length'],
            'pbs_gc': round(best_pbs['gc'], 3),
            'rt_sequence': best_rt['sequence'],
            'rt_length': best_rt['length'],
            'rt_gc': round(best_rt['gc'], 3),
            'extension_ss_dg': round(extension_ss_dg, 2),
            'gc': round(gc, 3),
            'delta_g': round(delta_g, 3),
            'mit_score': round(mit, 1),
            'cfd_score': round(cfd, 1),
            'coherence_R': round(R_bar, 4) if R_bar else None,
            'go_no_go': go_status,
            'combined_score': round(combined, 3),
        })

    if not candidates:
        return _empty_prime_edit_df()

    df = pd.DataFrame(candidates)
    df.sort_values(by='combined_score', ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)

    if verbose:
        print(f"Returning {min(len(df), config.top_n)} top candidates")

    return df.head(config.top_n)


def _empty_prime_edit_df() -> pd.DataFrame:
    """Return empty DataFrame with prime edit columns."""
    return pd.DataFrame(columns=[
        'spacer', 'pam', 'strand', 'nick_position', 'nick_to_edit_distance',
        'pbs_sequence', 'pbs_length', 'pbs_gc',
        'rt_sequence', 'rt_length', 'rt_gc',
        'extension_ss_dg', 'gc', 'delta_g', 'mit_score', 'cfd_score',
        'coherence_R', 'go_no_go', 'combined_score',
    ])


def _compute_guide_coherence(guide_seq: str) -> float:
    """Compute IR coherence for guide sequence."""
    try:
        H = build_grna_hamiltonian(guide_seq)
        terms = H.get_terms()
        if not terms:
            return 0.5

        energies = [abs(coeff) for coeff, _ in terms]
        mean_energy = np.mean(energies)
        std_energy = np.std(energies)

        if mean_energy > 0:
            V_phi = std_energy / mean_energy
            R_bar = np.exp(-V_phi / 2)
        else:
            R_bar = 0.5

        return float(np.clip(R_bar, 0, 1))
    except Exception:
        return 0.5


def _compute_prime_edit_score(
    pbs_score: float,
    rt_score: float,
    nick_to_edit: int,
    mit: float,
    cfd: float,
    R_bar: Optional[float],
    config: PrimeEditConfig,
) -> float:
    """Compute combined score for prime editing."""
    score = 0.0

    # PBS and RT template scores
    score += config.weight_pbs_score * pbs_score
    score += config.weight_rt_score * rt_score

    # Nick-to-edit distance (optimal 10-50bp)
    nick_score = 1.0
    if 10 <= nick_to_edit <= 50:
        nick_score = 1.2
    elif nick_to_edit < 0 or nick_to_edit > 80:
        nick_score = 0.6
    score += config.weight_nick_distance * nick_score

    # Specificity
    score += config.weight_mit * (mit / 100.0)
    score += config.weight_cfd * (cfd / 100.0)

    # Coherence
    if R_bar is not None:
        score += config.weight_coherence * R_bar

    return score


# Utility functions

def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(b, 'N') for b in reversed(seq.upper()))


def estimate_melting_temp(seq: str) -> float:
    """Estimate melting temperature using simplified formula."""
    seq = seq.upper()
    gc = sum(1 for b in seq if b in 'GC')
    at = sum(1 for b in seq if b in 'AT')

    if len(seq) < 14:
        # Short oligo formula
        tm = 2 * at + 4 * gc
    else:
        # Long oligo formula (simplified)
        tm = 64.9 + 41 * (gc - 16.4) / len(seq)

    return tm


def estimate_hairpin_dg(seq: str) -> float:
    """
    Estimate secondary structure ΔG (simplified).

    More negative = stronger secondary structure = worse.

    Returns:
        Estimated ΔG in kcal/mol.
    """
    seq = seq.upper()
    n = len(seq)

    if n < 8:
        return 0.0

    # Simple heuristic: check for complementary regions
    dg = 0.0

    # Look for potential hairpin stems (4-8bp)
    for stem_len in range(4, min(9, n // 2)):
        for i in range(n - 2 * stem_len - 3):
            stem1 = seq[i:i + stem_len]
            # Check region after potential loop (3-8bp loop)
            for loop_len in range(3, 9):
                j = i + stem_len + loop_len
                if j + stem_len > n:
                    continue
                stem2 = seq[j:j + stem_len]
                stem2_rc = reverse_complement(stem2)

                # Count matches
                matches = sum(1 for a, b in zip(stem1, stem2_rc) if a == b)
                if matches >= stem_len - 1:  # Allow 1 mismatch
                    # Penalize based on stem length
                    dg -= 1.5 * matches

    return max(-15.0, dg)  # Cap at reasonable minimum


def validate_prime_edit(
    spacer: str,
    pbs: str,
    rt_template: str,
    edit_type: str = "substitution",
) -> Dict[str, Any]:
    """
    Validate a pegRNA design.

    Args:
        spacer: 20bp spacer sequence.
        pbs: PBS sequence.
        rt_template: RT template sequence.
        edit_type: Type of edit.

    Returns:
        Validation results dictionary.
    """
    warnings = []

    # Check spacer
    spacer_gc = gc_content(spacer)
    if spacer_gc < 0.35 or spacer_gc > 0.75:
        warnings.append(f"Spacer GC {spacer_gc:.0%} out of optimal range")

    # Check PBS
    pbs_gc = gc_content(pbs)
    if pbs_gc < 0.35 or pbs_gc > 0.65:
        warnings.append(f"PBS GC {pbs_gc:.0%} may affect binding")
    if len(pbs) < PBS_MIN_LENGTH:
        warnings.append(f"PBS length {len(pbs)}bp too short")
    if len(pbs) > PBS_MAX_LENGTH:
        warnings.append(f"PBS length {len(pbs)}bp may reduce efficiency")

    # Check RT template
    rt_gc = gc_content(rt_template)
    if rt_gc < 0.35 or rt_gc > 0.65:
        warnings.append(f"RT template GC {rt_gc:.0%} may affect synthesis")
    if len(rt_template) < RT_MIN_LENGTH:
        warnings.append(f"RT template length {len(rt_template)}bp too short")
    if len(rt_template) > RT_MAX_LENGTH:
        warnings.append(f"RT template length {len(rt_template)}bp may reduce efficiency")

    # Check extension secondary structure
    extension = rt_template + pbs
    ss_dg = estimate_hairpin_dg(extension)
    if ss_dg < SECONDARY_STRUCTURE_FAIL:
        warnings.append(f"Strong secondary structure (ΔG={ss_dg:.1f} kcal/mol)")
    elif ss_dg < SECONDARY_STRUCTURE_WARNING:
        warnings.append(f"Moderate secondary structure (ΔG={ss_dg:.1f} kcal/mol)")

    R_bar = _compute_guide_coherence(spacer)

    return {
        'spacer': spacer,
        'spacer_gc': spacer_gc,
        'pbs': pbs,
        'pbs_gc': pbs_gc,
        'pbs_length': len(pbs),
        'rt_template': rt_template,
        'rt_gc': rt_gc,
        'rt_length': len(rt_template),
        'extension_ss_dg': ss_dg,
        'coherence_R': R_bar,
        'go_no_go': go_no_go(R_bar),
        'warnings': warnings,
        'valid': len(warnings) == 0,
    }
