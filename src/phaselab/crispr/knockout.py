"""
PhaseLab CRISPR Knockout: Cas9 cutting guide design for gene disruption.

Implements:
- Cut efficiency prediction (Rule Set 2 / Doench 2016)
- Repair pathway likelihood (NHEJ vs HDR)
- Off-target editing risk
- Frameshift probability
- IR coherence validation for simulation reliability
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass, field

from .pam_scan import find_pam_sites, filter_by_window, PAMHit
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


# Position-specific weights for cutting efficiency (Rule Set 2 approximation)
# PAM-distal (5') positions 1-10, PAM-proximal (3') positions 11-20
CUTTING_POSITION_WEIGHTS = {
    # Disfavored nucleotides at each position (penalty)
    1: {'T': -0.1},
    2: {'T': -0.1, 'C': -0.05},
    3: {},
    4: {'A': 0.05},  # A slightly favored
    5: {},
    6: {},
    7: {},
    8: {},
    9: {},
    10: {},
    # Seed region (positions 11-20) - more important
    11: {},
    12: {},
    13: {'G': 0.1},  # G favored
    14: {'G': 0.1, 'C': 0.05},
    15: {'G': 0.1},
    16: {'G': 0.15, 'C': 0.1},  # G strongly favored
    17: {'C': 0.1},
    18: {'G': 0.1, 'A': 0.05},
    19: {'G': 0.1, 'C': 0.1},
    20: {'G': 0.2, 'A': -0.1},  # G at position 20 (adjacent to PAM) very important
}


@dataclass
class KnockoutConfig:
    """Configuration for CRISPR knockout guide design."""

    # PAM settings
    pam: str = "NGG"
    guide_length: int = 20

    # Target window (relative to TSS or coding start)
    # For knockout, typically want early exons
    target_window: Tuple[int, int] = (0, 500)  # First 500bp of CDS

    # Filtering thresholds
    min_gc: float = 0.35
    max_gc: float = 0.75
    max_homopolymer: int = 4
    min_complexity: float = 0.5

    # Cut efficiency thresholds
    min_cut_efficiency: float = 0.3  # 30% minimum

    # Scoring weights
    weight_mit: float = 1.0
    weight_cfd: float = 1.0
    weight_cut_efficiency: float = 1.5  # Higher weight for cutting
    weight_delta_g: float = 0.3
    weight_coherence: float = 1.0

    # Coherence settings
    compute_coherence: bool = True
    coherence_shots: int = 2000

    # Output
    top_n: int = 10


def cut_efficiency_score(guide_seq: str) -> float:
    """
    Predict cutting efficiency using Rule Set 2-like model.

    Based on Doench et al. 2016 on-target scoring.
    Higher score = more efficient cutting.

    Args:
        guide_seq: 20bp guide sequence.

    Returns:
        Cut efficiency score (0.0 to 1.0).
    """
    guide_seq = guide_seq.upper()

    if len(guide_seq) != 20:
        return 0.5  # Default for non-standard length

    score = 0.5  # Base score

    # Position-specific preferences
    for i, base in enumerate(guide_seq):
        pos = i + 1  # 1-indexed
        if pos in CUTTING_POSITION_WEIGHTS:
            score += CUTTING_POSITION_WEIGHTS[pos].get(base, 0)

    # GC content preference (optimal around 50-60%)
    gc = gc_content(guide_seq)
    if 0.45 <= gc <= 0.65:
        score += 0.1
    elif gc < 0.30 or gc > 0.80:
        score -= 0.15

    # Penalize poly-T (termination signal)
    if 'TTTT' in guide_seq:
        score -= 0.2

    # Reward G at position 20 (adjacent to PAM)
    if guide_seq[-1] == 'G':
        score += 0.1

    # Penalize low complexity
    complexity = sequence_complexity(guide_seq)
    if complexity < 0.6:
        score -= 0.1

    return float(np.clip(score, 0, 1))


def frameshift_probability(
    guide_position: int,
    exon_length: int,
    cds_position: int,
) -> float:
    """
    Estimate probability of frameshift from indel at cut site.

    Cuts earlier in CDS have higher chance of causing functional knockouts.

    Args:
        guide_position: Position of guide in genomic coordinates.
        exon_length: Length of the exon.
        cds_position: Position within the coding sequence (0-based).

    Returns:
        Frameshift probability estimate (0.0 to 1.0).
    """
    # Cuts in first half of CDS more likely to cause functional KO
    relative_pos = cds_position / exon_length if exon_length > 0 else 0.5

    # Higher probability for early cuts
    if relative_pos < 0.25:
        base_prob = 0.90
    elif relative_pos < 0.50:
        base_prob = 0.75
    elif relative_pos < 0.75:
        base_prob = 0.60
    else:
        base_prob = 0.45

    # NHEJ typically produces ~66% frameshift (1bp or 2bp indels)
    return base_prob * 0.66


def repair_pathway_prediction(
    guide_seq: str,
    local_sequence: Optional[str] = None,
) -> Dict[str, float]:
    """
    Predict repair pathway preference (NHEJ vs HDR).

    For knockout, NHEJ is preferred (causes indels).

    Args:
        guide_seq: Guide RNA sequence.
        local_sequence: Sequence context around cut site.

    Returns:
        Dict with pathway probabilities.
    """
    # Base rates (cell-type dependent, these are defaults)
    nhej_base = 0.85  # NHEJ dominant in most cells
    hdr_base = 0.10
    mmej_base = 0.05  # Microhomology-mediated end joining

    # GC content affects repair
    gc = gc_content(guide_seq)

    # High GC slightly favors precise repair
    if gc > 0.65:
        nhej_base -= 0.05
        hdr_base += 0.03
        mmej_base += 0.02

    # Check for microhomologies if context provided
    if local_sequence:
        # Simplified: check for short repeats near cut site
        local_seq = local_sequence.upper()
        microhomology_score = 0
        for size in [2, 3, 4]:
            for i in range(len(local_seq) - size):
                pattern = local_seq[i:i+size]
                if pattern in local_seq[i+1:]:
                    microhomology_score += 1

        if microhomology_score > 3:
            mmej_base += 0.10
            nhej_base -= 0.10

    # Normalize
    total = nhej_base + hdr_base + mmej_base

    return {
        'NHEJ': nhej_base / total,
        'HDR': hdr_base / total,
        'MMEJ': mmej_base / total,
        'preferred': 'NHEJ' if nhej_base > max(hdr_base, mmej_base) else 'HDR/MMEJ',
    }


def design_knockout_guides(
    sequence: str,
    cds_start: int,
    config: Optional[KnockoutConfig] = None,
    exon_boundaries: Optional[List[Tuple[int, int]]] = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Design guide RNAs for CRISPR knockout.

    This pipeline is optimized for gene disruption through frameshift
    mutations caused by NHEJ repair of Cas9-induced DSBs.

    Args:
        sequence: Gene/exon DNA sequence.
        cds_start: Position of CDS start (ATG) in sequence.
        config: KnockoutConfig with parameters.
        exon_boundaries: Optional list of (start, end) for exon positions.
        verbose: Print progress messages.

    Returns:
        DataFrame with ranked knockout guide candidates.

    Example:
        >>> from phaselab.crispr import design_knockout_guides
        >>> guides = design_knockout_guides(
        ...     sequence=gene_sequence,
        ...     cds_start=200,
        ... )
        >>> print(guides[['sequence', 'cut_efficiency', 'frameshift_prob']].head())
    """
    if config is None:
        config = KnockoutConfig()

    sequence = sequence.upper()

    if verbose:
        print(f"Designing knockout guides for {len(sequence)} bp sequence...")

    # Find PAM sites
    all_hits = find_pam_sites(
        sequence,
        pam=config.pam,
        guide_length=config.guide_length,
        both_strands=True,
    )

    if verbose:
        print(f"Found {len(all_hits)} total PAM sites")

    # Filter to target window (relative to CDS start)
    window_hits = filter_by_window(
        all_hits,
        tss_position=cds_start,
        window=config.target_window,
    )

    if verbose:
        print(f"Filtered to {len(window_hits)} in target window")

    if not window_hits:
        return _empty_knockout_df()

    # Score candidates
    candidates = []
    exon_len = sequence[cds_start:cds_start + config.target_window[1]]
    exon_length = len(exon_len)

    for hit in window_hits:
        guide_seq = hit.guide

        # Quality filters
        gc = gc_content(guide_seq)
        if gc < config.min_gc or gc > config.max_gc:
            continue

        homo = max_homopolymer_run(guide_seq)
        if homo > config.max_homopolymer:
            continue

        complexity = sequence_complexity(guide_seq)
        if complexity < config.min_complexity:
            continue

        # Cutting efficiency
        cut_eff = cut_efficiency_score(guide_seq)
        if cut_eff < config.min_cut_efficiency:
            continue

        # Scores
        delta_g = delta_g_santalucia(guide_seq)
        mit = mit_specificity_score(guide_seq)
        cfd = cfd_score(guide_seq)

        # Position relative to CDS start
        cds_pos = ((hit.guide_start + hit.guide_end) // 2) - cds_start

        # Frameshift probability
        fs_prob = frameshift_probability(
            guide_position=hit.guide_start,
            exon_length=exon_length,
            cds_position=max(0, cds_pos),
        )

        # Repair pathway
        repair = repair_pathway_prediction(guide_seq)

        # IR coherence
        R_bar = None
        go_status = None
        if config.compute_coherence:
            R_bar = _compute_guide_coherence(guide_seq)
            go_status = go_no_go(R_bar)

        # Combined score (knockout-optimized)
        combined = _compute_knockout_score(
            cut_efficiency=cut_eff,
            frameshift_prob=fs_prob,
            mit=mit,
            cfd=cfd,
            delta_g=delta_g,
            R_bar=R_bar,
            config=config,
        )

        candidates.append({
            'sequence': guide_seq,
            'pam': hit.pam,
            'cds_position': cds_pos,
            'strand': hit.strand,
            'gc': round(gc, 3),
            'cut_efficiency': round(cut_eff, 3),
            'frameshift_prob': round(fs_prob, 3),
            'repair_pathway': repair['preferred'],
            'nhej_prob': round(repair['NHEJ'], 3),
            'delta_g': round(delta_g, 3),
            'mit_score': round(mit, 1),
            'cfd_score': round(cfd, 1),
            'coherence_R': round(R_bar, 4) if R_bar else None,
            'go_no_go': go_status,
            'combined_score': round(combined, 3),
        })

    if not candidates:
        return _empty_knockout_df()

    df = pd.DataFrame(candidates)
    df.sort_values(by='combined_score', ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)

    if verbose:
        print(f"Returning {min(len(df), config.top_n)} top candidates")

    return df.head(config.top_n)


def _empty_knockout_df() -> pd.DataFrame:
    """Return empty DataFrame with knockout columns."""
    return pd.DataFrame(columns=[
        'sequence', 'pam', 'cds_position', 'strand',
        'gc', 'cut_efficiency', 'frameshift_prob',
        'repair_pathway', 'nhej_prob', 'delta_g',
        'mit_score', 'cfd_score', 'coherence_R', 'go_no_go',
        'combined_score',
    ])


def _compute_guide_coherence(guide_seq: str) -> float:
    """Compute IR coherence using ATLAS-Q enhanced backend (v0.6.0+)."""
    from .coherence_utils import compute_guide_coherence
    return compute_guide_coherence(guide_seq, use_atlas_q=True)


def _compute_knockout_score(
    cut_efficiency: float,
    frameshift_prob: float,
    mit: float,
    cfd: float,
    delta_g: float,
    R_bar: Optional[float],
    config: KnockoutConfig,
) -> float:
    """Compute combined score optimized for knockout."""
    score = 0.0

    # Cutting efficiency is paramount
    score += config.weight_cut_efficiency * cut_efficiency

    # Frameshift probability
    score += 0.8 * frameshift_prob

    # Specificity
    score += config.weight_mit * (mit / 100.0)
    score += config.weight_cfd * (cfd / 100.0)

    # Binding energy
    dg_score = min(1.0, max(0, (-delta_g) / 25.0))
    score += config.weight_delta_g * dg_score

    # Coherence
    if R_bar is not None:
        score += config.weight_coherence * R_bar

    return score


def validate_knockout_guide(guide_seq: str) -> Dict[str, Any]:
    """
    Validate a single guide for knockout application.

    Args:
        guide_seq: Guide sequence to validate.

    Returns:
        Validation results dictionary.
    """
    guide_seq = guide_seq.upper()

    gc = gc_content(guide_seq)
    cut_eff = cut_efficiency_score(guide_seq)
    mit = mit_specificity_score(guide_seq)
    R_bar = _compute_guide_coherence(guide_seq)
    repair = repair_pathway_prediction(guide_seq)

    warnings = []
    if gc < 0.35 or gc > 0.75:
        warnings.append(f"Non-optimal GC content ({gc:.0%})")
    if cut_eff < 0.3:
        warnings.append(f"Low cutting efficiency ({cut_eff:.0%})")
    if mit < 50:
        warnings.append(f"Low specificity (MIT={mit:.0f})")
    if 'TTTT' in guide_seq:
        warnings.append("Contains poly-T (termination risk)")
    if repair['NHEJ'] < 0.7:
        warnings.append("Low NHEJ probability (knockout less efficient)")

    return {
        'sequence': guide_seq,
        'gc': gc,
        'cut_efficiency': cut_eff,
        'mit_score': mit,
        'repair_pathway': repair,
        'coherence_R': R_bar,
        'go_no_go': go_no_go(R_bar),
        'warnings': warnings,
        'suitable_for_knockout': len(warnings) == 0,
    }
