"""
PhaseLab CRISPRi: CRISPR interference for transcriptional repression.

Implements:
- TSS-proximal guide design for dCas9-KRAB repression
- Steric hindrance modeling
- Repression efficiency prediction
- Off-target risk assessment
- IR coherence validation
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass

from .pam_scan import find_pam_sites, filter_by_window, PAMHit
from .scoring import (
    gc_content,
    delta_g_santalucia,
    mit_specificity_score,
    cfd_score,
    max_homopolymer_run,
    sequence_complexity,
    chromatin_accessibility_score,
)
from ..core.coherence import coherence_score, go_no_go
from ..core.hamiltonians import build_grna_hamiltonian


# CRISPRi optimal position ranges (relative to TSS)
# Based on Gilbert et al. 2014 and Horlbeck et al. 2016
CRISPRI_OPTIMAL_RANGES = {
    'template_strand': (-50, +300),    # Optimal for blocking elongation
    'nontemplate_strand': (-50, +300), # Also effective
    'promoter': (-400, -50),           # For blocking initiation
}

# Position-specific repression efficiency weights
REPRESSION_POSITION_WEIGHTS = {
    # Template strand, relative to TSS
    'template': {
        (-50, 0): 0.6,     # Promoter proximal
        (0, 50): 0.9,      # TSS - highly effective
        (50, 150): 1.0,    # Most effective region
        (150, 300): 0.8,   # Still effective
        (300, 500): 0.5,   # Less effective
    },
    # Non-template strand
    'nontemplate': {
        (-50, 0): 0.5,
        (0, 50): 0.7,
        (50, 150): 0.85,
        (150, 300): 0.7,
        (300, 500): 0.4,
    },
}


@dataclass
class CRISPRiConfig:
    """Configuration for CRISPRi guide design."""

    # PAM settings
    pam: str = "NGG"
    guide_length: int = 20

    # CRISPRi window (relative to TSS)
    # Optimal: TSS to +150bp for steric blocking
    crispri_window: Tuple[int, int] = (-50, +300)

    # Filtering thresholds
    min_gc: float = 0.35
    max_gc: float = 0.75
    max_homopolymer: int = 4
    min_complexity: float = 0.5

    # Repression efficiency threshold
    min_repression_efficiency: float = 0.3

    # Scoring weights
    weight_mit: float = 1.0
    weight_cfd: float = 1.0
    weight_repression: float = 1.5
    weight_position: float = 1.2
    weight_chromatin: float = 0.8
    weight_delta_g: float = 0.3
    weight_coherence: float = 1.0

    # Coherence settings
    compute_coherence: bool = True
    coherence_shots: int = 2000

    # Repressor type affects optimal window
    repressor: str = "KRAB"  # Options: "KRAB", "MeCP2", "dCas9_only"

    # Output
    top_n: int = 10


def repression_efficiency_score(
    guide_seq: str,
    position: int,
    strand: str,
    repressor: str = "KRAB",
) -> float:
    """
    Predict repression efficiency for CRISPRi.

    Based on position relative to TSS and strand orientation.

    Args:
        guide_seq: Guide sequence.
        position: Position relative to TSS.
        strand: "+" (template) or "-" (non-template).
        repressor: Repressor domain type.

    Returns:
        Repression efficiency score (0.0 to 1.0).
    """
    guide_seq = guide_seq.upper()

    # Base score from position
    strand_key = 'template' if strand == '+' else 'nontemplate'
    position_score = 0.3  # Default

    for (start, end), score in REPRESSION_POSITION_WEIGHTS[strand_key].items():
        if start <= position < end:
            position_score = score
            break

    # Adjust for repressor type
    repressor_multiplier = {
        'KRAB': 1.0,           # Standard, effective
        'MeCP2': 0.9,          # Good for sustained repression
        'dCas9_only': 0.6,     # Steric only, less effective
        'SID4X': 0.95,         # Strong repressor
    }.get(repressor, 1.0)

    # GC content modifier
    gc = gc_content(guide_seq)
    gc_modifier = 1.0
    if gc < 0.35 or gc > 0.70:
        gc_modifier = 0.85

    # Complexity modifier
    complexity = sequence_complexity(guide_seq)
    complexity_modifier = min(1.0, complexity + 0.3)

    efficiency = position_score * repressor_multiplier * gc_modifier * complexity_modifier

    return float(np.clip(efficiency, 0, 1))


def steric_hindrance_score(
    position: int,
    strand: str,
) -> float:
    """
    Calculate steric hindrance potential.

    Guides that block RNA polymerase progression are most effective.

    Args:
        position: Position relative to TSS.
        strand: Strand orientation.

    Returns:
        Steric hindrance score (0.0 to 1.0).
    """
    # Maximum hindrance near TSS on template strand
    if strand == '+':  # Template strand
        if 0 <= position <= 100:
            # Peak at ~50bp downstream
            return 1.0 - abs(position - 50) / 100
        elif -50 <= position < 0:
            return 0.7 + position / 250  # Approaching TSS
        elif 100 < position <= 300:
            return max(0.3, 1.0 - (position - 100) / 400)
        else:
            return 0.2
    else:  # Non-template strand
        if 0 <= position <= 100:
            return 0.8 - abs(position - 50) / 125
        elif -50 <= position < 0:
            return 0.5 + position / 250
        elif 100 < position <= 300:
            return max(0.2, 0.8 - (position - 100) / 400)
        else:
            return 0.15


def design_crispri_guides(
    sequence: str,
    tss_index: int,
    config: Optional[CRISPRiConfig] = None,
    dnase_peaks: Optional[List[Tuple[int, int]]] = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Design guide RNAs for CRISPRi transcriptional repression.

    Optimized for dCas9-KRAB mediated gene silencing with
    guides targeting the TSS-proximal region.

    Args:
        sequence: Promoter/gene DNA sequence.
        tss_index: Position of TSS in sequence (0-based).
        config: CRISPRiConfig with parameters.
        dnase_peaks: Optional DNase hypersensitive sites.
        verbose: Print progress messages.

    Returns:
        DataFrame with ranked CRISPRi guide candidates.

    Example:
        >>> from phaselab.crispr import design_crispri_guides
        >>> guides = design_crispri_guides(
        ...     sequence=promoter_seq,
        ...     tss_index=500,
        ... )
        >>> print(guides[['sequence', 'position', 'repression_efficiency']].head())
    """
    if config is None:
        config = CRISPRiConfig()

    sequence = sequence.upper()

    if verbose:
        print(f"Designing CRISPRi guides for {len(sequence)} bp sequence...")
        print(f"Repressor: {config.repressor}")

    # Find PAM sites
    all_hits = find_pam_sites(
        sequence,
        pam=config.pam,
        guide_length=config.guide_length,
        both_strands=True,
    )

    if verbose:
        print(f"Found {len(all_hits)} total PAM sites")

    # Filter to CRISPRi window
    window_hits = filter_by_window(
        all_hits,
        tss_position=tss_index,
        window=config.crispri_window,
    )

    if verbose:
        print(f"Filtered to {len(window_hits)} in CRISPRi window {config.crispri_window}")

    if not window_hits:
        return _empty_crispri_df()

    # Score candidates
    candidates = []

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

        # Position relative to TSS
        rel_pos = ((hit.guide_start + hit.guide_end) // 2) - tss_index

        # Repression efficiency
        rep_eff = repression_efficiency_score(
            guide_seq=guide_seq,
            position=rel_pos,
            strand=hit.strand,
            repressor=config.repressor,
        )

        if rep_eff < config.min_repression_efficiency:
            continue

        # Steric hindrance
        steric = steric_hindrance_score(rel_pos, hit.strand)

        # Standard scores
        delta_g = delta_g_santalucia(guide_seq)
        mit = mit_specificity_score(guide_seq)
        cfd = cfd_score(guide_seq)

        # Chromatin accessibility
        chrom_state, chrom_access = chromatin_accessibility_score(
            position=hit.guide_start,
            tss_position=tss_index,
            dnase_peaks=dnase_peaks,
        )

        # IR coherence
        R_bar = None
        go_status = None
        if config.compute_coherence:
            R_bar = _compute_guide_coherence(guide_seq)
            go_status = go_no_go(R_bar)

        # Combined score
        combined = _compute_crispri_score(
            repression_efficiency=rep_eff,
            steric_hindrance=steric,
            mit=mit,
            cfd=cfd,
            chrom_access=chrom_access,
            delta_g=delta_g,
            R_bar=R_bar,
            config=config,
        )

        candidates.append({
            'sequence': guide_seq,
            'pam': hit.pam,
            'position': rel_pos,
            'strand': hit.strand,
            'gc': round(gc, 3),
            'repression_efficiency': round(rep_eff, 3),
            'steric_hindrance': round(steric, 3),
            'delta_g': round(delta_g, 3),
            'mit_score': round(mit, 1),
            'cfd_score': round(cfd, 1),
            'chromatin_state': chrom_state,
            'chromatin_accessibility': round(chrom_access, 3),
            'coherence_R': round(R_bar, 4) if R_bar else None,
            'go_no_go': go_status,
            'combined_score': round(combined, 3),
        })

    if not candidates:
        return _empty_crispri_df()

    df = pd.DataFrame(candidates)
    df.sort_values(by='combined_score', ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)

    if verbose:
        print(f"Returning {min(len(df), config.top_n)} top candidates")

    return df.head(config.top_n)


def _empty_crispri_df() -> pd.DataFrame:
    """Return empty DataFrame with CRISPRi columns."""
    return pd.DataFrame(columns=[
        'sequence', 'pam', 'position', 'strand',
        'gc', 'repression_efficiency', 'steric_hindrance',
        'delta_g', 'mit_score', 'cfd_score',
        'chromatin_state', 'chromatin_accessibility',
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


def _compute_crispri_score(
    repression_efficiency: float,
    steric_hindrance: float,
    mit: float,
    cfd: float,
    chrom_access: float,
    delta_g: float,
    R_bar: Optional[float],
    config: CRISPRiConfig,
) -> float:
    """Compute combined score for CRISPRi ranking."""
    score = 0.0

    # Repression efficiency is key
    score += config.weight_repression * repression_efficiency

    # Steric hindrance
    score += config.weight_position * steric_hindrance

    # Specificity
    score += config.weight_mit * (mit / 100.0)
    score += config.weight_cfd * (cfd / 100.0)

    # Chromatin accessibility
    score += config.weight_chromatin * chrom_access

    # Binding energy
    dg_score = min(1.0, max(0, (-delta_g) / 25.0))
    score += config.weight_delta_g * dg_score

    # Coherence
    if R_bar is not None:
        score += config.weight_coherence * R_bar

    return score


def validate_crispri_guide(
    guide_seq: str,
    position: int,
    strand: str = '+',
    repressor: str = "KRAB",
) -> Dict[str, Any]:
    """
    Validate a single guide for CRISPRi application.

    Args:
        guide_seq: Guide sequence to validate.
        position: Position relative to TSS.
        strand: Strand orientation.
        repressor: Repressor domain type.

    Returns:
        Validation results dictionary.
    """
    guide_seq = guide_seq.upper()

    gc = gc_content(guide_seq)
    rep_eff = repression_efficiency_score(guide_seq, position, strand, repressor)
    steric = steric_hindrance_score(position, strand)
    mit = mit_specificity_score(guide_seq)
    R_bar = _compute_guide_coherence(guide_seq)

    warnings = []
    if gc < 0.35 or gc > 0.75:
        warnings.append(f"Non-optimal GC content ({gc:.0%})")
    if rep_eff < 0.3:
        warnings.append(f"Low repression efficiency ({rep_eff:.0%})")
    if steric < 0.3:
        warnings.append(f"Low steric hindrance score ({steric:.2f})")
    if position > 300:
        warnings.append(f"Position {position}bp is downstream of optimal window")
    if position < -200:
        warnings.append(f"Position {position}bp is upstream of optimal window")
    if mit < 50:
        warnings.append(f"Low specificity (MIT={mit:.0f})")

    return {
        'sequence': guide_seq,
        'position': position,
        'strand': strand,
        'gc': gc,
        'repression_efficiency': rep_eff,
        'steric_hindrance': steric,
        'mit_score': mit,
        'coherence_R': R_bar,
        'go_no_go': go_no_go(R_bar),
        'warnings': warnings,
        'suitable_for_crispri': len(warnings) == 0,
    }


# Aliases for compatibility
design_interference_guides = design_crispri_guides
validate_interference_guide = validate_crispri_guide
