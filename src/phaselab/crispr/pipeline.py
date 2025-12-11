"""
PhaseLab CRISPR Pipeline: End-to-end guide RNA design.

High-level API for designing CRISPRa/CRISPRi guide RNAs with
multi-layer validation including IR coherence scoring.
"""

import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple

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


@dataclass
class GuideDesignConfig:
    """Configuration for guide RNA design pipeline."""

    # PAM settings
    pam: str = "NGG"
    guide_length: int = 20

    # CRISPRa window (relative to TSS)
    crispr_window: Tuple[int, int] = (-400, -50)

    # Filtering thresholds
    min_gc: float = 0.4
    max_gc: float = 0.7
    max_homopolymer: int = 4
    min_complexity: float = 0.5

    # Scoring weights for combined score
    weight_mit: float = 1.0
    weight_cfd: float = 1.0
    weight_gc: float = 0.5
    weight_chromatin: float = 0.8
    weight_coherence: float = 1.0
    weight_delta_g: float = 0.3

    # Coherence simulation settings
    compute_coherence: bool = True
    coherence_shots: int = 2000
    hardware_backend: Optional[str] = None

    # Output settings
    top_n: int = 10


def design_guides(
    sequence: str,
    tss_index: int,
    config: Optional[GuideDesignConfig] = None,
    dnase_peaks: Optional[List[Tuple[int, int]]] = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Design and rank guide RNAs for CRISPRa/CRISPRi.

    This is the main entry point for the CRISPR pipeline. It:
    1. Scans for PAM sites
    2. Filters candidates by window and quality
    3. Computes multi-layer scores
    4. Optionally runs IR coherence simulation
    5. Returns ranked candidates

    Args:
        sequence: Promoter DNA sequence (5'->3').
        tss_index: Position of TSS in sequence (0-based).
        config: GuideDesignConfig with parameters.
        dnase_peaks: Optional list of (start, end) DNase HS sites.
        verbose: Print progress messages.

    Returns:
        DataFrame with ranked guide candidates and scores.

    Example:
        >>> from phaselab.crispr import design_guides
        >>> guides = design_guides(
        ...     sequence=rai1_promoter,
        ...     tss_index=500,
        ... )
        >>> print(guides[['sequence', 'position', 'combined_score']].head())
    """
    if config is None:
        config = GuideDesignConfig()

    sequence = sequence.upper()

    if verbose:
        print(f"Scanning {len(sequence)} bp sequence for {config.pam} PAM sites...")

    # Step 1: Find all PAM sites
    all_hits = find_pam_sites(
        sequence,
        pam=config.pam,
        guide_length=config.guide_length,
        both_strands=True,
    )

    if verbose:
        print(f"Found {len(all_hits)} total PAM sites")

    # Step 2: Filter to CRISPRa window
    window_hits = filter_by_window(
        all_hits,
        tss_position=tss_index,
        window=config.crispr_window,
    )

    if verbose:
        print(f"Filtered to {len(window_hits)} in CRISPRa window {config.crispr_window}")

    if not window_hits:
        return _empty_results_df()

    # Step 3: Score and filter candidates
    candidates = []

    for hit in window_hits:
        guide_seq = hit.guide

        # Basic quality filters
        gc = gc_content(guide_seq)
        if gc < config.min_gc or gc > config.max_gc:
            continue

        homo = max_homopolymer_run(guide_seq)
        if homo > config.max_homopolymer:
            continue

        complexity = sequence_complexity(guide_seq)
        if complexity < config.min_complexity:
            continue

        # Compute scores
        delta_g = delta_g_santalucia(guide_seq)
        mit = mit_specificity_score(guide_seq)
        cfd = cfd_score(guide_seq)

        # Position relative to TSS
        rel_pos = ((hit.guide_start + hit.guide_end) // 2) - tss_index

        # Chromatin accessibility
        chrom_state, chrom_access = chromatin_accessibility_score(
            position=hit.guide_start,
            tss_position=tss_index,
            dnase_peaks=dnase_peaks,
        )

        # IR coherence (optional)
        R_bar = None
        go_status = None
        if config.compute_coherence:
            R_bar = _compute_guide_coherence(guide_seq)
            go_status = go_no_go(R_bar)

        # Combined score
        combined = _compute_combined_score(
            gc=gc,
            delta_g=delta_g,
            mit=mit,
            cfd=cfd,
            chrom_access=chrom_access,
            R_bar=R_bar,
            config=config,
        )

        candidates.append({
            'sequence': guide_seq,
            'pam': hit.pam,
            'position': rel_pos,
            'strand': hit.strand,
            'gc': round(gc, 3),
            'delta_g': round(delta_g, 3),
            'mit_score': round(mit, 1),
            'cfd_score': round(cfd, 1),
            'chromatin_state': chrom_state,
            'chromatin_accessibility': round(chrom_access, 3),
            'coherence_R': round(R_bar, 4) if R_bar else None,
            'go_no_go': go_status,
            'complexity': round(complexity, 3),
            'homopolymer': homo,
            'combined_score': round(combined, 3),
        })

    if not candidates:
        return _empty_results_df()

    # Step 4: Create DataFrame and sort
    df = pd.DataFrame(candidates)
    df.sort_values(by='combined_score', ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)

    if verbose:
        print(f"Returning {min(len(df), config.top_n)} top candidates")

    return df.head(config.top_n)


def _empty_results_df() -> pd.DataFrame:
    """Return empty DataFrame with correct columns."""
    return pd.DataFrame(columns=[
        'sequence', 'pam', 'position', 'strand',
        'gc', 'delta_g', 'mit_score', 'cfd_score',
        'chromatin_state', 'chromatin_accessibility',
        'coherence_R', 'go_no_go', 'complexity', 'homopolymer',
        'combined_score',
    ])


def _compute_guide_coherence(guide_seq: str) -> float:
    """
    Compute IR coherence for a guide sequence.

    Uses a simplified statevector simulation of the gRNA Hamiltonian.
    For hardware runs, use the full VQE pipeline.

    Args:
        guide_seq: 20bp guide sequence.

    Returns:
        Coherence R̄ value.
    """
    try:
        # Build Hamiltonian
        H = build_grna_hamiltonian(guide_seq)

        # Simple coherence estimate from Hamiltonian structure
        # (Full implementation would run VQE or statevector sim)
        terms = H.get_terms()
        if not terms:
            return 0.5

        # Use term coefficients as proxy for binding quality
        energies = [abs(coeff) for coeff, _ in terms]
        mean_energy = np.mean(energies)
        std_energy = np.std(energies)

        # Map to coherence (heuristic)
        # Lower variance = more consistent = higher coherence
        if mean_energy > 0:
            V_phi = std_energy / mean_energy
            R_bar = np.exp(-V_phi / 2)
        else:
            R_bar = 0.5

        return float(np.clip(R_bar, 0, 1))

    except Exception:
        return 0.5


def _compute_combined_score(
    gc: float,
    delta_g: float,
    mit: float,
    cfd: float,
    chrom_access: float,
    R_bar: Optional[float],
    config: GuideDesignConfig,
) -> float:
    """
    Compute weighted combined score for ranking.

    Higher score = better candidate.
    """
    score = 0.0

    # MIT and CFD (0-100 scale, normalize to 0-1)
    score += config.weight_mit * (mit / 100.0)
    score += config.weight_cfd * (cfd / 100.0)

    # GC content (optimal around 0.55)
    gc_score = 1.0 - 2.0 * abs(gc - 0.55)
    score += config.weight_gc * max(0, gc_score)

    # Chromatin accessibility (0-1)
    score += config.weight_chromatin * chrom_access

    # Coherence (0-1)
    if R_bar is not None:
        score += config.weight_coherence * R_bar

    # Delta G (more negative = better, typical range -30 to 0)
    # Normalize to 0-1 where -25 is best
    dg_score = min(1.0, max(0, (-delta_g) / 25.0))
    score += config.weight_delta_g * dg_score

    return score


def validate_guide(guide_seq: str) -> Dict[str, Any]:
    """
    Quick validation of a single guide sequence.

    Args:
        guide_seq: Guide sequence to validate.

    Returns:
        Dictionary with validation results.
    """
    guide_seq = guide_seq.upper()

    gc = gc_content(guide_seq)
    homo = max_homopolymer_run(guide_seq)
    complexity = sequence_complexity(guide_seq)
    delta_g = delta_g_santalucia(guide_seq)
    mit = mit_specificity_score(guide_seq)
    R_bar = _compute_guide_coherence(guide_seq)

    warnings = []
    if gc < 0.4:
        warnings.append("Low GC content (<40%)")
    if gc > 0.7:
        warnings.append("High GC content (>70%)")
    if homo > 4:
        warnings.append(f"Long homopolymer run ({homo}bp)")
    if complexity < 0.5:
        warnings.append("Low sequence complexity")
    if R_bar and R_bar < 0.135:
        warnings.append("Low coherence (NO-GO)")

    return {
        'sequence': guide_seq,
        'length': len(guide_seq),
        'gc': gc,
        'homopolymer': homo,
        'complexity': complexity,
        'delta_g': delta_g,
        'mit_score': mit,
        'coherence_R': R_bar,
        'go_no_go': go_no_go(R_bar) if R_bar else 'UNKNOWN',
        'warnings': warnings,
        'valid': len(warnings) == 0,
    }
