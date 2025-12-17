"""
PhaseLab CRISPR Pipeline: End-to-end guide RNA design.

High-level API for designing CRISPRa/CRISPRi guide RNAs with
multi-layer validation.

IMPORTANT (v1.0.0):
Guide-sequence coherence has been DEPRECATED based on E200-E211 experimental
validation showing r ≈ 0 correlation with outcomes. The validated approach
from E213-E216 uses SPATIAL COHERENCE of the response landscape, not probe
coherence.

For spatial coherence analysis, use:
- phaselab.spatial.regulatory for region classification
- phaselab.spatial.targeting for guide placement within stable regions

The guide-sequence coherence option is retained only for research comparison
and should NOT be used for primary guide ranking.
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
    poly_t_penalty,
    is_repeat_region,
    u6_compatibility_check,
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

    # U6/Pol III compatibility (v0.9.1+)
    filter_poly_t: bool = True  # Exclude guides with TTTT start
    filter_repeats: bool = True  # Exclude guides in repeat regions
    poly_t_threshold: int = 4  # TTTT or longer triggers filter

    # Scoring weights for combined score
    weight_mit: float = 1.0
    weight_cfd: float = 1.0
    weight_gc: float = 0.5
    weight_chromatin: float = 0.8
    weight_delta_g: float = 0.3

    # DEPRECATED: Guide-sequence coherence (v1.0.0)
    # E200-E211 showed r ≈ 0 correlation with outcomes.
    # Use spatial coherence from phaselab.spatial instead.
    # This option retained for research comparison only.
    compute_guide_coherence: bool = False  # DEPRECATED - default OFF
    weight_coherence: float = 0.0  # DEPRECATED - weight set to 0

    # Legacy compatibility (maps to compute_guide_coherence)
    compute_coherence: bool = False  # DEPRECATED alias
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
        warnings = []

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

        # U6/Pol III compatibility check (v0.9.1+)
        is_poly_t, poly_t_reason = poly_t_penalty(guide_seq, threshold=config.poly_t_threshold)
        if is_poly_t:
            if config.filter_poly_t:
                continue  # Hard filter
            else:
                warnings.append(poly_t_reason)

        # Repeat region check (v0.9.1+)
        is_repeat, repeat_reason = is_repeat_region(guide_seq)
        if is_repeat:
            if config.filter_repeats:
                continue  # Hard filter
            else:
                warnings.append(repeat_reason)

        # Full U6 compatibility check
        u6_compat, u6_warnings = u6_compatibility_check(guide_seq)
        warnings.extend(u6_warnings)

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

        # DEPRECATED: Guide-sequence coherence (v1.0.0)
        # E200-E211 showed this does NOT predict outcomes (r ≈ 0).
        # Use phaselab.spatial for region-based coherence instead.
        R_bar = None
        go_status = None
        if config.compute_guide_coherence or config.compute_coherence:
            import warnings
            warnings.warn(
                "Guide-sequence coherence is DEPRECATED (v1.0.0). "
                "E200-E211 showed r ≈ 0 correlation with outcomes. "
                "Use phaselab.spatial for region-based spatial coherence instead.",
                DeprecationWarning,
                stacklevel=2,
            )
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
            'u6_compatible': u6_compat,
            'warnings': warnings if warnings else None,
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
        'u6_compatible', 'warnings', 'combined_score',
    ])


def _compute_guide_coherence(guide_seq: str) -> float:
    """
    DEPRECATED: Compute IR coherence for a guide sequence.

    WARNING (v1.0.0): Guide-sequence coherence does NOT predict outcomes.
    E200-E211 experiments showed r ≈ 0 correlation between guide-sequence
    coherence and experimental results.

    The validated approach from E213-E216 uses SPATIAL COHERENCE of the
    response landscape. Use phaselab.spatial module instead:

        from phaselab.spatial import classify_regulatory_regions
        regions = classify_regulatory_regions(landscape)
        stable_regions = [r for r in regions if r.is_safe]

    This function is retained for research comparison only.

    Args:
        guide_seq: 20bp guide sequence.

    Returns:
        Coherence R̄ value (NOT predictive of outcomes).
    """
    from .coherence_utils import compute_guide_coherence
    return compute_guide_coherence(guide_seq, use_atlas_q=True)


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


def validate_guide(
    guide_seq: str,
    compute_deprecated_coherence: bool = False,
) -> Dict[str, Any]:
    """
    Quick validation of a single guide sequence.

    Includes U6/Pol III compatibility and repeat region checks (v0.9.1+).

    NOTE (v1.0.0): Guide-sequence coherence is now DEPRECATED.
    Use phaselab.spatial for region-based spatial coherence instead.
    The coherence_R field is only computed if compute_deprecated_coherence=True.

    Args:
        guide_seq: Guide sequence to validate.
        compute_deprecated_coherence: If True, compute guide-sequence coherence
            (DEPRECATED - does not predict outcomes).

    Returns:
        Dictionary with validation results.
    """
    guide_seq = guide_seq.upper()

    gc = gc_content(guide_seq)
    homo = max_homopolymer_run(guide_seq)
    complexity = sequence_complexity(guide_seq)
    delta_g = delta_g_santalucia(guide_seq)
    mit = mit_specificity_score(guide_seq)

    # DEPRECATED: Guide-sequence coherence
    R_bar = None
    if compute_deprecated_coherence:
        import warnings
        warnings.warn(
            "Guide-sequence coherence is DEPRECATED (v1.0.0). "
            "Use phaselab.spatial for validated region-based coherence.",
            DeprecationWarning,
            stacklevel=2,
        )
        R_bar = _compute_guide_coherence(guide_seq)

    # U6/Pol III compatibility (v0.9.1+)
    u6_compat, u6_warnings = u6_compatibility_check(guide_seq)
    is_repeat, repeat_reason = is_repeat_region(guide_seq)

    warnings = []

    # Sequence quality warnings
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

    # U6 compatibility warnings (critical for standard delivery)
    warnings.extend(u6_warnings)

    # Repeat region warning
    if is_repeat:
        warnings.append(f"REPEAT REGION: {repeat_reason}")

    # Determine exclusion status
    exclusions = []
    if not u6_compat:
        exclusions.append("U6_INCOMPATIBLE")
    if is_repeat:
        exclusions.append("REPEAT_REGION")

    return {
        'sequence': guide_seq,
        'length': len(guide_seq),
        'gc': gc,
        'homopolymer': homo,
        'complexity': complexity,
        'delta_g': delta_g,
        'mit_score': mit,
        'coherence_R': R_bar,  # DEPRECATED: None unless compute_deprecated_coherence=True
        'go_no_go': go_no_go(R_bar) if R_bar else 'N/A',  # DEPRECATED
        'coherence_note': (
            'Use phaselab.spatial for validated spatial coherence'
            if R_bar is None else 'DEPRECATED - does not predict outcomes'
        ),
        'u6_compatible': u6_compat,
        'is_repeat': is_repeat,
        'warnings': warnings,
        'exclusions': exclusions,
        'valid': len(warnings) == 0 and len(exclusions) == 0,
    }
