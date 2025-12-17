"""
PhaseLab Landscapes: Region classification.

Classifies landscape regions as STABLE, MIXED, AMPLIFYING, or IRRELEVANT
based on spatial coherence analysis.

This is the actionable output of landscape analysis:
- STABLE regions: Safe to perturb, predictable outcomes
- MIXED regions: Context-dependent, validate before use
- AMPLIFYING regions: Avoid - MYC-like behavior where perturbations amplify variance
- IRRELEVANT regions: No meaningful response detected
"""

import numpy as np
from typing import Optional, List, Tuple, Dict, Any
from scipy import ndimage

from .core import (
    ResponseLandscape,
    CoherenceProfile,
    StabilityClass,
    RegionClassification,
)
from .coherence import compute_spatial_coherence


def classify_regions(
    landscape: ResponseLandscape,
    profile: Optional[CoherenceProfile] = None,
    window: int = 50,
    stable_threshold: float = 0.7,
    mixed_threshold: float = 0.4,
    min_region_size: int = 20,
    signal_threshold: float = 0.5,
) -> RegionClassification:
    """
    Classify landscape regions by perturbation stability.

    This is the primary analysis function that produces actionable
    region classifications.

    Classification logic:
    1. Compute spatial coherence profile
    2. Identify regions by coherence level
    3. Check for amplifying behavior (positive coherence-variance correlation)
    4. Merge adjacent regions of same class
    5. Filter by minimum region size

    Args:
        landscape: Input ResponseLandscape.
        profile: Pre-computed CoherenceProfile (computed if not provided).
        window: Window size for coherence computation.
        stable_threshold: Coherence threshold for STABLE classification.
        mixed_threshold: Coherence threshold for MIXED (below this = IRRELEVANT).
        min_region_size: Minimum coordinate span for a valid region.
        signal_threshold: Minimum mean absolute response for relevance.

    Returns:
        RegionClassification with classified regions.

    Example:
        >>> regions = classify_regions(landscape, window=50)
        >>> print(regions.summary())
        >>> for start, end, score in regions.stable_regions:
        ...     print(f"Target region: {start}-{end}")
    """
    # Compute coherence profile if not provided
    if profile is None:
        profile = compute_spatial_coherence(landscape, window=window)

    # Check for amplifying behavior (positive correlation = MYC-like)
    is_amplifying_system = profile.correlation > 0.2 and profile.p_value < 0.05

    # Classify each position
    position_classes = []
    for i, (coord, coh, var) in enumerate(zip(
        profile.coords, profile.coherence, profile.local_variance
    )):
        # Check signal level
        if np.isnan(coh):
            cls = StabilityClass.IRRELEVANT
        elif is_amplifying_system and coh > stable_threshold:
            # In amplifying systems, high coherence = high risk
            cls = StabilityClass.AMPLIFYING
        elif coh >= stable_threshold:
            cls = StabilityClass.STABLE
        elif coh >= mixed_threshold:
            cls = StabilityClass.MIXED
        else:
            cls = StabilityClass.IRRELEVANT

        position_classes.append((coord, cls, coh))

    # Merge into contiguous regions
    regions = _merge_into_regions(position_classes, min_region_size)

    # Build classification params for reproducibility
    params = {
        'window': window,
        'stable_threshold': stable_threshold,
        'mixed_threshold': mixed_threshold,
        'min_region_size': min_region_size,
        'signal_threshold': signal_threshold,
        'is_amplifying_system': is_amplifying_system,
    }

    return RegionClassification(
        regions=regions,
        landscape=landscape,
        profile=profile,
        classification_params=params,
    )


def _merge_into_regions(
    position_classes: List[Tuple[float, StabilityClass, float]],
    min_region_size: int,
) -> List[Tuple[float, float, StabilityClass, float]]:
    """
    Merge adjacent positions of same class into contiguous regions.

    Args:
        position_classes: List of (coord, class, score) tuples.
        min_region_size: Minimum coordinate span for valid region.

    Returns:
        List of (start, end, class, mean_score) tuples.
    """
    if not position_classes:
        return []

    # Sort by coordinate
    position_classes = sorted(position_classes, key=lambda x: x[0])

    regions = []
    current_start = position_classes[0][0]
    current_class = position_classes[0][1]
    current_scores = [position_classes[0][2]]

    for coord, cls, score in position_classes[1:]:
        if cls == current_class:
            # Extend current region
            current_scores.append(score)
        else:
            # Close current region and start new one
            region_end = coord  # Use start of new region as end of current
            region_size = region_end - current_start

            if region_size >= min_region_size:
                mean_score = np.mean(current_scores)
                regions.append((current_start, region_end, current_class, mean_score))

            current_start = coord
            current_class = cls
            current_scores = [score]

    # Don't forget the last region
    if position_classes:
        region_end = position_classes[-1][0]
        region_size = region_end - current_start

        if region_size >= min_region_size:
            mean_score = np.mean(current_scores)
            regions.append((current_start, region_end, current_class, mean_score))

    return regions


def detect_amplifying_regions(
    landscape: ResponseLandscape,
    window: int = 50,
    amplification_threshold: float = 0.2,
) -> List[Tuple[float, float, float]]:
    """
    Specifically detect amplifying (MYC-like) regions.

    These are regions where perturbations amplify variance rather than
    having predictable effects. They should be AVOIDED for therapeutic
    targeting.

    Signs of amplification:
    1. Positive local coherence-variance correlation
    2. High variance despite apparent coherence
    3. Non-linear response patterns

    Args:
        landscape: Input ResponseLandscape.
        window: Window size for analysis.
        amplification_threshold: Positive correlation threshold for amplification.

    Returns:
        List of (start, end, amplification_score) tuples.
    """
    profile = compute_spatial_coherence(landscape, window=window)

    # Global amplification check
    if profile.correlation > amplification_threshold and profile.p_value < 0.05:
        # Entire landscape is amplifying
        coord_range = landscape.coord_range
        return [(coord_range[0], coord_range[1], profile.correlation)]

    # Local amplification detection using sliding correlation
    amplifying_regions = []
    local_window = max(20, len(profile.coords) // 10)

    for i in range(0, len(profile.coords) - local_window, local_window // 2):
        local_coh = profile.coherence[i:i + local_window]
        local_var = profile.local_variance[i:i + local_window]

        if len(local_coh) < 10:
            continue

        from scipy import stats
        local_corr, local_p = stats.pearsonr(local_coh, local_var)

        if local_corr > amplification_threshold and local_p < 0.1:
            start = profile.coords[i]
            end = profile.coords[min(i + local_window - 1, len(profile.coords) - 1)]
            amplifying_regions.append((start, end, local_corr))

    return amplifying_regions


def stability_boundaries(
    classification: RegionClassification,
) -> List[Tuple[float, StabilityClass, StabilityClass]]:
    """
    Find boundaries between stability classes.

    These boundaries are important because:
    - They mark transitions between safe and risky regions
    - Perturbations near boundaries may have unpredictable effects
    - Buffer zones around boundaries are recommended

    Args:
        classification: RegionClassification from classify_regions().

    Returns:
        List of (coordinate, class_before, class_after) tuples.
    """
    boundaries = []

    regions = sorted(classification.regions, key=lambda x: x[0])

    for i in range(len(regions) - 1):
        _, end1, class1, _ = regions[i]
        start2, _, class2, _ = regions[i + 1]

        if class1 != class2:
            boundary_coord = (end1 + start2) / 2
            boundaries.append((boundary_coord, class1, class2))

    return boundaries


def region_summary(
    classification: RegionClassification,
) -> Dict[str, Any]:
    """
    Generate detailed summary statistics for region classification.

    Args:
        classification: RegionClassification from classify_regions().

    Returns:
        Dictionary with summary statistics.
    """
    # Count regions by class
    class_counts = {cls: 0 for cls in StabilityClass}
    class_spans = {cls: 0.0 for cls in StabilityClass}

    for start, end, cls, _ in classification.regions:
        class_counts[cls] += 1
        class_spans[cls] += end - start

    total_span = sum(class_spans.values())

    # Compute fractions
    class_fractions = {
        cls: span / total_span if total_span > 0 else 0
        for cls, span in class_spans.items()
    }

    # Best stable region
    best_stable = None
    if classification.stable_regions:
        best_stable = max(classification.stable_regions, key=lambda x: x[2])

    return {
        'total_regions': classification.n_regions,
        'class_counts': {cls.value: count for cls, count in class_counts.items()},
        'class_spans': {cls.value: span for cls, span in class_spans.items()},
        'class_fractions': {cls.value: frac for cls, frac in class_fractions.items()},
        'total_span': total_span,
        'best_stable_region': best_stable,
        'profile_correlation': classification.profile.correlation,
        'variance_reduction': classification.profile.variance_reduction_estimate,
        'is_validated': classification.profile.is_validated,
    }


def recommend_targets(
    classification: RegionClassification,
    n_targets: int = 5,
    min_score: float = 0.5,
    avoid_boundaries: float = 20,
) -> List[Dict[str, Any]]:
    """
    Recommend specific target coordinates within stable regions.

    Args:
        classification: RegionClassification from classify_regions().
        n_targets: Maximum number of targets to recommend.
        min_score: Minimum coherence score for recommendation.
        avoid_boundaries: Coordinate distance to avoid near region boundaries.

    Returns:
        List of target recommendations with coordinates and scores.
    """
    recommendations = []

    for start, end, score in classification.stable_regions:
        if score < min_score:
            continue

        # Avoid boundaries
        safe_start = start + avoid_boundaries
        safe_end = end - avoid_boundaries

        if safe_end <= safe_start:
            continue

        # Center of region is safest
        center = (safe_start + safe_end) / 2

        recommendations.append({
            'coordinate': center,
            'region_start': start,
            'region_end': end,
            'coherence_score': score,
            'region_size': end - start,
            'recommendation': 'PRIMARY' if score > 0.7 else 'SECONDARY',
        })

    # Sort by score and return top n
    recommendations.sort(key=lambda x: x['coherence_score'], reverse=True)
    return recommendations[:n_targets]
