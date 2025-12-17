"""
PhaseLab Landscapes: Spatial coherence computation.

Implements the E213 methodology for computing spatial coherence
of response landscapes. This is the core algorithm that enables
prediction of perturbation stability.

Key insight:
    Spatial coherence of nearby responses predicts outcome variance.
    High coherence → low variance → reproducible results.
    Low coherence → high variance → unpredictable results.
"""

import numpy as np
from scipy import stats
from typing import Optional, Tuple
from dataclasses import dataclass

from .core import ResponseLandscape, CoherenceProfile


def compute_spatial_coherence(
    landscape: ResponseLandscape,
    window: int = 50,
    step: int = 1,
    min_points: int = 5,
) -> CoherenceProfile:
    """
    Compute spatial coherence profile using sliding window analysis.

    This is the core E213 algorithm. For each position, we compute the
    coherence (R̄) of responses in a local window, treating response
    values as phase-like observables.

    The key relationship validated in E213:
        Spatial coherence NEGATIVELY correlates with response variance
        (r = -0.24 to -0.50 across 4 datasets, 115,251 sgRNAs)

    Args:
        landscape: Input ResponseLandscape with perturbation-response data.
        window: Size of sliding window (in coordinate units).
        step: Step size for sliding (default 1 for dense profiles).
        min_points: Minimum points required in window for valid coherence.

    Returns:
        CoherenceProfile with spatial coherence at each position.

    Example:
        >>> landscape = ResponseLandscape(coords=positions, responses=logfc)
        >>> profile = compute_spatial_coherence(landscape, window=50)
        >>> print(f"Correlation: {profile.correlation:.3f}")
        >>> if profile.is_validated:
        ...     print("Spatial coherence predicts stability!")
    """
    # Sort landscape by coordinates
    landscape = landscape.sort_by_coords()
    coords = landscape.coords if landscape.coords.ndim == 1 else landscape.coords[:, 0]
    responses = landscape.mean_response

    # Determine window positions
    coord_min, coord_max = landscape.coord_range
    window_centers = np.arange(coord_min + window / 2, coord_max - window / 2, step)

    coherence_values = []
    variance_values = []
    valid_coords = []

    for center in window_centers:
        # Find points within window
        mask = np.abs(coords - center) <= window / 2
        window_responses = responses[mask]

        if len(window_responses) < min_points:
            continue

        # Remove NaN values
        window_responses = window_responses[~np.isnan(window_responses)]
        if len(window_responses) < min_points:
            continue

        # Compute local coherence (R̄)
        # Treat responses as phases by normalizing to [0, 2π]
        local_coherence = _response_coherence(window_responses)

        # Compute local variance
        local_variance = np.var(window_responses)

        coherence_values.append(local_coherence)
        variance_values.append(local_variance)
        valid_coords.append(center)

    coherence_values = np.array(coherence_values)
    variance_values = np.array(variance_values)
    valid_coords = np.array(valid_coords)

    # Compute correlation between coherence and variance
    if len(coherence_values) > 10:
        correlation, p_value = stats.pearsonr(coherence_values, variance_values)
    else:
        correlation, p_value = 0.0, 1.0

    # Estimate variance reduction in high-coherence regions
    variance_reduction = _estimate_variance_reduction(coherence_values, variance_values)

    return CoherenceProfile(
        coords=valid_coords,
        coherence=coherence_values,
        local_variance=variance_values,
        window_size=window,
        correlation=float(correlation),
        p_value=float(p_value),
        variance_reduction_estimate=variance_reduction,
    )


def _response_coherence(responses: np.ndarray) -> float:
    """
    Compute coherence (R̄) from response values.

    Maps responses to phases and computes the mean resultant length.
    This is the core IR metric applied to response variance.

    The mapping treats responses as phases by:
    1. Normalizing to zero mean
    2. Scaling to [-π, π] based on local range
    3. Computing R̄ = |mean(exp(i*φ))|

    Args:
        responses: Array of response values (e.g., logFC).

    Returns:
        R̄ coherence value in [0, 1].
    """
    if len(responses) < 2:
        return 1.0

    # Normalize responses to phases
    responses = responses - np.mean(responses)

    # Scale to [-π, π]
    resp_range = np.max(responses) - np.min(responses)
    if resp_range > 0:
        phases = responses / resp_range * 2 * np.pi
    else:
        # All same value = perfect coherence
        return 1.0

    # Compute R̄ = |mean(e^(iφ))|
    phasors = np.exp(1j * phases)
    R_bar = np.abs(np.mean(phasors))

    return float(R_bar)


def _estimate_variance_reduction(
    coherence: np.ndarray,
    variance: np.ndarray,
    high_coherence_threshold: float = 0.7,
) -> float:
    """
    Estimate variance reduction in high-coherence vs low-coherence regions.

    This quantifies the practical benefit of targeting high-coherence regions.
    E213 showed 32-49% variance reduction.

    Args:
        coherence: Array of coherence values.
        variance: Array of variance values.
        high_coherence_threshold: Threshold for "high" coherence.

    Returns:
        Estimated fraction variance reduction (0-1).
    """
    high_mask = coherence > high_coherence_threshold
    low_mask = coherence <= high_coherence_threshold

    if np.sum(high_mask) < 5 or np.sum(low_mask) < 5:
        return 0.0

    high_variance = np.mean(variance[high_mask])
    low_variance = np.mean(variance[low_mask])

    if low_variance > 0:
        reduction = 1 - (high_variance / low_variance)
        return max(0.0, min(1.0, float(reduction)))

    return 0.0


def compute_local_variance(
    landscape: ResponseLandscape,
    window: int = 50,
    step: int = 1,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute local variance profile across the landscape.

    Simpler than full coherence analysis - just computes variance
    in sliding windows.

    Args:
        landscape: Input ResponseLandscape.
        window: Size of sliding window.
        step: Step size for sliding.

    Returns:
        Tuple of (coords, variance_values).
    """
    landscape = landscape.sort_by_coords()
    coords = landscape.coords if landscape.coords.ndim == 1 else landscape.coords[:, 0]
    responses = landscape.mean_response

    coord_min, coord_max = landscape.coord_range
    window_centers = np.arange(coord_min + window / 2, coord_max - window / 2, step)

    variance_values = []
    valid_coords = []

    for center in window_centers:
        mask = np.abs(coords - center) <= window / 2
        window_responses = responses[mask]
        window_responses = window_responses[~np.isnan(window_responses)]

        if len(window_responses) < 3:
            continue

        variance_values.append(np.var(window_responses))
        valid_coords.append(center)

    return np.array(valid_coords), np.array(variance_values)


def coherence_variance_correlation(
    landscape: ResponseLandscape,
    window_sizes: Optional[list] = None,
) -> dict:
    """
    Compute coherence-variance correlation across multiple window sizes.

    E213 showed that window size affects correlation strength
    (r = -0.50 at window=100 vs -0.24 at window=50 for CD69).

    Args:
        landscape: Input ResponseLandscape.
        window_sizes: List of window sizes to test (default: [25, 50, 100, 200]).

    Returns:
        Dictionary mapping window_size → (correlation, p_value).
    """
    if window_sizes is None:
        # Adaptive window sizes based on landscape size
        coord_range = landscape.coord_range[1] - landscape.coord_range[0]
        window_sizes = [
            max(10, int(coord_range * f))
            for f in [0.02, 0.05, 0.1, 0.2]
        ]

    results = {}
    for window in window_sizes:
        profile = compute_spatial_coherence(landscape, window=window)
        results[window] = {
            'correlation': profile.correlation,
            'p_value': profile.p_value,
            'variance_reduction': profile.variance_reduction_estimate,
            'n_positions': profile.n_positions,
        }

    return results


def estimate_variance_reduction(
    landscape: ResponseLandscape,
    window: int = 50,
) -> dict:
    """
    Estimate variance reduction from targeting high-coherence regions.

    This is the key practical metric - how much more reproducible
    will results be if we target stable regions?

    Args:
        landscape: Input ResponseLandscape.
        window: Window size for coherence computation.

    Returns:
        Dictionary with variance reduction estimates and statistics.
    """
    profile = compute_spatial_coherence(landscape, window=window)

    # Quartile-based analysis
    q75 = np.percentile(profile.coherence, 75)
    q25 = np.percentile(profile.coherence, 25)

    high_mask = profile.coherence >= q75
    low_mask = profile.coherence <= q25

    high_variance = np.mean(profile.local_variance[high_mask]) if np.sum(high_mask) > 0 else 0
    low_variance = np.mean(profile.local_variance[low_mask]) if np.sum(low_mask) > 0 else 0

    if low_variance > 0:
        quartile_reduction = 1 - (high_variance / low_variance)
    else:
        quartile_reduction = 0.0

    return {
        'overall_reduction': profile.variance_reduction_estimate,
        'quartile_reduction': max(0.0, quartile_reduction),
        'high_coherence_mean_variance': high_variance,
        'low_coherence_mean_variance': low_variance,
        'correlation': profile.correlation,
        'p_value': profile.p_value,
        'is_validated': profile.is_validated,
        'interpretation': profile.interpretation,
    }


def optimal_window_size(
    landscape: ResponseLandscape,
    candidate_windows: Optional[list] = None,
) -> Tuple[int, dict]:
    """
    Find the optimal window size for coherence analysis.

    The optimal window maximizes the negative correlation between
    coherence and variance (strongest predictive relationship).

    Args:
        landscape: Input ResponseLandscape.
        candidate_windows: List of window sizes to test.

    Returns:
        Tuple of (optimal_window, full_results_dict).
    """
    results = coherence_variance_correlation(landscape, candidate_windows)

    # Find window with strongest negative correlation
    best_window = None
    best_correlation = 0

    for window, data in results.items():
        if data['correlation'] < best_correlation and data['p_value'] < 0.05:
            best_correlation = data['correlation']
            best_window = window

    # Default to middle window if no significant correlation found
    if best_window is None:
        windows = list(results.keys())
        best_window = windows[len(windows) // 2]

    return best_window, results
