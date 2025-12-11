"""
PhaseLab Coherence: IR coherence metrics for simulation reliability.

The core metrics from Informational Relativity:
- R̄ (R-bar): Order parameter / coherence score
- V_φ (V-phi): Phase variance
- Relationship: R̄ = exp(-V_φ/2)

These metrics assess whether a quantum or dynamical simulation is reliable.
"""

import numpy as np
from typing import Union, Optional, Tuple, List
from .constants import E_MINUS_2


def phase_variance(phases: np.ndarray) -> float:
    """
    Compute circular phase variance V_φ from an array of phases.

    V_φ = -2 * ln(R̄) where R̄ is the Kuramoto order parameter.

    Args:
        phases: Array of phase angles in radians.

    Returns:
        Phase variance V_φ (non-negative).
    """
    if len(phases) == 0:
        return np.inf

    # Compute complex order parameter
    z = np.mean(np.exp(1j * phases))
    R_bar = np.abs(z)

    # Avoid log(0)
    if R_bar < 1e-10:
        return np.inf

    V_phi = -2.0 * np.log(R_bar)
    return max(0.0, V_phi)


def coherence_score(
    data: Union[np.ndarray, List[float], float],
    mode: str = "auto"
) -> float:
    """
    Compute the coherence score R̄ from various input types.

    Args:
        data: Can be:
            - Array of phases → compute Kuramoto R̄
            - Array of expectation values → compute consistency R̄
            - Single V_φ value → compute R̄ = exp(-V_φ/2)
            - Statevector → compute from amplitudes
        mode: "phases", "expectations", "variance", or "auto"

    Returns:
        Coherence score R̄ in [0, 1].
    """
    data = np.asarray(data)

    # Auto-detect mode
    if mode == "auto":
        if data.ndim == 0:  # scalar
            mode = "variance"
        elif np.all(np.isreal(data)) and np.max(np.abs(data)) <= 2 * np.pi + 0.1:
            # Looks like phases
            mode = "phases"
        else:
            mode = "expectations"

    if mode == "variance":
        # data is V_φ
        V_phi = float(data)
        return np.exp(-V_phi / 2.0)

    elif mode == "phases":
        # Kuramoto order parameter
        z = np.mean(np.exp(1j * data))
        return float(np.abs(z))

    elif mode == "expectations":
        # For expectation values, use variance-based coherence
        # Higher consistency → higher coherence
        if len(data) < 2:
            return 1.0
        variance = np.var(data)
        # Map variance to coherence (heuristic)
        return float(np.exp(-variance / 2.0))

    else:
        raise ValueError(f"Unknown mode: {mode}")


def go_no_go(
    R_bar: float,
    threshold: float = E_MINUS_2
) -> str:
    """
    Determine GO/NO-GO classification based on coherence.

    The e^-2 threshold is the fundamental boundary from IR theory.
    Below this, simulations are considered unreliable.

    Args:
        R_bar: Coherence score.
        threshold: GO/NO-GO boundary (default: e^-2 ≈ 0.135).

    Returns:
        "GO" if R̄ > threshold, else "NO-GO".
    """
    if R_bar is None:
        return "NO-GO"
    return "GO" if R_bar > threshold else "NO-GO"


def classify_coherence(R_bar: float) -> str:
    """
    Classify coherence level into human-readable categories.

    Args:
        R_bar: Coherence score.

    Returns:
        Classification string.
    """
    if R_bar is None:
        return "UNKNOWN"
    if R_bar >= 0.9:
        return "EXCELLENT"
    elif R_bar >= 0.7:
        return "GOOD"
    elif R_bar >= E_MINUS_2:
        return "MARGINAL"
    else:
        return "UNRELIABLE"


def compare_sim_hardware(
    sim_R: float,
    hw_R: float,
    tolerance: float = 0.05
) -> Tuple[float, str]:
    """
    Compare simulator and hardware coherence values.

    Args:
        sim_R: Simulator coherence.
        hw_R: Hardware coherence.
        tolerance: Acceptable difference for "EXCELLENT" agreement.

    Returns:
        (difference, agreement_level)
    """
    diff = abs(hw_R - sim_R)

    if diff <= tolerance:
        agreement = "EXCELLENT"
    elif diff <= 2 * tolerance:
        agreement = "GOOD"
    elif diff <= 4 * tolerance:
        agreement = "MODERATE"
    else:
        agreement = "POOR"

    return diff, agreement


def ensemble_coherence(
    coherence_values: List[float],
    weights: Optional[List[float]] = None
) -> float:
    """
    Compute ensemble coherence from multiple measurements.

    Args:
        coherence_values: List of R̄ values.
        weights: Optional weights for weighted average.

    Returns:
        Ensemble coherence score.
    """
    if not coherence_values:
        return 0.0

    values = np.array(coherence_values)

    if weights is None:
        return float(np.mean(values))
    else:
        weights = np.array(weights)
        weights = weights / np.sum(weights)
        return float(np.sum(values * weights))
