"""
PhaseLab Quantum Coherence: Real circular statistics via ATLAS-Q.

Replaces the heuristic coherence calculation with proper circular
statistics from Pauli expectation values:

1. Map ⟨P⟩ ∈ [-1, 1] → φ = arccos(⟨P⟩) ∈ [0, π]
2. Compute mean phasor: ⟨e^(iφ)⟩
3. R̄ = |⟨e^(iφ)⟩| (mean resultant length)
4. V_φ = -2 ln(R̄) (circular variance)

This is the validated IR coherence metric from ATLAS-Q hardware runs.
"""

from dataclasses import dataclass
from typing import List, Optional, Tuple, Union

import numpy as np

from . import is_atlas_q_available
from ..core.constants import E_MINUS_2


@dataclass
class CoherenceResult:
    """
    Result from coherence calculation.

    Attributes:
        R_bar: Mean resultant length [0, 1]
        V_phi: Circular variance [0, ∞)
        is_go: Whether R̄ > e^-2 threshold
        n_measurements: Number of measurements used
        method: Calculation method used
    """
    R_bar: float
    V_phi: float
    is_go: bool
    n_measurements: int
    method: str

    def __str__(self) -> str:
        status = "GO" if self.is_go else "NO-GO"
        return f"Coherence(R̄={self.R_bar:.4f}, V_φ={self.V_phi:.4f}, {status})"


def compute_coherence_from_expectations(
    expectation_values: Union[np.ndarray, List[float]],
    e2_threshold: float = E_MINUS_2,
    use_atlas_q: bool = True,
) -> CoherenceResult:
    """
    Compute coherence from Pauli expectation values using circular statistics.

    This is the proper IR coherence calculation:
    1. Convert ⟨P⟩ → phases via arccos
    2. Compute mean resultant length R̄
    3. Derive circular variance V_φ = -2 ln(R̄)

    Args:
        expectation_values: Array of Pauli expectations ⟨P⟩ ∈ [-1, 1]
        e2_threshold: GO/NO-GO threshold (default: e^-2 ≈ 0.135)
        use_atlas_q: Use ATLAS-Q backend if available

    Returns:
        CoherenceResult with R̄, V_φ, and classification

    Example:
        >>> expectations = np.array([0.9, 0.85, 0.88, 0.92])
        >>> result = compute_coherence_from_expectations(expectations)
        >>> print(f"R̄ = {result.R_bar:.4f}, Status: {'GO' if result.is_go else 'NO-GO'}")
    """
    expectation_values = np.asarray(expectation_values)

    if expectation_values.size == 0:
        return CoherenceResult(
            R_bar=0.0,
            V_phi=np.inf,
            is_go=False,
            n_measurements=0,
            method="empty",
        )

    # Try ATLAS-Q backend
    if use_atlas_q and is_atlas_q_available():
        try:
            from atlas_q import compute_coherence as atlas_compute

            result = atlas_compute(expectation_values, e2_threshold)

            return CoherenceResult(
                R_bar=result.R_bar,
                V_phi=result.V_phi,
                is_go=result.is_above_e2_boundary,
                n_measurements=result.n_measurements,
                method="atlas_q",
            )
        except Exception:
            pass

    # Fallback: Native implementation
    return _compute_coherence_native(expectation_values, e2_threshold)


def _compute_coherence_native(
    expectation_values: np.ndarray,
    e2_threshold: float,
) -> CoherenceResult:
    """
    Native coherence calculation (fallback when ATLAS-Q unavailable).

    Uses same circular statistics algorithm as ATLAS-Q.
    """
    # Clip to valid range (handle numerical precision)
    expectation_values = np.clip(expectation_values, -1.0, 1.0)

    # Convert expectations to phases: ⟨P⟩ = cos(φ)
    phases = np.arccos(expectation_values)

    # Compute mean phasor
    phasors = np.exp(1j * phases)
    mean_phasor = np.mean(phasors)
    R_bar = float(np.abs(mean_phasor))

    # Compute circular variance
    if R_bar > 1e-10:
        V_phi = -2.0 * np.log(R_bar)
    else:
        V_phi = np.inf

    # Classify
    is_go = R_bar > e2_threshold

    return CoherenceResult(
        R_bar=R_bar,
        V_phi=V_phi,
        is_go=is_go,
        n_measurements=len(expectation_values),
        method="native",
    )


def compute_coherence_from_phases(
    phases: Union[np.ndarray, List[float]],
    e2_threshold: float = E_MINUS_2,
) -> CoherenceResult:
    """
    Compute coherence directly from phase angles.

    Uses the Kuramoto order parameter: R̄ = |⟨e^(iφ)⟩|

    Args:
        phases: Array of phase angles in radians
        e2_threshold: GO/NO-GO threshold

    Returns:
        CoherenceResult with R̄, V_φ, and classification
    """
    phases = np.asarray(phases)

    if phases.size == 0:
        return CoherenceResult(
            R_bar=0.0,
            V_phi=np.inf,
            is_go=False,
            n_measurements=0,
            method="phases_empty",
        )

    # Compute order parameter
    z = np.mean(np.exp(1j * phases))
    R_bar = float(np.abs(z))

    # Circular variance
    if R_bar > 1e-10:
        V_phi = -2.0 * np.log(R_bar)
    else:
        V_phi = np.inf

    is_go = R_bar > e2_threshold

    return CoherenceResult(
        R_bar=R_bar,
        V_phi=V_phi,
        is_go=is_go,
        n_measurements=len(phases),
        method="phases",
    )


def compute_coherence_from_hamiltonian(
    coefficients: np.ndarray,
    use_atlas_q: bool = True,
) -> CoherenceResult:
    """
    Estimate coherence from Hamiltonian structure (heuristic).

    This is a simplified method that estimates coherence from
    Hamiltonian coefficient variance without actual measurement.
    Use compute_coherence_from_expectations() for real coherence.

    Args:
        coefficients: Hamiltonian term coefficients
        use_atlas_q: Use ATLAS-Q for improved estimation

    Returns:
        CoherenceResult (heuristic estimate)
    """
    coefficients = np.asarray(coefficients)

    if coefficients.size == 0:
        return CoherenceResult(
            R_bar=0.5,
            V_phi=np.log(4),  # R̄ = 0.5 → V_φ ≈ 1.39
            is_go=True,
            n_measurements=0,
            method="heuristic_empty",
        )

    # Heuristic: Use coefficient variance as proxy
    energies = np.abs(coefficients)
    mean_energy = np.mean(energies)
    std_energy = np.std(energies)

    if mean_energy > 0:
        # Coefficient of variation as proxy for phase variance
        V_phi = std_energy / mean_energy
        R_bar = np.exp(-V_phi / 2)
    else:
        R_bar = 0.5
        V_phi = np.log(4)

    R_bar = float(np.clip(R_bar, 0, 1))
    is_go = R_bar > E_MINUS_2

    return CoherenceResult(
        R_bar=R_bar,
        V_phi=V_phi,
        is_go=is_go,
        n_measurements=len(coefficients),
        method="heuristic",
    )


def compare_coherence(
    sim_result: CoherenceResult,
    hw_result: CoherenceResult,
    tolerance: float = 0.05,
) -> Tuple[float, str]:
    """
    Compare simulator and hardware coherence values.

    Args:
        sim_result: Simulator coherence result
        hw_result: Hardware coherence result
        tolerance: Acceptable difference for "EXCELLENT"

    Returns:
        (absolute_difference, agreement_level)
    """
    diff = abs(hw_result.R_bar - sim_result.R_bar)

    if diff <= tolerance:
        agreement = "EXCELLENT"
    elif diff <= 2 * tolerance:
        agreement = "GOOD"
    elif diff <= 4 * tolerance:
        agreement = "MODERATE"
    else:
        agreement = "POOR"

    return diff, agreement
