"""
PhaseLab Protein Dynamics: IR coherence for molecular dynamics simulations.

Applies coherence metrics to assess MD simulation reliability:
- Trajectory convergence
- RMSD stability
- Ensemble sampling quality
"""

import numpy as np
from typing import List, Optional, Tuple, Dict
from dataclasses import dataclass

from ..core.constants import E_MINUS_2
from ..core.coherence import coherence_score, go_no_go


@dataclass
class MDCoherence:
    """Container for MD simulation coherence metrics."""

    R_bar: float
    go_no_go: str

    # Component scores
    rmsd_stability_R: float
    convergence_R: float
    sampling_R: float

    # Metadata
    n_frames: int
    simulation_time_ns: float
    temperature_K: float

    def to_dict(self) -> Dict:
        return {
            "R_bar": self.R_bar,
            "go_no_go": self.go_no_go,
            "rmsd_stability_R": self.rmsd_stability_R,
            "convergence_R": self.convergence_R,
            "sampling_R": self.sampling_R,
            "n_frames": self.n_frames,
            "simulation_time_ns": self.simulation_time_ns,
            "temperature_K": self.temperature_K,
        }


def rmsd_stability(
    rmsd_trajectory: np.ndarray,
    window_size: int = 100,
    equilibration_fraction: float = 0.2
) -> float:
    """
    Compute coherence from RMSD trajectory stability.

    A converged simulation has low RMSD variance after equilibration.
    High variance indicates the simulation hasn't stabilized.

    Args:
        rmsd_trajectory: RMSD values over time (Å)
        window_size: Rolling window for variance calculation
        equilibration_fraction: Fraction to discard as equilibration

    Returns:
        R̄ coherence for RMSD stability
    """
    if len(rmsd_trajectory) < window_size:
        return 0.0

    # Remove equilibration
    start_idx = int(len(rmsd_trajectory) * equilibration_fraction)
    production = rmsd_trajectory[start_idx:]

    if len(production) < window_size:
        return 0.0

    # Calculate rolling variance
    variances = []
    for i in range(len(production) - window_size):
        window = production[i:i + window_size]
        variances.append(np.var(window))

    mean_variance = np.mean(variances)

    # Convert variance to coherence
    # Typical stable RMSD variance ~ 0.1-0.5 Å²
    # Map to coherence scale
    R_rmsd = np.exp(-mean_variance / 0.5)

    return float(np.clip(R_rmsd, 0.0, 1.0))


def trajectory_coherence(
    coordinates: np.ndarray,
    reference_frame: int = -1
) -> float:
    """
    Compute coherence from coordinate trajectory.

    Uses per-atom position variance as a measure of sampling quality.
    High coherence = consistent positions, low = exploring conformations.

    Args:
        coordinates: Shape (n_frames, n_atoms, 3) trajectory
        reference_frame: Frame index for reference (-1 = last)

    Returns:
        R̄ coherence for trajectory
    """
    if coordinates.ndim != 3:
        raise ValueError("Expected shape (n_frames, n_atoms, 3)")

    n_frames, n_atoms, _ = coordinates.shape

    if n_frames < 2:
        return 1.0

    # Per-atom variance across frames
    per_atom_var = np.var(coordinates, axis=0)  # (n_atoms, 3)
    total_var = np.mean(per_atom_var)  # Average variance

    # Convert to coherence
    # Typical MD fluctuation ~ 0.5-2.0 Å² per coordinate
    R_traj = np.exp(-total_var / 2.0)

    return float(np.clip(R_traj, 0.0, 1.0))


def ensemble_convergence(
    observable_trajectory: np.ndarray,
    block_size: int = 100
) -> float:
    """
    Assess ensemble convergence using block averaging.

    Computes variance between block averages - converged simulations
    have low inter-block variance.

    Args:
        observable_trajectory: Time series of observable values
        block_size: Number of frames per block

    Returns:
        R̄ coherence for ensemble convergence
    """
    n = len(observable_trajectory)

    if n < 2 * block_size:
        return 0.0

    # Compute block averages
    n_blocks = n // block_size
    blocks = observable_trajectory[:n_blocks * block_size].reshape(n_blocks, block_size)
    block_means = np.mean(blocks, axis=1)

    # Variance of block means
    block_variance = np.var(block_means)

    # Expected variance decreases as 1/block_size for converged sampling
    overall_variance = np.var(observable_trajectory)

    if overall_variance < 1e-10:
        return 1.0

    # Ratio indicates convergence quality
    # For ideal sampling: block_variance ≈ overall_variance / n_blocks
    expected_block_var = overall_variance / n_blocks
    convergence_ratio = expected_block_var / (block_variance + 1e-10)

    # Map to coherence
    R_conv = np.clip(convergence_ratio, 0.0, 1.0)

    return float(R_conv)


def rmsf_coherence(
    rmsf_values: np.ndarray,
    expected_flexible: Optional[np.ndarray] = None
) -> float:
    """
    Compute coherence from RMSF (root mean square fluctuation).

    Args:
        rmsf_values: Per-residue RMSF values
        expected_flexible: Boolean mask of expected flexible regions

    Returns:
        R̄ coherence for RMSF distribution
    """
    if len(rmsf_values) == 0:
        return 0.0

    # Normalize RMSF
    rmsf_norm = rmsf_values / np.max(rmsf_values)

    if expected_flexible is not None:
        # Check if flexible regions match expectations
        predicted_flexible = rmsf_norm > 0.5
        agreement = np.mean(predicted_flexible == expected_flexible)
        return float(agreement)
    else:
        # Use variance of RMSF as coherence metric
        # Uniform RMSF = low coherence, structured = high
        variance = np.var(rmsf_norm)
        # Higher variance = more structured flexibility pattern
        R_rmsf = 1.0 - np.exp(-variance * 5)

        return float(np.clip(R_rmsf, 0.0, 1.0))


def compute_md_coherence(
    rmsd: Optional[np.ndarray] = None,
    coordinates: Optional[np.ndarray] = None,
    observable: Optional[np.ndarray] = None,
    simulation_time_ns: float = 0.0,
    temperature_K: float = 300.0
) -> MDCoherence:
    """
    Compute comprehensive MD simulation coherence.

    Args:
        rmsd: RMSD trajectory
        coordinates: Full coordinate trajectory
        observable: Any observable time series (energy, etc.)
        simulation_time_ns: Total simulation time
        temperature_K: Simulation temperature

    Returns:
        MDCoherence dataclass
    """
    scores = []

    # RMSD stability
    if rmsd is not None:
        rmsd_R = rmsd_stability(rmsd)
        scores.append(rmsd_R)
        n_frames = len(rmsd)
    else:
        rmsd_R = 0.0
        n_frames = 0

    # Trajectory coherence
    if coordinates is not None:
        traj_R = trajectory_coherence(coordinates)
        scores.append(traj_R)
        n_frames = coordinates.shape[0]
    else:
        traj_R = 0.0

    # Ensemble convergence
    if observable is not None:
        conv_R = ensemble_convergence(observable)
        scores.append(conv_R)
    else:
        conv_R = 0.0

    # Overall coherence
    if scores:
        R_bar = float(np.mean(scores))
    else:
        R_bar = 0.0

    return MDCoherence(
        R_bar=R_bar,
        go_no_go=go_no_go(R_bar),
        rmsd_stability_R=rmsd_R,
        convergence_R=conv_R,
        sampling_R=traj_R,
        n_frames=n_frames,
        simulation_time_ns=simulation_time_ns,
        temperature_K=temperature_K
    )
