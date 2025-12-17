"""
PhaseLab Transcriptional Phase Alignment: IR-native guide scoring.

Path B from the CRISPRa Breakthrough Design:
- Pure IR framework - no quantum chemistry needed
- Models how guide perturbations propagate through the promoter phase space
- Uses Kuramoto-style phase coupling to predict transcriptional outcome

Key Insight:
A CRISPRa guide doesn't just "bind" - it perturbs the local transcriptional
phase landscape. Effective guides create coherent phase perturbations that
propagate constructively to the TSS (transcription start site).

Pipeline:
    Guide position → Map to promoter phase space (θ_g)
    → Model local perturbation ξ(θ_g)
    → Simulate phase evolution to TSS (Kuramoto dynamics)
    → Compute final phase alignment at TSS
    → Return transcription enhancement factor

Theory:
    dθ_i/dt = ω_i + Σ_j K_{ij} sin(θ_j - θ_i) + ξ_i(t)

    Where:
    - θ_i: Phase at position i along the promoter
    - ω_i: Natural frequency (determined by sequence)
    - K_{ij}: Coupling strength (decays with distance)
    - ξ_i: External perturbation (the guide!)

Author: PhaseLab
Date: December 2025
"""

import numpy as np
from dataclasses import dataclass
from typing import List, Optional, Tuple, Dict, Any
from scipy.integrate import odeint
import logging

logger = logging.getLogger(__name__)

# Phase coupling parameters
DEFAULT_COUPLING_STRENGTH = 0.5  # K: base coupling
DEFAULT_PERTURBATION_STRENGTH = 1.0  # Guide perturbation amplitude
COUPLING_DECAY_LENGTH = 50  # bp: coupling falls off exponentially


@dataclass
class PhaseAlignmentResult:
    """
    Result from transcriptional phase alignment analysis.

    Attributes:
        phase_coherence: R̄ at TSS after perturbation (higher = better)
        enhancement_factor: Predicted fold-change in transcription
        critical_window: Whether guide is in optimal TSS window (-400 to -50)
        phase_velocity: Rate of phase propagation to TSS
        evidence: Evidence level
        details: Additional calculation details
    """
    phase_coherence: float
    enhancement_factor: float
    critical_window: bool
    phase_velocity: float
    evidence: str
    details: Dict[str, Any]

    def is_go(self) -> bool:
        """Check if phase coherence passes GO threshold (e^-2 ≈ 0.135)."""
        return self.phase_coherence > 0.135

    def __repr__(self):
        status = "GO" if self.is_go() else "NO-GO"
        return (
            f"PhaseAlignmentResult(R̄={self.phase_coherence:.3f} [{status}], "
            f"enhancement={self.enhancement_factor:.2f}×, "
            f"window={'YES' if self.critical_window else 'NO'})"
        )


def _sequence_to_frequencies(sequence: str) -> np.ndarray:
    """
    Convert DNA sequence to natural frequencies.

    GC-rich regions have higher "stiffness" (higher frequency)
    AT-rich regions are more flexible (lower frequency)
    """
    # Window size for frequency calculation
    window = 5
    n = len(sequence)
    sequence = sequence.upper()

    # Base frequency contributions
    freq_map = {
        'G': 1.2, 'C': 1.2,  # GC = higher frequency (stiffer)
        'A': 0.8, 'T': 0.8,  # AT = lower frequency (flexible)
        'U': 0.8, 'N': 1.0,
    }

    frequencies = np.ones(n)
    for i in range(n):
        # Average over local window
        start = max(0, i - window // 2)
        end = min(n, i + window // 2 + 1)
        local_freq = np.mean([freq_map.get(sequence[j], 1.0) for j in range(start, end)])
        frequencies[i] = local_freq

    return frequencies


def _build_coupling_matrix(
    n_sites: int,
    decay_length: float = COUPLING_DECAY_LENGTH,
    base_coupling: float = DEFAULT_COUPLING_STRENGTH,
) -> np.ndarray:
    """
    Build position-dependent coupling matrix.

    Coupling decays exponentially with distance:
    K_{ij} = K_0 * exp(-|i-j| / λ)
    """
    K = np.zeros((n_sites, n_sites))
    for i in range(n_sites):
        for j in range(n_sites):
            if i != j:
                distance = abs(i - j)
                K[i, j] = base_coupling * np.exp(-distance / decay_length)
    return K


def _kuramoto_dynamics(
    theta: np.ndarray,
    t: float,
    omega: np.ndarray,
    K: np.ndarray,
    perturbation: np.ndarray,
) -> np.ndarray:
    """
    Kuramoto model dynamics with external perturbation (vectorized).

    dθ_i/dt = ω_i + Σ_j K_{ij} sin(θ_j - θ_i) + ξ_i
    """
    # Vectorized: compute all phase differences at once
    # theta_j - theta_i as a matrix (broadcasting)
    phase_diff = theta[np.newaxis, :] - theta[:, np.newaxis]  # [n, n]

    # Coupling term: K[i,j] * sin(theta_j - theta_i)
    coupling = np.sum(K * np.sin(phase_diff), axis=1)

    return omega + coupling + perturbation


def _compute_order_parameter(theta: np.ndarray) -> float:
    """
    Compute Kuramoto order parameter R.

    R = |<e^{iθ}>| measures phase synchronization.
    R = 1: fully synchronized
    R = 0: fully desynchronized
    """
    phasors = np.exp(1j * theta)
    return float(np.abs(np.mean(phasors)))


def simulate_phase_evolution(
    promoter_sequence: str,
    tss_position: int,
    guide_position: int,
    guide_length: int = 20,
    perturbation_strength: float = DEFAULT_PERTURBATION_STRENGTH,
    t_max: float = 5.0,
    n_steps: int = 50,
    max_sites: int = 200,  # Reduce computation by coarse-graining
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Simulate phase evolution after guide perturbation.

    For long sequences (>max_sites), uses coarse-grained representation.

    Args:
        promoter_sequence: Full promoter DNA sequence
        tss_position: Position of TSS in sequence
        guide_position: Start position of guide binding
        guide_length: Length of guide
        perturbation_strength: Amplitude of guide perturbation
        t_max: Simulation time
        n_steps: Number of time steps

    Returns:
        (theta_history, t_history, final_coherence_at_tss)
    """
    n_original = len(promoter_sequence)

    # Coarse-grain if sequence too long
    if n_original > max_sites:
        # Downsample by factor
        factor = n_original // max_sites
        n = max_sites

        # Downsample sequence
        downsampled_seq = promoter_sequence[::factor][:n]

        # Rescale positions
        tss_pos_scaled = tss_position // factor
        guide_pos_scaled = guide_position // factor
        guide_len_scaled = max(1, guide_length // factor)
    else:
        n = n_original
        downsampled_seq = promoter_sequence
        tss_pos_scaled = tss_position
        guide_pos_scaled = guide_position
        guide_len_scaled = guide_length

    # Natural frequencies from sequence
    omega = _sequence_to_frequencies(downsampled_seq)

    # Coupling matrix
    K = _build_coupling_matrix(n)

    # Initial phases (random, representing unperturbed state)
    theta0 = np.random.uniform(0, 2 * np.pi, n)

    # External perturbation (localized at guide binding site)
    perturbation = np.zeros(n)
    for i in range(guide_pos_scaled, min(guide_pos_scaled + guide_len_scaled, n)):
        # Gaussian-like perturbation centered at guide
        rel_pos = (i - guide_pos_scaled) / max(1, guide_len_scaled)
        perturbation[i] = perturbation_strength * np.exp(-2 * (rel_pos - 0.5) ** 2)

    # Time evolution
    t = np.linspace(0, t_max, n_steps)

    # Solve ODE
    theta_history = odeint(
        _kuramoto_dynamics,
        theta0,
        t,
        args=(omega, K, perturbation),
    )

    # Compute coherence at TSS region (±5 sites around scaled TSS)
    tss_region = slice(max(0, tss_pos_scaled - 5), min(n, tss_pos_scaled + 5))
    final_theta_tss = theta_history[-1, tss_region]
    final_coherence = _compute_order_parameter(final_theta_tss)

    return theta_history, t, final_coherence


def compute_phase_alignment(
    guide_sequence: str,
    promoter_sequence: str,
    tss_position: int,
    guide_position: int,
    perturbation_strength: float = DEFAULT_PERTURBATION_STRENGTH,
) -> PhaseAlignmentResult:
    """
    Compute transcriptional phase alignment score for a guide.

    This is the main entry point for Path B: Transcriptional Phase Alignment.

    Args:
        guide_sequence: Guide RNA sequence (for verification)
        promoter_sequence: Full promoter DNA sequence
        tss_position: Position of TSS in promoter_sequence
        guide_position: Start position of guide binding

    Returns:
        PhaseAlignmentResult with coherence and enhancement prediction

    Example:
        >>> result = compute_phase_alignment(
        ...     guide_sequence="ATCGATCGATCGATCGATCG",
        ...     promoter_sequence=promoter_seq,
        ...     tss_position=500,
        ...     guide_position=200,
        ... )
        >>> print(f"Phase coherence: {result.phase_coherence:.3f}")
        >>> print(f"Enhancement: {result.enhancement_factor:.2f}×")
    """
    guide_length = len(guide_sequence)

    # Check if guide is in optimal CRISPRa window
    relative_position = guide_position - tss_position
    critical_window = -400 <= relative_position <= -50

    # Run phase simulation
    theta_history, t_history, final_coherence = simulate_phase_evolution(
        promoter_sequence=promoter_sequence,
        tss_position=tss_position,
        guide_position=guide_position,
        guide_length=guide_length,
        perturbation_strength=perturbation_strength,
    )

    # Also compute initial (unperturbed) coherence for comparison
    initial_coherence = _compute_order_parameter(theta_history[0, :])

    # Compute phase velocity (how fast perturbation reaches TSS)
    mid_idx = len(t_history) // 2
    mid_coherence = _compute_order_parameter(theta_history[mid_idx, :])
    phase_velocity = (final_coherence - initial_coherence) / t_history[-1]

    # Enhancement factor model:
    # Based on phase coherence increase and window position
    base_enhancement = 1 + 5 * final_coherence  # Up to 6× for perfect coherence

    # Bonus for optimal window position
    if critical_window:
        window_bonus = 1.5  # 50% bonus in optimal window
    else:
        # Penalty outside window
        window_bonus = 0.5 + 0.5 * np.exp(-abs(relative_position + 200) / 100)

    enhancement_factor = base_enhancement * window_bonus

    return PhaseAlignmentResult(
        phase_coherence=final_coherence,
        enhancement_factor=enhancement_factor,
        critical_window=critical_window,
        phase_velocity=phase_velocity,
        evidence="IR_DYNAMICS",
        details={
            'guide_position': guide_position,
            'tss_position': tss_position,
            'relative_position': relative_position,
            'initial_coherence': initial_coherence,
            'guide_length': guide_length,
            'simulation_time': t_history[-1],
        }
    )


def compute_phase_landscape(
    guides: List[Dict[str, Any]],
    promoter_sequence: str,
    tss_position: int,
) -> List[Dict[str, Any]]:
    """
    Compute phase alignment for multiple guides.

    Args:
        guides: List of guide dicts with 'sequence' and 'position' keys
        promoter_sequence: Full promoter DNA sequence
        tss_position: Position of TSS

    Returns:
        Guides with added 'phase_alignment' field

    Example:
        >>> guides = [{'sequence': 'ATG...', 'position': 200}, ...]
        >>> results = compute_phase_landscape(guides, promoter, tss=500)
    """
    for guide in guides:
        seq = guide.get('sequence', '')
        pos = guide.get('position', guide.get('tss_relative_position', 0))

        # If position is relative to TSS, convert to absolute
        if 'tss_relative_position' in guide:
            pos = tss_position + guide['tss_relative_position']

        try:
            result = compute_phase_alignment(
                guide_sequence=seq,
                promoter_sequence=promoter_sequence,
                tss_position=tss_position,
                guide_position=pos,
            )

            guide['phase_alignment'] = {
                'coherence': result.phase_coherence,
                'enhancement_factor': result.enhancement_factor,
                'critical_window': result.critical_window,
                'phase_velocity': result.phase_velocity,
                'is_go': result.is_go(),
            }
        except Exception as e:
            logger.warning(f"Phase alignment failed for {seq[:10]}...: {e}")
            guide['phase_alignment'] = {
                'coherence': 0.0,
                'enhancement_factor': 1.0,
                'critical_window': False,
                'phase_velocity': 0.0,
                'is_go': False,
                'error': str(e),
            }

    return guides


def rank_guides_by_phase(
    guides: List[Dict[str, Any]],
    promoter_sequence: str,
    tss_position: int,
) -> List[Dict[str, Any]]:
    """
    Rank guides by transcriptional phase alignment.

    This can be used as a post-processing step after standard guide design.

    Args:
        guides: List of guide dicts
        promoter_sequence: Full promoter DNA
        tss_position: TSS position

    Returns:
        Guides sorted by phase coherence (highest first)
    """
    guides = compute_phase_landscape(guides, promoter_sequence, tss_position)

    # Sort by phase coherence (descending)
    guides.sort(
        key=lambda g: g.get('phase_alignment', {}).get('coherence', 0),
        reverse=True
    )

    return guides


def optimal_guide_position(
    promoter_sequence: str,
    tss_position: int,
    search_window: Tuple[int, int] = (-400, -50),
    step_size: int = 20,
) -> Dict[str, Any]:
    """
    Find optimal guide position by scanning phase landscape.

    Args:
        promoter_sequence: Full promoter sequence
        tss_position: TSS position
        search_window: Window relative to TSS to search
        step_size: Position step size (bp)

    Returns:
        Dict with optimal position and predicted enhancement

    Example:
        >>> result = optimal_guide_position(promoter, tss=500)
        >>> print(f"Best position: {result['position']} ({result['enhancement']:.2f}×)")
    """
    positions = list(range(
        tss_position + search_window[0],
        tss_position + search_window[1],
        step_size
    ))

    best_result = None
    best_coherence = 0.0

    for pos in positions:
        # Use dummy 20bp guide for position scanning
        result = compute_phase_alignment(
            guide_sequence="N" * 20,  # Placeholder
            promoter_sequence=promoter_sequence,
            tss_position=tss_position,
            guide_position=pos,
        )

        if result.phase_coherence > best_coherence:
            best_coherence = result.phase_coherence
            best_result = {
                'position': pos,
                'relative_to_tss': pos - tss_position,
                'coherence': result.phase_coherence,
                'enhancement': result.enhancement_factor,
            }

    return best_result


__all__ = [
    'PhaseAlignmentResult',
    'compute_phase_alignment',
    'compute_phase_landscape',
    'rank_guides_by_phase',
    'optimal_guide_position',
]
