"""
PhaseLab Kuramoto: Base Kuramoto oscillator model.

The Kuramoto model describes coupled oscillators:
    dθ_i/dt = ω_i + Σ_j K_ij sin(θ_j - θ_i)

The order parameter R̄ measures synchronization:
    R̄ = |<e^{iθ}>| ∈ [0, 1]
    R̄ = 1: Perfect synchronization
    R̄ = 0: Complete desynchronization
"""

import numpy as np
from typing import Tuple, Optional, Callable
from scipy.integrate import solve_ivp


def kuramoto_order_parameter(phases: np.ndarray) -> Tuple[float, float]:
    """
    Compute Kuramoto order parameter from phase array.

    Args:
        phases: Array of oscillator phases (radians).

    Returns:
        (R, psi) where R is magnitude and psi is mean phase.
    """
    z = np.mean(np.exp(1j * phases))
    R = np.abs(z)
    psi = np.angle(z)
    return float(R), float(psi)


def kuramoto_ode(
    t: float,
    phases: np.ndarray,
    omegas: np.ndarray,
    K: np.ndarray,
) -> np.ndarray:
    """
    Kuramoto ODE right-hand side.

    dθ_i/dt = ω_i + Σ_j K_ij sin(θ_j - θ_i)

    Args:
        t: Time (not used, for ODE solver interface).
        phases: Current phases of all oscillators.
        omegas: Natural frequencies of oscillators.
        K: Coupling matrix (N x N).

    Returns:
        Phase derivatives dθ/dt.
    """
    N = len(phases)
    dphases = np.copy(omegas)

    for i in range(N):
        for j in range(N):
            if i != j:
                dphases[i] += K[i, j] * np.sin(phases[j] - phases[i])

    return dphases


def simulate_kuramoto(
    omegas: np.ndarray,
    K: np.ndarray,
    t_end: float = 240.0,
    dt: float = 0.1,
    phases_init: Optional[np.ndarray] = None,
) -> dict:
    """
    Simulate a Kuramoto oscillator network.

    Args:
        omegas: Natural frequencies (rad/time unit).
        K: Coupling matrix (N x N).
        t_end: Simulation duration.
        dt: Time step for output.
        phases_init: Initial phases (random if None).

    Returns:
        Dictionary with 't', 'phases', 'order_param' arrays.
    """
    N = len(omegas)
    t_eval = np.arange(0, t_end + dt, dt)

    if phases_init is None:
        phases_init = np.random.uniform(0, 2 * np.pi, N)

    sol = solve_ivp(
        fun=lambda t, y: kuramoto_ode(t, y, omegas, K),
        t_span=(0, t_end),
        y0=phases_init,
        t_eval=t_eval,
        method='RK45',
    )

    # Compute order parameter over time
    R_t = np.zeros(len(sol.t))
    for i, t in enumerate(sol.t):
        R_t[i], _ = kuramoto_order_parameter(sol.y[:, i])

    return {
        't': sol.t,
        'phases': np.mod(sol.y, 2 * np.pi),
        'order_param': R_t,
    }
