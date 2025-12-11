"""
PhaseLab SMS Model: Smith-Magenis Syndrome circadian clock simulation.

Extended Kuramoto model with:
- RAI1 dosage-dependent coupling
- PER gene delayed negative feedback
- REV-ERBα / RORα modulation of BMAL1
- Therapeutic window analysis

The model predicts how RAI1 haploinsufficiency affects circadian synchronization
and what boost level is needed to restore normal rhythms.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Any, List, Tuple, Optional
from scipy.integrate import solve_ivp

from ..core.coherence import coherence_score, go_no_go
from ..core.constants import E_MINUS_2, OMEGA_CIRCADIAN


@dataclass
class SMSClockParams:
    """
    Parameters for the Smith-Magenis Syndrome circadian model.

    All frequencies are in rad/hour, time is in hours.
    """
    # Natural frequencies (~2π/24 with small heterogeneity)
    omega_C: float = OMEGA_CIRCADIAN  # CLOCK
    omega_B: float = OMEGA_CIRCADIAN * 1.02  # BMAL1 (slight detuning)

    # Baseline coupling between CLOCK and BMAL1
    K_base: float = 0.6

    # === PER feedback parameters ===
    tau_P: float = 4.0      # hours, effective delay time constant
    alpha_P: float = 2.0    # how strongly PER suppresses coupling
    P_init: float = 0.0     # initial PER level

    # === RORα / REV-ERBα parameters ===
    tau_R: float = 12.0     # slow adjustment time for RORα
    tau_V: float = 12.0     # slow adjustment time for REV-ERBα
    R_target: float = 0.5   # steady-state RORα level
    V_target: float = 0.5   # steady-state REV-ERBα level
    R_init: float = 0.5     # initial RORα
    V_init: float = 0.5     # initial REV-ERBα

    # Sensitivity curve parameters (sigmoids)
    R_mid: float = 0.5      # RORα sigmoid midpoint
    R_k: float = 0.15       # RORα sigmoid steepness
    V_mid: float = 0.5      # REV-ERBα sigmoid midpoint
    V_k: float = 0.15       # REV-ERBα sigmoid steepness

    # Modulation strengths
    beta_R: float = 0.5     # RORα effect on BMAL1 coupling
    beta_V: float = 0.5     # REV-ERBα effect on BMAL1 coupling
    K0: float = 0.6         # Base BMAL1 coupling scale

    # RAI1 modulation
    rai1_scale_factor: float = 1.0  # multiplied by rai1_level at runtime

    # Additional oscillators (PER1/2, CRY1/2, REV-ERBα, RORα)
    include_per_cry: bool = False  # Extended model with PER/CRY oscillators

    # Noise
    noise_strength: float = 0.0  # Phase diffusion noise


def _sigmoid(x: float, x0: float, k: float) -> float:
    """Sigmoid function for sensitivity curves."""
    return 1.0 / (1.0 + np.exp(-(x - x0) / k))


def _sms_clock_ode(
    t: float,
    y: np.ndarray,
    params: SMSClockParams,
    rai1_level: float,
) -> np.ndarray:
    """
    SMS circadian clock ODE with PER delay and REV-ERBα/RORα modulation.

    State vector y = [theta_C, theta_B, P, R, V]
    - theta_C: CLOCK phase
    - theta_B: BMAL1 phase
    - P: PER inhibition level (delayed feedback)
    - R: RORα level (activates BMAL1)
    - V: REV-ERBα level (represses BMAL1)
    """
    theta_C, theta_B, P, R, V = y

    # Unpack parameters
    omega_C = params.omega_C
    omega_B0 = params.omega_B

    # RAI1 modulates base coupling
    K_base = params.K_base * params.rai1_scale_factor * rai1_level

    # === PER feedback ===
    # PER builds up when CLOCK/BMAL1 are in-phase (active)
    f_clock = 0.5 * (1.0 + np.cos(theta_C - theta_B))  # high when synchronized
    dP = (f_clock - P) / params.tau_P

    # === RORα / REV-ERBα slow dynamics ===
    # These regulate BMAL1 transcription on slower timescale
    # RORα target could depend on clock phase (optional enhancement)
    dR = (params.R_target - R) / params.tau_R
    dV = (params.V_target - V) / params.tau_V

    # === Sensitivity curves ===
    S_R = _sigmoid(R, params.R_mid, params.R_k)  # RORα activation
    S_V = _sigmoid(V, params.V_mid, params.V_k)  # REV-ERBα repression

    # === BMAL1 coupling modulation ===
    # RORα activates, REV-ERBα represses
    K_B = params.K0 * (1.0 + params.beta_R * S_R - params.beta_V * S_V)

    # === PER reduces effective coupling (delayed negative feedback) ===
    K_eff = (K_base + K_B) / (1.0 + params.alpha_P * P)

    # === Phase dynamics (Kuramoto-like) ===
    dtheta_C = omega_C + K_eff * np.sin(theta_B - theta_C)
    dtheta_B = omega_B0 + K_eff * np.sin(theta_C - theta_B)

    # Add noise if specified
    if params.noise_strength > 0:
        dtheta_C += params.noise_strength * np.random.randn()
        dtheta_B += params.noise_strength * np.random.randn()

    return np.array([dtheta_C, dtheta_B, dP, dR, dV])


def simulate_sms_clock(
    rai1_level: float = 0.5,
    t_end: float = 240.0,  # 10 days
    dt: float = 0.1,
    params: Optional[SMSClockParams] = None,
    y0: Optional[np.ndarray] = None,
    random_seed: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Simulate the SMS circadian clock under a given RAI1 level.

    Args:
        rai1_level: RAI1 expression level (0.0-1.5+).
            - 0.5 = typical SMS (haploinsufficient)
            - 1.0 = normal
            - 0.6-0.8 = therapeutic target
        t_end: Simulation duration in hours (default 240 = 10 days).
        dt: Output time step.
        params: SMSClockParams (uses defaults if None).
        y0: Initial state [theta_C, theta_B, P, R, V].
        random_seed: For reproducible initial conditions.

    Returns:
        Dictionary with:
            - t: time array
            - theta_C: CLOCK phase
            - theta_B: BMAL1 phase
            - P: PER feedback level
            - R: RORα level
            - V: REV-ERBα level
            - order_param: R̄(t) synchronization
            - final_R_bar: final synchronization score
            - classification: synchronization quality
            - go_no_go: GO/NO-GO status
    """
    if params is None:
        params = SMSClockParams()

    if random_seed is not None:
        np.random.seed(random_seed)

    t_eval = np.arange(0.0, t_end + dt, dt)

    if y0 is None:
        y0 = np.array([
            np.random.uniform(0, 2 * np.pi),  # theta_C
            np.random.uniform(0, 2 * np.pi),  # theta_B
            params.P_init,                     # P
            params.R_init,                     # R
            params.V_init,                     # V
        ])

    sol = solve_ivp(
        fun=lambda t, y: _sms_clock_ode(t, y, params, rai1_level),
        t_span=(0.0, t_end),
        y0=y0,
        t_eval=t_eval,
        method='RK45',
        rtol=1e-6,
        atol=1e-9,
    )

    theta_C = np.mod(sol.y[0], 2 * np.pi)
    theta_B = np.mod(sol.y[1], 2 * np.pi)
    P = sol.y[2]
    R = sol.y[3]
    V = sol.y[4]

    # Compute order parameter (synchronization) over time
    # Using just CLOCK and BMAL1 phases
    z = np.exp(1j * theta_C) + np.exp(1j * theta_B)
    R_bar = np.abs(z) / 2.0

    # Use last 25% of simulation for steady-state metrics
    steady_idx = int(0.75 * len(R_bar))
    final_R_bar = float(np.mean(R_bar[steady_idx:]))

    # Phase variance
    phase_diff = theta_C - theta_B
    V_phi = float(np.var(np.mod(phase_diff[steady_idx:], 2 * np.pi)))

    return {
        't': sol.t,
        'theta_C': theta_C,
        'theta_B': theta_B,
        'P': P,
        'R': R,
        'V': V,
        'order_param': R_bar,
        'final_R_bar': final_R_bar,
        'phase_variance': V_phi,
        'classification': classify_synchronization(final_R_bar),
        'go_no_go': go_no_go(final_R_bar),
        'rai1_level': rai1_level,
    }


def classify_synchronization(R_bar: float) -> str:
    """
    Classify circadian synchronization quality.

    Args:
        R_bar: Order parameter (0-1).

    Returns:
        Classification string.
    """
    if R_bar >= 0.9:
        return "SYNCHRONIZED"
    elif R_bar >= 0.7:
        return "PARTIALLY_SYNCHRONIZED"
    elif R_bar >= E_MINUS_2:
        return "WEAKLY_SYNCHRONIZED"
    else:
        return "DESYNCHRONIZED"


def therapeutic_scan(
    rai1_levels: Optional[List[float]] = None,
    params: Optional[SMSClockParams] = None,
    t_end: float = 240.0,
    n_trials: int = 3,
) -> Dict[str, Any]:
    """
    Scan RAI1 levels to find therapeutic window.

    Simulates clock at multiple RAI1 levels to identify
    the boost needed to restore synchronization.

    Args:
        rai1_levels: List of RAI1 levels to test.
        params: SMSClockParams.
        t_end: Simulation duration.
        n_trials: Number of trials per level (for averaging).

    Returns:
        Dictionary with:
            - levels: RAI1 levels tested
            - R_bars: mean synchronization at each level
            - R_bar_std: standard deviation
            - classifications: synchronization class
            - therapeutic_window: (min, max) for SYNCHRONIZED
            - optimal_level: RAI1 level with best sync
    """
    if rai1_levels is None:
        rai1_levels = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]

    results = {
        'levels': rai1_levels,
        'R_bars': [],
        'R_bar_std': [],
        'classifications': [],
    }

    for level in rai1_levels:
        trial_R_bars = []
        for trial in range(n_trials):
            res = simulate_sms_clock(
                rai1_level=level,
                t_end=t_end,
                params=params,
                random_seed=42 + trial,
            )
            trial_R_bars.append(res['final_R_bar'])

        mean_R = np.mean(trial_R_bars)
        std_R = np.std(trial_R_bars)

        results['R_bars'].append(mean_R)
        results['R_bar_std'].append(std_R)
        results['classifications'].append(classify_synchronization(mean_R))

    # Find therapeutic window
    R_bars = np.array(results['R_bars'])
    levels = np.array(rai1_levels)

    synchronized_mask = R_bars >= 0.9
    if np.any(synchronized_mask):
        sync_levels = levels[synchronized_mask]
        results['therapeutic_window'] = (float(sync_levels.min()), float(sync_levels.max()))
        results['optimal_level'] = float(levels[np.argmax(R_bars)])
    else:
        # Find best available
        results['therapeutic_window'] = None
        results['optimal_level'] = float(levels[np.argmax(R_bars)])

    # SMS baseline (50%) comparison
    sms_idx = np.argmin(np.abs(levels - 0.5))
    results['sms_baseline_R'] = float(R_bars[sms_idx])

    # Required boost
    if results['therapeutic_window']:
        min_therapeutic = results['therapeutic_window'][0]
        results['required_boost'] = min_therapeutic - 0.5  # from SMS baseline
    else:
        results['required_boost'] = None

    return results


def predict_sleep_quality(R_bar: float) -> str:
    """
    Predict sleep quality based on circadian synchronization.

    Based on SMS clinical observations and model predictions.

    Args:
        R_bar: Circadian order parameter.

    Returns:
        Sleep quality prediction string.
    """
    if R_bar >= 0.9:
        return "NORMAL - Regular sleep-wake cycles expected"
    elif R_bar >= 0.7:
        return "MILD_DISRUPTION - Some sleep fragmentation likely"
    elif R_bar >= 0.5:
        return "MODERATE_DISRUPTION - Significant sleep disturbances"
    elif R_bar >= E_MINUS_2:
        return "SEVERE_DISRUPTION - Inverted melatonin rhythm possible"
    else:
        return "CRITICAL - Circadian rhythm severely compromised"
