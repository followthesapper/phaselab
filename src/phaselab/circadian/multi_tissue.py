"""
PhaseLab Multi-Tissue Circadian: Inter-tissue clock coupling model.

Models the circadian clock network across multiple tissues, including:
- SCN (master pacemaker)
- Liver (metabolic clock)
- Muscle (peripheral clock)
- Heart (cardiovascular rhythms)
- Brain regions (cognitive performance)

The model captures how peripheral clocks synchronize to the SCN and
how this coupling is affected by disease or therapeutic interventions.

Version: 0.3.0
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any
from scipy.integrate import solve_ivp

from ..core.coherence import coherence_score, go_no_go
from ..core.constants import E_MINUS_2, OMEGA_CIRCADIAN


# Tissue-specific parameters
TISSUE_PARAMS = {
    "SCN": {
        "omega": OMEGA_CIRCADIAN,  # Master clock, defines rhythm
        "amplitude": 1.0,          # Strongest oscillator
        "coupling_to_scn": 0.0,    # Self-coupling N/A
        "intrinsic_period": 24.0,  # hours
    },
    "liver": {
        "omega": OMEGA_CIRCADIAN * 0.98,  # Slight period difference
        "amplitude": 0.8,
        "coupling_to_scn": 0.4,
        "intrinsic_period": 24.5,
    },
    "muscle": {
        "omega": OMEGA_CIRCADIAN * 1.01,
        "amplitude": 0.7,
        "coupling_to_scn": 0.3,
        "intrinsic_period": 23.8,
    },
    "heart": {
        "omega": OMEGA_CIRCADIAN * 0.99,
        "amplitude": 0.75,
        "coupling_to_scn": 0.35,
        "intrinsic_period": 24.2,
    },
    "brain_cortex": {
        "omega": OMEGA_CIRCADIAN * 1.0,
        "amplitude": 0.85,
        "coupling_to_scn": 0.45,
        "intrinsic_period": 24.0,
    },
    "adipose": {
        "omega": OMEGA_CIRCADIAN * 0.97,
        "amplitude": 0.6,
        "coupling_to_scn": 0.25,
        "intrinsic_period": 24.7,
    },
}


@dataclass
class MultiTissueParams:
    """Parameters for multi-tissue circadian model."""

    # Tissues to include
    tissues: List[str] = field(default_factory=lambda: ["SCN", "liver", "muscle"])

    # Global coupling strength (scaled by tissue-specific values)
    K_global: float = 0.5

    # Inter-peripheral coupling (tissues can influence each other)
    K_peripheral: float = 0.1

    # Light input to SCN
    light_amplitude: float = 0.3
    light_phase: float = 0.0  # Phase of light onset (0 = midnight)

    # Feeding input to peripheral clocks
    feeding_amplitude: float = 0.2
    feeding_phase: float = np.pi / 2  # Phase of feeding (π/2 = 6am)

    # Disease modulation
    disease_tissue: Optional[str] = None
    disease_severity: float = 0.0  # 0 = healthy, 1 = severe

    # Noise
    noise_strength: float = 0.01


@dataclass
class MultiTissueResult:
    """Results from multi-tissue simulation."""

    t: np.ndarray
    phases: Dict[str, np.ndarray]  # tissue -> phase trajectory
    amplitudes: Dict[str, np.ndarray]

    # Coherence metrics
    global_R_bar: float
    tissue_R_bars: Dict[str, float]
    scn_peripheral_coherence: float

    # Phase relationships
    phase_delays: Dict[str, float]  # Delay relative to SCN

    # Classification
    go_no_go: str
    classification: str

    def to_dict(self) -> Dict[str, Any]:
        return {
            "global_R_bar": self.global_R_bar,
            "tissue_R_bars": self.tissue_R_bars,
            "scn_peripheral_coherence": self.scn_peripheral_coherence,
            "phase_delays": self.phase_delays,
            "go_no_go": self.go_no_go,
            "classification": self.classification,
        }


def _multi_tissue_ode(
    t: float,
    y: np.ndarray,
    params: MultiTissueParams,
    tissue_list: List[str],
) -> np.ndarray:
    """
    Multi-tissue circadian ODE.

    State: [phase_tissue1, phase_tissue2, ..., amp_tissue1, amp_tissue2, ...]
    """
    n = len(tissue_list)
    phases = y[:n]
    amps = y[n:]

    dydt = np.zeros(2 * n)

    # Get SCN index if present
    scn_idx = tissue_list.index("SCN") if "SCN" in tissue_list else None
    scn_phase = phases[scn_idx] if scn_idx is not None else 0.0

    for i, tissue in enumerate(tissue_list):
        tp = TISSUE_PARAMS[tissue]
        omega = tp["omega"]
        coupling_to_scn = tp["coupling_to_scn"]

        # Base frequency
        dphase = omega

        # Coupling to SCN (for peripheral clocks)
        if tissue != "SCN" and scn_idx is not None:
            K_scn = params.K_global * coupling_to_scn
            dphase += K_scn * np.sin(scn_phase - phases[i])

        # Inter-peripheral coupling
        if tissue != "SCN":
            for j, other_tissue in enumerate(tissue_list):
                if i != j and other_tissue != "SCN":
                    dphase += params.K_peripheral * np.sin(phases[j] - phases[i])

        # Light input to SCN
        if tissue == "SCN":
            light_signal = params.light_amplitude * np.cos(
                OMEGA_CIRCADIAN * t - params.light_phase
            )
            dphase += light_signal

        # Feeding input to peripheral clocks
        if tissue in ["liver", "muscle", "adipose"]:
            feeding_signal = params.feeding_amplitude * np.cos(
                OMEGA_CIRCADIAN * t - params.feeding_phase
            )
            dphase += feeding_signal

        # Disease modulation
        if params.disease_tissue == tissue:
            # Reduce coupling strength
            dphase *= (1.0 - params.disease_severity * 0.5)

        # Add noise
        if params.noise_strength > 0:
            dphase += params.noise_strength * np.random.randn()

        dydt[i] = dphase

        # Amplitude dynamics (slow relaxation to intrinsic amplitude)
        target_amp = tp["amplitude"]
        if params.disease_tissue == tissue:
            target_amp *= (1.0 - params.disease_severity * 0.3)
        dydt[n + i] = 0.1 * (target_amp - amps[i])

    return dydt


def simulate_multi_tissue(
    params: Optional[MultiTissueParams] = None,
    t_end: float = 240.0,
    dt: float = 0.5,
    random_seed: Optional[int] = None,
) -> MultiTissueResult:
    """
    Simulate multi-tissue circadian clock.

    Args:
        params: MultiTissueParams (uses defaults if None)
        t_end: Simulation duration in hours
        dt: Output time step
        random_seed: For reproducible results

    Returns:
        MultiTissueResult with trajectories and coherence metrics
    """
    if params is None:
        params = MultiTissueParams()

    if random_seed is not None:
        np.random.seed(random_seed)

    tissues = params.tissues
    n = len(tissues)

    # Initial conditions
    y0 = np.zeros(2 * n)
    for i, tissue in enumerate(tissues):
        y0[i] = np.random.uniform(0, 2 * np.pi)  # Random initial phase
        y0[n + i] = TISSUE_PARAMS[tissue]["amplitude"]

    t_eval = np.arange(0, t_end + dt, dt)

    sol = solve_ivp(
        fun=lambda t, y: _multi_tissue_ode(t, y, params, tissues),
        t_span=(0, t_end),
        y0=y0,
        t_eval=t_eval,
        method='RK45',
    )

    # Extract results
    phases = {}
    amplitudes = {}
    for i, tissue in enumerate(tissues):
        phases[tissue] = np.mod(sol.y[i], 2 * np.pi)
        amplitudes[tissue] = sol.y[n + i]

    # Compute coherence metrics
    # Use last 25% for steady-state
    steady_idx = int(0.75 * len(sol.t))

    # Global coherence (all tissues)
    all_phases = np.array([phases[t][steady_idx:] for t in tissues])
    z_global = np.mean(np.exp(1j * all_phases), axis=0)
    global_R_bar = float(np.mean(np.abs(z_global)))

    # Per-tissue coherence (stability)
    tissue_R_bars = {}
    for tissue in tissues:
        ph = phases[tissue][steady_idx:]
        # Phase stability as coherence
        z = np.exp(1j * ph)
        tissue_R_bars[tissue] = float(np.abs(np.mean(z)))

    # SCN-peripheral coherence
    if "SCN" in tissues:
        peripheral_tissues = [t for t in tissues if t != "SCN"]
        if peripheral_tissues:
            scn_ph = phases["SCN"][steady_idx:]
            periph_phases = [phases[t][steady_idx:] for t in peripheral_tissues]
            # Mean phase relationship
            phase_diffs = [np.mean(np.exp(1j * (scn_ph - pp))) for pp in periph_phases]
            scn_peripheral_coherence = float(np.abs(np.mean(phase_diffs)))
        else:
            scn_peripheral_coherence = 1.0
    else:
        scn_peripheral_coherence = 0.0

    # Phase delays relative to SCN
    phase_delays = {}
    if "SCN" in tissues:
        scn_mean_phase = np.mean(phases["SCN"][steady_idx:])
        for tissue in tissues:
            if tissue != "SCN":
                tissue_mean_phase = np.mean(phases[tissue][steady_idx:])
                delay = tissue_mean_phase - scn_mean_phase
                delay_hours = (delay / OMEGA_CIRCADIAN) % 24
                phase_delays[tissue] = float(delay_hours)

    # Classification
    if global_R_bar >= 0.9:
        classification = "SYNCHRONIZED"
    elif global_R_bar >= 0.7:
        classification = "PARTIALLY_SYNCHRONIZED"
    elif global_R_bar >= E_MINUS_2:
        classification = "WEAKLY_SYNCHRONIZED"
    else:
        classification = "DESYNCHRONIZED"

    return MultiTissueResult(
        t=sol.t,
        phases=phases,
        amplitudes=amplitudes,
        global_R_bar=global_R_bar,
        tissue_R_bars=tissue_R_bars,
        scn_peripheral_coherence=scn_peripheral_coherence,
        phase_delays=phase_delays,
        go_no_go=go_no_go(global_R_bar),
        classification=classification,
    )


def jet_lag_simulation(
    time_shift: float = 8.0,  # hours
    direction: str = "east",  # "east" or "west"
    tissues: Optional[List[str]] = None,
    t_end: float = 336.0,  # 14 days
) -> Dict[str, Any]:
    """
    Simulate jet lag recovery across tissues.

    Args:
        time_shift: Number of hours shifted
        direction: "east" (phase advance) or "west" (phase delay)
        tissues: Tissues to simulate
        t_end: Simulation duration

    Returns:
        Dictionary with recovery trajectories and times
    """
    if tissues is None:
        tissues = ["SCN", "liver", "muscle", "heart"]

    # Initial simulation (pre-travel, synchronized)
    params_pre = MultiTissueParams(tissues=tissues, K_global=0.5)
    result_pre = simulate_multi_tissue(params_pre, t_end=72, random_seed=42)

    # Get final phases as new initial conditions with shift
    shift_radians = time_shift * OMEGA_CIRCADIAN
    if direction == "east":
        shift_radians = -shift_radians  # Phase advance

    # Shift peripheral clocks (SCN adapts faster)
    n = len(tissues)
    y0_shifted = np.zeros(2 * n)
    scn_idx = tissues.index("SCN") if "SCN" in tissues else None

    for i, tissue in enumerate(tissues):
        final_phase = result_pre.phases[tissue][-1]
        if tissue == "SCN":
            # SCN adapts faster
            y0_shifted[i] = final_phase + shift_radians * 0.3
        else:
            # Peripheral clocks shift fully
            y0_shifted[i] = final_phase + shift_radians
        y0_shifted[n + i] = result_pre.amplitudes[tissue][-1]

    # Post-travel simulation
    t_eval = np.arange(0, t_end, 0.5)

    params_post = MultiTissueParams(tissues=tissues, K_global=0.5)

    sol = solve_ivp(
        fun=lambda t, y: _multi_tissue_ode(t, y, params_post, tissues),
        t_span=(0, t_end),
        y0=y0_shifted,
        t_eval=t_eval,
        method='RK45',
    )

    # Track phase alignment over time
    phases = {}
    for i, tissue in enumerate(tissues):
        phases[tissue] = np.mod(sol.y[i], 2 * np.pi)

    # Calculate recovery time for each tissue
    recovery_times = {}
    if scn_idx is not None:
        scn_ph = phases["SCN"]
        for tissue in tissues:
            if tissue != "SCN":
                tissue_ph = phases[tissue]
                phase_diff = np.abs(np.mod(tissue_ph - scn_ph + np.pi, 2*np.pi) - np.pi)

                # Find when phase difference stabilizes < 0.5 rad (~2h)
                aligned = phase_diff < 0.5
                for j in range(len(aligned) - 24):  # Need 24 consecutive hours
                    if np.all(aligned[j:j+48]):  # 24 hours at dt=0.5
                        recovery_times[tissue] = float(sol.t[j])
                        break
                else:
                    recovery_times[tissue] = float(t_end)  # Didn't recover

    # Calculate global coherence over time
    coherence_trajectory = []
    for j in range(len(sol.t)):
        ph_j = [phases[t][j] for t in tissues]
        z = np.mean(np.exp(1j * np.array(ph_j)))
        coherence_trajectory.append(float(np.abs(z)))

    return {
        "t": sol.t,
        "phases": phases,
        "coherence_trajectory": np.array(coherence_trajectory),
        "recovery_times": recovery_times,
        "time_shift": time_shift,
        "direction": direction,
        "mean_recovery_time": np.mean(list(recovery_times.values())) if recovery_times else None,
    }


def shift_work_simulation(
    shift_schedule: str = "rotating",  # "rotating", "night", "day"
    shift_duration_days: int = 7,
    tissues: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Simulate circadian disruption from shift work.

    Args:
        shift_schedule: Type of shift work pattern
        shift_duration_days: Days per shift rotation
        tissues: Tissues to simulate

    Returns:
        Dictionary with chronic disruption metrics
    """
    if tissues is None:
        tissues = ["SCN", "liver", "muscle"]

    # Define light/feeding schedules for different shifts
    schedules = {
        "day": {
            "light_phase": 0.0,  # Normal
            "feeding_phase": np.pi / 2,
        },
        "night": {
            "light_phase": np.pi,  # Inverted
            "feeding_phase": 3 * np.pi / 2,
        },
        "rotating": None,  # Will alternate
    }

    results = []
    total_days = shift_duration_days * 4  # 4 rotations

    if shift_schedule == "rotating":
        # Alternate between day and night
        for rotation in range(4):
            is_night = rotation % 2 == 1
            schedule = schedules["night" if is_night else "day"]

            params = MultiTissueParams(
                tissues=tissues,
                light_phase=schedule["light_phase"],
                feeding_phase=schedule["feeding_phase"],
            )

            result = simulate_multi_tissue(
                params, t_end=shift_duration_days * 24, random_seed=42 + rotation
            )
            results.append({
                "rotation": rotation,
                "shift_type": "night" if is_night else "day",
                "global_R_bar": result.global_R_bar,
                "classification": result.classification,
            })
    else:
        schedule = schedules[shift_schedule]
        params = MultiTissueParams(
            tissues=tissues,
            light_phase=schedule["light_phase"],
            feeding_phase=schedule["feeding_phase"],
        )
        result = simulate_multi_tissue(params, t_end=total_days * 24)
        results.append({
            "shift_type": shift_schedule,
            "global_R_bar": result.global_R_bar,
            "classification": result.classification,
        })

    # Summary metrics
    R_bars = [r["global_R_bar"] for r in results]
    mean_R_bar = float(np.mean(R_bars))
    min_R_bar = float(np.min(R_bars))

    return {
        "schedule": shift_schedule,
        "rotations": results,
        "mean_R_bar": mean_R_bar,
        "min_R_bar": min_R_bar,
        "chronic_disruption": mean_R_bar < 0.7,
        "go_no_go": go_no_go(mean_R_bar),
    }
