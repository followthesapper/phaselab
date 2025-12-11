"""
PhaseLab Pharmacokinetics: Drug concentration modeling with circadian modulation.

Implements pharmacokinetic models that account for:
- Circadian variation in absorption, distribution, metabolism, elimination
- Time-of-day dependent drug efficacy
- Drug-clock interactions
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Any
from scipy.integrate import solve_ivp

from ..core.constants import OMEGA_CIRCADIAN


@dataclass
class PKParams:
    """
    Pharmacokinetic parameters with optional circadian modulation.

    Standard one-compartment model with first-order absorption and elimination.
    Circadian modulation can affect each parameter.
    """
    # Base parameters
    ka: float = 1.0          # Absorption rate constant (1/h)
    ke: float = 0.1          # Elimination rate constant (1/h)
    V: float = 70.0          # Volume of distribution (L)
    F: float = 1.0           # Bioavailability (0-1)

    # Circadian modulation amplitudes (fraction of base value)
    ka_circadian_amp: float = 0.0      # ka varies ± this fraction
    ke_circadian_amp: float = 0.0      # ke varies ± this fraction
    ka_circadian_phase: float = 0.0    # Phase of ka maximum
    ke_circadian_phase: float = np.pi  # Phase of ke maximum

    # Derived parameters
    @property
    def half_life(self) -> float:
        """Elimination half-life in hours."""
        return np.log(2) / self.ke

    @property
    def tmax_approx(self) -> float:
        """Approximate time to peak concentration."""
        return np.log(self.ka / self.ke) / (self.ka - self.ke)


@dataclass
class PKModel:
    """Pharmacokinetic model with circadian effects."""

    params: PKParams
    drug_name: str = "unnamed"

    def get_ka(self, t: float, circadian_phase: float = 0.0) -> float:
        """Get absorption rate at time t with circadian modulation."""
        if self.params.ka_circadian_amp > 0:
            modulation = 1.0 + self.params.ka_circadian_amp * np.cos(
                OMEGA_CIRCADIAN * t + circadian_phase - self.params.ka_circadian_phase
            )
            return self.params.ka * modulation
        return self.params.ka

    def get_ke(self, t: float, circadian_phase: float = 0.0) -> float:
        """Get elimination rate at time t with circadian modulation."""
        if self.params.ke_circadian_amp > 0:
            modulation = 1.0 + self.params.ke_circadian_amp * np.cos(
                OMEGA_CIRCADIAN * t + circadian_phase - self.params.ke_circadian_phase
            )
            return self.params.ke * modulation
        return self.params.ke


def _pk_ode(
    t: float,
    y: np.ndarray,
    model: PKModel,
    dose_times: List[float],
    doses: List[float],
    circadian_phase: float,
) -> np.ndarray:
    """
    PK ODE for one-compartment model.

    State: [A_gut, A_central]
    A_gut: Amount in absorption compartment
    A_central: Amount in central compartment
    """
    A_gut, A_central = y

    ka = model.get_ka(t, circadian_phase)
    ke = model.get_ke(t, circadian_phase)

    dA_gut = -ka * A_gut
    dA_central = ka * A_gut - ke * A_central

    return np.array([dA_gut, dA_central])


def simulate_pk(
    model: PKModel,
    doses: List[float],
    dose_times: List[float],
    t_end: float = 72.0,
    dt: float = 0.1,
    circadian_phase: float = 0.0,
) -> Dict[str, Any]:
    """
    Simulate pharmacokinetic profile.

    Args:
        model: PKModel with parameters
        doses: List of dose amounts (mg)
        dose_times: List of dosing times (hours)
        t_end: Simulation end time
        dt: Output time step
        circadian_phase: Phase of circadian clock at t=0

    Returns:
        Dictionary with concentration profile and metrics
    """
    t_eval = np.arange(0, t_end + dt, dt)

    # Initialize
    concentration = np.zeros_like(t_eval)
    A_central_total = np.zeros_like(t_eval)

    # Simulate each dose independently (superposition for linear PK)
    for dose, dose_time in zip(doses, dose_times):
        if dose_time >= t_end:
            continue

        # Initial condition: dose in gut
        y0 = np.array([dose * model.params.F, 0.0])

        # Time span for this dose
        t_span = (dose_time, t_end)
        t_eval_dose = t_eval[t_eval >= dose_time]

        sol = solve_ivp(
            fun=lambda t, y: _pk_ode(t, y, model, dose_times, doses, circadian_phase),
            t_span=t_span,
            y0=y0,
            t_eval=t_eval_dose,
            method='RK45',
        )

        # Add contribution to total
        for i, t in enumerate(sol.t):
            idx = np.argmin(np.abs(t_eval - t))
            A_central_total[idx] += sol.y[1, i]

    # Convert to concentration
    concentration = A_central_total / model.params.V

    # Compute metrics
    cmax = float(np.max(concentration))
    tmax = float(t_eval[np.argmax(concentration)])
    auc = float(np.trapz(concentration, t_eval))

    # Find time above threshold (if relevant)
    threshold = cmax * 0.5  # 50% of Cmax
    time_above_threshold = float(np.sum(concentration > threshold) * dt)

    return {
        "t": t_eval,
        "concentration": concentration,
        "cmax": cmax,
        "tmax": tmax,
        "auc": auc,
        "time_above_half_cmax": time_above_threshold,
        "doses": doses,
        "dose_times": dose_times,
        "circadian_phase": circadian_phase,
    }


def compute_auc(
    t: np.ndarray,
    concentration: np.ndarray,
    t_start: float = 0.0,
    t_end: Optional[float] = None,
) -> float:
    """
    Compute area under the concentration curve.

    Args:
        t: Time array
        concentration: Concentration array
        t_start: Start time for AUC calculation
        t_end: End time (default: last time point)

    Returns:
        AUC value
    """
    if t_end is None:
        t_end = t[-1]

    mask = (t >= t_start) & (t <= t_end)
    return float(np.trapz(concentration[mask], t[mask]))


def compute_cmax(concentration: np.ndarray) -> float:
    """Compute maximum concentration."""
    return float(np.max(concentration))


def compare_dosing_times(
    model: PKModel,
    dose: float,
    dosing_times: List[float],
    t_end: float = 48.0,
) -> Dict[str, Any]:
    """
    Compare PK profiles for different dosing times.

    Useful for identifying optimal time-of-day for dosing.

    Args:
        model: PKModel
        dose: Dose amount
        dosing_times: List of times to compare (hours from midnight)
        t_end: Simulation duration

    Returns:
        Comparison of PK metrics at each dosing time
    """
    results = []

    for dose_time in dosing_times:
        # Circadian phase at dosing time
        circadian_phase = OMEGA_CIRCADIAN * dose_time

        pk_result = simulate_pk(
            model=model,
            doses=[dose],
            dose_times=[dose_time],
            t_end=t_end,
            circadian_phase=circadian_phase,
        )

        results.append({
            "dose_time": dose_time,
            "dose_time_clock": f"{int(dose_time % 24):02d}:{int((dose_time % 1) * 60):02d}",
            "cmax": pk_result["cmax"],
            "tmax": pk_result["tmax"],
            "auc": pk_result["auc"],
        })

    # Find optimal
    cmax_values = [r["cmax"] for r in results]
    auc_values = [r["auc"] for r in results]

    optimal_cmax_time = dosing_times[np.argmax(cmax_values)]
    optimal_auc_time = dosing_times[np.argmax(auc_values)]

    return {
        "results": results,
        "optimal_cmax_time": optimal_cmax_time,
        "optimal_auc_time": optimal_auc_time,
        "cmax_variation": (max(cmax_values) - min(cmax_values)) / np.mean(cmax_values),
        "auc_variation": (max(auc_values) - min(auc_values)) / np.mean(auc_values),
    }


# Pre-defined drug models


def get_drug_model(drug_name: str) -> PKModel:
    """
    Get pre-configured PK model for common drugs.

    Args:
        drug_name: Name of drug

    Returns:
        PKModel configured for the drug
    """
    models = {
        "melatonin": PKModel(
            params=PKParams(
                ka=2.0,      # Fast absorption
                ke=0.5,      # Short half-life (~1.5h)
                V=35.0,
                F=0.15,      # Low oral bioavailability
                ka_circadian_amp=0.2,
                ke_circadian_amp=0.1,
            ),
            drug_name="melatonin",
        ),
        "tacrolimus": PKModel(
            params=PKParams(
                ka=1.5,
                ke=0.07,     # ~10h half-life
                V=1000.0,    # Large Vd
                F=0.25,
                ka_circadian_amp=0.3,  # Strong circadian variation
                ke_circadian_amp=0.2,
            ),
            drug_name="tacrolimus",
        ),
        "atorvastatin": PKModel(
            params=PKParams(
                ka=1.0,
                ke=0.05,     # ~14h half-life
                V=380.0,
                F=0.14,
                ka_circadian_amp=0.15,
            ),
            drug_name="atorvastatin",
        ),
        "metformin": PKModel(
            params=PKParams(
                ka=0.5,      # Slower absorption
                ke=0.12,     # ~6h half-life
                V=650.0,
                F=0.55,
                ka_circadian_amp=0.1,
            ),
            drug_name="metformin",
        ),
    }

    if drug_name.lower() not in models:
        raise ValueError(f"Unknown drug: {drug_name}. Available: {list(models.keys())}")

    return models[drug_name.lower()]
