"""
PhaseLab Population Variability: Individual variation in drug response.

Models inter-individual variability in:
- Pharmacokinetic parameters
- Circadian rhythm characteristics
- Drug response
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Optional, Any

from .pharmacokinetics import PKModel, PKParams, simulate_pk
from .chronotherapy import ChronotherapyOptimizer, DoseSchedule
from ..core.coherence import coherence_score, go_no_go


@dataclass
class PopulationPK:
    """
    Population pharmacokinetic model with inter-individual variability.

    Models variability using log-normal distributions for PK parameters.
    """

    # Base model
    base_model: PKModel

    # Inter-individual variability (coefficient of variation)
    cv_ka: float = 0.3      # 30% CV in absorption
    cv_ke: float = 0.25     # 25% CV in elimination
    cv_V: float = 0.2       # 20% CV in volume
    cv_F: float = 0.3       # 30% CV in bioavailability

    # Circadian variability
    cv_circadian_amp: float = 0.4  # Variability in circadian amplitude
    cv_circadian_phase: float = 0.2  # Variability in circadian phase (fraction of 24h)

    def sample_individual(self, random_state: Optional[int] = None) -> PKModel:
        """
        Sample an individual's PK parameters from population distribution.

        Args:
            random_state: Random seed for reproducibility

        Returns:
            PKModel with sampled parameters
        """
        if random_state is not None:
            np.random.seed(random_state)

        # Sample from log-normal distributions
        base = self.base_model.params

        ka = base.ka * np.exp(np.random.normal(0, self.cv_ka))
        ke = base.ke * np.exp(np.random.normal(0, self.cv_ke))
        V = base.V * np.exp(np.random.normal(0, self.cv_V))
        F = np.clip(base.F * np.exp(np.random.normal(0, self.cv_F)), 0, 1)

        # Sample circadian variability
        ka_circ_amp = base.ka_circadian_amp * np.exp(
            np.random.normal(0, self.cv_circadian_amp)
        )
        ke_circ_amp = base.ke_circadian_amp * np.exp(
            np.random.normal(0, self.cv_circadian_amp)
        )

        # Phase variability
        phase_shift = np.random.normal(0, self.cv_circadian_phase * 2 * np.pi)

        params = PKParams(
            ka=ka,
            ke=ke,
            V=V,
            F=F,
            ka_circadian_amp=ka_circ_amp,
            ke_circadian_amp=ke_circ_amp,
            ka_circadian_phase=base.ka_circadian_phase + phase_shift,
            ke_circadian_phase=base.ke_circadian_phase + phase_shift,
        )

        return PKModel(params=params, drug_name=self.base_model.drug_name)


def simulate_population(
    pop_model: PopulationPK,
    doses: List[float],
    dose_times: List[float],
    n_individuals: int = 100,
    t_end: float = 48.0,
    random_seed: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Simulate drug response across a population.

    Args:
        pop_model: PopulationPK model
        doses: Dose amounts
        dose_times: Dosing times
        n_individuals: Number of individuals to simulate
        t_end: Simulation duration
        random_seed: Random seed

    Returns:
        Population simulation results
    """
    if random_seed is not None:
        np.random.seed(random_seed)

    individual_results = []

    for i in range(n_individuals):
        # Sample individual
        ind_model = pop_model.sample_individual(random_state=random_seed + i if random_seed else None)

        # Simulate PK
        result = simulate_pk(
            model=ind_model,
            doses=doses,
            dose_times=dose_times,
            t_end=t_end,
        )

        individual_results.append({
            "individual": i,
            "cmax": result["cmax"],
            "tmax": result["tmax"],
            "auc": result["auc"],
            "ka": ind_model.params.ka,
            "ke": ind_model.params.ke,
            "V": ind_model.params.V,
            "F": ind_model.params.F,
        })

    # Population statistics
    cmax_values = [r["cmax"] for r in individual_results]
    auc_values = [r["auc"] for r in individual_results]
    tmax_values = [r["tmax"] for r in individual_results]

    # Coherence of population response
    # Lower variance = higher coherence
    cmax_cv = np.std(cmax_values) / np.mean(cmax_values)
    auc_cv = np.std(auc_values) / np.mean(auc_values)

    # Map CV to coherence (lower CV = higher coherence)
    pop_coherence = np.exp(-np.mean([cmax_cv, auc_cv]))

    return {
        "individuals": individual_results,
        "population_stats": {
            "cmax_mean": np.mean(cmax_values),
            "cmax_std": np.std(cmax_values),
            "cmax_cv": cmax_cv,
            "auc_mean": np.mean(auc_values),
            "auc_std": np.std(auc_values),
            "auc_cv": auc_cv,
            "tmax_mean": np.mean(tmax_values),
            "tmax_std": np.std(tmax_values),
        },
        "population_coherence": pop_coherence,
        "go_no_go": go_no_go(pop_coherence),
        "n_individuals": n_individuals,
    }


def individual_optimization(
    pop_model: PopulationPK,
    individual_params: Dict[str, float],
    total_dose: float,
    frequency: str = "once_daily",
) -> DoseSchedule:
    """
    Optimize dosing for a specific individual.

    Args:
        pop_model: PopulationPK for drug
        individual_params: Individual's measured parameters (ka, ke, etc.)
        total_dose: Total daily dose
        frequency: Dosing frequency

    Returns:
        Personalized DoseSchedule
    """
    # Create individual model
    base = pop_model.base_model.params

    ind_params = PKParams(
        ka=individual_params.get("ka", base.ka),
        ke=individual_params.get("ke", base.ke),
        V=individual_params.get("V", base.V),
        F=individual_params.get("F", base.F),
        ka_circadian_amp=individual_params.get("ka_circadian_amp", base.ka_circadian_amp),
        ke_circadian_amp=individual_params.get("ke_circadian_amp", base.ke_circadian_amp),
    )

    ind_model = PKModel(params=ind_params, drug_name=pop_model.base_model.drug_name)

    # Optimize for this individual
    optimizer = ChronotherapyOptimizer(ind_model)

    if frequency == "once_daily":
        return optimizer.optimize_single_dose(total_dose)
    elif frequency == "twice_daily":
        return optimizer.optimize_twice_daily(total_dose)
    else:
        raise ValueError(f"Unknown frequency: {frequency}")


def chronotype_adjustment(
    base_schedule: DoseSchedule,
    chronotype: str = "intermediate",
) -> DoseSchedule:
    """
    Adjust dosing schedule based on chronotype.

    Args:
        base_schedule: Base optimized schedule
        chronotype: "early", "intermediate", or "late"

    Returns:
        Adjusted DoseSchedule
    """
    # Chronotype phase shifts (hours)
    shifts = {
        "early": -2.0,      # Early birds: shift earlier
        "intermediate": 0.0,
        "late": 2.0,        # Night owls: shift later
    }

    shift = shifts.get(chronotype, 0.0)

    # Adjust times
    adjusted_times = [(t + shift) % 24 for t in base_schedule.times]

    return DoseSchedule(
        drug_name=base_schedule.drug_name,
        doses=base_schedule.doses,
        times=adjusted_times,
        frequency=base_schedule.frequency,
        predicted_cmax=base_schedule.predicted_cmax,
        predicted_auc=base_schedule.predicted_auc,
        predicted_efficacy=base_schedule.predicted_efficacy,
        chronotherapy_R=base_schedule.chronotherapy_R,
        go_no_go=base_schedule.go_no_go,
    )


def therapeutic_drug_monitoring(
    pop_model: PopulationPK,
    measured_concentrations: List[float],
    measurement_times: List[float],
    dose_history: List[Tuple[float, float]],  # (dose, time) pairs
) -> Dict[str, Any]:
    """
    Bayesian estimation of individual parameters from TDM data.

    Simplified implementation using least squares fitting.

    Args:
        pop_model: Population PK model
        measured_concentrations: Measured drug concentrations
        measurement_times: Times of measurements
        dose_history: List of (dose, time) pairs

    Returns:
        Estimated individual parameters and predictions
    """
    from scipy.optimize import minimize

    doses = [d[0] for d in dose_history]
    dose_times = [d[1] for d in dose_history]

    def objective(params):
        """Objective function: sum of squared errors."""
        ka, ke = params

        # Create model with these parameters
        test_params = PKParams(
            ka=ka,
            ke=ke,
            V=pop_model.base_model.params.V,
            F=pop_model.base_model.params.F,
        )
        test_model = PKModel(params=test_params)

        # Simulate
        result = simulate_pk(
            model=test_model,
            doses=doses,
            dose_times=dose_times,
            t_end=max(measurement_times) + 12,
        )

        # Interpolate to measurement times
        predicted = np.interp(measurement_times, result["t"], result["concentration"])

        # Sum of squared errors
        sse = np.sum((predicted - measured_concentrations) ** 2)

        return sse

    # Initial guess from population mean
    x0 = [pop_model.base_model.params.ka, pop_model.base_model.params.ke]

    # Bounds
    bounds = [(0.1, 10.0), (0.01, 1.0)]

    # Optimize
    result = minimize(objective, x0, bounds=bounds, method='L-BFGS-B')

    estimated_ka, estimated_ke = result.x

    # Predict future concentrations
    pred_params = PKParams(
        ka=estimated_ka,
        ke=estimated_ke,
        V=pop_model.base_model.params.V,
        F=pop_model.base_model.params.F,
    )
    pred_model = PKModel(params=pred_params)

    prediction = simulate_pk(
        model=pred_model,
        doses=doses,
        dose_times=dose_times,
        t_end=72,
    )

    return {
        "estimated_ka": estimated_ka,
        "estimated_ke": estimated_ke,
        "estimated_half_life": np.log(2) / estimated_ke,
        "fit_error": result.fun,
        "predicted_profile": {
            "t": prediction["t"],
            "concentration": prediction["concentration"],
        },
        "predicted_cmax": prediction["cmax"],
        "predicted_auc": prediction["auc"],
    }
