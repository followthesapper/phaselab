"""
PhaseLab Chronotherapy: Optimal timing of drug administration.

Provides tools for:
- Finding optimal dosing times based on circadian rhythms
- Accounting for drug-clock interactions
- Personalizing chronotherapy schedules
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Any
from scipy.optimize import minimize_scalar

from .pharmacokinetics import PKModel, simulate_pk
from ..core.constants import OMEGA_CIRCADIAN, E_MINUS_2
from ..core.coherence import go_no_go


@dataclass
class DoseSchedule:
    """Optimized dosing schedule."""

    drug_name: str
    doses: List[float]
    times: List[float]  # Hours from midnight
    frequency: str  # "once_daily", "twice_daily", etc.

    # Predicted outcomes
    predicted_cmax: float
    predicted_auc: float
    predicted_efficacy: float

    # Coherence with circadian rhythm
    chronotherapy_R: float
    go_no_go: str

    def to_dict(self) -> Dict[str, Any]:
        return {
            "drug_name": self.drug_name,
            "doses": self.doses,
            "times": self.times,
            "times_formatted": [f"{int(t):02d}:{int((t % 1) * 60):02d}" for t in self.times],
            "frequency": self.frequency,
            "predicted_cmax": self.predicted_cmax,
            "predicted_auc": self.predicted_auc,
            "predicted_efficacy": self.predicted_efficacy,
            "chronotherapy_R": self.chronotherapy_R,
            "go_no_go": self.go_no_go,
        }


class ChronotherapyOptimizer:
    """
    Optimizer for chronotherapy dosing schedules.

    Finds optimal dosing times that:
    1. Maximize drug efficacy
    2. Minimize side effects
    3. Align with patient's circadian rhythm
    """

    def __init__(
        self,
        model: PKModel,
        target_tissue: str = "generic",
        efficacy_window: Tuple[float, float] = (8, 20),  # Hours when efficacy matters
        side_effect_window: Tuple[float, float] = (0, 8),  # Hours to minimize exposure
    ):
        """
        Initialize optimizer.

        Args:
            model: PKModel for the drug
            target_tissue: Target tissue for drug action
            efficacy_window: Time window for desired drug effect
            side_effect_window: Time window to minimize exposure
        """
        self.model = model
        self.target_tissue = target_tissue
        self.efficacy_window = efficacy_window
        self.side_effect_window = side_effect_window

    def optimize_single_dose(
        self,
        dose: float,
        t_end: float = 48.0,
    ) -> DoseSchedule:
        """
        Find optimal time for single daily dose.

        Args:
            dose: Dose amount
            t_end: Simulation duration

        Returns:
            DoseSchedule with optimal timing
        """
        def objective(dose_time):
            # Simulate PK
            result = simulate_pk(
                model=self.model,
                doses=[dose],
                dose_times=[dose_time],
                t_end=t_end,
                circadian_phase=OMEGA_CIRCADIAN * dose_time,
            )

            # Score based on:
            # 1. AUC in efficacy window (maximize)
            # 2. AUC in side effect window (minimize)
            # 3. Cmax (don't exceed limits)

            t = result["t"]
            conc = result["concentration"]

            # AUC in efficacy window
            eff_mask = (t >= self.efficacy_window[0]) & (t <= self.efficacy_window[1])
            auc_efficacy = np.trapz(conc[eff_mask], t[eff_mask]) if np.any(eff_mask) else 0

            # AUC in side effect window
            se_mask = (t >= self.side_effect_window[0]) & (t <= self.side_effect_window[1])
            auc_side_effect = np.trapz(conc[se_mask], t[se_mask]) if np.any(se_mask) else 0

            # Objective: maximize efficacy, minimize side effects
            # Negative because minimize_scalar minimizes
            score = -(auc_efficacy - 0.5 * auc_side_effect)

            return score

        # Search over 24 hours
        result = minimize_scalar(objective, bounds=(0, 24), method='bounded')
        optimal_time = result.x

        # Get final PK profile at optimal time
        pk_result = simulate_pk(
            model=self.model,
            doses=[dose],
            dose_times=[optimal_time],
            t_end=t_end,
            circadian_phase=OMEGA_CIRCADIAN * optimal_time,
        )

        # Calculate chronotherapy coherence
        # Higher coherence if dosing aligns with circadian-optimal window
        chrono_R = self._calculate_chrono_coherence(optimal_time, pk_result)

        return DoseSchedule(
            drug_name=self.model.drug_name,
            doses=[dose],
            times=[optimal_time],
            frequency="once_daily",
            predicted_cmax=pk_result["cmax"],
            predicted_auc=pk_result["auc"],
            predicted_efficacy=self._estimate_efficacy(pk_result),
            chronotherapy_R=chrono_R,
            go_no_go=go_no_go(chrono_R),
        )

    def optimize_twice_daily(
        self,
        total_dose: float,
        dose_split: float = 0.5,
        t_end: float = 48.0,
    ) -> DoseSchedule:
        """
        Find optimal times for twice-daily dosing.

        Args:
            total_dose: Total daily dose
            dose_split: Fraction of dose in first administration
            t_end: Simulation duration

        Returns:
            DoseSchedule with optimal timing
        """
        dose1 = total_dose * dose_split
        dose2 = total_dose * (1 - dose_split)

        best_score = float('inf')
        best_times = (8, 20)

        # Grid search over possible time combinations
        for t1 in range(6, 14):  # Morning dose 6am-2pm
            for t2 in range(14, 24):  # Evening dose 2pm-midnight
                result = simulate_pk(
                    model=self.model,
                    doses=[dose1, dose2],
                    dose_times=[t1, t2],
                    t_end=t_end,
                    circadian_phase=0,
                )

                # Score: stable levels + efficacy
                conc = result["concentration"]
                variance = np.var(conc[conc > 0])
                auc = result["auc"]

                score = variance - 0.1 * auc  # Minimize variance, maximize AUC

                if score < best_score:
                    best_score = score
                    best_times = (t1, t2)

        # Final simulation at optimal times
        pk_result = simulate_pk(
            model=self.model,
            doses=[dose1, dose2],
            dose_times=list(best_times),
            t_end=t_end,
        )

        chrono_R = self._calculate_chrono_coherence(best_times[0], pk_result)

        return DoseSchedule(
            drug_name=self.model.drug_name,
            doses=[dose1, dose2],
            times=list(best_times),
            frequency="twice_daily",
            predicted_cmax=pk_result["cmax"],
            predicted_auc=pk_result["auc"],
            predicted_efficacy=self._estimate_efficacy(pk_result),
            chronotherapy_R=chrono_R,
            go_no_go=go_no_go(chrono_R),
        )

    def _calculate_chrono_coherence(
        self,
        dose_time: float,
        pk_result: Dict[str, Any],
    ) -> float:
        """Calculate coherence of dosing with circadian rhythm."""
        # Phase of dosing relative to circadian cycle
        dose_phase = (dose_time / 24.0) * 2 * np.pi

        # Ideal phase for this drug's target tissue
        # (simplified - real implementation would use tissue-specific data)
        ideal_phase = np.pi / 2  # Default: morning dosing

        # Phase coherence
        phase_diff = np.abs(dose_phase - ideal_phase)
        R_phase = np.cos(phase_diff / 2) ** 2

        # Exposure stability coherence
        conc = pk_result["concentration"]
        if np.max(conc) > 0:
            conc_norm = conc / np.max(conc)
            R_stability = 1.0 - np.std(conc_norm[conc_norm > 0.1])
        else:
            R_stability = 0

        # Combined coherence
        R_chrono = (R_phase + R_stability) / 2

        return float(np.clip(R_chrono, 0, 1))

    def _estimate_efficacy(self, pk_result: Dict[str, Any]) -> float:
        """Estimate drug efficacy based on PK profile."""
        # Simplified efficacy model based on AUC
        auc = pk_result["auc"]
        cmax = pk_result["cmax"]

        # Assume efficacy scales with exposure up to a point
        # Sigmoid relationship
        ec50 = cmax * 0.5  # 50% of max concentration
        efficacy = auc / (auc + ec50 * 24)  # Normalize by time

        return float(np.clip(efficacy, 0, 1))


def optimize_dosing_time(
    model: PKModel,
    dose: float,
    frequency: str = "once_daily",
    target_tissue: str = "generic",
) -> DoseSchedule:
    """
    Convenience function for dosing time optimization.

    Args:
        model: PKModel
        dose: Total daily dose
        frequency: "once_daily" or "twice_daily"
        target_tissue: Target tissue

    Returns:
        Optimized DoseSchedule
    """
    optimizer = ChronotherapyOptimizer(model, target_tissue)

    if frequency == "once_daily":
        return optimizer.optimize_single_dose(dose)
    elif frequency == "twice_daily":
        return optimizer.optimize_twice_daily(dose)
    else:
        raise ValueError(f"Unknown frequency: {frequency}")


def circadian_drug_interaction(
    drug_model: PKModel,
    circadian_amplitude: float = 1.0,
    circadian_phase: float = 0.0,
) -> Dict[str, Any]:
    """
    Analyze interaction between drug and circadian rhythm.

    Args:
        drug_model: PKModel
        circadian_amplitude: Strength of circadian rhythm (0-1)
        circadian_phase: Phase offset of circadian rhythm

    Returns:
        Analysis of drug-clock interaction
    """
    # Simulate at different times of day
    times = np.arange(0, 24, 2)
    results = []

    for t in times:
        pk = simulate_pk(
            model=drug_model,
            doses=[100],  # Standard dose
            dose_times=[t],
            t_end=48,
            circadian_phase=circadian_phase,
        )
        results.append({
            "dose_time": t,
            "cmax": pk["cmax"],
            "auc": pk["auc"],
        })

    # Analyze variation
    cmax_values = [r["cmax"] for r in results]
    auc_values = [r["auc"] for r in results]

    cmax_variation = (max(cmax_values) - min(cmax_values)) / np.mean(cmax_values)
    auc_variation = (max(auc_values) - min(auc_values)) / np.mean(auc_values)

    # Determine if chronotherapy is warranted
    chronotherapy_benefit = cmax_variation > 0.2 or auc_variation > 0.2

    return {
        "time_profiles": results,
        "cmax_variation": cmax_variation,
        "auc_variation": auc_variation,
        "best_time_cmax": times[np.argmax(cmax_values)],
        "best_time_auc": times[np.argmax(auc_values)],
        "worst_time_cmax": times[np.argmin(cmax_values)],
        "chronotherapy_recommended": chronotherapy_benefit,
    }
