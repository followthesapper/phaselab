"""
Circadian Rescue Simulation for Smith-Magenis Syndrome.

This trial simulates how partial RAI1 rescue restores circadian
coherence using the SMS clock model.

Key question: Given a predicted expression boost from CRISPRa,
does the circadian system recover to a therapeutic state?

The simulation:
1. Takes expected RAI1 expression level from CRISPRa trial
2. Runs the SMS clock model at that expression level
3. Evaluates synchronization (R̄), sleep quality prediction
4. Reports whether circadian rescue is predicted
"""

import numpy as np
from typing import Dict, Any, List, Optional
import logging

from .core import (
    SMSTrialResult,
    SMSTrialConfig,
    TrialType,
    TrialStatus,
)

logger = logging.getLogger(__name__)


def run_circadian_rescue_simulation(
    predicted_rai1_expression: float = 0.80,
    config: Optional[SMSTrialConfig] = None,
) -> SMSTrialResult:
    """
    Simulate circadian rescue with predicted RAI1 expression boost.

    Uses the SMS clock model to predict whether the expected
    expression change from CRISPRa treatment will restore
    circadian synchronization.

    Args:
        predicted_rai1_expression: Expected RAI1 expression level
            as fraction of normal (0.5 = SMS baseline, 1.0 = normal).
        config: Trial configuration.

    Returns:
        SMSTrialResult with circadian rescue predictions.

    Example:
        >>> from phaselab.trials.sms import run_circadian_rescue_simulation
        >>> result = run_circadian_rescue_simulation(predicted_rai1_expression=0.80)
        >>> print(f"Final R̄: {result.metrics['final_R_bar']:.3f}")
        >>> print(f"Sleep quality: {result.metrics['sleep_quality_prediction']}")
    """
    if config is None:
        config = SMSTrialConfig()

    if config.verbose:
        print("=" * 60)
        print("SMS Circadian Rescue Simulation")
        print("=" * 60)
        print(f"Baseline RAI1: {config.baseline_expression:.0%}")
        print(f"Predicted RAI1: {predicted_rai1_expression:.0%}")
        print(f"Simulation duration: {config.simulation_hours}h ({config.simulation_hours/24:.0f} days)")
        print()

    warnings = []
    errors = []

    try:
        # Import SMS clock model
        from phaselab.circadian.sms_model import (
            simulate_sms_clock,
            therapeutic_scan,
            predict_sleep_quality,
            SMSClockParams,
        )
        from phaselab.core.constants import E_MINUS_2
        from phaselab.fusion import ClaimLevel

        # Run baseline simulation (SMS state)
        if config.verbose:
            print("Running baseline simulation (SMS state)...")

        baseline_results = []
        for trial in range(config.n_circadian_trials):
            result = simulate_sms_clock(
                rai1_level=config.baseline_expression,
                t_end=config.simulation_hours,
                random_seed=42 + trial,
            )
            baseline_results.append(result)

        baseline_R_bar = np.mean([r['final_R_bar'] for r in baseline_results])
        baseline_R_bar_std = np.std([r['final_R_bar'] for r in baseline_results])

        if config.verbose:
            print(f"  Baseline R̄: {baseline_R_bar:.3f} ± {baseline_R_bar_std:.3f}")

        # Run treated simulation
        if config.verbose:
            print(f"Running treated simulation (RAI1 = {predicted_rai1_expression:.0%})...")

        treated_results = []
        for trial in range(config.n_circadian_trials):
            result = simulate_sms_clock(
                rai1_level=predicted_rai1_expression,
                t_end=config.simulation_hours,
                random_seed=42 + trial,
            )
            treated_results.append(result)

        treated_R_bar = np.mean([r['final_R_bar'] for r in treated_results])
        treated_R_bar_std = np.std([r['final_R_bar'] for r in treated_results])

        if config.verbose:
            print(f"  Treated R̄: {treated_R_bar:.3f} ± {treated_R_bar_std:.3f}")

        # Calculate improvement
        R_bar_improvement = treated_R_bar - baseline_R_bar
        relative_improvement = R_bar_improvement / baseline_R_bar if baseline_R_bar > 0 else 0

        # Therapeutic scan to find optimal RAI1 level
        if config.verbose:
            print("Running therapeutic scan...")

        scan_result = therapeutic_scan(
            rai1_levels=np.linspace(0.3, 1.2, 10).tolist(),
            t_end=config.simulation_hours,
            n_trials=3,
        )

        # Classifications
        baseline_classification = baseline_results[0]['classification']
        treated_classification = treated_results[0]['classification']

        # Sleep quality predictions
        baseline_sleep = predict_sleep_quality(baseline_R_bar)
        treated_sleep = predict_sleep_quality(treated_R_bar)

        # Determine rescue status
        if treated_R_bar >= 0.9:
            rescue_status = "FULL_RESCUE"
            rescue_description = "Full circadian rescue predicted. Synchronization restored to near-normal levels."
        elif treated_R_bar >= 0.7:
            rescue_status = "PARTIAL_RESCUE"
            rescue_description = "Partial circadian rescue predicted. Significant improvement but not fully normalized."
        elif treated_R_bar > baseline_R_bar + 0.1:
            rescue_status = "MODEST_IMPROVEMENT"
            rescue_description = "Modest improvement predicted. Some benefit but circadian disruption persists."
        else:
            rescue_status = "NO_RESCUE"
            rescue_description = "No significant circadian rescue predicted at this expression level."

        if config.verbose:
            print(f"\nRescue status: {rescue_status}")
            print(f"  {rescue_description}")
            print(f"  Baseline sleep: {baseline_sleep}")
            print(f"  Treated sleep: {treated_sleep}")

        # Build candidates (time series snapshots)
        candidates = []

        # Candidate 1: Predicted treatment outcome
        candidates.append({
            'scenario': 'predicted_treatment',
            'rai1_level': predicted_rai1_expression,
            'final_R_bar': treated_R_bar,
            'final_R_bar_std': treated_R_bar_std,
            'classification': treated_classification,
            'sleep_quality': treated_sleep,
            'rescue_status': rescue_status,
        })

        # Candidate 2: Optimal treatment from scan
        if scan_result.get('optimal_level'):
            optimal_level = scan_result['optimal_level']
            optimal_idx = scan_result['levels'].index(
                min(scan_result['levels'], key=lambda x: abs(x - optimal_level))
            )
            optimal_R_bar = scan_result['R_bars'][optimal_idx]

            candidates.append({
                'scenario': 'optimal_expression',
                'rai1_level': optimal_level,
                'final_R_bar': optimal_R_bar,
                'classification': scan_result['classifications'][optimal_idx],
                'sleep_quality': predict_sleep_quality(optimal_R_bar),
                'rescue_status': 'FULL_RESCUE' if optimal_R_bar >= 0.9 else 'PARTIAL_RESCUE',
            })

        # Candidate 3: Therapeutic window bounds
        if scan_result.get('therapeutic_window'):
            tw = scan_result['therapeutic_window']
            candidates.append({
                'scenario': 'therapeutic_window',
                'rai1_level_min': tw[0],
                'rai1_level_max': tw[1],
                'required_boost': scan_result.get('required_boost', 0),
                'rescue_status': 'THERAPEUTIC_WINDOW_FOUND',
            })

        # Determine overall claim level
        if rescue_status in ["FULL_RESCUE", "PARTIAL_RESCUE"]:
            if treated_R_bar_std < 0.1:
                overall_claim = ClaimLevel.CONTEXT_DEPENDENT
                claim_desc = (
                    f"Context-dependent evidence for circadian rescue. "
                    f"Simulations consistently predict {rescue_status.lower().replace('_', ' ')} "
                    f"with R̄ = {treated_R_bar:.3f} at {predicted_rai1_expression:.0%} RAI1."
                )
            else:
                overall_claim = ClaimLevel.EXPLORATORY
                claim_desc = (
                    f"Exploratory evidence for circadian improvement. "
                    f"Simulation variability is moderate (±{treated_R_bar_std:.3f})."
                )
        else:
            overall_claim = ClaimLevel.EXPLORATORY
            claim_desc = (
                f"Exploratory evidence. Model predicts {rescue_status.lower().replace('_', ' ')}. "
                f"Higher RAI1 boost may be needed for circadian rescue."
            )

        # Warnings
        if predicted_rai1_expression > 1.1:
            warnings.append(
                f"Predicted expression ({predicted_rai1_expression:.0%}) exceeds 110% of normal. "
                f"Risk of overexpression effects not captured by circadian model."
            )

        if treated_R_bar_std > 0.15:
            warnings.append(
                f"High simulation variability (±{treated_R_bar_std:.3f}). "
                f"Circadian outcome may be sensitive to initial conditions."
            )

        # Standard model limitation warning
        warnings.append(
            "Circadian model is a simplified 2-oscillator Kuramoto system. "
            "Real SMS pathophysiology involves additional factors not captured here."
        )

        # Build metrics
        metrics = {
            # Baseline
            'baseline_rai1': config.baseline_expression,
            'baseline_R_bar': baseline_R_bar,
            'baseline_R_bar_std': baseline_R_bar_std,
            'baseline_classification': baseline_classification,
            'baseline_sleep_quality': baseline_sleep,

            # Treated
            'predicted_rai1': predicted_rai1_expression,
            'final_R_bar': treated_R_bar,
            'final_R_bar_std': treated_R_bar_std,
            'treated_classification': treated_classification,
            'sleep_quality_prediction': treated_sleep,

            # Improvement
            'R_bar_improvement': R_bar_improvement,
            'relative_improvement': relative_improvement,
            'rescue_status': rescue_status,

            # Therapeutic scan
            'therapeutic_window': scan_result.get('therapeutic_window'),
            'optimal_level': scan_result.get('optimal_level'),
            'required_boost': scan_result.get('required_boost'),

            # Simulation params
            'simulation_hours': config.simulation_hours,
            'n_trials': config.n_circadian_trials,
        }

        # Build summary
        summary = (
            f"Circadian Rescue Simulation: {rescue_status}. "
            f"RAI1 boost from {config.baseline_expression:.0%} to {predicted_rai1_expression:.0%} "
            f"predicts R̄ improvement from {baseline_R_bar:.3f} to {treated_R_bar:.3f}. "
            f"Sleep quality: '{treated_sleep.split(' - ')[0]}'. "
            f"Claim level: {overall_claim.value}."
        )

        return SMSTrialResult(
            trial_type=TrialType.CIRCADIAN_SIMULATION,
            status=TrialStatus.COMPLETED,
            summary=summary,
            candidates=candidates,
            best_candidate=candidates[0] if candidates else None,
            metrics=metrics,
            claim_level=overall_claim.value,
            claim_description=claim_desc,
            warnings=warnings,
            errors=errors,
            metadata={
                'model': 'SMS_Kuramoto',
                'modality': 'CircadianSimulation',
                'baseline_results': [r['final_R_bar'] for r in baseline_results],
                'treated_results': [r['final_R_bar'] for r in treated_results],
            },
        )

    except Exception as e:
        logger.error(f"Circadian simulation failed: {e}")
        return SMSTrialResult(
            trial_type=TrialType.CIRCADIAN_SIMULATION,
            status=TrialStatus.FAILED,
            summary=f"Circadian Rescue Simulation failed: {str(e)}",
            errors=[str(e)],
            claim_level="unknown",
            claim_description="Simulation execution failed.",
        )
