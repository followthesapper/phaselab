"""
PhaseLab Dosage: Therapeutic window and expression level optimization.

Implements:
- Expression level estimation from guide binding quality
- Therapeutic window calculation for haploinsufficiency disorders
- Dosage-response modeling (Hill equation)
- Toxicity threshold estimation
- Optimal intervention timing
- IR coherence integration for reliability scoring
"""

import numpy as np
from typing import Dict, Any, Optional, List, Tuple
from dataclasses import dataclass, field
from scipy.optimize import minimize_scalar, brentq
from scipy.integrate import odeint

from ..core.coherence import coherence_score, go_no_go
from ..core.constants import E_MINUS_2


@dataclass
class TherapeuticWindow:
    """
    Defines the therapeutic window for a gene therapy target.

    For haploinsufficiency disorders:
    - Baseline: ~50% expression (one functional allele)
    - Therapeutic minimum: Expression level that provides benefit
    - Therapeutic maximum: Expression level before toxicity
    - Optimal: Target expression for maximum benefit
    """
    baseline_expression: float = 0.50      # Baseline as fraction of normal
    therapeutic_min: float = 0.70          # Minimum therapeutic benefit
    therapeutic_max: float = 1.20          # Maximum before overexpression toxicity
    optimal_expression: float = 0.85       # Optimal target
    normal_expression: float = 1.00        # Wild-type level

    # Safety margins
    lower_safety_margin: float = 0.10      # Below therapeutic_min
    upper_safety_margin: float = 0.15      # Above optimal

    # Disease-specific parameters
    disease_name: Optional[str] = None
    gene_symbol: Optional[str] = None

    @property
    def therapeutic_range(self) -> Tuple[float, float]:
        """Return the therapeutic range as a tuple."""
        return (self.therapeutic_min, self.therapeutic_max)

    @property
    def required_upregulation(self) -> float:
        """Required fold-change from baseline to reach therapeutic minimum."""
        return self.therapeutic_min / self.baseline_expression

    @property
    def optimal_upregulation(self) -> float:
        """Required fold-change from baseline to reach optimal."""
        return self.optimal_expression / self.baseline_expression


# Pre-defined therapeutic windows for known disorders
THERAPEUTIC_WINDOWS = {
    'SMS': TherapeuticWindow(
        baseline_expression=0.50,
        therapeutic_min=0.70,
        therapeutic_max=1.10,
        optimal_expression=0.80,
        disease_name="Smith-Magenis Syndrome",
        gene_symbol="RAI1",
    ),
    'SCN2A_NDD': TherapeuticWindow(
        baseline_expression=0.50,
        therapeutic_min=0.65,
        therapeutic_max=1.15,
        optimal_expression=0.85,
        disease_name="SCN2A Neurodevelopmental Disorder",
        gene_symbol="SCN2A",
    ),
    'SHANK3': TherapeuticWindow(
        baseline_expression=0.50,
        therapeutic_min=0.60,
        therapeutic_max=1.10,
        optimal_expression=0.80,
        disease_name="Phelan-McDermid Syndrome",
        gene_symbol="SHANK3",
    ),
    'CHD8': TherapeuticWindow(
        baseline_expression=0.50,
        therapeutic_min=0.65,
        therapeutic_max=1.15,
        optimal_expression=0.85,
        disease_name="CHD8-related ASD",
        gene_symbol="CHD8",
    ),
}


@dataclass
class DosageConfig:
    """Configuration for dosage modeling."""

    # Hill equation parameters
    hill_coefficient: float = 2.0          # Cooperativity
    ec50: float = 0.5                      # Half-maximal concentration
    max_effect: float = 1.0                # Maximum possible effect

    # CRISPRa parameters
    crispra_fold_change_min: float = 1.5   # Minimum expected upregulation
    crispra_fold_change_max: float = 3.0   # Maximum expected upregulation
    crispra_fold_change_typical: float = 2.0

    # Timing parameters
    onset_hours: float = 24.0              # Time to first effect
    peak_hours: float = 72.0               # Time to peak effect
    duration_hours: float = 168.0          # Duration of effect

    # Safety
    max_safe_overexpression: float = 1.5   # 150% of normal

    # Coherence integration
    coherence_weight: float = 0.3          # How much coherence affects prediction


def estimate_expression_change(
    guide_coherence: float,
    binding_energy: float,
    config: Optional[DosageConfig] = None,
) -> Dict[str, float]:
    """
    Estimate expression change from CRISPRa based on guide quality.

    Higher coherence and stronger binding → more reliable upregulation.

    Args:
        guide_coherence: IR coherence R̄ value (0-1).
        binding_energy: ΔG of binding (kcal/mol, negative = stronger).
        config: DosageConfig parameters.

    Returns:
        Dictionary with estimated expression changes.
    """
    if config is None:
        config = DosageConfig()

    # Base fold change depends on binding strength
    # Stronger binding (more negative ΔG) → better upregulation
    binding_factor = min(1.0, max(0, (-binding_energy - 5) / 20))

    base_fold = (
        config.crispra_fold_change_min +
        (config.crispra_fold_change_max - config.crispra_fold_change_min) * binding_factor
    )

    # Coherence modulates reliability of prediction
    # High coherence → prediction more reliable
    reliability = guide_coherence if guide_coherence > E_MINUS_2 else 0.5

    # Final estimated fold change
    estimated_fold = base_fold * (0.7 + 0.3 * reliability)

    # Confidence intervals based on coherence
    uncertainty = 0.3 * (1 - reliability)
    fold_low = estimated_fold * (1 - uncertainty)
    fold_high = estimated_fold * (1 + uncertainty)

    return {
        'estimated_fold_change': estimated_fold,
        'fold_change_low': fold_low,
        'fold_change_high': fold_high,
        'base_fold_change': base_fold,
        'reliability': reliability,
        'binding_factor': binding_factor,
        'confidence': reliability,
    }


def calculate_therapeutic_window(
    disease: str = None,
    baseline: float = 0.50,
    therapeutic_min: float = 0.70,
    therapeutic_max: float = 1.20,
    optimal: float = 0.85,
) -> TherapeuticWindow:
    """
    Calculate therapeutic window for a target.

    Args:
        disease: Pre-defined disease name (e.g., "SMS", "SCN2A_NDD").
        baseline: Baseline expression if custom.
        therapeutic_min: Minimum therapeutic level if custom.
        therapeutic_max: Maximum safe level if custom.
        optimal: Optimal target level if custom.

    Returns:
        TherapeuticWindow object.

    Example:
        >>> window = calculate_therapeutic_window(disease="SMS")
        >>> print(f"Required upregulation: {window.required_upregulation:.1f}x")
    """
    if disease and disease in THERAPEUTIC_WINDOWS:
        return THERAPEUTIC_WINDOWS[disease]

    return TherapeuticWindow(
        baseline_expression=baseline,
        therapeutic_min=therapeutic_min,
        therapeutic_max=therapeutic_max,
        optimal_expression=optimal,
    )


def dosage_response_curve(
    concentrations: np.ndarray,
    ec50: float = 0.5,
    hill_n: float = 2.0,
    max_effect: float = 1.0,
) -> np.ndarray:
    """
    Calculate Hill equation dose-response curve.

    Response = max_effect * C^n / (EC50^n + C^n)

    Args:
        concentrations: Array of concentrations/doses.
        ec50: Half-maximal effective concentration.
        hill_n: Hill coefficient (cooperativity).
        max_effect: Maximum response.

    Returns:
        Array of responses.
    """
    c = np.asarray(concentrations)
    return max_effect * np.power(c, hill_n) / (np.power(ec50, hill_n) + np.power(c, hill_n))


def expression_over_time(
    t: np.ndarray,
    peak_expression: float,
    onset_hours: float = 24.0,
    peak_hours: float = 72.0,
    decay_hours: float = 168.0,
) -> np.ndarray:
    """
    Model expression change over time after CRISPRa intervention.

    Uses a rise-and-decay model.

    Args:
        t: Time points (hours).
        peak_expression: Maximum expression fold-change.
        onset_hours: Time to significant effect.
        peak_hours: Time to peak.
        decay_hours: Time for effect to decay.

    Returns:
        Expression fold-change at each time point.
    """
    t = np.asarray(t)

    # Parameters for rise and decay
    k_rise = 2.0 / onset_hours
    k_decay = 2.0 / (decay_hours - peak_hours)

    expression = np.zeros_like(t, dtype=float)

    # Rising phase
    mask_rise = t < peak_hours
    expression[mask_rise] = 1.0 + (peak_expression - 1.0) * (1 - np.exp(-k_rise * t[mask_rise]))

    # Decay phase
    mask_decay = t >= peak_hours
    t_decay = t[mask_decay] - peak_hours
    expression[mask_decay] = 1.0 + (peak_expression - 1.0) * np.exp(-k_decay * t_decay)

    return expression


def optimize_dosage(
    therapeutic_window: TherapeuticWindow,
    guide_coherence: float,
    binding_energy: float,
    config: Optional[DosageConfig] = None,
) -> Dict[str, Any]:
    """
    Optimize dosage to achieve therapeutic expression levels.

    Args:
        therapeutic_window: Target therapeutic window.
        guide_coherence: IR coherence of selected guide.
        binding_energy: Binding energy of selected guide.
        config: Dosage configuration.

    Returns:
        Optimization results with recommended dosage parameters.
    """
    if config is None:
        config = DosageConfig()

    # Estimate achievable expression change
    expression = estimate_expression_change(guide_coherence, binding_energy, config)
    estimated_fold = expression['estimated_fold_change']

    # Calculate resulting expression level
    baseline = therapeutic_window.baseline_expression
    achieved_expression = baseline * estimated_fold

    # Check if within therapeutic window
    in_window = (
        therapeutic_window.therapeutic_min <=
        achieved_expression <=
        therapeutic_window.therapeutic_max
    )

    # Calculate required dosage adjustment
    target_expression = therapeutic_window.optimal_expression
    required_fold = target_expression / baseline

    # Dosage factor (relative to standard)
    if estimated_fold > 0:
        dosage_factor = required_fold / estimated_fold
    else:
        dosage_factor = 1.0

    # Safety assessment
    safety_margin = therapeutic_window.therapeutic_max - achieved_expression
    risk_of_overexpression = achieved_expression > therapeutic_window.therapeutic_max

    # Time to therapeutic level
    time_profile = expression_over_time(
        np.linspace(0, 168, 100),
        peak_expression=estimated_fold,
        onset_hours=config.onset_hours,
        peak_hours=config.peak_hours,
        decay_hours=config.duration_hours,
    )

    # Find time within therapeutic window
    expression_levels = baseline * time_profile
    in_therapeutic = (
        (expression_levels >= therapeutic_window.therapeutic_min) &
        (expression_levels <= therapeutic_window.therapeutic_max)
    )
    therapeutic_duration = np.sum(in_therapeutic) * 168 / 100  # hours

    return {
        'estimated_fold_change': estimated_fold,
        'achieved_expression': achieved_expression,
        'target_expression': target_expression,
        'in_therapeutic_window': in_window,
        'dosage_factor': dosage_factor,
        'safety_margin': safety_margin,
        'risk_of_overexpression': risk_of_overexpression,
        'therapeutic_duration_hours': therapeutic_duration,
        'reliability': expression['reliability'],
        'recommendation': _get_dosage_recommendation(
            achieved_expression,
            therapeutic_window,
            expression['reliability'],
        ),
    }


def _get_dosage_recommendation(
    achieved: float,
    window: TherapeuticWindow,
    reliability: float,
) -> str:
    """Generate dosage recommendation based on results."""
    if achieved < window.therapeutic_min:
        deficit = window.therapeutic_min - achieved
        return f"INSUFFICIENT: Expression {deficit:.0%} below therapeutic minimum. Consider higher dosage or alternative guide."
    elif achieved > window.therapeutic_max:
        excess = achieved - window.therapeutic_max
        return f"WARNING: Expression {excess:.0%} above safe maximum. Reduce dosage or use weaker guide."
    elif achieved > window.optimal_expression:
        return f"GOOD: Expression within window but {(achieved - window.optimal_expression):.0%} above optimal. Monitor for side effects."
    elif reliability < 0.5:
        return f"UNCERTAIN: Expression likely therapeutic but prediction reliability is low ({reliability:.0%}). Verify experimentally."
    else:
        return f"OPTIMAL: Predicted expression {achieved:.0%} is within therapeutic window with good reliability."


def validate_therapeutic_level(
    expression_level: float,
    therapeutic_window: TherapeuticWindow,
) -> Dict[str, Any]:
    """
    Validate whether an expression level is therapeutic.

    Args:
        expression_level: Expression as fraction of normal (1.0 = normal).
        therapeutic_window: Target therapeutic window.

    Returns:
        Validation results.
    """
    tw = therapeutic_window

    # Classifications
    below_baseline = expression_level < tw.baseline_expression
    below_therapeutic = expression_level < tw.therapeutic_min
    in_therapeutic = tw.therapeutic_min <= expression_level <= tw.therapeutic_max
    above_optimal = expression_level > tw.optimal_expression
    above_max = expression_level > tw.therapeutic_max

    # Calculate distances
    distance_to_optimal = abs(expression_level - tw.optimal_expression)
    distance_to_min = max(0, tw.therapeutic_min - expression_level)
    distance_to_max = max(0, expression_level - tw.therapeutic_max)

    # Therapeutic score (1.0 at optimal, decreasing away)
    if in_therapeutic:
        therapeutic_score = 1.0 - distance_to_optimal / (tw.therapeutic_max - tw.therapeutic_min)
    else:
        therapeutic_score = 0.0

    # Classification
    if above_max:
        classification = "OVEREXPRESSION_RISK"
        status = "NO-GO"
    elif in_therapeutic:
        if distance_to_optimal < 0.05:
            classification = "OPTIMAL"
        else:
            classification = "THERAPEUTIC"
        status = "GO"
    elif below_therapeutic:
        if below_baseline:
            classification = "BELOW_BASELINE"
        else:
            classification = "SUBTHERAPEUTIC"
        status = "NO-GO"
    else:
        classification = "UNKNOWN"
        status = "NO-GO"

    return {
        'expression_level': expression_level,
        'classification': classification,
        'status': status,
        'in_therapeutic_window': in_therapeutic,
        'therapeutic_score': therapeutic_score,
        'distance_to_optimal': distance_to_optimal,
        'distance_to_min': distance_to_min,
        'distance_to_max': distance_to_max,
        'required_change': tw.therapeutic_min - expression_level if below_therapeutic else 0,
        'safety_margin': tw.therapeutic_max - expression_level,
    }


def compare_guides_therapeutic(
    guides: List[Dict[str, Any]],
    therapeutic_window: TherapeuticWindow,
    config: Optional[DosageConfig] = None,
) -> List[Dict[str, Any]]:
    """
    Compare multiple guides for therapeutic potential.

    Args:
        guides: List of guide dictionaries with 'coherence_R' and 'delta_g'.
        therapeutic_window: Target therapeutic window.
        config: Dosage configuration.

    Returns:
        List of guides with therapeutic assessment added.
    """
    if config is None:
        config = DosageConfig()

    results = []

    for guide in guides:
        coherence = guide.get('coherence_R', 0.5)
        delta_g = guide.get('delta_g', -15)

        optimization = optimize_dosage(
            therapeutic_window=therapeutic_window,
            guide_coherence=coherence,
            binding_energy=delta_g,
            config=config,
        )

        result = {**guide}
        result['therapeutic_assessment'] = optimization
        result['achieved_expression'] = optimization['achieved_expression']
        result['in_therapeutic_window'] = optimization['in_therapeutic_window']
        result['therapeutic_recommendation'] = optimization['recommendation']

        results.append(result)

    # Sort by therapeutic score (in window first, then by closeness to optimal)
    results.sort(
        key=lambda x: (
            x['in_therapeutic_window'],
            -abs(x['achieved_expression'] - therapeutic_window.optimal_expression)
        ),
        reverse=True,
    )

    return results
