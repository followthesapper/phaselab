"""
PhaseLab Efficacy: Therapeutic efficacy prediction and comparison.

Implements:
- Multi-factor efficacy prediction
- Guide comparison for therapeutic applications
- Risk-benefit analysis
- IR coherence integration
"""

import numpy as np
from typing import Dict, Any, Optional, List
from dataclasses import dataclass

from ..core.coherence import go_no_go
from ..core.constants import E_MINUS_2


@dataclass
class EfficacyPredictor:
    """
    Predicts therapeutic efficacy from guide and target parameters.

    Integrates:
    - Guide binding quality (ΔG, MIT, CFD)
    - IR coherence (simulation reliability)
    - Chromatin accessibility
    - Expression target match
    """

    # Weights for different factors
    weight_binding: float = 0.25
    weight_specificity: float = 0.20
    weight_coherence: float = 0.20
    weight_chromatin: float = 0.15
    weight_expression_match: float = 0.20

    # Thresholds
    min_coherence: float = E_MINUS_2
    min_mit_score: float = 50
    min_chromatin_accessibility: float = 0.3

    def predict(
        self,
        delta_g: float,
        mit_score: float,
        cfd_score: float,
        coherence_R: float,
        chromatin_accessibility: float,
        expression_match: float = 1.0,
    ) -> Dict[str, Any]:
        """
        Predict therapeutic efficacy.

        Args:
            delta_g: Binding energy (kcal/mol).
            mit_score: MIT specificity score (0-100).
            cfd_score: CFD cutting score (0-100).
            coherence_R: IR coherence value.
            chromatin_accessibility: Chromatin accessibility (0-1).
            expression_match: How well achieved expression matches target (0-1).

        Returns:
            Efficacy prediction with component scores.
        """
        # Normalize binding energy to 0-1 (more negative = better)
        binding_score = min(1.0, max(0, (-delta_g - 5) / 20))

        # Normalize specificity (average of MIT and CFD)
        specificity_score = (mit_score / 100 + cfd_score / 100) / 2

        # Coherence score (already 0-1)
        coherence_score = coherence_R if coherence_R else 0.5

        # Chromatin score (already 0-1)
        chromatin_score = chromatin_accessibility

        # Expression match (already 0-1)
        expression_score = expression_match

        # Weighted combination
        efficacy = (
            self.weight_binding * binding_score +
            self.weight_specificity * specificity_score +
            self.weight_coherence * coherence_score +
            self.weight_chromatin * chromatin_score +
            self.weight_expression_match * expression_score
        )

        # Flags for potential issues
        warnings = []
        if coherence_R and coherence_R < self.min_coherence:
            warnings.append(f"Low coherence ({coherence_R:.3f}) indicates unreliable simulation")
        if mit_score < self.min_mit_score:
            warnings.append(f"Low specificity (MIT={mit_score:.0f}) may cause off-target effects")
        if chromatin_accessibility < self.min_chromatin_accessibility:
            warnings.append(f"Low chromatin accessibility ({chromatin_accessibility:.2f}) may reduce efficacy")

        # Classification
        if efficacy >= 0.8:
            classification = "EXCELLENT"
        elif efficacy >= 0.6:
            classification = "GOOD"
        elif efficacy >= 0.4:
            classification = "MODERATE"
        else:
            classification = "LOW"

        return {
            'efficacy_score': efficacy,
            'classification': classification,
            'component_scores': {
                'binding': binding_score,
                'specificity': specificity_score,
                'coherence': coherence_score,
                'chromatin': chromatin_score,
                'expression_match': expression_score,
            },
            'warnings': warnings,
            'go_no_go': go_no_go(coherence_R) if coherence_R else "UNKNOWN",
        }


def predict_therapeutic_efficacy(
    guide_data: Dict[str, Any],
    expression_target: float = 0.85,
    achieved_expression: float = None,
) -> Dict[str, Any]:
    """
    Predict therapeutic efficacy from guide data.

    Args:
        guide_data: Dictionary with guide metrics.
        expression_target: Target expression level.
        achieved_expression: Predicted achieved expression.

    Returns:
        Efficacy prediction.
    """
    predictor = EfficacyPredictor()

    # Extract values with defaults
    delta_g = guide_data.get('delta_g', -15)
    mit_score = guide_data.get('mit_score', 70)
    cfd_score = guide_data.get('cfd_score', 70)
    coherence_R = guide_data.get('coherence_R', 0.5)
    chromatin = guide_data.get('chromatin_accessibility', 0.5)

    # Calculate expression match
    if achieved_expression is not None and expression_target > 0:
        expression_match = 1.0 - min(1.0, abs(achieved_expression - expression_target) / expression_target)
    else:
        expression_match = 0.7  # Default moderate match

    return predictor.predict(
        delta_g=delta_g,
        mit_score=mit_score,
        cfd_score=cfd_score,
        coherence_R=coherence_R,
        chromatin_accessibility=chromatin,
        expression_match=expression_match,
    )


def compare_interventions(
    interventions: List[Dict[str, Any]],
    therapeutic_window: Optional[Any] = None,
) -> Dict[str, Any]:
    """
    Compare multiple intervention strategies.

    Args:
        interventions: List of intervention dictionaries with:
            - name: Intervention name
            - type: "CRISPRa", "CRISPRi", "knockout", etc.
            - guide_data: Guide metrics
            - achieved_expression: Predicted expression level
        therapeutic_window: Target therapeutic window (optional).

    Returns:
        Comparison results with rankings.
    """
    results = []

    for intervention in interventions:
        name = intervention.get('name', 'Unknown')
        int_type = intervention.get('type', 'Unknown')
        guide_data = intervention.get('guide_data', {})
        achieved = intervention.get('achieved_expression')

        # Get efficacy prediction
        target_expression = 0.85
        if therapeutic_window:
            target_expression = therapeutic_window.optimal_expression

        efficacy = predict_therapeutic_efficacy(
            guide_data=guide_data,
            expression_target=target_expression,
            achieved_expression=achieved,
        )

        # Check therapeutic window if provided
        in_window = None
        if therapeutic_window and achieved is not None:
            in_window = (
                therapeutic_window.therapeutic_min <=
                achieved <=
                therapeutic_window.therapeutic_max
            )

        results.append({
            'name': name,
            'type': int_type,
            'efficacy_score': efficacy['efficacy_score'],
            'classification': efficacy['classification'],
            'achieved_expression': achieved,
            'in_therapeutic_window': in_window,
            'warnings': efficacy['warnings'],
            'go_no_go': efficacy['go_no_go'],
        })

    # Sort by efficacy
    results.sort(key=lambda x: x['efficacy_score'], reverse=True)

    # Add rankings
    for i, r in enumerate(results):
        r['rank'] = i + 1

    # Summary
    best = results[0] if results else None
    in_window_count = sum(1 for r in results if r['in_therapeutic_window'])

    return {
        'interventions': results,
        'best_intervention': best['name'] if best else None,
        'best_efficacy': best['efficacy_score'] if best else None,
        'n_in_therapeutic_window': in_window_count,
        'recommendation': _generate_recommendation(results, therapeutic_window),
    }


def _generate_recommendation(
    results: List[Dict[str, Any]],
    therapeutic_window: Optional[Any],
) -> str:
    """Generate a recommendation based on comparison results."""
    if not results:
        return "No interventions to compare."

    best = results[0]

    if best['go_no_go'] == "NO-GO":
        return (
            f"CAUTION: Best option ({best['name']}) has low coherence. "
            "Consider alternative guides or validate experimentally."
        )

    if therapeutic_window and not best['in_therapeutic_window']:
        return (
            f"WARNING: Best option ({best['name']}) may not achieve therapeutic levels. "
            "Consider dose adjustment or alternative strategy."
        )

    if best['efficacy_score'] >= 0.7:
        return (
            f"RECOMMENDED: {best['name']} shows high predicted efficacy "
            f"({best['efficacy_score']:.2f}). Proceed with experimental validation."
        )
    elif best['efficacy_score'] >= 0.5:
        return (
            f"ACCEPTABLE: {best['name']} shows moderate efficacy "
            f"({best['efficacy_score']:.2f}). May require optimization."
        )
    else:
        return (
            f"SUBOPTIMAL: Best option ({best['name']}) has low predicted efficacy "
            f"({best['efficacy_score']:.2f}). Recommend screening additional guides."
        )


def risk_benefit_analysis(
    guide_data: Dict[str, Any],
    therapeutic_window: Any,
    achieved_expression: float,
) -> Dict[str, Any]:
    """
    Perform risk-benefit analysis for a therapeutic intervention.

    Args:
        guide_data: Guide metrics dictionary.
        therapeutic_window: Target therapeutic window.
        achieved_expression: Predicted achieved expression.

    Returns:
        Risk-benefit analysis.
    """
    # Benefits
    efficacy = predict_therapeutic_efficacy(
        guide_data=guide_data,
        expression_target=therapeutic_window.optimal_expression,
        achieved_expression=achieved_expression,
    )

    in_window = (
        therapeutic_window.therapeutic_min <=
        achieved_expression <=
        therapeutic_window.therapeutic_max
    )

    expression_deficit_resolved = (
        achieved_expression >= therapeutic_window.therapeutic_min
        and therapeutic_window.baseline_expression < therapeutic_window.therapeutic_min
    )

    # Risks
    mit_score = guide_data.get('mit_score', 70)
    coherence = guide_data.get('coherence_R', 0.5)

    off_target_risk = "LOW" if mit_score >= 70 else ("MODERATE" if mit_score >= 50 else "HIGH")
    overexpression_risk = "HIGH" if achieved_expression > therapeutic_window.therapeutic_max else "LOW"
    reliability_risk = "LOW" if coherence > 0.5 else ("MODERATE" if coherence > E_MINUS_2 else "HIGH")

    # Overall assessment
    benefit_score = (
        0.4 * efficacy['efficacy_score'] +
        0.3 * (1.0 if in_window else 0.0) +
        0.3 * (1.0 if expression_deficit_resolved else 0.0)
    )

    risk_score = (
        0.4 * (1.0 if off_target_risk == "HIGH" else (0.5 if off_target_risk == "MODERATE" else 0.0)) +
        0.3 * (1.0 if overexpression_risk == "HIGH" else 0.0) +
        0.3 * (1.0 if reliability_risk == "HIGH" else (0.5 if reliability_risk == "MODERATE" else 0.0))
    )

    ratio = benefit_score / (risk_score + 0.01)  # Avoid division by zero

    if ratio >= 3.0:
        recommendation = "STRONGLY_FAVORABLE"
    elif ratio >= 1.5:
        recommendation = "FAVORABLE"
    elif ratio >= 0.8:
        recommendation = "NEUTRAL"
    else:
        recommendation = "UNFAVORABLE"

    return {
        'benefit_score': benefit_score,
        'risk_score': risk_score,
        'benefit_risk_ratio': ratio,
        'recommendation': recommendation,
        'benefits': {
            'efficacy_score': efficacy['efficacy_score'],
            'in_therapeutic_window': in_window,
            'expression_deficit_resolved': expression_deficit_resolved,
        },
        'risks': {
            'off_target_risk': off_target_risk,
            'overexpression_risk': overexpression_risk,
            'reliability_risk': reliability_risk,
        },
        'efficacy_details': efficacy,
    }
