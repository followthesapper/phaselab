"""
PhaseLab Therapy: Therapeutic window and dosage optimization.

Provides tools for:
- Expression level estimation from CRISPRa/CRISPRi
- Therapeutic window identification
- Dosage-response modeling
- Toxicity threshold estimation
- IR coherence-validated efficacy prediction
"""

from .dosage import (
    TherapeuticWindow,
    DosageConfig,
    estimate_expression_change,
    calculate_therapeutic_window,
    dosage_response_curve,
    optimize_dosage,
    validate_therapeutic_level,
)

from .efficacy import (
    EfficacyPredictor,
    predict_therapeutic_efficacy,
    compare_interventions,
)

__all__ = [
    # Dosage module
    "TherapeuticWindow",
    "DosageConfig",
    "estimate_expression_change",
    "calculate_therapeutic_window",
    "dosage_response_curve",
    "optimize_dosage",
    "validate_therapeutic_level",
    # Efficacy module
    "EfficacyPredictor",
    "predict_therapeutic_efficacy",
    "compare_interventions",
]
