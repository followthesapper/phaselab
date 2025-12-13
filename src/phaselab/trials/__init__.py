"""
PhaseLab Trials: Preclinical decision engine for gene therapy development.

Provides comprehensive trial runners for disease-specific therapeutic
strategies including CRISPRa, CRISPRi, knockout, base editing, and
prime editing modalities.

NEW in v0.9.0:
- SMS (Smith-Magenis Syndrome) trial pipeline
- CRISPRa RAI1 activation trial with therapeutic window validation
- CRISPRi modifier gene suppression trials
- Circadian rescue simulation and prediction
- Delivery feasibility assessment for CNS-targeted therapies
"""

from phaselab.trials.sms import (
    SMSTrial,
    SMSTrialResult,
    SMSTrialConfig,
    SMSPipeline,
    run_sms_crispra_trial,
    run_sms_crispri_trial,
    run_sms_knockout_trial,
    run_sms_base_editing_trial,
    run_sms_prime_editing_trial,
    run_circadian_rescue_simulation,
    run_delivery_assessment,
)

__all__ = [
    # SMS Trials
    "SMSTrial",
    "SMSTrialResult",
    "SMSTrialConfig",
    "SMSPipeline",
    "run_sms_crispra_trial",
    "run_sms_crispri_trial",
    "run_sms_knockout_trial",
    "run_sms_base_editing_trial",
    "run_sms_prime_editing_trial",
    "run_circadian_rescue_simulation",
    "run_delivery_assessment",
]
