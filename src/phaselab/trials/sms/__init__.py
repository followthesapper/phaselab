"""
PhaseLab SMS Trials: Smith-Magenis Syndrome preclinical decision engine.

This module provides comprehensive trial runners for SMS gene therapy
development, integrating:

1. CRISPRa RAI1 activation (core therapeutic strategy)
2. CRISPRi modifier gene suppression
3. Knockout guides for model validation
4. Base editing for hypomorphic variant rescue
5. Prime editing for regulatory motif correction
6. Circadian rescue simulation
7. AAV delivery feasibility assessment

All trials use the PhaseLab Virtual Assay Stack (v0.7.0+) with
claim levels (v0.8.0) for honest uncertainty quantification.

SMS Background
--------------
Smith-Magenis Syndrome is caused by RAI1 haploinsufficiency,
typically from a 17p11.2 deletion or RAI1 point mutations.

Key features:
- Circadian rhythm inversion (melatonin peaks during day)
- Sleep disturbances
- Behavioral challenges
- Intellectual disability

Therapeutic goal: Boost RAI1 expression from ~50% to 70-110%
(therapeutic window) without overexpression toxicity.

Version: 0.9.0
"""

from phaselab.trials.sms.core import (
    SMSTrial,
    SMSTrialResult,
    SMSTrialConfig,
    TrialStatus,
)
from phaselab.trials.sms.crispra_trial import run_sms_crispra_trial
from phaselab.trials.sms.crispri_trial import run_sms_crispri_trial
from phaselab.trials.sms.knockout_trial import run_sms_knockout_trial
from phaselab.trials.sms.base_editing_trial import run_sms_base_editing_trial
from phaselab.trials.sms.prime_editing_trial import run_sms_prime_editing_trial
from phaselab.trials.sms.circadian_simulation import run_circadian_rescue_simulation
from phaselab.trials.sms.delivery_assessment import run_delivery_assessment
from phaselab.trials.sms.pipeline import SMSPipeline

__all__ = [
    # Core
    "SMSTrial",
    "SMSTrialResult",
    "SMSTrialConfig",
    "TrialStatus",
    # Trial runners
    "run_sms_crispra_trial",
    "run_sms_crispri_trial",
    "run_sms_knockout_trial",
    "run_sms_base_editing_trial",
    "run_sms_prime_editing_trial",
    "run_circadian_rescue_simulation",
    "run_delivery_assessment",
    # Pipeline
    "SMSPipeline",
]
