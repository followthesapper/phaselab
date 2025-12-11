"""
PhaseLab Circadian: Clock gene network models with IR coherence.

Provides:
- Kuramoto-based circadian oscillator models
- SMS-specific RAI1 dosage models
- PER gene delay dynamics
- REV-ERBα / RORα modulation
- Therapeutic window analysis
- Multi-tissue circadian models (v0.3.0)
- Jet lag and shift work simulations (v0.3.0)
"""

from .sms_model import (
    simulate_sms_clock,
    SMSClockParams,
    therapeutic_scan,
    classify_synchronization,
)
from .kuramoto import (
    kuramoto_order_parameter,
    kuramoto_ode,
)
from .multi_tissue import (
    MultiTissueParams,
    MultiTissueResult,
    simulate_multi_tissue,
    jet_lag_simulation,
    shift_work_simulation,
    TISSUE_PARAMS,
)

__all__ = [
    # SMS model
    "simulate_sms_clock",
    "SMSClockParams",
    "therapeutic_scan",
    "classify_synchronization",
    # Kuramoto
    "kuramoto_order_parameter",
    "kuramoto_ode",
    # Multi-tissue (v0.3.0)
    "MultiTissueParams",
    "MultiTissueResult",
    "simulate_multi_tissue",
    "jet_lag_simulation",
    "shift_work_simulation",
    "TISSUE_PARAMS",
]
