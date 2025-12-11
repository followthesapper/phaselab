"""
PhaseLab Drug: Drug response modeling with circadian considerations.

Provides:
- Pharmacokinetic phase dynamics
- Drug-clock interactions
- Chronotherapy optimization
- Dosing schedule optimization
- Individual variability modeling

Version: 0.3.0
"""

from .pharmacokinetics import (
    PKModel,
    PKParams,
    simulate_pk,
    compute_auc,
    compute_cmax,
)

from .chronotherapy import (
    ChronotherapyOptimizer,
    DoseSchedule,
    optimize_dosing_time,
    circadian_drug_interaction,
)

from .variability import (
    PopulationPK,
    simulate_population,
    individual_optimization,
)

__all__ = [
    # PK
    "PKModel",
    "PKParams",
    "simulate_pk",
    "compute_auc",
    "compute_cmax",
    # Chronotherapy
    "ChronotherapyOptimizer",
    "DoseSchedule",
    "optimize_dosing_time",
    "circadian_drug_interaction",
    # Variability
    "PopulationPK",
    "simulate_population",
    "individual_optimization",
]
