"""
PhaseLab: A framework for phase-coherence analysis across quantum, biological, and dynamical systems.

Built on the Informational Relativity (IR) framework, PhaseLab provides:
- Quantum coherence metrics (R̄, V_φ) for simulation reliability
- CRISPR/CRISPRa guide RNA design pipelines
- Circadian clock modeling for gene therapy dosage optimization
- Protein folding and dynamics coherence assessment (v0.2.0)
- Tissue-specific chromatin accessibility models (v0.2.0)

Author: Dylan Vaca
License: MIT
"""

__version__ = "0.2.0-dev"
__author__ = "Dylan Vaca"

from .core.coherence import coherence_score, go_no_go, phase_variance
from .core.constants import E_MINUS_2, FOUR_PI_SQUARED

__all__ = [
    "coherence_score",
    "go_no_go",
    "phase_variance",
    "E_MINUS_2",
    "FOUR_PI_SQUARED",
    "__version__",
]
