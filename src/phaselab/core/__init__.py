"""
PhaseLab Core: IR coherence metrics and quantum utilities.
"""

from .coherence import coherence_score, go_no_go, phase_variance
from .constants import E_MINUS_2, FOUR_PI_SQUARED
from .hamiltonians import build_pauli_hamiltonian

__all__ = [
    "coherence_score",
    "go_no_go",
    "phase_variance",
    "E_MINUS_2",
    "FOUR_PI_SQUARED",
    "build_pauli_hamiltonian",
]
