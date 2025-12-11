"""
PhaseLab Protein: IR coherence metrics for protein folding simulations.

This module applies the Informational Relativity framework to assess:
- Protein structure prediction reliability
- Molecular dynamics simulation coherence
- Binding affinity prediction confidence

Version: 0.2.0
"""

from .folding import (
    FoldingCoherence,
    compute_structure_coherence,
    ramachandran_coherence,
    contact_map_coherence,
    go_no_go_structure,
)

from .dynamics import (
    MDCoherence,
    trajectory_coherence,
    rmsd_stability,
    ensemble_convergence,
)

from .binding import (
    BindingCoherence,
    docking_coherence,
    affinity_confidence,
)

__all__ = [
    # Folding
    "FoldingCoherence",
    "compute_structure_coherence",
    "ramachandran_coherence",
    "contact_map_coherence",
    "go_no_go_structure",
    # Dynamics
    "MDCoherence",
    "trajectory_coherence",
    "rmsd_stability",
    "ensemble_convergence",
    # Binding
    "BindingCoherence",
    "docking_coherence",
    "affinity_confidence",
]
