"""
PhaseLab Protein: IR coherence metrics for protein folding and mutation analysis.

This module applies the Informational Relativity framework to assess:
- Protein structure prediction reliability
- Molecular dynamics simulation coherence
- Binding affinity prediction confidence
- Mutational scanning landscape analysis (E213 methodology)

The mutscan submodule applies the validated E213 spatial coherence methodology
to deep mutational scanning data, identifying functional domains with coherent
mutational effects.

Version: 0.3.0
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

from .mutscan import (
    MutScanLandscape,
    MutScanResult,
    FunctionalDomain,
    load_mutscan_data,
    analyze_mutscan_coherence,
    local_coherence_profile,
    map_coherence_to_structure,
    compare_mutscan_datasets,
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
    # Mutational Scanning (E213 methodology)
    "MutScanLandscape",
    "MutScanResult",
    "FunctionalDomain",
    "load_mutscan_data",
    "analyze_mutscan_coherence",
    "local_coherence_profile",
    "map_coherence_to_structure",
    "compare_mutscan_datasets",
]
