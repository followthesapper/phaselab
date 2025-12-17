"""
PhaseLab Chem: Spatial coherence for chemical and biochemical systems.

This module applies the E213-validated spatial coherence methodology to:
- Binding landscapes (protein-ligand, protein-protein)
- Reaction optimization (enzymatic, synthetic)
- Chemical space exploration

Key insight:
- Chemical perturbations (mutations, modifications, conditions) create response landscapes
- Binding affinity, reaction rate, selectivity = response signals
- Spatial coherence identifies reliable regions for optimization

Example use cases:
1. Binding: Which binding site mutations reliably affect affinity?
2. Reactions: Which condition ranges give reproducible yields?
3. Optimization: Where in chemical space is behavior predictable?

Version: 0.1.0
"""

from .binding import (
    BindingLandscape,
    BindingResult,
    load_binding_data,
    analyze_binding_coherence,
    identify_stable_binding_modes,
    hot_spot_analysis,
)

from .reaction import (
    ReactionLandscape,
    ReactionResult,
    load_reaction_data,
    analyze_reaction_coherence,
    identify_stable_conditions,
    optimize_with_stability,
)

from .screening import (
    ScreeningLandscape,
    ScreeningResult,
    load_screening_data,
    analyze_screening_coherence,
    identify_reliable_hits,
)

__all__ = [
    # Binding
    "BindingLandscape",
    "BindingResult",
    "load_binding_data",
    "analyze_binding_coherence",
    "identify_stable_binding_modes",
    "hot_spot_analysis",
    # Reaction
    "ReactionLandscape",
    "ReactionResult",
    "load_reaction_data",
    "analyze_reaction_coherence",
    "identify_stable_conditions",
    "optimize_with_stability",
    # Screening
    "ScreeningLandscape",
    "ScreeningResult",
    "load_screening_data",
    "analyze_screening_coherence",
    "identify_reliable_hits",
]
