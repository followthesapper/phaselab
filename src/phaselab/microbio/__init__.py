"""
PhaseLab Microbio: Spatial coherence for microbial perturbation experiments.

This module applies the E213-validated spatial coherence methodology to:
- TnSeq (transposon sequencing) fitness landscapes
- RB-TnSeq (random barcoded TnSeq)
- CRISPRi essential gene screens
- Drug response profiling

Key insight:
- Microbial genomes have perturbation landscapes too
- Transposon insertion sites = perturbation positions
- Fitness scores = response signal
- Spatial coherence identifies operationally stable gene regions

Example use cases:
1. TnSeq: Which insertions in a gene reliably affect fitness?
2. CRISPRi: Which guide positions give reproducible knockdown?
3. Drug response: Which concentrations show stable pharmacology?

Version: 0.1.0
"""

from .tnseq import (
    TnSeqLandscape,
    TnSeqResult,
    load_tnseq_data,
    analyze_tnseq_coherence,
    identify_essential_domains,
    compare_conditions,
)

from .crispri import (
    CRISPRiLandscape,
    CRISPRiResult,
    load_crispri_screen,
    analyze_crispri_coherence,
    rank_crispri_guides,
)

from .drug import (
    DrugResponseLandscape,
    DrugResponseResult,
    load_dose_response,
    analyze_drug_coherence,
    identify_stable_dosing_window,
    pharmacodynamic_stability,
)

from .fitness import (
    FitnessLandscape,
    compute_fitness_coherence,
    identify_functional_domains,
    essential_gene_analysis,
)

__all__ = [
    # TnSeq
    "TnSeqLandscape",
    "TnSeqResult",
    "load_tnseq_data",
    "analyze_tnseq_coherence",
    "identify_essential_domains",
    "compare_conditions",
    # CRISPRi
    "CRISPRiLandscape",
    "CRISPRiResult",
    "load_crispri_screen",
    "analyze_crispri_coherence",
    "rank_crispri_guides",
    # Drug response
    "DrugResponseLandscape",
    "DrugResponseResult",
    "load_dose_response",
    "analyze_drug_coherence",
    "identify_stable_dosing_window",
    "pharmacodynamic_stability",
    # Fitness
    "FitnessLandscape",
    "compute_fitness_coherence",
    "identify_functional_domains",
    "essential_gene_analysis",
]
