"""
PhaseLab Omics: Spatial coherence for genomic and epigenomic data.

This module applies the E213-validated spatial coherence methodology to:
- ATAC-seq (chromatin accessibility)
- ChIP-seq (histone marks, TF binding)
- RNA-seq (expression profiles)
- Multi-omics integration

Key insight:
- Genomic position = perturbation coordinate
- Signal intensity = response value
- Spatial coherence identifies reliable signal regions

Example use cases:
1. ATAC: Which accessible regions are reliably detected?
2. ChIP: Which peaks represent stable binding?
3. RNA-seq: Which expression changes are reproducible?

Version: 0.1.0
"""

from .atac import (
    ATACLandscape,
    ATACResult,
    load_atac_data,
    analyze_atac_coherence,
    identify_stable_peaks,
    compare_atac_conditions,
)

from .chip import (
    ChIPLandscape,
    ChIPResult,
    load_chip_data,
    analyze_chip_coherence,
    identify_stable_binding,
)

from .expression import (
    ExpressionLandscape,
    ExpressionResult,
    load_expression_data,
    analyze_expression_coherence,
    identify_reliable_changes,
)

__all__ = [
    # ATAC-seq
    "ATACLandscape",
    "ATACResult",
    "load_atac_data",
    "analyze_atac_coherence",
    "identify_stable_peaks",
    "compare_atac_conditions",
    # ChIP-seq
    "ChIPLandscape",
    "ChIPResult",
    "load_chip_data",
    "analyze_chip_coherence",
    "identify_stable_binding",
    # Expression
    "ExpressionLandscape",
    "ExpressionResult",
    "load_expression_data",
    "analyze_expression_coherence",
    "identify_reliable_changes",
]
