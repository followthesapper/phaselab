"""
Biological Context Layer for PhaseLab Virtual Assay Stack.

This module provides integration with epigenomic databases (ENCODE, Roadmap)
for context-aware CRISPR guide scoring.

Components:
- ChromatinContext: ATAC-seq, DNase-seq accessibility
- MethylationContext: CpG methylation from WGBS/RRBS
- HistoneContext: Histone modification marks (H3K4me3, H3K27ac, etc.)
- ContextStack: Unified interface combining all context sources

Example
-------
>>> from phaselab.context import ContextStack, CellType
>>>
>>> # Create context for K562 cells
>>> ctx = ContextStack(cell_type=CellType.K562)
>>>
>>> # Get accessibility at a genomic position
>>> accessibility = ctx.get_accessibility("chr4", 55095264)
>>>
>>> # Get full context for guide scoring
>>> context = ctx.get_context("chr4", 55095264, 55095284)
"""

from phaselab.context.chromatin import (
    ChromatinContext,
    AccessibilityScore,
)
from phaselab.context.methylation import (
    MethylationContext,
    MethylationScore,
)
from phaselab.context.histone import (
    HistoneContext,
    HistoneMarks,
)
from phaselab.context.stack import (
    ContextStack,
    BiologicalContext,
    CellType,
)

__all__ = [
    # Chromatin
    "ChromatinContext",
    "AccessibilityScore",
    # Methylation
    "MethylationContext",
    "MethylationScore",
    # Histone
    "HistoneContext",
    "HistoneMarks",
    # Unified
    "ContextStack",
    "BiologicalContext",
    "CellType",
]
