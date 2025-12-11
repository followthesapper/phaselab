"""
PhaseLab CRISPR: Guide RNA design pipeline with IR coherence validation.

Provides:
- PAM site scanning (NGG, NNGRRT, etc.)
- CRISPRa window filtering
- Off-target scoring (MIT, CFD algorithms)
- Thermodynamic binding energy (SantaLucia)
- Chromatin accessibility modeling
- IR coherence-based reliability scoring
"""

from .pipeline import design_guides, GuideDesignConfig
from .pam_scan import find_pam_sites, PAM_PATTERNS
from .scoring import (
    gc_content,
    delta_g_santalucia,
    mit_specificity_score,
    cfd_score,
    max_homopolymer_run,
)

__all__ = [
    "design_guides",
    "GuideDesignConfig",
    "find_pam_sites",
    "PAM_PATTERNS",
    "gc_content",
    "delta_g_santalucia",
    "mit_specificity_score",
    "cfd_score",
    "max_homopolymer_run",
]
