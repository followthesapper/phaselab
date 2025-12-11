"""
PhaseLab Chromatin: Expanded chromatin accessibility models for CRISPRa/i.

Provides tissue-specific chromatin accessibility scoring using:
- ENCODE DNase-seq data
- PsychENCODE brain-specific profiles
- ATAC-seq integration
- Cell-type specific models

Version: 0.2.0
"""

from .accessibility import (
    ChromatinModel,
    ChromatinScore,
    compute_accessibility,
    load_encode_model,
    load_psychencode_model,
)

from .tissue_models import (
    AVAILABLE_TISSUES,
    get_tissue_model,
    brain_accessibility,
    liver_accessibility,
    blood_accessibility,
)

from .integration import (
    combine_with_crispr,
    accessibility_weighted_coherence,
)

__all__ = [
    # Core
    "ChromatinModel",
    "ChromatinScore",
    "compute_accessibility",
    # Data sources
    "load_encode_model",
    "load_psychencode_model",
    # Tissue models
    "AVAILABLE_TISSUES",
    "get_tissue_model",
    "brain_accessibility",
    "liver_accessibility",
    "blood_accessibility",
    # Integration
    "combine_with_crispr",
    "accessibility_weighted_coherence",
]
