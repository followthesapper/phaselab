"""
PhaseLab Chromatin: Expanded chromatin accessibility models for CRISPRa/i.

Provides tissue-specific chromatin accessibility scoring using:
- ENCODE DNase-seq data
- PsychENCODE brain-specific profiles
- ATAC-seq integration (BigWig/BED)
- CpG methylation modeling
- Nucleosome occupancy prediction
- Cell-type specific models

Version: 0.5.0
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

# v0.5.0: ATAC-seq integration
from .atac_integration import (
    ATACProfile,
    ATACPeak,
    load_encode_atac,
    create_tissue_profile,
    score_guide_accessibility,
    compare_tissue_accessibility,
)

# v0.5.0: CpG methylation
from .methylation import (
    MethylationScore,
    MethylationConfig,
    MethylationProfile,
    score_methylation,
    score_guide_methylation,
    adjust_crispra_efficiency,
    find_cpg_islands,
    CpGIsland,
)

# v0.5.0: Nucleosome positioning
from .nucleosome import (
    NucleosomeProfile,
    NucleosomeConfig,
    NucleosomeCall,
    predict_nucleosome_profile,
    predict_nucleosome_affinity,
    predict_nucleosome_occupancy,
    score_guide_nucleosome,
    find_nucleosome_free_regions,
)

__all__ = [
    # Core accessibility
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
    # ATAC-seq (v0.5.0)
    "ATACProfile",
    "ATACPeak",
    "load_encode_atac",
    "create_tissue_profile",
    "score_guide_accessibility",
    "compare_tissue_accessibility",
    # Methylation (v0.5.0)
    "MethylationScore",
    "MethylationConfig",
    "MethylationProfile",
    "score_methylation",
    "score_guide_methylation",
    "adjust_crispra_efficiency",
    "find_cpg_islands",
    "CpGIsland",
    # Nucleosome (v0.5.0)
    "NucleosomeProfile",
    "NucleosomeConfig",
    "NucleosomeCall",
    "predict_nucleosome_profile",
    "predict_nucleosome_affinity",
    "predict_nucleosome_occupancy",
    "score_guide_nucleosome",
    "find_nucleosome_free_regions",
]
