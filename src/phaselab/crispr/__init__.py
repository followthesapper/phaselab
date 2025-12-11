"""
PhaseLab CRISPR: Comprehensive guide RNA design pipeline with IR coherence validation.

Provides:
- CRISPRa guide design for transcriptional activation
- CRISPRi guide design for transcriptional interference/repression
- CRISPR knockout guide design for gene disruption
- Prime editing pegRNA design for precise edits
- Base editing guide design (ABE/CBE) for single-nucleotide changes
- PAM site scanning (NGG, NNGRRT, etc.)
- Off-target scoring (MIT, CFD algorithms)
- Thermodynamic binding energy (SantaLucia)
- Chromatin accessibility modeling
- IR coherence-based reliability scoring
"""

# CRISPRa (activation)
from .pipeline import design_guides, GuideDesignConfig, validate_guide

# CRISPRi (interference/repression)
from .interference import (
    design_crispri_guides,
    CRISPRiConfig,
    validate_crispri_guide,
    repression_efficiency_score,
    steric_hindrance_score,
)

# CRISPR knockout
from .knockout import (
    design_knockout_guides,
    KnockoutConfig,
    validate_knockout_guide,
    cut_efficiency_score,
    frameshift_probability,
    repair_pathway_prediction,
)

# Prime editing
from .prime_editing import (
    design_prime_edit,
    PrimeEditConfig,
    validate_prime_edit,
    design_pbs,
    design_rt_template,
    pbs_score,
    rt_template_score,
    reverse_complement,
    estimate_hairpin_dg,
)

# Base editing
from .base_editing import (
    design_base_edit_guides,
    BaseEditConfig,
    validate_base_edit,
    design_abe_guides,
    design_cbe_guides,
    editing_efficiency_at_position,
    find_bystanders,
    get_activity_window,
)

# PAM scanning
from .pam_scan import find_pam_sites, PAM_PATTERNS

# Scoring
from .scoring import (
    gc_content,
    delta_g_santalucia,
    mit_specificity_score,
    cfd_score,
    max_homopolymer_run,
    chromatin_accessibility_score,
)

__all__ = [
    # CRISPRa (activation)
    "design_guides",
    "GuideDesignConfig",
    "validate_guide",

    # CRISPRi (interference)
    "design_crispri_guides",
    "CRISPRiConfig",
    "validate_crispri_guide",
    "repression_efficiency_score",
    "steric_hindrance_score",

    # Knockout
    "design_knockout_guides",
    "KnockoutConfig",
    "validate_knockout_guide",
    "cut_efficiency_score",
    "frameshift_probability",
    "repair_pathway_prediction",

    # Prime editing
    "design_prime_edit",
    "PrimeEditConfig",
    "validate_prime_edit",
    "design_pbs",
    "design_rt_template",
    "pbs_score",
    "rt_template_score",
    "reverse_complement",
    "estimate_hairpin_dg",

    # Base editing
    "design_base_edit_guides",
    "BaseEditConfig",
    "validate_base_edit",
    "design_abe_guides",
    "design_cbe_guides",
    "editing_efficiency_at_position",
    "find_bystanders",
    "get_activity_window",

    # PAM scanning
    "find_pam_sites",
    "PAM_PATTERNS",

    # Scoring
    "gc_content",
    "delta_g_santalucia",
    "mit_specificity_score",
    "cfd_score",
    "max_homopolymer_run",
    "chromatin_accessibility_score",
]
