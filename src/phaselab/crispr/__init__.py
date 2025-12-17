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

NEW in v0.9.3 - CRISPRa Binding Register Model:
- NucleaseRole: Explicit BINDING vs CUTTING mode
- Relaxed PAM patterns for dCas9 binding (e.g., SaCas9 NNGRRN)
- Sliding binding register (±2bp) for GC-dense promoters
- Configurable guide length for literature reproduction
- Validated against Chang et al. 2022 sg2 winner

NEW in v0.9.2 - Guide Enumeration & Policy System:
- Region declaration with multi-TSS support
- PAM scanning for SpCas9, SaCas9, Cas12a
- Policy-based dominance ranking
- Benchmark mode for validation against published guides
- Reproducibility manifests

NEW in v0.7.0 - Enhanced Pipeline:
- Full Virtual Assay Stack integration
- Biological context from ENCODE (ATAC-seq, methylation, histones)
- ML efficiency predictions (DeepCRISPR, DeepSpCas9 adapters)
- Evidence fusion with uncertainty quantification
- Context-aware guide ranking
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
    # v0.9.1: U6/Pol III compatibility and repeat detection
    poly_t_penalty,
    is_repeat_region,
    u6_compatibility_check,
    # v0.9.1: CRISPOR composite scoring (legacy)
    OFFTARGET_MISMATCH_WEIGHTS,
    CRISPORMetrics,
    crispor_composite_score,
    rank_guides_crispor_style,
    validate_and_rerank_with_crispor,
    # v0.9.2: Policy-based dominance ranking
    RankingPolicy,
    PolicyConfig,
    POLICY_CONFIGS,
    GateResult,
    GuideTier,
    apply_hard_gates,
    rank_guides,
    emit_manifest,
    print_ranking_report,
)

# v0.6.0: Unified ATLAS-Q enhanced coherence
# v0.6.1: Added CoherenceMode for heuristic vs quantum selection
from .coherence_utils import (
    CoherenceMode,
    compute_guide_coherence,
    compute_guide_coherence_with_details,
    compute_coherence_batch,
    compute_coherence_with_zscore,
    is_guide_coherent,
    get_coherence_method,
    get_coherence_eligibility_info,
)

# v0.5.0: Multi-guide synergy
from .multiguide import (
    GuideCandidate,
    GuidePair,
    MultiGuideSet,
    MultiGuideConfig,
    design_multiguide_set,
    analyze_guide_pair,
    predict_pairwise_synergy,
    validate_guide_set,
    optimize_guide_spacing,
)

# v0.5.0: Enhancer targeting
from .enhancer import (
    Enhancer,
    EnhancerConfig,
    EnhancerGuideResult,
    design_enhancer_guides,
    identify_target_enhancers,
    score_enhancer_for_activation,
    predict_enhancer_activation_effect,
    compare_promoter_vs_enhancer,
    get_known_enhancers,
)

# v0.7.0: Enhanced Pipeline with Virtual Assay Stack
from .enhanced_pipeline import (
    Modality,
    EnhancedGuideConfig,
    EnhancedGuide,
    EnhancedDesignResult,
    design_enhanced_guides,
    compare_guides_with_without_context,
)

# v0.9.1: Local CRISPOR integration for genome-wide off-target validation
# Note: Requires separate CRISPOR installation due to large genome files (~6GB)
from .crispor_integration import (
    CrisporConfig,
    CrisporValidator,
    setup_crispor,
    validate_with_crispor,
)

# v0.9.2: Region declaration, enumeration, and full design pipeline
from .region import (
    GenomeBuild,
    TSSSource,
    TSSAnnotation,
    Window,
    Region,
    RegionSet,
    RegionBuilder,
    build_regions_for_gene,
    get_tss_for_gene,
)
from .enumerate import (
    Nuclease,
    NucleaseRole,  # v0.9.3: BINDING vs CUTTING mode
    NucleaseConfig,
    NUCLEASE_CONFIGS,
    CandidateGuide,
    PAMScanner,
    EnumerationResult,
    enumerate_guides,
    enumerate_from_sequence,
)
from .design import (
    DesignResult,
    BenchmarkResult,
    design_crispra_guides,  # v0.9.3: Primary CRISPRa API with binding mode
    design_crispra_guides as design_crispra_guides_v2,  # Alias for compatibility
    design_knockout_guides as design_knockout_guides_v2,
    evaluate_candidate,
    evaluate_candidates,
    benchmark_against_published,
)

# v0.9.4: Three Breakthrough Paths for CRISPRa Guide Ranking
# Path A: Binding Energy Landscape (quantum chemistry)
from .binding_landscape import (
    BindingEnergyResult,
    compute_binding_energy,
    compute_binding_landscape,
    rank_guides_by_binding,
)
# Path B: Transcriptional Phase Alignment (IR dynamics)
from .transcriptional_phase import (
    PhaseAlignmentResult,
    compute_phase_alignment,
    compute_phase_landscape,
    rank_guides_by_phase,
    optimal_guide_position,
)
# Path C: Off-Target Landscape Geometry (coherence contrast)
from .offtarget_geometry import (
    OffTargetGeometryResult,
    compute_offtarget_geometry,
    compute_geometry_landscape,
    rank_guides_by_geometry,
    integrate_crispor_offtargets,
)
# Unified multi-evidence scorer
from .multi_evidence import (
    MultiEvidenceResult,
    compute_multi_evidence_score,
    rank_guides_multi_evidence,
)

# v0.9.5: Quantum Discriminator for Late-Stage Guide Selection
from .quantum_discriminator import (
    DiscriminatorStatus,
    QuantumGuideResult,
    DiscriminatorResult,
    run_quantum_discriminator,
    design_guides_with_quantum_discriminator,
    DISCRIMINATOR_GATES,
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
    # v0.9.1: U6/Pol III compatibility and repeat detection
    "poly_t_penalty",
    "is_repeat_region",
    "u6_compatibility_check",
    # v0.9.1: CRISPOR composite scoring (legacy)
    "OFFTARGET_MISMATCH_WEIGHTS",
    "CRISPORMetrics",
    "crispor_composite_score",
    "rank_guides_crispor_style",
    "validate_and_rerank_with_crispor",
    # v0.9.2: Policy-based dominance ranking
    "RankingPolicy",
    "PolicyConfig",
    "POLICY_CONFIGS",
    "GateResult",
    "GuideTier",
    "apply_hard_gates",
    "rank_guides",
    "emit_manifest",
    "print_ranking_report",

    # ATLAS-Q Enhanced Coherence (v0.6.0+)
    # v0.6.1: CoherenceMode for heuristic vs quantum selection
    "CoherenceMode",
    "compute_guide_coherence",
    "compute_guide_coherence_with_details",
    "compute_coherence_batch",
    "compute_coherence_with_zscore",
    "is_guide_coherent",
    "get_coherence_method",
    "get_coherence_eligibility_info",

    # Multi-guide (v0.5.0)
    "GuideCandidate",
    "GuidePair",
    "MultiGuideSet",
    "MultiGuideConfig",
    "design_multiguide_set",
    "analyze_guide_pair",
    "predict_pairwise_synergy",
    "validate_guide_set",
    "optimize_guide_spacing",

    # Enhancer targeting (v0.5.0)
    "Enhancer",
    "EnhancerConfig",
    "EnhancerGuideResult",
    "design_enhancer_guides",
    "identify_target_enhancers",
    "score_enhancer_for_activation",
    "predict_enhancer_activation_effect",
    "compare_promoter_vs_enhancer",
    "get_known_enhancers",

    # Enhanced Pipeline with Virtual Assay Stack (v0.7.0)
    "Modality",
    "EnhancedGuideConfig",
    "EnhancedGuide",
    "EnhancedDesignResult",
    "design_enhanced_guides",
    "compare_guides_with_without_context",

    # CRISPOR Integration (v0.9.1)
    "CrisporConfig",
    "CrisporValidator",
    "setup_crispor",
    "validate_with_crispor",

    # v0.9.2: Region declaration
    "GenomeBuild",
    "TSSSource",
    "TSSAnnotation",
    "Window",
    "Region",
    "RegionSet",
    "RegionBuilder",
    "build_regions_for_gene",
    "get_tss_for_gene",

    # v0.9.2: Guide enumeration
    "Nuclease",
    "NucleaseConfig",
    "NUCLEASE_CONFIGS",
    "CandidateGuide",
    "PAMScanner",
    "EnumerationResult",
    "enumerate_guides",
    "enumerate_from_sequence",

    # v0.9.3: NucleaseRole for binding vs cutting mode
    "NucleaseRole",

    # v0.9.2/v0.9.3: Full design pipeline
    "DesignResult",
    "BenchmarkResult",
    "design_crispra_guides",  # v0.9.3: Primary CRISPRa API
    "design_crispra_guides_v2",  # Alias for compatibility
    "design_knockout_guides_v2",
    "evaluate_candidate",
    "evaluate_candidates",
    "benchmark_against_published",

    # v0.9.4: Three Breakthrough Paths for CRISPRa
    # Path A: Binding Energy Landscape (quantum chemistry)
    "BindingEnergyResult",
    "compute_binding_energy",
    "compute_binding_landscape",
    "rank_guides_by_binding",
    # Path B: Transcriptional Phase Alignment (IR dynamics)
    "PhaseAlignmentResult",
    "compute_phase_alignment",
    "compute_phase_landscape",
    "rank_guides_by_phase",
    "optimal_guide_position",
    # Path C: Off-Target Landscape Geometry (coherence contrast)
    "OffTargetGeometryResult",
    "compute_offtarget_geometry",
    "compute_geometry_landscape",
    "rank_guides_by_geometry",
    "integrate_crispor_offtargets",
    # Unified multi-evidence scorer
    "MultiEvidenceResult",
    "compute_multi_evidence_score",
    "rank_guides_multi_evidence",

    # v0.9.5: Quantum Discriminator for Late-Stage Guide Selection
    "DiscriminatorStatus",
    "QuantumGuideResult",
    "DiscriminatorResult",
    "run_quantum_discriminator",
    "design_guides_with_quantum_discriminator",
    "DISCRIMINATOR_GATES",
]
