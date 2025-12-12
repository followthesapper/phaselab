"""
PhaseLab CRISPOR Integration.

Provides integration with CRISPOR for comprehensive guide RNA design:
- Off-target scoring (MIT, CFD)
- On-target activity prediction (Doench 2016)
- Genome-wide off-target enumeration
- Combined IR coherence + specificity scoring
- IR-enhanced off-target analysis (v0.6.0+):
  - Off-target entropy (risk distribution)
  - Coherence contrast (ΔR̄)
  - Energy spectrum analysis
  - Exonic risk flagging

Usage:
    from phaselab.integrations.crispor import CrisporClient, design_guides_with_crispor

    # Initialize client
    client = CrisporClient(crispor_path="/path/to/crispor")

    # Design guides with combined scoring
    results = design_guides_with_crispor(
        sequence="ATCG...",
        crispor_client=client,
    )

    # For IR-enhanced off-target analysis
    from phaselab.integrations.crispor import analyze_offtarget_landscape
    ir_analysis = analyze_offtarget_landscape(guide_seq, r_bar, offtargets)
"""

from .client import CrisporClient, CrisporConfig, CrisporOutput
from .parsers import (
    CrisporGuideRow,
    CrisporOffTarget,
    parse_guides_tsv,
    parse_offtargets_tsv,
    index_guides_by_sequence,
)
from .pipeline import (
    EvidenceLevel,
    CoherenceSource,
    GuideCandidate,
    ScoringWeights,
    check_safety_gate,
    design_guides_with_crispor,
    compute_final_score,
    merge_crispor_data,
    rank_guides,
    generate_report,
    save_results_json,
)
from .offtarget_ir import (
    OffTargetSite,
    OffTargetIRAnalysis,
    analyze_offtarget_landscape,
    compute_offtarget_entropy,
    compute_coherence_contrast,
    compute_energy_spectrum,
    compute_ir_enhanced_score,
    compute_offtarget_clustering,
    compute_region_difficulty,
    compute_risk_mass,
    compute_tail_risk,
    compute_concentration_measures,
)

__all__ = [
    # Client
    "CrisporClient",
    "CrisporConfig",
    "CrisporOutput",
    # Parsers
    "CrisporGuideRow",
    "CrisporOffTarget",
    "parse_guides_tsv",
    "parse_offtargets_tsv",
    "index_guides_by_sequence",
    # Pipeline - Evidence & Enums (v0.6.1)
    "EvidenceLevel",
    "CoherenceSource",
    # Pipeline - Core
    "GuideCandidate",
    "ScoringWeights",
    "check_safety_gate",
    "design_guides_with_crispor",
    "compute_final_score",
    "merge_crispor_data",
    "rank_guides",
    "generate_report",
    "save_results_json",
    # IR-Enhanced Off-Target Analysis (v0.6.0+)
    "OffTargetSite",
    "OffTargetIRAnalysis",
    "analyze_offtarget_landscape",
    "compute_offtarget_entropy",
    "compute_coherence_contrast",
    "compute_energy_spectrum",
    "compute_ir_enhanced_score",
    "compute_offtarget_clustering",
    "compute_region_difficulty",
    # Risk Mass & Concentration (v0.6.1)
    "compute_risk_mass",
    "compute_tail_risk",
    "compute_concentration_measures",
]
