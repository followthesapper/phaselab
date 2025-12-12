"""
PhaseLab: A framework for phase-coherence analysis across quantum, biological, and dynamical systems.

Built on the Informational Relativity (IR) framework, PhaseLab provides:
- Quantum coherence metrics (R̄, V_φ) for simulation reliability
- CRISPR/CRISPRa guide RNA design pipelines
- Circadian clock modeling for gene therapy dosage optimization
- Protein folding and dynamics coherence assessment (v0.2.0)
- Tissue-specific chromatin accessibility models (v0.2.0)
- Multi-tissue circadian models with inter-tissue coupling (v0.3.0)
- Drug response modeling and chronotherapy optimization (v0.3.0)
- Expanded CRISPR editors: base editing (ABE/CBE) and prime editing (v0.3.0)
- Complete CRISPR toolkit with knockout and CRISPRi (v0.4.0)
- Therapeutic dosage optimization with haploinsufficiency models (v0.4.0)

NEW in v0.5.0:
- Real ATAC-seq BigWig integration for tissue-specific accessibility
- CpG methylation modeling for CRISPRa efficiency prediction
- Nucleosome occupancy prediction (NuPoP-like algorithm)
- Multi-guide synergy modeling for combinatorial CRISPR
- Enhancer targeting for CRISPRa
- AAV serotype selection and delivery modeling
- Immunogenicity prediction for Cas9 and guide RNA
- Comprehensive validation report generation

NEW in v0.6.0:
- ATLAS-Q integration for advanced quantum simulation
- IR measurement grouping (5× variance reduction)
- Real circular statistics coherence (replaces heuristic)
- Coherence-aware VQE optimization
- Optional GPU acceleration via Triton kernels
- Rust backend support for 30-77× faster simulation
- IR-enhanced off-target analysis for CRISPOR integration
- Unified ATLAS-Q coherence across all CRISPR modalities
- Off-target entropy and coherence contrast metrics

NEW in v0.6.1:
- Coherence mode parameter: mode="heuristic" (fast) vs mode="quantum" (VQE)
- Honest coherence weighting: heuristic demoted to tie-breaker (0.05 vs 0.30)
- Two-stage scoring: hard safety gates + soft ranking
- Risk mass metrics: risk_mass_close, risk_mass_exonic, tail_risk_score
- Evidence levels: A/B/C classification for validation status
- Score capping: unvalidated guides capped to prevent misleading rankings

Author: Dylan Vaca
License: MIT
"""

__version__ = "0.6.1"
__author__ = "Dylan Vaca"

from .core.coherence import coherence_score, go_no_go, phase_variance
from .core.constants import E_MINUS_2, FOUR_PI_SQUARED

# Quantum coherence (v0.6.0 - ATLAS-Q enhanced)
from .quantum.coherence import (
    CoherenceResult,
    compute_coherence_from_expectations,
    compute_coherence_from_phases,
    compute_coherence_from_hamiltonian,
)
from .quantum import is_atlas_q_available

# CRISPR coherence utilities (v0.6.0, v0.6.1 adds CoherenceMode)
from .crispr.coherence_utils import (
    CoherenceMode,
    compute_guide_coherence,
    compute_guide_coherence_with_details,
    is_guide_coherent,
    get_coherence_method,
    get_coherence_eligibility_info,
)

# CRISPOR Integration IR-Enhanced Analysis (v0.6.0)
from .integrations.crispor.offtarget_ir import (
    OffTargetSite,
    OffTargetIRAnalysis,
    analyze_offtarget_landscape,
    compute_ir_enhanced_score,
)

__all__ = [
    # Core coherence
    "coherence_score",
    "go_no_go",
    "phase_variance",
    "E_MINUS_2",
    "FOUR_PI_SQUARED",
    "__version__",

    # Quantum coherence (v0.6.0)
    "CoherenceResult",
    "compute_coherence_from_expectations",
    "compute_coherence_from_phases",
    "compute_coherence_from_hamiltonian",
    "is_atlas_q_available",

    # CRISPR coherence (v0.6.0, v0.6.1)
    "CoherenceMode",
    "compute_guide_coherence",
    "compute_guide_coherence_with_details",
    "is_guide_coherent",
    "get_coherence_method",
    "get_coherence_eligibility_info",

    # CRISPOR IR-Enhanced (v0.6.0)
    "OffTargetSite",
    "OffTargetIRAnalysis",
    "analyze_offtarget_landscape",
    "compute_ir_enhanced_score",
]
