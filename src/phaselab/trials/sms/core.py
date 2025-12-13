"""
Core data structures for SMS trials.

Defines trial types, results, and configuration shared across all
SMS trial runners.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, Any, List, Optional
from datetime import datetime


class TrialStatus(Enum):
    """Trial execution status."""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    BLOCKED = "blocked"  # Prerequisites not met


class TrialType(Enum):
    """Types of SMS trials."""
    CRISPRA_RAI1 = "crispra_rai1"
    CRISPRI_MODIFIER = "crispri_modifier"
    KNOCKOUT_VALIDATION = "knockout_validation"
    BASE_EDITING = "base_editing"
    PRIME_EDITING = "prime_editing"
    CIRCADIAN_SIMULATION = "circadian_simulation"
    DELIVERY_ASSESSMENT = "delivery_assessment"


@dataclass
class SMSTrialConfig:
    """
    Configuration for SMS trial execution.

    Attributes:
        therapeutic_window: Target expression range (70-110% of normal).
        target_tissues: Primary tissues for delivery.
        use_virtual_assay: Whether to use v0.7.0+ Virtual Assay Stack.
        coherence_mode: "heuristic" or "quantum" for IR coherence.
        require_go_status: Require GO status for candidates.
        min_claim_level: Minimum claim level for results.
        verbose: Print progress messages.
    """
    # Therapeutic parameters
    therapeutic_window: tuple = (0.70, 1.10)  # 70-110% of normal
    optimal_expression: float = 0.80  # 80% target for safety margin
    baseline_expression: float = 0.50  # SMS haploinsufficient baseline

    # Target tissues
    target_tissues: List[str] = field(default_factory=lambda: ["brain", "hypothalamus"])

    # Analysis options
    use_virtual_assay: bool = True
    coherence_mode: str = "heuristic"  # "heuristic" or "quantum"
    require_go_status: bool = True
    min_claim_level: str = "exploratory"  # "unknown", "exploratory", "context_dependent", "strong_computational"

    # Guide design
    top_n_guides: int = 10

    # Delivery
    delivery_route: str = "ICV"  # Intracerebroventricular for CNS
    preferred_serotype: str = "AAV9"  # FDA-approved, crosses BBB

    # Circadian simulation
    simulation_hours: float = 240.0  # 10 days
    n_circadian_trials: int = 5

    # Output
    verbose: bool = False
    include_raw_data: bool = False


@dataclass
class SMSTrialResult:
    """
    Result from an SMS trial.

    Attributes:
        trial_type: Type of trial executed.
        status: Execution status.
        summary: Human-readable summary.
        candidates: List of candidate guides/designs.
        best_candidate: Top recommended candidate.
        metrics: Key metrics dictionary.
        claim_level: v0.8.0 claim level for results.
        warnings: List of warning messages.
        errors: List of error messages.
        metadata: Additional trial metadata.
    """
    trial_type: TrialType
    status: TrialStatus
    summary: str

    # Results
    candidates: List[Dict[str, Any]] = field(default_factory=list)
    best_candidate: Optional[Dict[str, Any]] = None

    # Metrics
    metrics: Dict[str, Any] = field(default_factory=dict)

    # Claim level (v0.8.0)
    claim_level: str = "unknown"
    claim_description: str = ""

    # Diagnostics
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)

    # Metadata
    metadata: Dict[str, Any] = field(default_factory=dict)
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())

    @property
    def is_successful(self) -> bool:
        """Check if trial completed successfully."""
        return self.status == TrialStatus.COMPLETED

    @property
    def has_viable_candidates(self) -> bool:
        """Check if any viable candidates were found."""
        return len(self.candidates) > 0 and self.best_candidate is not None

    @property
    def n_candidates(self) -> int:
        """Number of candidates found."""
        return len(self.candidates)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "trial_type": self.trial_type.value,
            "status": self.status.value,
            "summary": self.summary,
            "n_candidates": self.n_candidates,
            "best_candidate": self.best_candidate,
            "metrics": self.metrics,
            "claim_level": self.claim_level,
            "claim_description": self.claim_description,
            "warnings": self.warnings,
            "errors": self.errors,
            "timestamp": self.timestamp,
        }


@dataclass
class SMSTrial:
    """
    Container for an SMS trial specification.

    Attributes:
        trial_type: Type of trial.
        name: Human-readable trial name.
        description: Trial description.
        config: Trial configuration.
        prerequisites: List of prerequisite trial types.
        sequence_data: DNA sequence data for the trial.
    """
    trial_type: TrialType
    name: str
    description: str
    config: SMSTrialConfig = field(default_factory=SMSTrialConfig)
    prerequisites: List[TrialType] = field(default_factory=list)
    sequence_data: Dict[str, str] = field(default_factory=dict)

    def check_prerequisites(self, completed_trials: List[TrialType]) -> bool:
        """Check if all prerequisites are met."""
        return all(p in completed_trials for p in self.prerequisites)


# RAI1 gene information
RAI1_INFO = {
    "gene_symbol": "RAI1",
    "gene_name": "Retinoic Acid Induced 1",
    "chromosome": "17",
    "band": "17p11.2",
    "strand": "+",
    "ncbi_gene_id": "10743",
    "ensembl_id": "ENSG00000108557",
    "uniprot_id": "Q7Z5J4",

    # Key regulatory regions (GRCh38 coordinates)
    "tss": 17584792,  # Transcription start site
    "promoter_start": 17583792,  # -1000bp from TSS
    "promoter_end": 17585292,    # +500bp from TSS
    "first_exon_start": 17584792,
    "first_exon_end": 17585100,

    # CRISPRa optimal window (relative to TSS)
    "crispra_window": (-400, +100),

    # Known enhancers
    "enhancers": [
        {"name": "ENH1", "start": 17581000, "end": 17582000, "activity": 0.8},
        {"name": "ENH2", "start": 17578500, "end": 17579500, "activity": 0.6},
    ],

    # Expression data
    "normal_expression_tpm": {
        "brain": 45.2,
        "hypothalamus": 62.1,
        "liver": 12.3,
        "muscle": 8.7,
    },
}

# SMS modifier genes that could be targeted with CRISPRi
SMS_MODIFIER_GENES = {
    "PER1": {
        "role": "Circadian clock component, period gene 1",
        "rationale": "Suppression may help adjust circadian phase",
        "caution": "Critical clock gene, careful dosing needed",
    },
    "PER2": {
        "role": "Circadian clock component, period gene 2",
        "rationale": "May compensate for RAI1-mediated PER dysregulation",
        "caution": "Essential for normal circadian function",
    },
    "CRY1": {
        "role": "Cryptochrome 1, negative feedback arm",
        "rationale": "Modulating CRY1 can adjust clock period",
        "caution": "Interacts with multiple pathways",
    },
    "CLOCK": {
        "role": "Core clock transcription factor",
        "rationale": "Reducing CLOCK may normalize BMAL1:CLOCK ratio",
        "caution": "Master regulator, avoid strong suppression",
    },
    "BMAL1": {
        "role": "Core clock transcription factor partner",
        "rationale": "Not recommended - already dysregulated in SMS",
        "caution": "Avoid suppression",
    },
}

# Example sequences (truncated for demonstration)
# In production, these would be fetched from genome databases
EXAMPLE_RAI1_PROMOTER = (
    "GCGCGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGC"
    "TGCGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCT"
    "GCGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTG"
    "CGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTGC"
    "GGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTGCG"
    "GCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTGCGGCGGCGGCTGCGG"
    "ATGCGATCGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGATCGATCG"
    "CGATCGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGATCGATCGCG"
)
