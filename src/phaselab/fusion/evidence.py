"""
Evidence data structures for the fusion system.

Defines standardized evidence objects that can be combined
from multiple sources.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional, Dict, Any, List, Tuple
import logging

import numpy as np

logger = logging.getLogger(__name__)


class EvidenceType(Enum):
    """Type of evidence."""
    HARD_GATE = "hard_gate"      # Binary pass/fail
    SOFT_SCORE = "soft_score"    # Continuous [0, 1]
    COHERENCE = "coherence"      # IR coherence metric


class ClaimLevel(Enum):
    """
    Claim level classification for guide predictions.

    Indicates the strength of computational support for the prediction.
    This helps users understand what the evidence actually supports.

    Levels
    ------
    STRONG_COMPUTATIONAL : str
        Multiple layers agree (ML, context, sequence)
        High confidence (>0.8)
        Validated against wet-lab data
    CONTEXT_DEPENDENT : str
        Prediction depends on biological context
        May vary by cell type, tissue, or conditions
        Moderate confidence (0.5-0.8)
    EXPLORATORY : str
        Limited evidence, single layer agreement
        Useful for hypothesis generation
        Lower confidence (0.3-0.5)
    UNKNOWN : str
        Insufficient evidence to make a call
        Confidence below threshold (<0.3)
        Explicit "we don't know" bucket
    """
    STRONG_COMPUTATIONAL = "strong_computational"
    CONTEXT_DEPENDENT = "context_dependent"
    EXPLORATORY = "exploratory"
    UNKNOWN = "unknown"

    @classmethod
    def from_confidence(
        cls,
        confidence: float,
        n_agreeing_layers: int = 1,
        has_wet_lab_validation: bool = False,
    ) -> "ClaimLevel":
        """
        Determine claim level from confidence and layer agreement.

        Parameters
        ----------
        confidence : float
            Overall confidence in prediction [0, 1]
        n_agreeing_layers : int
            Number of evidence layers that agree
        has_wet_lab_validation : bool
            Whether prediction is validated against wet-lab data

        Returns
        -------
        ClaimLevel
            Appropriate claim level
        """
        # Unknown: insufficient evidence
        if confidence < 0.3 or n_agreeing_layers == 0:
            return cls.UNKNOWN

        # Exploratory: limited evidence
        if confidence < 0.5 or n_agreeing_layers == 1:
            return cls.EXPLORATORY

        # Context-dependent: moderate evidence
        if confidence < 0.8 or n_agreeing_layers < 3:
            return cls.CONTEXT_DEPENDENT

        # Strong computational: high evidence + agreement
        return cls.STRONG_COMPUTATIONAL

    @property
    def description(self) -> str:
        """Human-readable description of claim level."""
        descriptions = {
            ClaimLevel.STRONG_COMPUTATIONAL: (
                "Strong computational support: Multiple evidence layers agree "
                "with high confidence. Suitable for prioritization."
            ),
            ClaimLevel.CONTEXT_DEPENDENT: (
                "Context-dependent prediction: May vary by cell type or conditions. "
                "Validate in target context before proceeding."
            ),
            ClaimLevel.EXPLORATORY: (
                "Exploratory only: Limited evidence. Use for hypothesis generation, "
                "not definitive ranking."
            ),
            ClaimLevel.UNKNOWN: (
                "Unknown: Insufficient evidence to make a reliable prediction. "
                "Additional data needed."
            ),
        }
        return descriptions.get(self, "Unknown claim level")

    @property
    def color_code(self) -> str:
        """Color code for visualization (green/yellow/orange/gray)."""
        colors = {
            ClaimLevel.STRONG_COMPUTATIONAL: "green",
            ClaimLevel.CONTEXT_DEPENDENT: "yellow",
            ClaimLevel.EXPLORATORY: "orange",
            ClaimLevel.UNKNOWN: "gray",
        }
        return colors.get(self, "gray")


class EvidenceSource(Enum):
    """Source of evidence."""
    # Sequence-based
    GC_CONTENT = "gc_content"
    THERMODYNAMICS = "thermodynamics"
    POSITION = "position"
    HOMOPOLYMER = "homopolymer"

    # Context-based
    ACCESSIBILITY = "accessibility"
    METHYLATION = "methylation"
    HISTONE = "histone"
    CHROMATIN_STATE = "chromatin_state"

    # ML predictions
    DEEPCRISPR = "deepcrispr"
    DEEPSPCAS9 = "deepspcas9"
    ENFORMER = "enformer"
    RULE_BASED = "rule_based"

    # Specificity
    OFF_TARGET_COUNT = "off_target_count"
    OFF_TARGET_CFD = "off_target_cfd"
    RISK_MASS = "risk_mass"

    # Coherence
    COHERENCE_HEURISTIC = "coherence_heuristic"
    COHERENCE_QUANTUM = "coherence_quantum"
    COHERENCE_HARDWARE = "coherence_hardware"

    # Custom
    CUSTOM = "custom"


@dataclass
class Evidence:
    """
    Single piece of evidence for guide scoring.

    Attributes
    ----------
    source : EvidenceSource
        Where this evidence comes from
    evidence_type : EvidenceType
        Type of evidence (gate, score, coherence)
    value : float
        The evidence value
        - For HARD_GATE: 1.0 = pass, 0.0 = fail
        - For SOFT_SCORE: [0, 1] continuous
        - For COHERENCE: [0, 1] R-bar value
    confidence : float
        Confidence in this evidence [0, 1]
    weight : float
        Relative weight in fusion (before normalization)
    metadata : Dict[str, Any]
        Additional information
    timestamp : float, optional
        When evidence was generated (for caching)
    """
    source: EvidenceSource
    evidence_type: EvidenceType
    value: float
    confidence: float = 1.0
    weight: float = 1.0
    metadata: Dict[str, Any] = field(default_factory=dict)
    timestamp: Optional[float] = None

    def __post_init__(self):
        # Validate value based on type
        if self.evidence_type == EvidenceType.HARD_GATE:
            self.value = 1.0 if self.value > 0.5 else 0.0
        else:
            self.value = float(np.clip(self.value, 0.0, 1.0))

        self.confidence = float(np.clip(self.confidence, 0.0, 1.0))
        self.weight = max(0.0, float(self.weight))

    @property
    def is_gate(self) -> bool:
        """Whether this is a hard gate."""
        return self.evidence_type == EvidenceType.HARD_GATE

    @property
    def passes_gate(self) -> bool:
        """Whether hard gate passes (True if not a gate)."""
        if not self.is_gate:
            return True
        return self.value > 0.5

    @property
    def effective_weight(self) -> float:
        """Weight adjusted by confidence."""
        return self.weight * self.confidence

    @property
    def is_coherence(self) -> bool:
        """Whether this is coherence evidence."""
        return self.evidence_type == EvidenceType.COHERENCE

    @classmethod
    def hard_gate(
        cls,
        source: EvidenceSource,
        passes: bool,
        confidence: float = 1.0,
        **metadata,
    ) -> "Evidence":
        """Create a hard gate evidence."""
        return cls(
            source=source,
            evidence_type=EvidenceType.HARD_GATE,
            value=1.0 if passes else 0.0,
            confidence=confidence,
            weight=1.0,  # Gates always have weight 1
            metadata=metadata,
        )

    @classmethod
    def soft_score(
        cls,
        source: EvidenceSource,
        score: float,
        confidence: float = 1.0,
        weight: float = 1.0,
        **metadata,
    ) -> "Evidence":
        """Create a soft score evidence."""
        return cls(
            source=source,
            evidence_type=EvidenceType.SOFT_SCORE,
            value=score,
            confidence=confidence,
            weight=weight,
            metadata=metadata,
        )

    @classmethod
    def coherence(
        cls,
        r_bar: float,
        mode: str = "heuristic",
        confidence: float = 1.0,
        weight: float = 1.0,
        **metadata,
    ) -> "Evidence":
        """Create coherence evidence."""
        source_map = {
            "heuristic": EvidenceSource.COHERENCE_HEURISTIC,
            "quantum": EvidenceSource.COHERENCE_QUANTUM,
            "hardware": EvidenceSource.COHERENCE_HARDWARE,
        }
        return cls(
            source=source_map.get(mode, EvidenceSource.COHERENCE_HEURISTIC),
            evidence_type=EvidenceType.COHERENCE,
            value=r_bar,
            confidence=confidence,
            weight=weight,
            metadata={"mode": mode, **metadata},
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "source": self.source.value,
            "type": self.evidence_type.value,
            "value": self.value,
            "confidence": self.confidence,
            "weight": self.weight,
            "effective_weight": self.effective_weight,
            "metadata": self.metadata,
        }


@dataclass
class EvidenceCollection:
    """
    Collection of evidence for a single guide.

    Organizes evidence by type and provides aggregation helpers.
    """
    guide_id: str
    evidences: List[Evidence] = field(default_factory=list)

    def add(self, evidence: Evidence) -> None:
        """Add evidence to collection."""
        self.evidences.append(evidence)

    @property
    def gates(self) -> List[Evidence]:
        """All hard gate evidence."""
        return [e for e in self.evidences if e.is_gate]

    @property
    def scores(self) -> List[Evidence]:
        """All soft score evidence."""
        return [
            e for e in self.evidences
            if e.evidence_type == EvidenceType.SOFT_SCORE
        ]

    @property
    def coherence(self) -> Optional[Evidence]:
        """Coherence evidence (highest confidence if multiple)."""
        coh = [e for e in self.evidences if e.is_coherence]
        if not coh:
            return None
        return max(coh, key=lambda e: e.confidence)

    @property
    def passes_all_gates(self) -> bool:
        """Whether all hard gates pass."""
        return all(e.passes_gate for e in self.gates)

    @property
    def failed_gates(self) -> List[Evidence]:
        """Gates that failed."""
        return [e for e in self.gates if not e.passes_gate]

    def by_source(self, source: EvidenceSource) -> Optional[Evidence]:
        """Get evidence by source."""
        for e in self.evidences:
            if e.source == source:
                return e
        return None

    def sources_present(self) -> List[EvidenceSource]:
        """List of evidence sources present."""
        return list(set(e.source for e in self.evidences))

    def total_effective_weight(self) -> float:
        """Sum of effective weights for soft scores."""
        return sum(e.effective_weight for e in self.scores)

    def summary(self) -> Dict[str, Any]:
        """Get summary of evidence collection."""
        return {
            "guide_id": self.guide_id,
            "n_gates": len(self.gates),
            "n_scores": len(self.scores),
            "has_coherence": self.coherence is not None,
            "passes_all_gates": self.passes_all_gates,
            "failed_gates": [e.source.value for e in self.failed_gates],
            "sources": [s.value for s in self.sources_present()],
        }
