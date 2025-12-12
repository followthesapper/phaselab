"""
Evidence fusion engine.

Combines multiple evidence sources into a unified guide score
with proper uncertainty propagation.
"""

from dataclasses import dataclass, field
from typing import Optional, Dict, List, Any, Tuple
import logging

import numpy as np

from phaselab.core.constants import E_MINUS_2
from phaselab.fusion.evidence import (
    Evidence,
    EvidenceSource,
    EvidenceType,
    EvidenceCollection,
)

logger = logging.getLogger(__name__)


@dataclass
class FusionConfig:
    """
    Configuration for evidence fusion.

    Attributes
    ----------
    sequence_weight : float
        Weight for sequence-based scores
    context_weight : float
        Weight for biological context
    ml_weight : float
        Weight for ML predictions
    specificity_weight : float
        Weight for off-target specificity
    coherence_weight_heuristic : float
        Weight for heuristic coherence
    coherence_weight_quantum : float
        Weight for quantum coherence
    hard_gate_sources : List[EvidenceSource]
        Sources that act as hard gates
    go_threshold : float
        Coherence threshold for GO status
    """
    sequence_weight: float = 0.25
    context_weight: float = 0.15
    ml_weight: float = 0.30
    specificity_weight: float = 0.25
    coherence_weight_heuristic: float = 0.05
    coherence_weight_quantum: float = 0.30
    hard_gate_sources: List[EvidenceSource] = field(
        default_factory=lambda: [
            EvidenceSource.OFF_TARGET_COUNT,
            EvidenceSource.GC_CONTENT,
        ]
    )
    go_threshold: float = E_MINUS_2

    @classmethod
    def default(cls) -> "FusionConfig":
        """Default configuration."""
        return cls()

    @classmethod
    def conservative(cls) -> "FusionConfig":
        """Conservative configuration with higher specificity weight."""
        return cls(
            sequence_weight=0.20,
            context_weight=0.10,
            ml_weight=0.25,
            specificity_weight=0.40,
            coherence_weight_heuristic=0.05,
            coherence_weight_quantum=0.25,
        )

    @classmethod
    def ml_focused(cls) -> "FusionConfig":
        """ML-focused configuration."""
        return cls(
            sequence_weight=0.15,
            context_weight=0.10,
            ml_weight=0.45,
            specificity_weight=0.25,
            coherence_weight_heuristic=0.05,
            coherence_weight_quantum=0.30,
        )


@dataclass
class FusedResult:
    """
    Result of evidence fusion.

    Attributes
    ----------
    score : float
        Final combined score [0, 1]
    confidence : float
        Overall confidence in the score
    confidence_interval : Tuple[float, float]
        95% CI for the score
    go_status : str
        GO or NO-GO based on coherence
    go_reason : str
        Explanation for GO/NO-GO
    passes_gates : bool
        Whether all hard gates passed
    failed_gates : List[str]
        Names of failed gates
    component_scores : Dict[str, float]
        Individual component scores
    coherence : Optional[float]
        Coherence value if available
    evidence_count : int
        Number of evidence sources used
    warnings : List[str]
        Fusion warnings
    """
    score: float
    confidence: float
    confidence_interval: Tuple[float, float]
    go_status: str
    go_reason: str
    passes_gates: bool
    failed_gates: List[str]
    component_scores: Dict[str, float]
    coherence: Optional[float]
    evidence_count: int
    warnings: List[str] = field(default_factory=list)

    @property
    def is_go(self) -> bool:
        """Whether status is GO."""
        return self.go_status == "GO"

    @property
    def is_viable(self) -> bool:
        """Whether guide is viable (passes gates and GO)."""
        return self.passes_gates and self.is_go

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "score": self.score,
            "confidence": self.confidence,
            "confidence_interval": self.confidence_interval,
            "go_status": self.go_status,
            "go_reason": self.go_reason,
            "passes_gates": self.passes_gates,
            "failed_gates": self.failed_gates,
            "component_scores": self.component_scores,
            "coherence": self.coherence,
            "evidence_count": self.evidence_count,
            "warnings": self.warnings,
        }


class EvidenceFusion:
    """
    Evidence fusion engine for CRISPR guide scoring.

    Combines evidence from multiple sources into a unified score
    using a hierarchical approach:

    1. Hard gates (must pass all)
    2. Weighted soft scores
    3. Coherence overlay

    Parameters
    ----------
    config : FusionConfig
        Fusion configuration

    Examples
    --------
    >>> fusion = EvidenceFusion()
    >>>
    >>> # Add sequence-based evidence
    >>> fusion.add_evidence(Evidence.soft_score(
    ...     EvidenceSource.GC_CONTENT, 0.75, confidence=0.95
    ... ))
    >>>
    >>> # Add ML prediction
    >>> fusion.add_evidence(Evidence.soft_score(
    ...     EvidenceSource.DEEPCRISPR, 0.82, confidence=0.85
    ... ))
    >>>
    >>> # Add coherence
    >>> fusion.add_evidence(Evidence.coherence(0.45, mode="heuristic"))
    >>>
    >>> # Fuse all evidence
    >>> result = fusion.fuse()
    >>> print(f"Score: {result.score:.2f}, Status: {result.go_status}")
    """

    # Evidence source categories
    SEQUENCE_SOURCES = {
        EvidenceSource.GC_CONTENT,
        EvidenceSource.THERMODYNAMICS,
        EvidenceSource.POSITION,
        EvidenceSource.HOMOPOLYMER,
    }

    CONTEXT_SOURCES = {
        EvidenceSource.ACCESSIBILITY,
        EvidenceSource.METHYLATION,
        EvidenceSource.HISTONE,
        EvidenceSource.CHROMATIN_STATE,
    }

    ML_SOURCES = {
        EvidenceSource.DEEPCRISPR,
        EvidenceSource.DEEPSPCAS9,
        EvidenceSource.ENFORMER,
        EvidenceSource.RULE_BASED,
    }

    SPECIFICITY_SOURCES = {
        EvidenceSource.OFF_TARGET_COUNT,
        EvidenceSource.OFF_TARGET_CFD,
        EvidenceSource.RISK_MASS,
    }

    COHERENCE_SOURCES = {
        EvidenceSource.COHERENCE_HEURISTIC,
        EvidenceSource.COHERENCE_QUANTUM,
        EvidenceSource.COHERENCE_HARDWARE,
    }

    def __init__(self, config: Optional[FusionConfig] = None):
        self.config = config or FusionConfig.default()
        self._evidence = EvidenceCollection(guide_id="")
        self._warnings: List[str] = []

    def reset(self, guide_id: str = "") -> None:
        """Reset fusion state for new guide."""
        self._evidence = EvidenceCollection(guide_id=guide_id)
        self._warnings = []

    def add_evidence(self, evidence: Evidence) -> None:
        """Add evidence to the fusion."""
        self._evidence.add(evidence)

    def add_sequence_score(
        self,
        score: float,
        source: EvidenceSource = EvidenceSource.GC_CONTENT,
        confidence: float = 1.0,
        **metadata,
    ) -> None:
        """Add sequence-based score."""
        self.add_evidence(Evidence.soft_score(
            source=source,
            score=score,
            confidence=confidence,
            weight=self.config.sequence_weight,
            **metadata,
        ))

    def add_context_score(
        self,
        score: float,
        source: EvidenceSource = EvidenceSource.ACCESSIBILITY,
        confidence: float = 1.0,
        **metadata,
    ) -> None:
        """Add biological context score."""
        self.add_evidence(Evidence.soft_score(
            source=source,
            score=score,
            confidence=confidence,
            weight=self.config.context_weight,
            **metadata,
        ))

    def add_ml_prediction(
        self,
        score: float,
        source: EvidenceSource = EvidenceSource.DEEPCRISPR,
        confidence: float = 1.0,
        **metadata,
    ) -> None:
        """Add ML prediction score."""
        self.add_evidence(Evidence.soft_score(
            source=source,
            score=score,
            confidence=confidence,
            weight=self.config.ml_weight,
            **metadata,
        ))

    def add_specificity_score(
        self,
        score: float,
        source: EvidenceSource = EvidenceSource.OFF_TARGET_CFD,
        confidence: float = 1.0,
        is_gate: bool = False,
        **metadata,
    ) -> None:
        """Add specificity/off-target score."""
        if is_gate:
            # Convert to binary gate
            passes = score > 0.5
            self.add_evidence(Evidence.hard_gate(
                source=source,
                passes=passes,
                confidence=confidence,
                **metadata,
            ))
        else:
            self.add_evidence(Evidence.soft_score(
                source=source,
                score=score,
                confidence=confidence,
                weight=self.config.specificity_weight,
                **metadata,
            ))

    def add_coherence(
        self,
        r_bar: float,
        mode: str = "heuristic",
        confidence: float = 1.0,
        **metadata,
    ) -> None:
        """Add coherence evidence."""
        weight = (
            self.config.coherence_weight_quantum
            if mode in ("quantum", "hardware")
            else self.config.coherence_weight_heuristic
        )
        self.add_evidence(Evidence.coherence(
            r_bar=r_bar,
            mode=mode,
            confidence=confidence,
            weight=weight,
            **metadata,
        ))

    def add_hard_gate(
        self,
        source: EvidenceSource,
        passes: bool,
        confidence: float = 1.0,
        **metadata,
    ) -> None:
        """Add hard gate evidence."""
        self.add_evidence(Evidence.hard_gate(
            source=source,
            passes=passes,
            confidence=confidence,
            **metadata,
        ))

    def fuse(self) -> FusedResult:
        """
        Fuse all evidence into final result.

        Returns
        -------
        FusedResult
            Combined score with metadata
        """
        warnings = list(self._warnings)

        # Stage 1: Check hard gates
        passes_gates = self._evidence.passes_all_gates
        failed_gates = [e.source.value for e in self._evidence.failed_gates]

        if not passes_gates:
            warnings.append(f"Failed gates: {failed_gates}")

        # Stage 2: Calculate component scores
        component_scores = self._calculate_component_scores()

        # Stage 3: Weight and combine
        combined_score = self._combine_scores(component_scores)

        # Stage 4: Apply coherence overlay
        coherence_evidence = self._evidence.coherence
        coherence_value = coherence_evidence.value if coherence_evidence else None

        go_status, go_reason = self._determine_go_status(coherence_value)

        # Apply NO-GO penalty
        if go_status == "NO-GO":
            combined_score *= 0.5
            warnings.append("Score penalized due to NO-GO status")

        # Calculate confidence
        confidence = self._calculate_confidence()

        # Calculate confidence interval
        ci = self._calculate_confidence_interval(combined_score, confidence)

        return FusedResult(
            score=float(combined_score),
            confidence=float(confidence),
            confidence_interval=ci,
            go_status=go_status,
            go_reason=go_reason,
            passes_gates=passes_gates,
            failed_gates=failed_gates,
            component_scores=component_scores,
            coherence=coherence_value,
            evidence_count=len(self._evidence.evidences),
            warnings=warnings,
        )

    def _calculate_component_scores(self) -> Dict[str, float]:
        """Calculate scores for each evidence category."""
        components = {}

        # Sequence scores
        seq_scores = [
            e for e in self._evidence.scores
            if e.source in self.SEQUENCE_SOURCES
        ]
        if seq_scores:
            components["sequence"] = self._weighted_mean(seq_scores)

        # Context scores
        ctx_scores = [
            e for e in self._evidence.scores
            if e.source in self.CONTEXT_SOURCES
        ]
        if ctx_scores:
            components["context"] = self._weighted_mean(ctx_scores)

        # ML scores
        ml_scores = [
            e for e in self._evidence.scores
            if e.source in self.ML_SOURCES
        ]
        if ml_scores:
            components["ml"] = self._weighted_mean(ml_scores)

        # Specificity scores
        spec_scores = [
            e for e in self._evidence.scores
            if e.source in self.SPECIFICITY_SOURCES
        ]
        if spec_scores:
            components["specificity"] = self._weighted_mean(spec_scores)

        return components

    def _weighted_mean(self, evidences: List[Evidence]) -> float:
        """Calculate confidence-weighted mean of evidences."""
        if not evidences:
            return 0.5

        values = np.array([e.value for e in evidences])
        weights = np.array([e.effective_weight for e in evidences])

        if weights.sum() == 0:
            return float(np.mean(values))

        return float(np.average(values, weights=weights))

    def _combine_scores(self, components: Dict[str, float]) -> float:
        """Combine component scores into final score."""
        if not components:
            return 0.5

        # Get weights for present components
        weight_map = {
            "sequence": self.config.sequence_weight,
            "context": self.config.context_weight,
            "ml": self.config.ml_weight,
            "specificity": self.config.specificity_weight,
        }

        total_weight = 0.0
        weighted_sum = 0.0

        for name, score in components.items():
            weight = weight_map.get(name, 0.0)
            weighted_sum += score * weight
            total_weight += weight

        if total_weight == 0:
            return 0.5

        return weighted_sum / total_weight

    def _determine_go_status(
        self,
        coherence: Optional[float],
    ) -> Tuple[str, str]:
        """Determine GO/NO-GO status from coherence."""
        if coherence is None:
            return "GO", "No coherence data (assuming GO)"

        threshold = self.config.go_threshold

        if coherence >= threshold:
            margin = coherence - threshold
            if coherence >= 0.5:
                return "GO", f"Strong coherence (R={coherence:.3f}, +{margin:.3f} above threshold)"
            else:
                return "GO", f"Adequate coherence (R={coherence:.3f}, +{margin:.3f} above threshold)"
        else:
            deficit = threshold - coherence
            if coherence >= 0.05:
                return "NO-GO", f"Below threshold (R={coherence:.3f}, -{deficit:.3f} below e^-2)"
            else:
                return "NO-GO", f"Critical coherence deficit (R={coherence:.3f})"

    def _calculate_confidence(self) -> float:
        """Calculate overall confidence in fusion."""
        if not self._evidence.scores:
            return 0.0

        confidences = [e.confidence for e in self._evidence.scores]

        # Weight by source importance
        weights = [e.weight for e in self._evidence.scores]

        if sum(weights) == 0:
            return float(np.mean(confidences))

        return float(np.average(confidences, weights=weights))

    def _calculate_confidence_interval(
        self,
        score: float,
        confidence: float,
    ) -> Tuple[float, float]:
        """Calculate 95% confidence interval."""
        # Width inversely proportional to confidence
        # and number of evidence sources
        n_sources = len(self._evidence.scores)
        base_width = 0.3 * (1 - confidence)

        # Shrink with more sources
        if n_sources > 1:
            shrink = 1.0 / np.sqrt(n_sources)
            width = base_width * shrink
        else:
            width = base_width

        low = max(0.0, score - width)
        high = min(1.0, score + width)

        return (low, high)

    @classmethod
    def fuse_guide(
        cls,
        sequence_scores: Dict[str, float],
        context_scores: Optional[Dict[str, float]] = None,
        ml_scores: Optional[Dict[str, float]] = None,
        specificity_scores: Optional[Dict[str, float]] = None,
        coherence: Optional[Tuple[float, str]] = None,
        config: Optional[FusionConfig] = None,
    ) -> FusedResult:
        """
        Convenience method to fuse guide evidence in one call.

        Parameters
        ----------
        sequence_scores : Dict[str, float]
            Sequence-based scores {source_name: score}
        context_scores : Dict[str, float], optional
            Biological context scores
        ml_scores : Dict[str, float], optional
            ML prediction scores
        specificity_scores : Dict[str, float], optional
            Off-target specificity scores
        coherence : Tuple[float, str], optional
            (r_bar, mode) tuple
        config : FusionConfig, optional
            Fusion configuration

        Returns
        -------
        FusedResult
            Fused score
        """
        fusion = cls(config=config)

        # Add sequence scores
        source_map = {
            "gc": EvidenceSource.GC_CONTENT,
            "gc_content": EvidenceSource.GC_CONTENT,
            "thermo": EvidenceSource.THERMODYNAMICS,
            "thermodynamics": EvidenceSource.THERMODYNAMICS,
            "position": EvidenceSource.POSITION,
            "homopolymer": EvidenceSource.HOMOPOLYMER,
        }
        for name, score in sequence_scores.items():
            source = source_map.get(name.lower(), EvidenceSource.GC_CONTENT)
            fusion.add_sequence_score(score, source=source)

        # Add context scores
        if context_scores:
            ctx_map = {
                "accessibility": EvidenceSource.ACCESSIBILITY,
                "methylation": EvidenceSource.METHYLATION,
                "histone": EvidenceSource.HISTONE,
            }
            for name, score in context_scores.items():
                source = ctx_map.get(name.lower(), EvidenceSource.ACCESSIBILITY)
                fusion.add_context_score(score, source=source)

        # Add ML scores
        if ml_scores:
            ml_map = {
                "deepcrispr": EvidenceSource.DEEPCRISPR,
                "deepspcas9": EvidenceSource.DEEPSPCAS9,
                "enformer": EvidenceSource.ENFORMER,
                "rule_based": EvidenceSource.RULE_BASED,
            }
            for name, score in ml_scores.items():
                source = ml_map.get(name.lower(), EvidenceSource.RULE_BASED)
                fusion.add_ml_prediction(score, source=source)

        # Add specificity scores
        if specificity_scores:
            spec_map = {
                "off_target_count": EvidenceSource.OFF_TARGET_COUNT,
                "off_target_cfd": EvidenceSource.OFF_TARGET_CFD,
                "risk_mass": EvidenceSource.RISK_MASS,
            }
            for name, score in specificity_scores.items():
                source = spec_map.get(name.lower(), EvidenceSource.OFF_TARGET_CFD)
                fusion.add_specificity_score(score, source=source)

        # Add coherence
        if coherence:
            r_bar, mode = coherence
            fusion.add_coherence(r_bar, mode=mode)

        return fusion.fuse()
