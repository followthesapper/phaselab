"""
Predictor stack for combining multiple ML predictors.

Provides ensemble predictions and manages predictor registration,
fallbacks, and confidence-weighted aggregation.
"""

from dataclasses import dataclass, field
from typing import Optional, Dict, List, Any, Callable
import logging

import numpy as np

from phaselab.predictors.protocol import (
    OutcomePredictor,
    PredictionResult,
    PredictorType,
    EvidenceLevel,
)

logger = logging.getLogger(__name__)


@dataclass
class EnsemblePrediction:
    """
    Combined prediction from multiple predictors.

    Attributes
    ----------
    score : float
        Aggregated score [0, 1]
    confidence : float
        Overall confidence
    confidence_interval : tuple
        95% CI for aggregated score
    predictions : Dict[str, PredictionResult]
        Individual predictions by predictor name
    agreement : float
        Agreement between predictors [0, 1]
    dominant_type : PredictorType
        Most represented prediction type
    warnings : List[str]
        Aggregated warnings
    """
    score: float
    confidence: float
    confidence_interval: tuple
    predictions: Dict[str, PredictionResult]
    agreement: float
    dominant_type: PredictorType
    warnings: List[str] = field(default_factory=list)

    @property
    def n_predictors(self) -> int:
        """Number of contributing predictors."""
        return len(self.predictions)

    @property
    def is_reliable(self) -> bool:
        """Whether ensemble prediction is reliable."""
        return self.confidence >= 0.5 and self.n_predictors >= 2

    def get_by_type(self, pred_type: PredictorType) -> List[PredictionResult]:
        """Get predictions of a specific type."""
        return [
            p for p in self.predictions.values()
            if p.predictor_type == pred_type
        ]

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "score": self.score,
            "confidence": self.confidence,
            "confidence_interval": self.confidence_interval,
            "agreement": self.agreement,
            "n_predictors": self.n_predictors,
            "dominant_type": self.dominant_type.value,
            "predictions": {
                name: {
                    "score": p.score,
                    "confidence": p.confidence,
                    "type": p.predictor_type.value,
                    "evidence": p.evidence_level.value,
                }
                for name, p in self.predictions.items()
            },
            "warnings": self.warnings,
        }


class PredictorStack:
    """
    Stack of ML predictors for ensemble predictions.

    Manages multiple predictors, handles unavailable dependencies gracefully,
    and combines predictions using confidence-weighted aggregation.

    Parameters
    ----------
    aggregation : str
        How to combine predictions: "weighted", "mean", "max", "min"
    min_confidence : float
        Minimum confidence to include prediction in ensemble
    fallback_predictor : OutcomePredictor, optional
        Predictor to use when others fail

    Examples
    --------
    >>> from phaselab.predictors import PredictorStack, DeepCRISPRAdapter, RuleBasedPredictor
    >>>
    >>> # Create stack with predictors
    >>> stack = PredictorStack(aggregation="weighted")
    >>> stack.register(DeepCRISPRAdapter())
    >>> stack.register(RuleBasedPredictor())  # Fallback
    >>>
    >>> # Get ensemble prediction
    >>> result = stack.predict("ATGCGATCGATCGATCGATCNGG")
    >>> print(f"Score: {result.score:.2f} (n={result.n_predictors})")
    """

    def __init__(
        self,
        aggregation: str = "weighted",
        min_confidence: float = 0.0,
        fallback_predictor: Optional[OutcomePredictor] = None,
    ):
        self.aggregation = aggregation
        self.min_confidence = min_confidence
        self.fallback_predictor = fallback_predictor
        self._predictors: Dict[str, OutcomePredictor] = {}
        self._weights: Dict[str, float] = {}

    def register(
        self,
        predictor: OutcomePredictor,
        weight: float = 1.0,
        name: Optional[str] = None,
    ) -> None:
        """
        Register a predictor with the stack.

        Parameters
        ----------
        predictor : OutcomePredictor
            Predictor to register
        weight : float
            Weight for aggregation (higher = more influence)
        name : str, optional
            Override predictor name
        """
        key = name or predictor.name
        self._predictors[key] = predictor
        self._weights[key] = weight
        logger.info(f"Registered predictor: {key} (available={predictor.is_available})")

    def unregister(self, name: str) -> bool:
        """Remove a predictor from the stack."""
        if name in self._predictors:
            del self._predictors[name]
            del self._weights[name]
            return True
        return False

    @property
    def predictors(self) -> Dict[str, OutcomePredictor]:
        """Registered predictors."""
        return self._predictors.copy()

    @property
    def available_predictors(self) -> List[str]:
        """Names of predictors with available dependencies."""
        return [
            name for name, pred in self._predictors.items()
            if pred.is_available
        ]

    def predict(
        self,
        sequence: str,
        target_gene: Optional[str] = None,
        cell_type: Optional[str] = None,
        context_5p: Optional[str] = None,
        context_3p: Optional[str] = None,
        **kwargs,
    ) -> EnsemblePrediction:
        """
        Get ensemble prediction from all registered predictors.

        Parameters
        ----------
        sequence : str
            Guide sequence
        target_gene : str, optional
            Target gene
        cell_type : str, optional
            Cell type
        context_5p : str, optional
            5' context
        context_3p : str, optional
            3' context
        **kwargs
            Additional predictor parameters

        Returns
        -------
        EnsemblePrediction
            Combined prediction from all available predictors
        """
        predictions: Dict[str, PredictionResult] = {}
        warnings: List[str] = []

        # Collect predictions from all available predictors
        for name, predictor in self._predictors.items():
            if not predictor.is_available:
                warnings.append(f"{name}: unavailable (missing dependencies)")
                continue

            try:
                result = predictor.predict(
                    sequence=sequence,
                    target_gene=target_gene,
                    cell_type=cell_type,
                    context_5p=context_5p,
                    context_3p=context_3p,
                    **kwargs,
                )

                if result.confidence >= self.min_confidence:
                    predictions[name] = result
                else:
                    warnings.append(f"{name}: low confidence ({result.confidence:.2f})")

            except Exception as e:
                warnings.append(f"{name}: error ({e})")
                logger.warning(f"Predictor {name} failed: {e}")

        # Use fallback if no predictions
        if not predictions and self.fallback_predictor:
            try:
                result = self.fallback_predictor.predict(
                    sequence=sequence,
                    target_gene=target_gene,
                    cell_type=cell_type,
                    **kwargs,
                )
                predictions[self.fallback_predictor.name] = result
                warnings.append("Using fallback predictor")
            except Exception as e:
                warnings.append(f"Fallback failed: {e}")

        # Return neutral if still no predictions
        if not predictions:
            return EnsemblePrediction(
                score=0.5,
                confidence=0.0,
                confidence_interval=(0.0, 1.0),
                predictions={},
                agreement=0.0,
                dominant_type=PredictorType.EFFICIENCY,
                warnings=warnings + ["No predictions available"],
            )

        # Aggregate predictions
        return self._aggregate(predictions, warnings)

    def _aggregate(
        self,
        predictions: Dict[str, PredictionResult],
        warnings: List[str],
    ) -> EnsemblePrediction:
        """Aggregate multiple predictions into ensemble."""
        scores = []
        weights = []
        confidences = []

        for name, pred in predictions.items():
            scores.append(pred.score)
            confidences.append(pred.confidence)
            weights.append(self._weights.get(name, 1.0) * pred.confidence)

        scores = np.array(scores)
        weights = np.array(weights)
        confidences = np.array(confidences)

        # Normalize weights
        if weights.sum() > 0:
            weights = weights / weights.sum()
        else:
            weights = np.ones_like(weights) / len(weights)

        # Aggregate based on method
        if self.aggregation == "weighted":
            agg_score = float(np.dot(scores, weights))
        elif self.aggregation == "mean":
            agg_score = float(np.mean(scores))
        elif self.aggregation == "max":
            agg_score = float(np.max(scores))
        elif self.aggregation == "min":
            agg_score = float(np.min(scores))
        else:
            agg_score = float(np.mean(scores))

        # Calculate agreement (inverse of std dev)
        if len(scores) > 1:
            agreement = float(1.0 - np.std(scores))
        else:
            agreement = 0.5

        # Aggregate confidence
        agg_confidence = float(np.mean(confidences) * agreement)

        # Calculate confidence interval
        all_lows = [p.confidence_interval[0] for p in predictions.values()]
        all_highs = [p.confidence_interval[1] for p in predictions.values()]
        ci = (
            float(np.percentile(all_lows, 10)),
            float(np.percentile(all_highs, 90)),
        )

        # Determine dominant prediction type
        type_counts: Dict[PredictorType, int] = {}
        for pred in predictions.values():
            type_counts[pred.predictor_type] = type_counts.get(pred.predictor_type, 0) + 1
        dominant_type = max(type_counts.keys(), key=lambda t: type_counts[t])

        return EnsemblePrediction(
            score=agg_score,
            confidence=agg_confidence,
            confidence_interval=ci,
            predictions=predictions,
            agreement=agreement,
            dominant_type=dominant_type,
            warnings=warnings,
        )

    def predict_batch(
        self,
        sequences: List[str],
        target_genes: Optional[List[str]] = None,
        cell_type: Optional[str] = None,
        **kwargs,
    ) -> List[EnsemblePrediction]:
        """
        Get ensemble predictions for multiple sequences.

        Parameters
        ----------
        sequences : List[str]
            Guide sequences
        target_genes : List[str], optional
            Target genes
        cell_type : str, optional
            Cell type
        **kwargs
            Additional parameters

        Returns
        -------
        List[EnsemblePrediction]
            Ensemble predictions for each sequence
        """
        genes = target_genes or [None] * len(sequences)
        if len(genes) == 1:
            genes = genes * len(sequences)

        return [
            self.predict(
                sequence=seq,
                target_gene=gene,
                cell_type=cell_type,
                **kwargs,
            )
            for seq, gene in zip(sequences, genes)
        ]

    @classmethod
    def default_stack(cls) -> "PredictorStack":
        """
        Create a default predictor stack.

        Includes all available predictors with sensible weights.
        """
        from phaselab.predictors.adapters import (
            DeepCRISPRAdapter,
            DeepSpCas9Adapter,
            RuleBasedPredictor,
        )

        stack = cls(aggregation="weighted", min_confidence=0.3)

        # Register ML predictors with higher weight
        stack.register(DeepCRISPRAdapter(), weight=1.5)
        stack.register(DeepSpCas9Adapter(), weight=1.5)

        # Rule-based as fallback with lower weight
        stack.register(RuleBasedPredictor(), weight=0.5)

        return stack

    def summary(self) -> str:
        """Get summary of registered predictors."""
        lines = ["PredictorStack Summary:", "=" * 40]
        for name, pred in self._predictors.items():
            status = "available" if pred.is_available else "unavailable"
            weight = self._weights[name]
            lines.append(f"  {name}: {status} (weight={weight:.2f})")
        lines.append(f"  Aggregation: {self.aggregation}")
        lines.append(f"  Min confidence: {self.min_confidence}")
        return "\n".join(lines)
