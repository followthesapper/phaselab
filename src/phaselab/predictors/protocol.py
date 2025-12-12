"""
Protocol definition for ML outcome predictors.

Defines the interface that all prediction adapters must implement
for integration with the PhaseLab scoring pipeline.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional, Dict, Any, List, Tuple
import logging

import numpy as np

logger = logging.getLogger(__name__)


class PredictorType(Enum):
    """Types of ML predictors."""
    EFFICIENCY = "efficiency"          # On-target editing efficiency
    SPECIFICITY = "specificity"        # Off-target propensity
    EXPRESSION = "expression"          # Gene expression change
    INDEL = "indel"                    # Indel outcome prediction
    REPAIR = "repair"                  # Repair pathway prediction
    COMBINED = "combined"              # Multi-output predictor


class EvidenceLevel(Enum):
    """Evidence level for predictions."""
    VALIDATED = "A"      # Experimentally validated model
    PUBLISHED = "B"      # Peer-reviewed but not independently validated
    PRETRAINED = "C"     # Pre-trained, not validated on target cell type
    HEURISTIC = "D"      # Rule-based, not ML


@dataclass
class PredictorMetadata:
    """
    Metadata describing a predictor.

    Attributes
    ----------
    name : str
        Predictor identifier
    version : str
        Model version
    predictor_type : PredictorType
        What the predictor outputs
    evidence_level : EvidenceLevel
        Validation status
    citation : str, optional
        Publication reference
    model_url : str, optional
        URL to model weights/repository
    cell_types : List[str]
        Cell types model was trained on
    pam_sequences : List[str]
        PAM sequences supported (e.g., ["NGG", "NAG"])
    """
    name: str
    version: str
    predictor_type: PredictorType
    evidence_level: EvidenceLevel
    citation: Optional[str] = None
    model_url: Optional[str] = None
    cell_types: List[str] = field(default_factory=list)
    pam_sequences: List[str] = field(default_factory=lambda: ["NGG"])


@dataclass
class PredictionResult:
    """
    Standardized prediction output.

    All ML predictors return this structure regardless of their
    internal implementation.

    Attributes
    ----------
    score : float
        Primary prediction score [0, 1]
        For efficiency: 1 = high efficiency
        For specificity: 1 = high specificity (few off-targets)
    confidence : float
        Prediction confidence [0, 1]
    confidence_interval : Tuple[float, float]
        95% confidence interval for score
    predictor_name : str
        Name of predictor that generated this
    predictor_type : PredictorType
        Type of prediction
    evidence_level : EvidenceLevel
        Evidence strength
    raw_output : Dict[str, Any]
        Raw predictor output for debugging
    warnings : List[str]
        Any warnings generated during prediction
    """
    score: float
    confidence: float
    confidence_interval: Tuple[float, float]
    predictor_name: str
    predictor_type: PredictorType
    evidence_level: EvidenceLevel
    raw_output: Dict[str, Any] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)

    def __post_init__(self):
        # Ensure score is in valid range
        self.score = float(np.clip(self.score, 0.0, 1.0))
        self.confidence = float(np.clip(self.confidence, 0.0, 1.0))

        # Validate confidence interval
        low, high = self.confidence_interval
        self.confidence_interval = (
            float(np.clip(low, 0.0, 1.0)),
            float(np.clip(high, 0.0, 1.0)),
        )

    @property
    def is_reliable(self) -> bool:
        """Whether prediction meets minimum reliability threshold."""
        return self.confidence >= 0.5

    @property
    def interval_width(self) -> float:
        """Width of confidence interval."""
        return self.confidence_interval[1] - self.confidence_interval[0]

    @classmethod
    def neutral(
        cls,
        predictor_name: str,
        predictor_type: PredictorType = PredictorType.EFFICIENCY,
        reason: str = "No prediction available",
    ) -> "PredictionResult":
        """Create a neutral/uninformative prediction."""
        return cls(
            score=0.5,
            confidence=0.0,
            confidence_interval=(0.0, 1.0),
            predictor_name=predictor_name,
            predictor_type=predictor_type,
            evidence_level=EvidenceLevel.HEURISTIC,
            warnings=[reason],
        )


class OutcomePredictor(ABC):
    """
    Abstract base class for ML outcome predictors.

    All prediction adapters must implement this interface to integrate
    with the PhaseLab scoring pipeline.

    Subclasses should:
    1. Implement metadata property with model information
    2. Implement predict() method for single sequences
    3. Optionally override predict_batch() for efficiency
    4. Handle missing dependencies gracefully
    """

    @property
    @abstractmethod
    def metadata(self) -> PredictorMetadata:
        """Return predictor metadata."""
        pass

    @property
    def name(self) -> str:
        """Predictor name."""
        return self.metadata.name

    @property
    def is_available(self) -> bool:
        """
        Check if predictor dependencies are available.

        Returns False if required packages are not installed.
        """
        return True

    @abstractmethod
    def predict(
        self,
        sequence: str,
        target_gene: Optional[str] = None,
        cell_type: Optional[str] = None,
        context_5p: Optional[str] = None,
        context_3p: Optional[str] = None,
        **kwargs,
    ) -> PredictionResult:
        """
        Generate prediction for a single guide sequence.

        Parameters
        ----------
        sequence : str
            Guide sequence (20-23bp typically, including PAM)
        target_gene : str, optional
            Target gene symbol for expression predictors
        cell_type : str, optional
            Cell type for cell-type-specific predictions
        context_5p : str, optional
            5' flanking sequence
        context_3p : str, optional
            3' flanking sequence
        **kwargs
            Additional predictor-specific parameters

        Returns
        -------
        PredictionResult
            Standardized prediction output
        """
        pass

    def predict_batch(
        self,
        sequences: List[str],
        target_genes: Optional[List[str]] = None,
        cell_type: Optional[str] = None,
        **kwargs,
    ) -> List[PredictionResult]:
        """
        Generate predictions for multiple sequences.

        Default implementation calls predict() in a loop.
        Override for batched efficiency.

        Parameters
        ----------
        sequences : List[str]
            List of guide sequences
        target_genes : List[str], optional
            Target genes (same length as sequences, or single value)
        cell_type : str, optional
            Cell type
        **kwargs
            Additional parameters

        Returns
        -------
        List[PredictionResult]
            Predictions for each sequence
        """
        results = []
        genes = target_genes or [None] * len(sequences)
        if len(genes) == 1:
            genes = genes * len(sequences)

        for seq, gene in zip(sequences, genes):
            results.append(
                self.predict(
                    sequence=seq,
                    target_gene=gene,
                    cell_type=cell_type,
                    **kwargs,
                )
            )
        return results

    def validate_sequence(self, sequence: str) -> Tuple[bool, str]:
        """
        Validate input sequence.

        Parameters
        ----------
        sequence : str
            Sequence to validate

        Returns
        -------
        Tuple[bool, str]
            (is_valid, error_message)
        """
        if not sequence:
            return False, "Empty sequence"

        seq_upper = sequence.upper()

        # Check for valid nucleotides
        valid_chars = set("ACGTN")
        if not all(c in valid_chars for c in seq_upper):
            invalid = set(seq_upper) - valid_chars
            return False, f"Invalid characters: {invalid}"

        # Check length
        if len(sequence) < 20:
            return False, f"Sequence too short: {len(sequence)} < 20"

        if len(sequence) > 30:
            return False, f"Sequence too long: {len(sequence)} > 30"

        # Check PAM
        pam = seq_upper[-3:]
        supported_pams = self.metadata.pam_sequences
        pam_valid = any(
            self._matches_pam(pam, supported)
            for supported in supported_pams
        )
        if not pam_valid:
            return False, f"Unsupported PAM: {pam}"

        return True, ""

    @staticmethod
    def _matches_pam(pam: str, pattern: str) -> bool:
        """Check if PAM matches pattern (N = any)."""
        if len(pam) != len(pattern):
            return False
        for p, pat in zip(pam, pattern):
            if pat == "N":
                continue
            if p != pat:
                return False
        return True

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(name={self.name!r}, available={self.is_available})"
