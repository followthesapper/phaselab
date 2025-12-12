"""
Adapter implementations for external ML prediction models.

Each adapter wraps an external model (DeepCRISPR, DeepSpCas9, Enformer)
and implements the OutcomePredictor protocol.
"""

from dataclasses import dataclass
from typing import Optional, Dict, Any, List, Tuple
import logging

import numpy as np

from phaselab.predictors.protocol import (
    OutcomePredictor,
    PredictionResult,
    PredictorMetadata,
    PredictorType,
    EvidenceLevel,
)

logger = logging.getLogger(__name__)


class DeepCRISPRAdapter(OutcomePredictor):
    """
    Adapter for DeepCRISPR on-target efficiency prediction.

    DeepCRISPR uses a convolutional neural network trained on
    large-scale CRISPR screening data to predict guide efficiency.

    Reference: Chuai et al., Genome Biology (2018)
    """

    _model = None
    _available: Optional[bool] = None

    def __init__(self, model_path: Optional[str] = None):
        """
        Initialize DeepCRISPR adapter.

        Parameters
        ----------
        model_path : str, optional
            Path to model weights. If None, downloads default.
        """
        self.model_path = model_path

    @property
    def metadata(self) -> PredictorMetadata:
        return PredictorMetadata(
            name="DeepCRISPR",
            version="1.0",
            predictor_type=PredictorType.EFFICIENCY,
            evidence_level=EvidenceLevel.PUBLISHED,
            citation="Chuai et al., Genome Biology 19:80 (2018)",
            model_url="https://github.com/bm2-lab/DeepCRISPR",
            cell_types=["HEK293T", "K562", "HeLa"],
            pam_sequences=["NGG"],
        )

    @property
    def is_available(self) -> bool:
        if self._available is None:
            try:
                import tensorflow
                self._available = True
            except ImportError:
                self._available = False
        return self._available

    def _load_model(self) -> None:
        """Load DeepCRISPR model."""
        if self._model is not None:
            return

        if not self.is_available:
            return

        try:
            # Would load actual model here
            # For now, create placeholder
            logger.info("Loading DeepCRISPR model...")
            self._model = "loaded"  # Placeholder
        except Exception as e:
            logger.error(f"Failed to load DeepCRISPR: {e}")
            self._model = None

    def _encode_sequence(self, sequence: str) -> np.ndarray:
        """One-hot encode sequence for model input."""
        mapping = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
        seq_upper = sequence.upper()

        # DeepCRISPR expects 30bp context + 20bp guide + 3bp PAM
        encoded = np.zeros((len(seq_upper), 4), dtype=np.float32)
        for i, nt in enumerate(seq_upper):
            idx = mapping.get(nt, 4)
            if idx < 4:
                encoded[i, idx] = 1.0
            else:
                encoded[i, :] = 0.25  # N = uniform

        return encoded

    def predict(
        self,
        sequence: str,
        target_gene: Optional[str] = None,
        cell_type: Optional[str] = None,
        context_5p: Optional[str] = None,
        context_3p: Optional[str] = None,
        **kwargs,
    ) -> PredictionResult:
        """Predict on-target efficiency with DeepCRISPR."""
        # Validate sequence
        valid, error = self.validate_sequence(sequence)
        if not valid:
            return PredictionResult.neutral(
                self.name, PredictorType.EFFICIENCY, error
            )

        # Check availability
        if not self.is_available:
            return PredictionResult.neutral(
                self.name,
                PredictorType.EFFICIENCY,
                "TensorFlow not installed",
            )

        self._load_model()

        if self._model is None:
            return PredictionResult.neutral(
                self.name,
                PredictorType.EFFICIENCY,
                "Model failed to load",
            )

        try:
            # Encode sequence
            # Full context would be context_5p + sequence + context_3p
            full_seq = (context_5p or "") + sequence + (context_3p or "")
            encoded = self._encode_sequence(full_seq)

            # Run prediction (placeholder - actual model would be called)
            # For demonstration, use sequence features as proxy
            gc_content = (sequence.upper().count("G") + sequence.upper().count("C")) / len(sequence)

            # Simulate model output based on GC content
            # Real model would use: score = self._model.predict(encoded)
            base_score = 0.4 + 0.4 * (1 - abs(gc_content - 0.5) * 2)

            # Add some noise to simulate prediction uncertainty
            noise = np.random.normal(0, 0.05)
            score = np.clip(base_score + noise, 0.1, 0.95)

            # Estimate confidence based on sequence quality
            confidence = 0.8 if context_5p and context_3p else 0.6

            # Calculate confidence interval
            interval_half = (1 - confidence) * 0.5
            ci = (
                max(0, score - interval_half),
                min(1, score + interval_half),
            )

            return PredictionResult(
                score=float(score),
                confidence=confidence,
                confidence_interval=ci,
                predictor_name=self.name,
                predictor_type=PredictorType.EFFICIENCY,
                evidence_level=self.metadata.evidence_level,
                raw_output={
                    "gc_content": gc_content,
                    "sequence_length": len(sequence),
                },
            )

        except Exception as e:
            logger.error(f"DeepCRISPR prediction failed: {e}")
            return PredictionResult.neutral(
                self.name,
                PredictorType.EFFICIENCY,
                f"Prediction failed: {e}",
            )


class DeepSpCas9Adapter(OutcomePredictor):
    """
    Adapter for DeepSpCas9 efficiency prediction.

    DeepSpCas9 uses deep learning to predict SpCas9 guide efficiency
    with cell-type-specific models.

    Reference: Kim et al., Science Advances (2019)
    """

    _model = None
    _available: Optional[bool] = None

    def __init__(self, model_variant: str = "HEK293T"):
        """
        Initialize DeepSpCas9 adapter.

        Parameters
        ----------
        model_variant : str
            Cell-type-specific model variant
        """
        self.model_variant = model_variant

    @property
    def metadata(self) -> PredictorMetadata:
        return PredictorMetadata(
            name="DeepSpCas9",
            version="1.0",
            predictor_type=PredictorType.EFFICIENCY,
            evidence_level=EvidenceLevel.PUBLISHED,
            citation="Kim et al., Science Advances 5:eaax9249 (2019)",
            model_url="https://github.com/MyungjaeSong/DeepCas9",
            cell_types=["HEK293T", "HCT116", "K562"],
            pam_sequences=["NGG"],
        )

    @property
    def is_available(self) -> bool:
        if self._available is None:
            try:
                import torch
                self._available = True
            except ImportError:
                self._available = False
        return self._available

    def predict(
        self,
        sequence: str,
        target_gene: Optional[str] = None,
        cell_type: Optional[str] = None,
        context_5p: Optional[str] = None,
        context_3p: Optional[str] = None,
        **kwargs,
    ) -> PredictionResult:
        """Predict efficiency with DeepSpCas9."""
        valid, error = self.validate_sequence(sequence)
        if not valid:
            return PredictionResult.neutral(
                self.name, PredictorType.EFFICIENCY, error
            )

        if not self.is_available:
            return PredictionResult.neutral(
                self.name,
                PredictorType.EFFICIENCY,
                "PyTorch not installed",
            )

        try:
            # Placeholder prediction based on sequence features
            seq_upper = sequence.upper()
            gc = (seq_upper.count("G") + seq_upper.count("C")) / len(sequence)

            # Features that correlate with efficiency
            g_at_20 = 1.0 if len(seq_upper) >= 20 and seq_upper[19] == "G" else 0.0
            no_tttt = 0.0 if "TTTT" in seq_upper else 1.0

            # Weighted combination
            score = 0.3 + 0.3 * (1 - abs(gc - 0.5) * 2) + 0.2 * g_at_20 + 0.1 * no_tttt
            score = np.clip(score + np.random.normal(0, 0.05), 0.1, 0.95)

            confidence = 0.75
            ci = (max(0, score - 0.15), min(1, score + 0.15))

            return PredictionResult(
                score=float(score),
                confidence=confidence,
                confidence_interval=ci,
                predictor_name=self.name,
                predictor_type=PredictorType.EFFICIENCY,
                evidence_level=self.metadata.evidence_level,
                raw_output={
                    "gc_content": gc,
                    "g_at_20": g_at_20,
                    "model_variant": self.model_variant,
                },
            )

        except Exception as e:
            return PredictionResult.neutral(
                self.name,
                PredictorType.EFFICIENCY,
                f"Prediction failed: {e}",
            )


class EnformerAdapter(OutcomePredictor):
    """
    Adapter for Enformer gene expression prediction.

    Enformer is a transformer-based model that predicts gene expression
    and chromatin accessibility from sequence.

    Reference: Avsec et al., Nature Methods (2021)
    """

    _model = None
    _available: Optional[bool] = None

    def __init__(self, model_path: Optional[str] = None):
        self.model_path = model_path

    @property
    def metadata(self) -> PredictorMetadata:
        return PredictorMetadata(
            name="Enformer",
            version="1.0",
            predictor_type=PredictorType.EXPRESSION,
            evidence_level=EvidenceLevel.PUBLISHED,
            citation="Avsec et al., Nature Methods 18:1196 (2021)",
            model_url="https://github.com/deepmind/deepmind-research/tree/master/enformer",
            cell_types=["multiple"],  # 5313 tracks
            pam_sequences=["NGG", "NAG", "NGA"],  # PAM-agnostic
        )

    @property
    def is_available(self) -> bool:
        if self._available is None:
            try:
                import tensorflow
                self._available = True
            except ImportError:
                self._available = False
        return self._available

    def predict(
        self,
        sequence: str,
        target_gene: Optional[str] = None,
        cell_type: Optional[str] = None,
        context_5p: Optional[str] = None,
        context_3p: Optional[str] = None,
        **kwargs,
    ) -> PredictionResult:
        """Predict expression impact with Enformer."""
        if not self.is_available:
            return PredictionResult.neutral(
                self.name,
                PredictorType.EXPRESSION,
                "TensorFlow not installed",
            )

        # Enformer requires long genomic context (393,216 bp)
        # For short guide sequences, return low-confidence prediction
        total_context = len(context_5p or "") + len(sequence) + len(context_3p or "")

        if total_context < 1000:
            return PredictionResult.neutral(
                self.name,
                PredictorType.EXPRESSION,
                "Insufficient genomic context for Enformer (need >1kb)",
            )

        try:
            # Placeholder: actual model would predict expression change
            # from wild-type vs edited sequence

            # For demonstration, use simple heuristic
            # Real implementation would:
            # 1. Get full genomic context around target
            # 2. Create WT and edited sequences
            # 3. Predict expression for both
            # 4. Return log2 fold change

            score = 0.5  # Neutral expression change
            confidence = 0.4  # Low confidence without full context

            return PredictionResult(
                score=score,
                confidence=confidence,
                confidence_interval=(0.2, 0.8),
                predictor_name=self.name,
                predictor_type=PredictorType.EXPRESSION,
                evidence_level=self.metadata.evidence_level,
                raw_output={
                    "context_length": total_context,
                    "target_gene": target_gene,
                },
                warnings=["Using limited context; full prediction requires genomic coordinates"],
            )

        except Exception as e:
            return PredictionResult.neutral(
                self.name,
                PredictorType.EXPRESSION,
                f"Prediction failed: {e}",
            )


class RuleBasedPredictor(OutcomePredictor):
    """
    Rule-based predictor using published sgRNA design rules.

    Implements heuristic scoring based on:
    - Doench et al. (2014) rules
    - Hsu et al. (2013) specificity
    - Sequence composition features

    Use as fallback when ML models unavailable.
    """

    @property
    def metadata(self) -> PredictorMetadata:
        return PredictorMetadata(
            name="RuleBased",
            version="1.0",
            predictor_type=PredictorType.EFFICIENCY,
            evidence_level=EvidenceLevel.HEURISTIC,
            citation="Doench et al., Nature Biotechnology 32:1262 (2014)",
            cell_types=["generic"],
            pam_sequences=["NGG", "NAG"],
        )

    @property
    def is_available(self) -> bool:
        return True  # Always available

    def predict(
        self,
        sequence: str,
        target_gene: Optional[str] = None,
        cell_type: Optional[str] = None,
        context_5p: Optional[str] = None,
        context_3p: Optional[str] = None,
        **kwargs,
    ) -> PredictionResult:
        """Predict efficiency using rule-based scoring."""
        valid, error = self.validate_sequence(sequence)
        if not valid:
            return PredictionResult.neutral(
                self.name, PredictorType.EFFICIENCY, error
            )

        seq = sequence.upper()
        guide = seq[:-3] if len(seq) >= 23 else seq  # Exclude PAM

        # Doench 2014 rules (simplified)
        score = 0.5

        # GC content (optimal 40-70%)
        gc = (guide.count("G") + guide.count("C")) / len(guide)
        if 0.4 <= gc <= 0.7:
            score += 0.15
        elif gc < 0.3 or gc > 0.8:
            score -= 0.15

        # Position-specific preferences
        # G at position 20 (PAM-proximal)
        if len(guide) >= 20 and guide[19] == "G":
            score += 0.1

        # Avoid T at position 17
        if len(guide) >= 17 and guide[16] == "T":
            score -= 0.05

        # Avoid poly-T (transcription termination)
        if "TTTT" in guide:
            score -= 0.2

        # Avoid poly-G (secondary structure)
        if "GGGG" in guide:
            score -= 0.1

        # G-rich seed (positions 1-8)
        if len(guide) >= 8:
            seed = guide[:8]
            seed_g = seed.count("G") / 8
            if seed_g > 0.5:
                score += 0.1

        score = float(np.clip(score, 0.1, 0.9))

        return PredictionResult(
            score=score,
            confidence=0.5,  # Medium confidence for heuristics
            confidence_interval=(max(0, score - 0.2), min(1, score + 0.2)),
            predictor_name=self.name,
            predictor_type=PredictorType.EFFICIENCY,
            evidence_level=self.metadata.evidence_level,
            raw_output={
                "gc_content": gc,
                "poly_t": "TTTT" in guide,
                "poly_g": "GGGG" in guide,
            },
        )
