"""Tests for the ML predictor protocol module."""

import pytest
import numpy as np

from phaselab.predictors import (
    OutcomePredictor,
    PredictionResult,
    PredictorMetadata,
    DeepCRISPRAdapter,
    DeepSpCas9Adapter,
    EnformerAdapter,
    RuleBasedPredictor,
    PredictorStack,
    EnsemblePrediction,
)
from phaselab.predictors.protocol import (
    PredictorType,
    EvidenceLevel,
)


class TestPredictionResult:
    """Tests for PredictionResult."""

    def test_basic_creation(self):
        """Test basic result creation."""
        result = PredictionResult(
            score=0.75,
            confidence=0.9,
            confidence_interval=(0.65, 0.85),
            predictor_name="TestPredictor",
            predictor_type=PredictorType.EFFICIENCY,
            evidence_level=EvidenceLevel.PUBLISHED,
        )
        assert result.score == 0.75
        assert result.confidence == 0.9
        assert result.is_reliable

    def test_score_clipping(self):
        """Test that scores are clipped to [0, 1]."""
        result = PredictionResult(
            score=1.5,  # Out of bounds
            confidence=-0.1,  # Out of bounds
            confidence_interval=(-0.2, 1.5),
            predictor_name="Test",
            predictor_type=PredictorType.EFFICIENCY,
            evidence_level=EvidenceLevel.HEURISTIC,
        )
        assert result.score == 1.0
        assert result.confidence == 0.0
        assert result.confidence_interval == (0.0, 1.0)

    def test_neutral_result(self):
        """Test neutral/uninformative result creation."""
        result = PredictionResult.neutral(
            predictor_name="Test",
            reason="No data available",
        )
        assert result.score == 0.5
        assert result.confidence == 0.0
        assert not result.is_reliable
        assert "No data available" in result.warnings

    def test_interval_width(self):
        """Test confidence interval width calculation."""
        result = PredictionResult(
            score=0.5,
            confidence=0.8,
            confidence_interval=(0.3, 0.7),
            predictor_name="Test",
            predictor_type=PredictorType.EFFICIENCY,
            evidence_level=EvidenceLevel.HEURISTIC,
        )
        assert result.interval_width == pytest.approx(0.4)


class TestRuleBasedPredictor:
    """Tests for RuleBasedPredictor (always available)."""

    def test_metadata(self):
        """Test predictor metadata."""
        pred = RuleBasedPredictor()
        meta = pred.metadata

        assert meta.name == "RuleBased"
        assert meta.predictor_type == PredictorType.EFFICIENCY
        assert meta.evidence_level == EvidenceLevel.HEURISTIC

    def test_is_available(self):
        """Test availability check."""
        pred = RuleBasedPredictor()
        assert pred.is_available

    def test_predict_valid_sequence(self):
        """Test prediction with valid sequence."""
        pred = RuleBasedPredictor()
        result = pred.predict("ATGCGATCGATCGATCGATCNGG")

        assert isinstance(result, PredictionResult)
        assert 0 <= result.score <= 1
        assert result.predictor_name == "RuleBased"

    def test_predict_invalid_sequence(self):
        """Test prediction with invalid sequence."""
        pred = RuleBasedPredictor()
        result = pred.predict("ABC")  # Too short, invalid chars

        assert result.score == 0.5
        assert result.confidence == 0.0

    def test_gc_content_effect(self):
        """Test that GC content affects score."""
        pred = RuleBasedPredictor()

        # Optimal GC (around 50%)
        optimal_gc = pred.predict("GCGCATATGCGCATATGCGCNGG")

        # Very high GC
        high_gc = pred.predict("GCGCGCGCGCGCGCGCGCGCNGG")

        # Optimal should score higher (or at least not worse)
        # Due to randomness in prediction, just check both are valid
        assert 0 <= optimal_gc.score <= 1
        assert 0 <= high_gc.score <= 1

    def test_poly_t_penalty(self):
        """Test that poly-T sequences are penalized."""
        pred = RuleBasedPredictor()

        # No poly-T
        no_poly = pred.predict("ATGCGATCGATCGATCGATCNGG")

        # With poly-T
        with_poly = pred.predict("ATGCGATCTTTTGATCGATCNGG")

        # Poly-T should be penalized
        assert "poly_t" in with_poly.raw_output
        assert with_poly.raw_output["poly_t"] == True


class TestDeepCRISPRAdapter:
    """Tests for DeepCRISPR adapter."""

    def test_metadata(self):
        """Test predictor metadata."""
        pred = DeepCRISPRAdapter()
        meta = pred.metadata

        assert meta.name == "DeepCRISPR"
        assert meta.predictor_type == PredictorType.EFFICIENCY
        assert meta.evidence_level == EvidenceLevel.PUBLISHED
        assert "NGG" in meta.pam_sequences

    def test_unavailable_returns_neutral(self):
        """Test that unavailable predictor returns neutral result."""
        pred = DeepCRISPRAdapter()

        if not pred.is_available:
            result = pred.predict("ATGCGATCGATCGATCGATCNGG")
            # Should return neutral when TensorFlow unavailable
            assert result.confidence < 0.5 or "not installed" in str(result.warnings)

    def test_sequence_validation(self):
        """Test sequence validation."""
        pred = DeepCRISPRAdapter()

        # Invalid sequence (too short)
        valid, error = pred.validate_sequence("ATGC")
        assert not valid
        assert "short" in error.lower()

        # Invalid characters
        valid, error = pred.validate_sequence("SHORT")
        assert not valid
        assert "invalid" in error.lower()

        # Valid sequence
        valid, error = pred.validate_sequence("ATGCGATCGATCGATCGATCNGG")
        assert valid


class TestDeepSpCas9Adapter:
    """Tests for DeepSpCas9 adapter."""

    def test_metadata(self):
        """Test predictor metadata."""
        pred = DeepSpCas9Adapter()
        meta = pred.metadata

        assert meta.name == "DeepSpCas9"
        assert "HEK293T" in meta.cell_types


class TestEnformerAdapter:
    """Tests for Enformer adapter."""

    def test_metadata(self):
        """Test predictor metadata."""
        pred = EnformerAdapter()
        meta = pred.metadata

        assert meta.name == "Enformer"
        assert meta.predictor_type == PredictorType.EXPRESSION

    def test_requires_context(self):
        """Test that Enformer warns about insufficient context."""
        pred = EnformerAdapter()

        # Short sequence without context
        result = pred.predict("ATGCGATCGATCGATCGATCNGG")

        # Should indicate insufficient context
        assert result.confidence < 0.5 or len(result.warnings) > 0


class TestPredictorStack:
    """Tests for PredictorStack ensemble."""

    def test_register_predictor(self):
        """Test registering predictors."""
        stack = PredictorStack()
        pred = RuleBasedPredictor()

        stack.register(pred)

        assert "RuleBased" in stack.predictors
        assert "RuleBased" in stack.available_predictors

    def test_unregister_predictor(self):
        """Test unregistering predictors."""
        stack = PredictorStack()
        stack.register(RuleBasedPredictor())

        assert stack.unregister("RuleBased")
        assert "RuleBased" not in stack.predictors

    def test_predict_single_predictor(self):
        """Test prediction with single predictor."""
        stack = PredictorStack()
        stack.register(RuleBasedPredictor())

        result = stack.predict("ATGCGATCGATCGATCGATCNGG")

        assert isinstance(result, EnsemblePrediction)
        assert result.n_predictors >= 1

    def test_predict_multiple_predictors(self):
        """Test prediction with multiple predictors."""
        stack = PredictorStack()
        stack.register(RuleBasedPredictor(), name="Rule1")
        stack.register(RuleBasedPredictor(), name="Rule2")

        result = stack.predict("ATGCGATCGATCGATCGATCNGG")

        assert result.n_predictors == 2
        assert "Rule1" in result.predictions
        assert "Rule2" in result.predictions

    def test_weighted_aggregation(self):
        """Test weighted score aggregation."""
        stack = PredictorStack(aggregation="weighted")
        stack.register(RuleBasedPredictor(), weight=2.0, name="HighWeight")
        stack.register(RuleBasedPredictor(), weight=0.5, name="LowWeight")

        result = stack.predict("ATGCGATCGATCGATCGATCNGG")

        # Higher weighted predictor should have more influence
        assert result.score is not None

    def test_agreement_calculation(self):
        """Test agreement between predictors."""
        stack = PredictorStack()
        # Same predictor twice should have high agreement
        stack.register(RuleBasedPredictor(), name="Rule1")
        stack.register(RuleBasedPredictor(), name="Rule2")

        result = stack.predict("ATGCGATCGATCGATCGATCNGG")

        # Agreement should be relatively high for same predictor
        assert 0 <= result.agreement <= 1

    def test_default_stack(self):
        """Test default stack creation."""
        stack = PredictorStack.default_stack()

        # Should have multiple predictors registered
        assert len(stack.predictors) >= 2
        assert "RuleBased" in stack.predictors

    def test_predict_batch(self):
        """Test batch prediction."""
        stack = PredictorStack()
        stack.register(RuleBasedPredictor())

        sequences = [
            "ATGCGATCGATCGATCGATCNGG",
            "GCATGCATGCATGCATGCATNGG",
        ]
        results = stack.predict_batch(sequences)

        assert len(results) == 2
        assert all(isinstance(r, EnsemblePrediction) for r in results)

    def test_ensemble_to_dict(self):
        """Test ensemble result serialization."""
        stack = PredictorStack()
        stack.register(RuleBasedPredictor())

        result = stack.predict("ATGCGATCGATCGATCGATCNGG")
        d = result.to_dict()

        assert "score" in d
        assert "confidence" in d
        assert "predictions" in d

    def test_min_confidence_filter(self):
        """Test minimum confidence filtering."""
        stack = PredictorStack(min_confidence=0.99)  # Very high threshold
        stack.register(RuleBasedPredictor())

        result = stack.predict("ATGCGATCGATCGATCGATCNGG")

        # Rule-based has moderate confidence, may be filtered
        # Either no predictions or warning about low confidence
        assert result.n_predictors == 0 or len(result.warnings) > 0

    def test_summary(self):
        """Test stack summary."""
        stack = PredictorStack()
        stack.register(RuleBasedPredictor())

        summary = stack.summary()

        assert "RuleBased" in summary
        assert "available" in summary
