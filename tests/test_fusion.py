"""Tests for the evidence fusion module."""

import pytest
import numpy as np

from phaselab.fusion import (
    Evidence,
    EvidenceSource,
    EvidenceType,
    EvidenceFusion,
    FusedResult,
    FusionConfig,
    Calibrator,
    CalibrationCurve,
)
from phaselab.fusion.evidence import EvidenceCollection
from phaselab.core.constants import E_MINUS_2


class TestEvidence:
    """Tests for Evidence dataclass."""

    def test_soft_score_creation(self):
        """Test creating soft score evidence."""
        ev = Evidence.soft_score(
            source=EvidenceSource.GC_CONTENT,
            score=0.75,
            confidence=0.9,
        )
        assert ev.value == 0.75
        assert ev.confidence == 0.9
        assert ev.evidence_type == EvidenceType.SOFT_SCORE
        assert not ev.is_gate

    def test_hard_gate_creation(self):
        """Test creating hard gate evidence."""
        ev = Evidence.hard_gate(
            source=EvidenceSource.OFF_TARGET_COUNT,
            passes=True,
            confidence=1.0,
        )
        assert ev.value == 1.0
        assert ev.is_gate
        assert ev.passes_gate

    def test_gate_fail(self):
        """Test failed hard gate."""
        ev = Evidence.hard_gate(
            source=EvidenceSource.OFF_TARGET_COUNT,
            passes=False,
        )
        assert ev.value == 0.0
        assert not ev.passes_gate

    def test_coherence_creation(self):
        """Test creating coherence evidence."""
        ev = Evidence.coherence(
            r_bar=0.45,
            mode="heuristic",
        )
        assert ev.value == 0.45
        assert ev.is_coherence
        assert ev.source == EvidenceSource.COHERENCE_HEURISTIC

    def test_quantum_coherence(self):
        """Test quantum coherence evidence."""
        ev = Evidence.coherence(
            r_bar=0.65,
            mode="quantum",
        )
        assert ev.source == EvidenceSource.COHERENCE_QUANTUM

    def test_value_clipping(self):
        """Test that values are clipped to valid range."""
        ev = Evidence.soft_score(
            source=EvidenceSource.GC_CONTENT,
            score=1.5,  # Out of bounds
            confidence=-0.1,  # Out of bounds
        )
        assert ev.value == 1.0
        assert ev.confidence == 0.0

    def test_effective_weight(self):
        """Test effective weight calculation."""
        ev = Evidence.soft_score(
            source=EvidenceSource.GC_CONTENT,
            score=0.5,
            confidence=0.8,
            weight=2.0,
        )
        assert ev.effective_weight == 1.6  # 0.8 * 2.0

    def test_to_dict(self):
        """Test evidence serialization."""
        ev = Evidence.soft_score(
            source=EvidenceSource.GC_CONTENT,
            score=0.75,
        )
        d = ev.to_dict()

        assert d["source"] == "gc_content"
        assert d["value"] == 0.75


class TestEvidenceCollection:
    """Tests for EvidenceCollection."""

    def test_add_evidence(self):
        """Test adding evidence to collection."""
        coll = EvidenceCollection(guide_id="test")
        coll.add(Evidence.soft_score(EvidenceSource.GC_CONTENT, 0.5))
        coll.add(Evidence.hard_gate(EvidenceSource.OFF_TARGET_COUNT, True))

        assert len(coll.evidences) == 2

    def test_gates_property(self):
        """Test filtering gates."""
        coll = EvidenceCollection(guide_id="test")
        coll.add(Evidence.soft_score(EvidenceSource.GC_CONTENT, 0.5))
        coll.add(Evidence.hard_gate(EvidenceSource.OFF_TARGET_COUNT, True))
        coll.add(Evidence.hard_gate(EvidenceSource.GC_CONTENT, False))

        assert len(coll.gates) == 2

    def test_scores_property(self):
        """Test filtering soft scores."""
        coll = EvidenceCollection(guide_id="test")
        coll.add(Evidence.soft_score(EvidenceSource.GC_CONTENT, 0.5))
        coll.add(Evidence.hard_gate(EvidenceSource.OFF_TARGET_COUNT, True))
        coll.add(Evidence.coherence(0.45, "heuristic"))

        assert len(coll.scores) == 1

    def test_coherence_property(self):
        """Test getting coherence evidence."""
        coll = EvidenceCollection(guide_id="test")
        coll.add(Evidence.soft_score(EvidenceSource.GC_CONTENT, 0.5))
        coll.add(Evidence.coherence(0.45, "heuristic"))

        assert coll.coherence is not None
        assert coll.coherence.value == 0.45

    def test_passes_all_gates(self):
        """Test gate checking."""
        coll = EvidenceCollection(guide_id="test")
        coll.add(Evidence.hard_gate(EvidenceSource.OFF_TARGET_COUNT, True))
        coll.add(Evidence.hard_gate(EvidenceSource.GC_CONTENT, True))

        assert coll.passes_all_gates

        # Add failing gate
        coll.add(Evidence.hard_gate(EvidenceSource.HOMOPOLYMER, False))
        assert not coll.passes_all_gates

    def test_failed_gates(self):
        """Test getting failed gates."""
        coll = EvidenceCollection(guide_id="test")
        coll.add(Evidence.hard_gate(EvidenceSource.OFF_TARGET_COUNT, True))
        coll.add(Evidence.hard_gate(EvidenceSource.GC_CONTENT, False))

        failed = coll.failed_gates
        assert len(failed) == 1
        assert failed[0].source == EvidenceSource.GC_CONTENT


class TestFusionConfig:
    """Tests for FusionConfig."""

    def test_default_config(self):
        """Test default configuration."""
        config = FusionConfig.default()

        assert config.sequence_weight == 0.25
        assert config.go_threshold == E_MINUS_2

    def test_conservative_config(self):
        """Test conservative configuration."""
        config = FusionConfig.conservative()

        assert config.specificity_weight > config.sequence_weight

    def test_ml_focused_config(self):
        """Test ML-focused configuration."""
        config = FusionConfig.ml_focused()

        assert config.ml_weight > config.sequence_weight


class TestEvidenceFusion:
    """Tests for EvidenceFusion engine."""

    def test_basic_fusion(self):
        """Test basic evidence fusion."""
        fusion = EvidenceFusion()

        fusion.add_sequence_score(0.75)
        fusion.add_coherence(0.45, mode="heuristic")

        result = fusion.fuse()

        assert isinstance(result, FusedResult)
        assert 0 <= result.score <= 1

    def test_go_status_above_threshold(self):
        """Test GO status when above threshold."""
        fusion = EvidenceFusion()

        fusion.add_sequence_score(0.75)
        fusion.add_coherence(0.5, mode="heuristic")  # Above e^-2

        result = fusion.fuse()

        assert result.go_status == "GO"
        assert result.is_go

    def test_nogo_status_below_threshold(self):
        """Test NO-GO status when below threshold."""
        fusion = EvidenceFusion()

        fusion.add_sequence_score(0.75)
        fusion.add_coherence(0.1, mode="heuristic")  # Below e^-2

        result = fusion.fuse()

        assert result.go_status == "NO-GO"
        assert not result.is_go

    def test_nogo_penalty(self):
        """Test that NO-GO penalizes score."""
        fusion_go = EvidenceFusion()
        fusion_go.add_sequence_score(0.80)
        fusion_go.add_coherence(0.5, mode="heuristic")
        result_go = fusion_go.fuse()

        fusion_nogo = EvidenceFusion()
        fusion_nogo.add_sequence_score(0.80)
        fusion_nogo.add_coherence(0.1, mode="heuristic")
        result_nogo = fusion_nogo.fuse()

        # NO-GO should have lower score due to penalty
        assert result_nogo.score < result_go.score

    def test_hard_gate_failure(self):
        """Test hard gate failure."""
        fusion = EvidenceFusion()

        fusion.add_sequence_score(0.90)
        fusion.add_hard_gate(EvidenceSource.OFF_TARGET_COUNT, passes=False)
        fusion.add_coherence(0.5, mode="heuristic")

        result = fusion.fuse()

        assert not result.passes_gates
        assert len(result.failed_gates) > 0

    def test_multiple_evidence_sources(self):
        """Test fusion with multiple evidence sources."""
        fusion = EvidenceFusion()

        fusion.add_sequence_score(0.75, source=EvidenceSource.GC_CONTENT)
        fusion.add_sequence_score(0.80, source=EvidenceSource.THERMODYNAMICS)
        fusion.add_context_score(0.70, source=EvidenceSource.ACCESSIBILITY)
        fusion.add_ml_prediction(0.85, source=EvidenceSource.DEEPCRISPR)
        fusion.add_specificity_score(0.90, source=EvidenceSource.OFF_TARGET_CFD)
        fusion.add_coherence(0.45, mode="heuristic")

        result = fusion.fuse()

        assert result.evidence_count == 6
        assert "sequence" in result.component_scores
        assert "context" in result.component_scores
        assert "ml" in result.component_scores
        assert "specificity" in result.component_scores

    def test_reset(self):
        """Test fusion reset."""
        fusion = EvidenceFusion()

        fusion.add_sequence_score(0.75)
        fusion.reset()

        result = fusion.fuse()

        # After reset, should have no evidence
        assert result.evidence_count == 0

    def test_fuse_guide_convenience(self):
        """Test convenience method for guide fusion."""
        result = EvidenceFusion.fuse_guide(
            sequence_scores={"gc": 0.75, "thermo": 0.80},
            context_scores={"accessibility": 0.70},
            ml_scores={"deepcrispr": 0.85},
            coherence=(0.45, "heuristic"),
        )

        assert isinstance(result, FusedResult)
        assert result.coherence == 0.45

    def test_confidence_interval(self):
        """Test confidence interval calculation."""
        fusion = EvidenceFusion()

        fusion.add_sequence_score(0.75, confidence=0.9)
        fusion.add_ml_prediction(0.80, confidence=0.85)
        fusion.add_coherence(0.5, mode="heuristic")

        result = fusion.fuse()

        low, high = result.confidence_interval
        assert low <= result.score <= high
        assert low >= 0.0
        assert high <= 1.0

    def test_result_to_dict(self):
        """Test result serialization."""
        fusion = EvidenceFusion()

        fusion.add_sequence_score(0.75)
        fusion.add_coherence(0.45, mode="heuristic")

        result = fusion.fuse()
        d = result.to_dict()

        assert "score" in d
        assert "go_status" in d
        assert "confidence_interval" in d

    def test_is_viable(self):
        """Test viability check."""
        # Viable: passes gates and GO
        fusion1 = EvidenceFusion()
        fusion1.add_sequence_score(0.75)
        fusion1.add_hard_gate(EvidenceSource.OFF_TARGET_COUNT, True)
        fusion1.add_coherence(0.5, "heuristic")
        result1 = fusion1.fuse()
        assert result1.is_viable

        # Not viable: fails gate
        fusion2 = EvidenceFusion()
        fusion2.add_sequence_score(0.75)
        fusion2.add_hard_gate(EvidenceSource.OFF_TARGET_COUNT, False)
        fusion2.add_coherence(0.5, "heuristic")
        result2 = fusion2.fuse()
        assert not result2.is_viable

        # Not viable: NO-GO
        fusion3 = EvidenceFusion()
        fusion3.add_sequence_score(0.75)
        fusion3.add_hard_gate(EvidenceSource.OFF_TARGET_COUNT, True)
        fusion3.add_coherence(0.1, "heuristic")
        result3 = fusion3.fuse()
        assert not result3.is_viable


class TestCalibrationCurve:
    """Tests for CalibrationCurve."""

    def test_identity_curve(self):
        """Test identity calibration curve."""
        curve = CalibrationCurve.identity("test")

        assert curve.calibrate(0.5) == pytest.approx(0.5)
        assert curve.calibrate(0.0) == pytest.approx(0.0)
        assert curve.calibrate(1.0) == pytest.approx(1.0)

    def test_from_bins(self):
        """Test creating curve from bins."""
        curve = CalibrationCurve.from_bins(
            name="test",
            bin_edges=[0.0, 0.5, 1.0],
            observed_rates=[0.3, 0.7],
        )

        # Should calibrate according to observed rates
        assert curve.calibrate(0.25) < 0.5
        assert curve.calibrate(0.75) > 0.5

    def test_calibrate_batch(self):
        """Test batch calibration."""
        curve = CalibrationCurve.identity("test")

        scores = np.array([0.0, 0.5, 1.0])
        calibrated = curve.calibrate_batch(scores)

        np.testing.assert_array_almost_equal(scores, calibrated)


class TestCalibrator:
    """Tests for Calibrator."""

    def test_default_curves(self):
        """Test default calibration curves are loaded."""
        cal = Calibrator()

        assert "DeepCRISPR" in cal.available_predictors
        assert "RuleBased" in cal.available_predictors

    def test_calibrate_known_predictor(self):
        """Test calibrating known predictor."""
        cal = Calibrator()

        calibrated = cal.calibrate("DeepCRISPR", 0.5)

        assert 0 <= calibrated <= 1

    def test_calibrate_unknown_predictor(self):
        """Test calibrating unknown predictor returns raw score."""
        cal = Calibrator()

        calibrated = cal.calibrate("UnknownPredictor", 0.75)

        assert calibrated == 0.75

    def test_fit_binned(self):
        """Test fitting binned calibration curve."""
        cal = Calibrator(method="binned")

        # Generate synthetic calibration data
        np.random.seed(42)
        n = 100
        predicted = np.random.uniform(0, 1, n)
        # Simulated outcomes correlated with predictions
        outcomes = (np.random.uniform(0, 1, n) < predicted).astype(float)

        curve = cal.fit(
            predictor_name="TestPredictor",
            predicted_scores=predicted,
            actual_outcomes=outcomes,
            n_bins=5,
        )

        assert curve.name == "TestPredictor"
        assert "TestPredictor" in cal.available_predictors

    def test_reliability_diagram(self):
        """Test reliability diagram data generation."""
        cal = Calibrator()

        np.random.seed(42)
        predicted = np.random.uniform(0, 1, 50)
        outcomes = (np.random.uniform(0, 1, 50) < predicted).astype(float)

        data = cal.reliability_diagram(
            predictor_name="Test",
            predicted_scores=predicted,
            actual_outcomes=outcomes,
            n_bins=5,
        )

        assert "bin_centers" in data
        assert "observed_rates" in data
        assert "predicted_rates" in data
        assert "counts" in data
        assert len(data["bin_centers"]) == 5
