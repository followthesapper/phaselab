"""Tests for v0.8.0 claim levels, layer disagreement, and unknown bucket features."""

import pytest
import numpy as np

from phaselab.fusion import (
    ClaimLevel,
    LayerDisagreement,
    EvidenceFusion,
    FusedResult,
    FusionConfig,
    EvidenceSource,
)


class TestClaimLevel:
    """Tests for ClaimLevel enum and classification."""

    def test_claim_level_values(self):
        """Test claim level enum values."""
        assert ClaimLevel.STRONG_COMPUTATIONAL.value == "strong_computational"
        assert ClaimLevel.CONTEXT_DEPENDENT.value == "context_dependent"
        assert ClaimLevel.EXPLORATORY.value == "exploratory"
        assert ClaimLevel.UNKNOWN.value == "unknown"

    def test_from_confidence_unknown(self):
        """Test claim level classification for unknown cases."""
        # Low confidence -> UNKNOWN
        assert ClaimLevel.from_confidence(0.2, 2) == ClaimLevel.UNKNOWN
        assert ClaimLevel.from_confidence(0.1, 3) == ClaimLevel.UNKNOWN

        # Zero agreeing layers -> UNKNOWN
        assert ClaimLevel.from_confidence(0.9, 0) == ClaimLevel.UNKNOWN

    def test_from_confidence_exploratory(self):
        """Test claim level classification for exploratory cases."""
        # Single layer agreement -> EXPLORATORY
        assert ClaimLevel.from_confidence(0.6, 1) == ClaimLevel.EXPLORATORY

        # Low-moderate confidence -> EXPLORATORY
        assert ClaimLevel.from_confidence(0.4, 2) == ClaimLevel.EXPLORATORY

    def test_from_confidence_context_dependent(self):
        """Test claim level classification for context-dependent cases."""
        # Moderate confidence, 2 layers
        assert ClaimLevel.from_confidence(0.6, 2) == ClaimLevel.CONTEXT_DEPENDENT
        assert ClaimLevel.from_confidence(0.7, 2) == ClaimLevel.CONTEXT_DEPENDENT

    def test_from_confidence_strong_computational(self):
        """Test claim level classification for strong computational support."""
        # High confidence, 3+ layers
        assert ClaimLevel.from_confidence(0.85, 3) == ClaimLevel.STRONG_COMPUTATIONAL
        assert ClaimLevel.from_confidence(0.9, 4) == ClaimLevel.STRONG_COMPUTATIONAL

    def test_claim_level_description(self):
        """Test claim level descriptions."""
        assert "Strong computational" in ClaimLevel.STRONG_COMPUTATIONAL.description
        assert "Context-dependent" in ClaimLevel.CONTEXT_DEPENDENT.description
        assert "Exploratory" in ClaimLevel.EXPLORATORY.description
        assert "Unknown" in ClaimLevel.UNKNOWN.description

    def test_claim_level_color_code(self):
        """Test claim level color codes."""
        assert ClaimLevel.STRONG_COMPUTATIONAL.color_code == "green"
        assert ClaimLevel.CONTEXT_DEPENDENT.color_code == "yellow"
        assert ClaimLevel.EXPLORATORY.color_code == "orange"
        assert ClaimLevel.UNKNOWN.color_code == "gray"


class TestLayerDisagreement:
    """Tests for LayerDisagreement detection."""

    def test_no_disagreement_small_delta(self):
        """Test no disagreement for small score differences."""
        result = LayerDisagreement.check_pair("ml", 0.7, "context", 0.75)
        assert result is None  # Delta < 0.2

    def test_no_disagreement_none_values(self):
        """Test no disagreement when values are None."""
        result = LayerDisagreement.check_pair("ml", None, "context", 0.7)
        assert result is None

        result = LayerDisagreement.check_pair("ml", 0.7, "context", None)
        assert result is None

    def test_moderate_disagreement(self):
        """Test moderate disagreement detection."""
        result = LayerDisagreement.check_pair("ml", 0.8, "context", 0.55)

        assert result is not None
        assert result.delta == pytest.approx(0.25, abs=0.01)
        assert not result.is_critical
        assert "Moderate disagreement" in result.message

    def test_critical_disagreement(self):
        """Test critical disagreement detection."""
        result = LayerDisagreement.check_pair("ml", 0.9, "context", 0.4)

        assert result is not None
        assert result.delta == pytest.approx(0.5, abs=0.01)
        assert result.is_critical
        assert "CRITICAL" in result.message

    def test_disagreement_message_content(self):
        """Test disagreement message contains relevant info."""
        result = LayerDisagreement.check_pair("ml", 0.9, "context", 0.4)

        assert "ml" in result.message
        assert "context" in result.message
        assert "0.9" in result.message or "0.90" in result.message
        assert "0.4" in result.message or "0.40" in result.message

    def test_custom_critical_threshold(self):
        """Test custom critical threshold."""
        # With default threshold (0.4), this should be moderate
        result1 = LayerDisagreement.check_pair("ml", 0.8, "context", 0.5)
        assert not result1.is_critical

        # With lower threshold (0.2), this should be critical
        result2 = LayerDisagreement.check_pair(
            "ml", 0.8, "context", 0.5, critical_threshold=0.2
        )
        assert result2.is_critical


class TestFusionClaimLevel:
    """Tests for claim level integration in EvidenceFusion."""

    def test_fusion_unknown_low_confidence(self):
        """Test fusion returns UNKNOWN for low confidence."""
        fusion = EvidenceFusion()

        # Add only weak evidence
        fusion.add_sequence_score(0.5, confidence=0.2)

        result = fusion.fuse()

        assert result.claim_level == ClaimLevel.UNKNOWN
        assert result.is_unknown

    def test_fusion_unknown_no_signal(self):
        """Test fusion returns UNKNOWN when all scores near 0.5."""
        fusion = EvidenceFusion()

        # All scores near 0.5 - no signal
        fusion.add_sequence_score(0.5, confidence=0.8)
        fusion.add_context_score(0.5, confidence=0.8)
        fusion.add_ml_prediction(0.5, confidence=0.8)

        result = fusion.fuse()

        assert result.is_unknown

    def test_fusion_exploratory_single_layer(self):
        """Test fusion returns EXPLORATORY for single layer evidence."""
        fusion = EvidenceFusion()

        # Only ML prediction
        fusion.add_ml_prediction(0.8, confidence=0.7)

        result = fusion.fuse()

        # Should be exploratory (single source)
        assert result.claim_level in [ClaimLevel.EXPLORATORY, ClaimLevel.CONTEXT_DEPENDENT]

    def test_fusion_strong_computational_multi_layer(self):
        """Test fusion returns STRONG_COMPUTATIONAL for multi-layer agreement."""
        fusion = EvidenceFusion()

        # Multiple agreeing layers with high confidence
        fusion.add_sequence_score(0.85, confidence=0.9)
        fusion.add_context_score(0.82, confidence=0.85)
        fusion.add_ml_prediction(0.88, confidence=0.9)
        fusion.add_specificity_score(0.80, confidence=0.9)
        fusion.add_coherence(0.5, mode="quantum")

        result = fusion.fuse()

        assert result.claim_level == ClaimLevel.STRONG_COMPUTATIONAL
        assert not result.is_unknown

    def test_fusion_detects_disagreements(self):
        """Test fusion detects layer disagreements."""
        fusion = EvidenceFusion()

        # ML says high, context says low
        fusion.add_ml_prediction(0.9, confidence=0.9)
        fusion.add_context_score(0.3, confidence=0.9)

        result = fusion.fuse()

        assert len(result.disagreements) > 0
        assert result.has_critical_disagreement

    def test_fusion_to_dict_includes_claim_level(self):
        """Test FusedResult.to_dict includes v0.8.0 fields."""
        fusion = EvidenceFusion()
        fusion.add_sequence_score(0.8, confidence=0.9)
        fusion.add_coherence(0.5)

        result = fusion.fuse()
        d = result.to_dict()

        assert "claim_level" in d
        assert "claim_description" in d
        assert "disagreements" in d
        assert "is_unknown" in d
        assert "n_agreeing_layers" in d
        assert "has_critical_disagreement" in d


class TestFusionDisagreementTracking:
    """Tests for layer disagreement tracking in fusion."""

    def test_ml_vs_context_disagreement(self):
        """Test ML vs context disagreement detection."""
        fusion = EvidenceFusion()

        fusion.add_ml_prediction(0.9, confidence=0.9)
        fusion.add_context_score(0.3, confidence=0.9)
        fusion.add_sequence_score(0.6, confidence=0.9)

        result = fusion.fuse()

        # Should detect ML vs context disagreement
        ml_ctx = [d for d in result.disagreements
                  if "ml" in d.layer1 or "ml" in d.layer2]
        assert len(ml_ctx) > 0

    def test_multiple_critical_disagreements_trigger_unknown(self):
        """Test multiple critical disagreements trigger UNKNOWN."""
        fusion = EvidenceFusion()

        # Create scenario where multiple layers critically disagree
        fusion.add_ml_prediction(0.9, confidence=0.9)
        fusion.add_context_score(0.2, confidence=0.9)
        fusion.add_sequence_score(0.1, confidence=0.9)
        fusion.add_specificity_score(0.95, confidence=0.9)

        result = fusion.fuse()

        # Multiple critical disagreements should make it UNKNOWN
        n_critical = sum(1 for d in result.disagreements if d.is_critical)
        if n_critical >= 2:
            assert result.is_unknown

    def test_no_disagreement_when_layers_agree(self):
        """Test no disagreements when all layers agree."""
        fusion = EvidenceFusion()

        # All layers agree
        fusion.add_ml_prediction(0.8, confidence=0.9)
        fusion.add_context_score(0.75, confidence=0.9)
        fusion.add_sequence_score(0.78, confidence=0.9)
        fusion.add_specificity_score(0.82, confidence=0.9)

        result = fusion.fuse()

        # Should have no critical disagreements
        assert not result.has_critical_disagreement


class TestUnknownBucket:
    """Tests for explicit unknown bucket handling."""

    def test_unknown_no_evidence(self):
        """Test UNKNOWN when no evidence provided."""
        fusion = EvidenceFusion()
        # Add only coherence (doesn't count as component)
        fusion.add_coherence(0.5)

        result = fusion.fuse()

        assert result.is_unknown

    def test_unknown_warning_message(self):
        """Test warning message added for UNKNOWN status."""
        fusion = EvidenceFusion()
        fusion.add_sequence_score(0.5, confidence=0.2)

        result = fusion.fuse()

        assert any("UNKNOWN" in w for w in result.warnings)

    def test_not_unknown_with_good_evidence(self):
        """Test not UNKNOWN with sufficient evidence."""
        fusion = EvidenceFusion()

        fusion.add_ml_prediction(0.8, confidence=0.9)
        fusion.add_context_score(0.75, confidence=0.85)
        fusion.add_coherence(0.5)

        result = fusion.fuse()

        assert not result.is_unknown


class TestFusedResultProperties:
    """Tests for FusedResult properties."""

    def test_n_agreeing_layers(self):
        """Test n_agreeing_layers property."""
        fusion = EvidenceFusion()

        # All three scores agree (within 0.2 of mean)
        fusion.add_ml_prediction(0.8, confidence=0.9)
        fusion.add_context_score(0.75, confidence=0.9)
        fusion.add_sequence_score(0.7, confidence=0.9)

        result = fusion.fuse()

        # All layers should agree (mean ~0.75, all within 0.2 of mean)
        assert result.n_agreeing_layers >= 2

    def test_claim_description_property(self):
        """Test claim_description property returns readable text."""
        fusion = EvidenceFusion()
        fusion.add_sequence_score(0.8, confidence=0.9)

        result = fusion.fuse()

        assert len(result.claim_description) > 20  # Meaningful description
