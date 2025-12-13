"""Tests for v0.9.0 SMS trials module."""

import pytest
import numpy as np

from phaselab.trials.sms import (
    SMSTrial,
    SMSTrialResult,
    SMSTrialConfig,
    SMSPipeline,
    TrialStatus,
    run_sms_crispra_trial,
    run_sms_crispri_trial,
    run_sms_knockout_trial,
    run_sms_base_editing_trial,
    run_sms_prime_editing_trial,
    run_circadian_rescue_simulation,
    run_delivery_assessment,
)
from phaselab.trials.sms.core import (
    TrialType,
    RAI1_INFO,
    SMS_MODIFIER_GENES,
    EXAMPLE_RAI1_PROMOTER,
)


class TestSMSTrialConfig:
    """Tests for SMSTrialConfig."""

    def test_default_config(self):
        """Test default configuration values."""
        config = SMSTrialConfig()

        assert config.therapeutic_window == (0.70, 1.10)
        assert config.optimal_expression == 0.80
        assert config.baseline_expression == 0.50
        assert config.coherence_mode in ["heuristic", "quantum"]

    def test_custom_config(self):
        """Test custom configuration."""
        config = SMSTrialConfig(
            therapeutic_window=(0.60, 1.20),
            optimal_expression=0.90,
            top_n_guides=20,
        )

        assert config.therapeutic_window == (0.60, 1.20)
        assert config.optimal_expression == 0.90
        assert config.top_n_guides == 20


class TestSMSTrialResult:
    """Tests for SMSTrialResult."""

    def test_result_properties(self):
        """Test result properties."""
        result = SMSTrialResult(
            trial_type=TrialType.CRISPRA_RAI1,
            status=TrialStatus.COMPLETED,
            summary="Test summary",
            candidates=[{"sequence": "ATCG" * 5}],
            best_candidate={"sequence": "ATCG" * 5},
        )

        assert result.is_successful
        assert result.has_viable_candidates
        assert result.n_candidates == 1

    def test_failed_result(self):
        """Test failed result properties."""
        result = SMSTrialResult(
            trial_type=TrialType.CRISPRA_RAI1,
            status=TrialStatus.FAILED,
            summary="Failed",
            errors=["Error message"],
        )

        assert not result.is_successful
        assert not result.has_viable_candidates

    def test_to_dict(self):
        """Test result serialization."""
        result = SMSTrialResult(
            trial_type=TrialType.CRISPRA_RAI1,
            status=TrialStatus.COMPLETED,
            summary="Test",
            claim_level="context_dependent",
        )

        d = result.to_dict()

        assert d["trial_type"] == "crispra_rai1"
        assert d["status"] == "completed"
        assert d["claim_level"] == "context_dependent"


class TestCRISPRaTrial:
    """Tests for CRISPRa RAI1 activation trial."""

    def test_crispra_trial_runs(self):
        """Test CRISPRa trial executes."""
        config = SMSTrialConfig(
            verbose=False,
            use_virtual_assay=False,  # Use basic pipeline for speed
            require_go_status=False,
        )

        result = run_sms_crispra_trial(config=config)

        assert result.status in [TrialStatus.COMPLETED, TrialStatus.FAILED]
        assert result.trial_type == TrialType.CRISPRA_RAI1

    def test_crispra_trial_with_sequence(self):
        """Test CRISPRa trial with custom sequence."""
        config = SMSTrialConfig(
            verbose=False,
            use_virtual_assay=False,
            require_go_status=False,
        )

        result = run_sms_crispra_trial(
            promoter_sequence=EXAMPLE_RAI1_PROMOTER,
            config=config,
        )

        assert result.status in [TrialStatus.COMPLETED, TrialStatus.FAILED]
        assert "RAI1" in result.summary

    def test_crispra_trial_metrics(self):
        """Test CRISPRa trial produces expected metrics."""
        config = SMSTrialConfig(
            verbose=False,
            use_virtual_assay=False,
            require_go_status=False,
        )

        result = run_sms_crispra_trial(config=config)

        if result.is_successful:
            assert "therapeutic_window" in result.metrics
            assert "target_gene" in result.metrics
            assert result.metrics["target_gene"] == "RAI1"


class TestCRISPRiTrial:
    """Tests for CRISPRi modifier suppression trial."""

    def test_crispri_trial_per1(self):
        """Test CRISPRi trial for PER1."""
        config = SMSTrialConfig(
            verbose=False,
            use_virtual_assay=False,
            require_go_status=False,
        )

        result = run_sms_crispri_trial(
            target_gene="PER1",
            config=config,
        )

        assert result.status in [TrialStatus.COMPLETED, TrialStatus.FAILED]
        assert result.trial_type == TrialType.CRISPRI_MODIFIER

    def test_crispri_trial_bmal1_warning(self):
        """Test CRISPRi trial warns against BMAL1."""
        config = SMSTrialConfig(
            verbose=False,
            use_virtual_assay=False,
            require_go_status=False,
        )

        result = run_sms_crispri_trial(
            target_gene="BMAL1",
            config=config,
        )

        # Should have warning about BMAL1
        assert any("BMAL1" in w and "NOT recommended" in w for w in result.warnings)


class TestKnockoutTrial:
    """Tests for knockout validation trial."""

    def test_knockout_trial_runs(self):
        """Test knockout trial executes."""
        config = SMSTrialConfig(
            verbose=False,
            require_go_status=False,
        )

        result = run_sms_knockout_trial(
            target_gene="PER1",
            config=config,
        )

        assert result.status in [TrialStatus.COMPLETED, TrialStatus.FAILED]
        assert result.trial_type == TrialType.KNOCKOUT_VALIDATION

    def test_knockout_trial_research_warning(self):
        """Test knockout trial includes research-only warning."""
        config = SMSTrialConfig(verbose=False, require_go_status=False)

        result = run_sms_knockout_trial(config=config)

        # Should warn this is research use only
        assert any("RESEARCH USE" in w for w in result.warnings)


class TestBaseEditingTrial:
    """Tests for base editing trial."""

    def test_base_editing_trial_runs(self):
        """Test base editing trial executes."""
        config = SMSTrialConfig(
            verbose=False,
            require_go_status=False,
        )

        result = run_sms_base_editing_trial(
            variant_id="p.R1217Q",
            config=config,
        )

        assert result.status in [TrialStatus.COMPLETED, TrialStatus.FAILED]
        assert result.trial_type == TrialType.BASE_EDITING

    def test_base_editing_applicability_warning(self):
        """Test base editing includes applicability warning."""
        config = SMSTrialConfig(verbose=False, require_go_status=False)

        result = run_sms_base_editing_trial(config=config)

        # Should warn about limited applicability
        assert any("only applicable" in w.lower() for w in result.warnings)


class TestPrimeEditingTrial:
    """Tests for prime editing trial."""

    def test_prime_editing_trial_runs(self):
        """Test prime editing trial executes."""
        config = SMSTrialConfig(
            verbose=False,
            require_go_status=False,
        )

        result = run_sms_prime_editing_trial(
            edit_target="codon_correction",
            config=config,
        )

        assert result.status in [TrialStatus.COMPLETED, TrialStatus.FAILED]
        assert result.trial_type == TrialType.PRIME_EDITING

    def test_prime_editing_delivery_warning(self):
        """Test prime editing includes delivery warning."""
        config = SMSTrialConfig(verbose=False, require_go_status=False)

        result = run_sms_prime_editing_trial(config=config)

        # Should warn about payload size
        assert any("dual-AAV" in w.lower() or "payload" in w.lower() for w in result.warnings)


class TestCircadianSimulation:
    """Tests for circadian rescue simulation."""

    def test_circadian_simulation_runs(self):
        """Test circadian simulation executes."""
        config = SMSTrialConfig(
            verbose=False,
            simulation_hours=48.0,  # Shorter for testing
            n_circadian_trials=2,
        )

        result = run_circadian_rescue_simulation(
            predicted_rai1_expression=0.80,
            config=config,
        )

        assert result.status == TrialStatus.COMPLETED
        assert result.trial_type == TrialType.CIRCADIAN_SIMULATION

    def test_circadian_simulation_metrics(self):
        """Test circadian simulation produces expected metrics."""
        config = SMSTrialConfig(
            verbose=False,
            simulation_hours=48.0,
            n_circadian_trials=2,
        )

        result = run_circadian_rescue_simulation(
            predicted_rai1_expression=0.80,
            config=config,
        )

        assert "baseline_R_bar" in result.metrics
        assert "final_R_bar" in result.metrics
        assert "rescue_status" in result.metrics
        assert "sleep_quality_prediction" in result.metrics

    def test_circadian_rescue_at_normal_expression(self):
        """Test that normal RAI1 level predicts rescue."""
        config = SMSTrialConfig(
            verbose=False,
            simulation_hours=120.0,
            n_circadian_trials=3,
        )

        result = run_circadian_rescue_simulation(
            predicted_rai1_expression=1.0,  # Normal level
            config=config,
        )

        # Should predict rescue at normal expression
        rescue_status = result.metrics.get("rescue_status", "")
        assert rescue_status in ["FULL_RESCUE", "PARTIAL_RESCUE"]


class TestDeliveryAssessment:
    """Tests for delivery feasibility assessment."""

    def test_delivery_assessment_runs(self):
        """Test delivery assessment executes."""
        config = SMSTrialConfig(verbose=False)

        result = run_delivery_assessment(
            modality="CRISPRa_VP64",
            target_tissue="brain",
            config=config,
        )

        assert result.status == TrialStatus.COMPLETED
        assert result.trial_type == TrialType.DELIVERY_ASSESSMENT

    def test_delivery_assessment_metrics(self):
        """Test delivery assessment produces expected metrics."""
        config = SMSTrialConfig(verbose=False)

        result = run_delivery_assessment(config=config)

        assert "payload_size" in result.metrics
        assert "delivery_feasibility" in result.metrics
        assert "packaging_status" in result.metrics

    def test_delivery_assessment_oversized_payload(self):
        """Test delivery assessment handles oversized payloads."""
        config = SMSTrialConfig(verbose=False)

        result = run_delivery_assessment(
            modality="prime_editing",  # Very large payload
            config=config,
        )

        # Should identify as not feasible or challenging
        feasibility = result.metrics.get("delivery_feasibility", "")
        assert feasibility in ["NOT_FEASIBLE", "CHALLENGING"]


class TestSMSPipeline:
    """Tests for SMS pipeline orchestrator."""

    def test_pipeline_initialization(self):
        """Test pipeline can be initialized."""
        pipeline = SMSPipeline()

        assert pipeline.config is not None

    def test_pipeline_crispra_only(self):
        """Test pipeline CRISPRa-only mode."""
        config = SMSTrialConfig(
            verbose=False,
            use_virtual_assay=False,
            require_go_status=False,
        )
        pipeline = SMSPipeline(config=config)

        result = pipeline.run_crispra_only()

        assert result.trial_type == TrialType.CRISPRA_RAI1

    def test_pipeline_circadian_only(self):
        """Test pipeline circadian-only mode."""
        config = SMSTrialConfig(
            verbose=False,
            simulation_hours=48.0,
            n_circadian_trials=2,
        )
        pipeline = SMSPipeline(config=config)

        result = pipeline.run_circadian_only(predicted_expression=0.80)

        assert result.trial_type == TrialType.CIRCADIAN_SIMULATION

    def test_pipeline_delivery_only(self):
        """Test pipeline delivery-only mode."""
        config = SMSTrialConfig(verbose=False)
        pipeline = SMSPipeline(config=config)

        result = pipeline.run_delivery_only()

        assert result.trial_type == TrialType.DELIVERY_ASSESSMENT

    def test_full_pipeline_basic(self):
        """Test full pipeline execution (basic mode)."""
        config = SMSTrialConfig(
            verbose=False,
            use_virtual_assay=False,
            require_go_status=False,
            simulation_hours=48.0,
            n_circadian_trials=2,
        )
        pipeline = SMSPipeline(config=config)

        result = pipeline.run_full_pipeline(
            include_modifiers=False,
            include_editing=False,
        )

        assert result.status == "complete"
        assert result.crispra_result is not None
        assert result.circadian_result is not None
        assert result.delivery_result is not None
        assert result.overall_go_nogo in ["GO", "CONDITIONAL", "NO-GO", "UNKNOWN"]

    def test_pipeline_produces_falsification_tests(self):
        """Test pipeline generates falsification tests."""
        config = SMSTrialConfig(
            verbose=False,
            use_virtual_assay=False,
            require_go_status=False,
            simulation_hours=48.0,
            n_circadian_trials=2,
        )
        pipeline = SMSPipeline(config=config)

        result = pipeline.run_full_pipeline()

        # Should have falsification tests
        assert len(result.falsification_tests) > 0

        # Should have test A (ranking validity)
        test_ids = [t["id"] for t in result.falsification_tests]
        assert "A" in test_ids

    def test_pipeline_produces_validation_priorities(self):
        """Test pipeline generates validation priorities."""
        config = SMSTrialConfig(
            verbose=False,
            use_virtual_assay=False,
            require_go_status=False,
            simulation_hours=48.0,
            n_circadian_trials=2,
        )
        pipeline = SMSPipeline(config=config)

        result = pipeline.run_full_pipeline()

        # Should have validation priorities
        assert len(result.validation_priorities) > 0

    def test_pipeline_produces_wet_lab_recommendations(self):
        """Test pipeline generates wet lab recommendations."""
        config = SMSTrialConfig(
            verbose=False,
            use_virtual_assay=False,
            require_go_status=False,
            simulation_hours=48.0,
            n_circadian_trials=2,
        )
        pipeline = SMSPipeline(config=config)

        result = pipeline.run_full_pipeline()

        # Should have wet lab recommendations
        assert len(result.wet_lab_recommendations) > 0

    def test_pipeline_result_to_dict(self):
        """Test pipeline result serialization."""
        config = SMSTrialConfig(
            verbose=False,
            use_virtual_assay=False,
            require_go_status=False,
            simulation_hours=48.0,
            n_circadian_trials=2,
        )
        pipeline = SMSPipeline(config=config)

        result = pipeline.run_full_pipeline()
        d = result.to_dict()

        assert "status" in d
        assert "overall_go_nogo" in d
        assert "overall_claim_level" in d
        assert "falsification_tests" in d


class TestCoreData:
    """Tests for core SMS data structures."""

    def test_rai1_info(self):
        """Test RAI1 info is complete."""
        assert "gene_symbol" in RAI1_INFO
        assert RAI1_INFO["gene_symbol"] == "RAI1"
        assert "tss" in RAI1_INFO
        assert "crispra_window" in RAI1_INFO

    def test_modifier_genes(self):
        """Test modifier genes data."""
        assert "PER1" in SMS_MODIFIER_GENES
        assert "PER2" in SMS_MODIFIER_GENES
        assert "CRY1" in SMS_MODIFIER_GENES

        for gene, info in SMS_MODIFIER_GENES.items():
            assert "role" in info
            assert "rationale" in info
            assert "caution" in info

    def test_example_sequence(self):
        """Test example promoter sequence is valid."""
        assert len(EXAMPLE_RAI1_PROMOTER) > 100
        assert all(c in "ATCGN" for c in EXAMPLE_RAI1_PROMOTER)
