"""Tests for the enhanced CRISPR pipeline with Virtual Assay Stack."""

import pytest
import numpy as np

from phaselab.crispr import (
    Modality,
    EnhancedGuideConfig,
    EnhancedGuide,
    EnhancedDesignResult,
    design_enhanced_guides,
    compare_guides_with_without_context,
)
from phaselab.crispr.coherence_utils import CoherenceMode
from phaselab.context import CellType
from phaselab.fusion import FusionConfig
from phaselab.core.constants import E_MINUS_2


# Test promoter sequence (synthetic, ~1000bp around TSS)
TEST_PROMOTER = (
    "ATGCGATCGATCGATCGATCG" * 10 +  # Upstream region
    "GGGGGGGGGGGGGGGGGGGG" +  # GC-rich region
    "ATATATATATATATATATAT" +  # AT-rich region
    "CGATCGATCGATCGATCGAT" +  # Mixed region with PAM sites
    "AGGTCGATCGATCGATCGAT" +  # NGG PAM site
    "TGGTCGATCGATCGATCGAT" +  # Another NGG PAM
    "CGGTCGATCGATCGATCGAT" +  # Another NGG PAM
    "ATGCGATCGATCGATCGATCG" * 10 +  # More upstream
    "ATGCATGCATGCATGCATGC" * 5 +  # TSS region
    "GCGATCGATCGATCGATCGA" * 10  # Downstream
)
TEST_TSS = 400  # TSS position in sequence


class TestEnhancedGuideConfig:
    """Tests for EnhancedGuideConfig."""

    def test_default_config(self):
        """Test default configuration."""
        config = EnhancedGuideConfig()

        assert config.pam == "NGG"
        assert config.guide_length == 20
        assert config.modality == Modality.CRISPRA
        assert config.use_biological_context == True
        assert config.use_ml_predictors == True

    def test_window_by_modality(self):
        """Test window selection by modality."""
        config_a = EnhancedGuideConfig(modality=Modality.CRISPRA)
        assert config_a.window == (-400, -50)

        config_i = EnhancedGuideConfig(modality=Modality.CRISPRI)
        assert config_i.window == (-50, 300)

        config_ko = EnhancedGuideConfig(modality=Modality.KNOCKOUT)
        assert config_ko.window == (0, 500)

    def test_coherence_mode(self):
        """Test coherence mode configuration."""
        config = EnhancedGuideConfig(coherence_mode=CoherenceMode.QUANTUM)
        assert config.coherence_mode == CoherenceMode.QUANTUM

    def test_custom_fusion_config(self):
        """Test custom fusion configuration."""
        fusion = FusionConfig.conservative()
        config = EnhancedGuideConfig(fusion_config=fusion)
        assert config.fusion_config.specificity_weight > 0.3


class TestEnhancedGuide:
    """Tests for EnhancedGuide dataclass."""

    def test_guide_creation(self):
        """Test creating an enhanced guide."""
        guide = EnhancedGuide(
            sequence="ATGCGATCGATCGATCGATC",
            pam="NGG",
            position=-200,
            strand="+",
            gc=0.55,
            coherence=0.45,
            fused_score=0.75,
            go_status="GO",
            passes_gates=True,
        )

        assert guide.sequence == "ATGCGATCGATCGATCGATC"
        assert guide.is_go
        assert guide.is_viable

    def test_guide_not_viable_failed_gates(self):
        """Test guide not viable when gates fail."""
        guide = EnhancedGuide(
            sequence="ATGCGATCGATCGATCGATC",
            pam="NGG",
            position=-200,
            strand="+",
            go_status="GO",
            passes_gates=False,
            failed_gates=["gc_content"],
        )

        assert guide.is_go
        assert not guide.is_viable

    def test_guide_not_viable_nogo(self):
        """Test guide not viable when NO-GO."""
        guide = EnhancedGuide(
            sequence="ATGCGATCGATCGATCGATC",
            pam="NGG",
            position=-200,
            strand="+",
            go_status="NO-GO",
            passes_gates=True,
        )

        assert not guide.is_go
        assert not guide.is_viable

    def test_guide_to_dict(self):
        """Test guide serialization."""
        guide = EnhancedGuide(
            sequence="ATGCGATCGATCGATCGATC",
            pam="NGG",
            position=-200,
            strand="+",
            fused_score=0.75,
        )

        d = guide.to_dict()
        assert "sequence" in d
        assert "fused_score" in d
        assert "is_viable" in d


class TestEnhancedDesignResult:
    """Tests for EnhancedDesignResult."""

    def test_result_creation(self):
        """Test creating a design result."""
        config = EnhancedGuideConfig()
        guides = [
            EnhancedGuide(
                sequence="ATGCGATCGATCGATCGATC",
                pam="NGG",
                position=-200,
                strand="+",
                go_status="GO",
                passes_gates=True,
                fused_score=0.8,
            ),
            EnhancedGuide(
                sequence="GCATGCATGCATGCATGCAT",
                pam="NGG",
                position=-150,
                strand="+",
                go_status="NO-GO",
                passes_gates=True,
                fused_score=0.6,
            ),
        ]

        result = EnhancedDesignResult(
            guides=guides,
            config=config,
            target_gene="TEST",
            total_candidates=10,
        )

        assert len(result.guides) == 2
        assert len(result.go_guides) == 1
        assert len(result.viable_guides) == 1

    def test_result_to_dataframe(self):
        """Test converting result to DataFrame."""
        config = EnhancedGuideConfig()
        guides = [
            EnhancedGuide(
                sequence="ATGCGATCGATCGATCGATC",
                pam="NGG",
                position=-200,
                strand="+",
            ),
        ]

        result = EnhancedDesignResult(guides=guides, config=config)
        df = result.to_dataframe()

        assert len(df) == 1
        assert "sequence" in df.columns
        assert "fused_score" in df.columns

    def test_result_summary(self):
        """Test result summary."""
        config = EnhancedGuideConfig()
        result = EnhancedDesignResult(
            guides=[],
            config=config,
            total_candidates=100,
            filtered_by_gates=30,
            filtered_by_coherence=20,
        )

        summary = result.summary()
        assert summary["total_candidates"] == 100
        assert summary["filtered_by_gates"] == 30
        assert summary["filtered_by_coherence"] == 20


class TestDesignEnhancedGuides:
    """Tests for the enhanced guide design function."""

    def test_basic_design(self):
        """Test basic guide design without context."""
        config = EnhancedGuideConfig(
            use_biological_context=False,
            use_ml_predictors=False,
            top_n=5,
        )

        result = design_enhanced_guides(
            sequence=TEST_PROMOTER,
            tss_index=TEST_TSS,
            config=config,
        )

        assert isinstance(result, EnhancedDesignResult)
        # May find guides depending on PAM sites in test sequence

    def test_design_with_ml_predictors(self):
        """Test design with ML predictors enabled."""
        config = EnhancedGuideConfig(
            use_biological_context=False,
            use_ml_predictors=True,
            top_n=5,
        )

        result = design_enhanced_guides(
            sequence=TEST_PROMOTER,
            tss_index=TEST_TSS,
            config=config,
            target_gene="TEST",
        )

        # Should run without error
        assert isinstance(result, EnhancedDesignResult)

    def test_design_with_context(self):
        """Test design with biological context."""
        config = EnhancedGuideConfig(
            use_biological_context=True,
            use_ml_predictors=False,
            cell_type=CellType.K562,
            top_n=5,
        )

        result = design_enhanced_guides(
            sequence=TEST_PROMOTER,
            tss_index=TEST_TSS,
            config=config,
            chrom="chr1",
            chrom_offset=1000000,
        )

        # Should run without error (context may not have data)
        assert isinstance(result, EnhancedDesignResult)

    def test_design_full_stack(self):
        """Test design with full Virtual Assay Stack."""
        config = EnhancedGuideConfig(
            use_biological_context=True,
            use_ml_predictors=True,
            cell_type=CellType.K562,
            coherence_mode=CoherenceMode.HEURISTIC,
            top_n=10,
        )

        result = design_enhanced_guides(
            sequence=TEST_PROMOTER,
            tss_index=TEST_TSS,
            config=config,
            target_gene="SCN2A",
            chrom="chr2",
            chrom_offset=165000000,
        )

        assert isinstance(result, EnhancedDesignResult)

    def test_design_different_modalities(self):
        """Test design for different CRISPR modalities."""
        for modality in [Modality.CRISPRA, Modality.CRISPRI, Modality.KNOCKOUT]:
            config = EnhancedGuideConfig(
                modality=modality,
                use_biological_context=False,
                use_ml_predictors=False,
                top_n=3,
            )

            result = design_enhanced_guides(
                sequence=TEST_PROMOTER,
                tss_index=TEST_TSS,
                config=config,
            )

            assert result.config.modality == modality

    def test_include_nogo_guides(self):
        """Test including NO-GO guides in output."""
        config_exclude = EnhancedGuideConfig(
            use_biological_context=False,
            use_ml_predictors=False,
            include_nogo=False,
            top_n=20,
        )

        config_include = EnhancedGuideConfig(
            use_biological_context=False,
            use_ml_predictors=False,
            include_nogo=True,
            top_n=20,
        )

        result_exclude = design_enhanced_guides(
            sequence=TEST_PROMOTER,
            tss_index=TEST_TSS,
            config=config_exclude,
        )

        result_include = design_enhanced_guides(
            sequence=TEST_PROMOTER,
            tss_index=TEST_TSS,
            config=config_include,
        )

        # Including NO-GO should give >= results
        assert len(result_include.guides) >= len(result_exclude.guides)

    def test_empty_sequence(self):
        """Test with empty sequence."""
        config = EnhancedGuideConfig()

        result = design_enhanced_guides(
            sequence="",
            tss_index=0,
            config=config,
        )

        assert len(result.guides) == 0

    def test_no_pam_sites(self):
        """Test sequence with no PAM sites."""
        no_pam_seq = "ATATATAT" * 100  # No NGG

        config = EnhancedGuideConfig()

        result = design_enhanced_guides(
            sequence=no_pam_seq,
            tss_index=400,
            config=config,
        )

        assert len(result.guides) == 0


class TestGuideScoring:
    """Tests for guide scoring components."""

    def test_evidence_levels(self):
        """Test evidence level assignment."""
        # Heuristic mode should give level C
        config_h = EnhancedGuideConfig(
            coherence_mode=CoherenceMode.HEURISTIC,
            use_biological_context=False,
            use_ml_predictors=False,
        )

        result_h = design_enhanced_guides(
            sequence=TEST_PROMOTER,
            tss_index=TEST_TSS,
            config=config_h,
        )

        for guide in result_h.guides:
            assert guide.evidence_level == "C"

    def test_go_threshold(self):
        """Test GO/NO-GO threshold application."""
        config = EnhancedGuideConfig(
            go_threshold=E_MINUS_2,
            use_biological_context=False,
            use_ml_predictors=False,
            include_nogo=True,
        )

        result = design_enhanced_guides(
            sequence=TEST_PROMOTER,
            tss_index=TEST_TSS,
            config=config,
        )

        for guide in result.guides:
            if guide.coherence >= E_MINUS_2:
                assert guide.go_status == "GO"
            else:
                assert guide.go_status == "NO-GO"


class TestCompareGuidesWithWithoutContext:
    """Tests for the comparison function."""

    def test_comparison_runs(self):
        """Test that comparison function runs."""
        result = compare_guides_with_without_context(
            sequence=TEST_PROMOTER,
            tss_index=TEST_TSS,
            chrom="chr1",
            chrom_offset=1000000,
            cell_type=CellType.K562,
        )

        assert "basic" in result
        assert "enhanced" in result
        assert "rank_changes" in result
        assert "basic_summary" in result
        assert "enhanced_summary" in result

    def test_comparison_with_target_gene(self):
        """Test comparison with target gene."""
        result = compare_guides_with_without_context(
            sequence=TEST_PROMOTER,
            tss_index=TEST_TSS,
            chrom="chr2",
            chrom_offset=165000000,
            cell_type=CellType.K562,
            target_gene="SCN2A",
        )

        assert "basic" in result
        assert "enhanced" in result
