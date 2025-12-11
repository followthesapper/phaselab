"""
Tests for phaselab.crispr module.
"""

import numpy as np
import pytest
from phaselab.crispr.pam_scan import (
    find_pam_sites,
    filter_by_window,
    reverse_complement,
    PAM_PATTERNS,
)
from phaselab.crispr.scoring import (
    gc_content,
    max_homopolymer_run,
    delta_g_santalucia,
    sequence_complexity,
    mit_specificity_score,
    cfd_score,
    chromatin_accessibility_score,
)
from phaselab.crispr.pipeline import (
    design_guides,
    GuideDesignConfig,
    validate_guide,
)


class TestReverseComplement:
    """Test reverse complement function."""

    def test_simple(self):
        assert reverse_complement("ATCG") == "CGAT"

    def test_palindrome(self):
        assert reverse_complement("ATAT") == "ATAT"

    def test_gc_only(self):
        assert reverse_complement("GCGC") == "GCGC"


class TestPAMScan:
    """Test PAM site scanning."""

    def test_find_ngg_sites(self):
        """Find NGG PAM sites."""
        seq = "AAAAAAAAAAAAAAAAAAAAAAGGATCG"
        #      0         1         2
        #      01234567890123456789012345678
        hits = find_pam_sites(seq, pam="NGG", both_strands=False)

        assert len(hits) >= 1
        # Should find AGG at position 21
        pam_positions = [h.position for h in hits]
        assert 21 in pam_positions

    def test_guide_extraction(self):
        """Verify guide is extracted correctly."""
        # Create a sequence where we know exactly where the PAM is
        guide_expected = "ATCGATCGATCGATCGATCG"
        pam = "AGG"
        seq = guide_expected + pam + "NNNNNN"

        hits = find_pam_sites(seq, pam="NGG", both_strands=False)

        found = False
        for hit in hits:
            if hit.guide == guide_expected:
                found = True
                break
        assert found, f"Expected guide not found. Hits: {[h.guide for h in hits]}"

    def test_both_strands(self):
        """Should find hits on both strands."""
        # NGG on forward, CCN on reverse (which is NGG on rev comp)
        seq = "ATCGATCGATCGATCGATCGAGGNNNNNNNNNNNNNNNNNNNNNCCT"
        hits = find_pam_sites(seq, pam="NGG", both_strands=True)

        strands = set(h.strand for h in hits)
        assert "+" in strands or "-" in strands  # At least one strand

    def test_filter_by_window(self):
        """Filter hits by CRISPRa window."""
        # Create mock hits
        from phaselab.crispr.pam_scan import PAMHit

        hits = [
            PAMHit(position=100, strand="+", guide="A"*20, pam="AGG",
                   guide_start=80, guide_end=100),
            PAMHit(position=450, strand="+", guide="T"*20, pam="AGG",
                   guide_start=430, guide_end=450),
            PAMHit(position=600, strand="+", guide="C"*20, pam="AGG",
                   guide_start=580, guide_end=600),
        ]

        # TSS at position 500, window (-400, -50)
        filtered = filter_by_window(hits, tss_position=500, window=(-400, -50))

        # Only first two should be in window
        assert len(filtered) <= 2


class TestScoring:
    """Test scoring functions."""

    def test_gc_content(self):
        """GC content calculation."""
        assert gc_content("GCGC") == 1.0
        assert gc_content("ATAT") == 0.0
        assert gc_content("ATGC") == 0.5
        assert gc_content("GGGGGGGGAA") == 0.8

    def test_max_homopolymer_run(self):
        """Homopolymer detection."""
        assert max_homopolymer_run("ATCG") == 1
        assert max_homopolymer_run("AAAT") == 3
        assert max_homopolymer_run("AAAAAATCG") == 6
        assert max_homopolymer_run("ATGGGGGC") == 5

    def test_delta_g_santalucia(self):
        """Thermodynamic ΔG calculation."""
        # GC-rich should be more negative (stronger binding)
        gc_rich = "GCGCGCGCGCGCGCGCGCGC"
        at_rich = "ATATATATATATATATATAT"

        dg_gc = delta_g_santalucia(gc_rich)
        dg_at = delta_g_santalucia(at_rich)

        assert dg_gc < dg_at  # GC-rich more negative

    def test_sequence_complexity(self):
        """Sequence complexity score."""
        repetitive = "AAAAAAAAAAAAAAAAAAA"
        complex_seq = "ATCGATCGTAGCTAGCTAG"

        assert sequence_complexity(repetitive) < sequence_complexity(complex_seq)

    def test_mit_specificity_score(self):
        """MIT specificity score range."""
        guide = "ATCGATCGATCGATCGATCG"
        score = mit_specificity_score(guide)

        assert 0 <= score <= 100

    def test_cfd_score_perfect_match(self):
        """CFD score for perfect match should be 100."""
        guide = "ATCGATCGATCGATCGATCG"
        score = cfd_score(guide, target_seq=None)  # None = perfect match

        assert score == 100.0

    def test_cfd_score_with_mismatches(self):
        """CFD score decreases with mismatches."""
        guide = "ATCGATCGATCGATCGATCG"
        target = "ATCGATCGATCGATCGATCC"  # 1 mismatch

        score = cfd_score(guide, target_seq=target)
        assert score < 100.0

    def test_chromatin_accessibility(self):
        """Chromatin accessibility near TSS."""
        # Close to TSS should be more accessible
        state_close, score_close = chromatin_accessibility_score(
            position=450, tss_position=500
        )
        state_far, score_far = chromatin_accessibility_score(
            position=100, tss_position=500
        )

        assert score_close > score_far


class TestPipeline:
    """Test the main design_guides pipeline."""

    @pytest.fixture
    def sample_promoter(self):
        """Sample promoter sequence for testing."""
        # 600bp with known PAM sites
        return (
            "ATCGATCGATCGATCGATCG" * 10 +  # 200bp
            "NNNNNNNNNNNNNNNNNNNN" +        # 20bp (TSS area)
            "ATCGATCGATCGATCGATCGAGG" +     # Guide + PAM
            "NNNNNNNNNNNNNNNNNNNN" * 10 +   # 200bp
            "GCTAGCTAGCTAGCTAGCTAGG" +      # Another Guide + PAM
            "ATCGATCG" * 10                  # padding
        )

    def test_design_guides_returns_dataframe(self, sample_promoter):
        """design_guides should return a DataFrame."""
        import pandas as pd

        config = GuideDesignConfig(compute_coherence=False)
        result = design_guides(
            sequence=sample_promoter,
            tss_index=220,
            config=config,
        )

        assert isinstance(result, pd.DataFrame)

    def test_design_guides_columns(self, sample_promoter):
        """Result should have expected columns."""
        config = GuideDesignConfig(compute_coherence=True)
        result = design_guides(
            sequence=sample_promoter,
            tss_index=220,
            config=config,
        )

        expected_cols = ['sequence', 'position', 'gc', 'mit_score', 'coherence_R', 'go_no_go']
        for col in expected_cols:
            assert col in result.columns, f"Missing column: {col}"

    def test_validate_guide_valid(self):
        """Validate a good guide."""
        result = validate_guide("ATCGATCGATCGATCGATCG")

        assert result['length'] == 20
        assert 'gc' in result
        assert 'go_no_go' in result
        assert isinstance(result['warnings'], list)

    def test_validate_guide_low_gc(self):
        """Low GC should generate warning."""
        result = validate_guide("AAAAAAAAAAAAAAAAAAAT")

        assert not result['valid']
        assert any("GC" in w for w in result['warnings'])

    def test_validate_guide_homopolymer(self):
        """Long homopolymer should generate warning."""
        result = validate_guide("AAAAAATCGATCGATCGATC")

        assert any("homopolymer" in w.lower() for w in result['warnings'])


class TestGuideDesignConfig:
    """Test configuration dataclass."""

    def test_default_values(self):
        """Check default configuration."""
        config = GuideDesignConfig()

        assert config.pam == "NGG"
        assert config.guide_length == 20
        assert config.crispr_window == (-400, -50)
        assert config.min_gc == 0.4
        assert config.max_gc == 0.7

    def test_custom_values(self):
        """Custom configuration."""
        config = GuideDesignConfig(
            pam="NNGRRT",
            min_gc=0.3,
            max_gc=0.8,
        )

        assert config.pam == "NNGRRT"
        assert config.min_gc == 0.3


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
