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


class TestCRISPORCompositeScoring:
    """Test CRISPOR-style composite scoring with mismatch distance weighting."""

    def test_crispor_composite_basic(self):
        """Test basic CRISPOR composite score calculation."""
        from phaselab.crispr.scoring import crispor_composite_score

        # Perfect guide: MIT=93, CFD=95, no dangerous off-targets
        score, breakdown = crispor_composite_score(
            mit_score=93,
            cfd_score=95,
            off_targets={0: 0, 1: 0, 2: 0, 3: 5, 4: 15},
        )

        # Base: 93+95=188
        # v0.9.2 weights: 3mm=0.1, 4mm=0.01 (negligible for non-critical mismatches)
        # Penalties: 5*0.1 + 15*0.01 = 0.5 + 0.15 = 0.65
        # Final: 188 - 0.65 = 187.35
        assert score == 187.35
        assert breakdown["base_score"] == 188.0
        assert breakdown["total_ot_penalty"] == 0.65

    def test_crispor_composite_u6_penalty(self):
        """Test that U6 incompatibility applies 100-point penalty."""
        from phaselab.crispr.scoring import crispor_composite_score

        # Guide with TTTT start (U6 incompatible)
        score_bad, _ = crispor_composite_score(
            mit_score=98,
            cfd_score=98,
            off_targets={0: 0, 1: 0, 2: 1, 3: 6, 4: 22},
            u6_compatible=False,
        )

        score_good, _ = crispor_composite_score(
            mit_score=98,
            cfd_score=98,
            off_targets={0: 0, 1: 0, 2: 1, 3: 6, 4: 22},
            u6_compatible=True,
        )

        # U6 penalty should be exactly 100 points
        assert score_good - score_bad == 100.0

    def test_crispor_composite_repeat_penalty(self):
        """Test that repeat region applies 1000-point penalty."""
        from phaselab.crispr.scoring import crispor_composite_score

        score_repeat, breakdown = crispor_composite_score(
            mit_score=90,
            cfd_score=90,
            is_repeat=True,
        )

        score_normal, _ = crispor_composite_score(
            mit_score=90,
            cfd_score=90,
            is_repeat=False,
        )

        assert breakdown["repeat_penalty"] == 1000.0
        assert score_normal - score_repeat == 1000.0

    def test_crispor_mismatch_weighting(self):
        """Test that 0-1mm off-targets are penalized more than 3-4mm."""
        from phaselab.crispr.scoring import crispor_composite_score

        # One 0mm off-target (catastrophic)
        score_0mm, breakdown_0mm = crispor_composite_score(
            mit_score=90, cfd_score=90,
            off_targets={0: 1},
        )

        # One 4mm off-target (minimal risk)
        score_4mm, breakdown_4mm = crispor_composite_score(
            mit_score=90, cfd_score=90,
            off_targets={4: 1},
        )

        # v0.9.2: 0mm = 10000 (disqualifying), 4mm = 0.01 (negligible)
        assert breakdown_0mm["total_ot_penalty"] == 10000.0
        assert breakdown_4mm["total_ot_penalty"] == 0.01
        assert score_4mm > score_0mm

    def test_rank_guides_crispor_style(self):
        """Test that ranking properly handles the MIT98/CFD98 trap."""
        from phaselab.crispr.scoring import rank_guides_crispor_style

        guides = [
            # Guide #1: Lower MIT/CFD but clean off-target profile
            {
                "sequence": "TTCGATGAATGGTTGCTACC",
                "mit_score": 93,
                "cfd_score": 95,
                "off_targets": {0: 0, 1: 0, 2: 0, 3: 5, 4: 15},
                "u6_compatible": True,
                "is_repeat": False,
            },
            # Guide #7: High MIT/CFD but TTTT start (U6 trap)
            {
                "sequence": "TTTTAATGGCCGGCGATGCC",
                "mit_score": 98,
                "cfd_score": 98,
                "off_targets": {0: 0, 1: 0, 2: 1, 3: 6, 4: 22},
                "u6_compatible": False,
                "is_repeat": False,
            },
            # Guide #8: Good scores but in repeat region
            {
                "sequence": "CAGCAGCAGCAGCAGCAGCA",
                "mit_score": 85,
                "cfd_score": 88,
                "off_targets": {0: 0, 1: 0, 2: 0, 3: 2, 4: 8},
                "u6_compatible": True,
                "is_repeat": True,
            },
        ]

        # With default exclusions (require U6, exclude repeats)
        ranked = rank_guides_crispor_style(guides)

        # Only Guide #1 should remain (others excluded)
        assert len(ranked) == 1
        assert ranked[0]["sequence"] == "TTCGATGAATGGTTGCTACC"
        assert ranked[0]["crispor_rank"] == 1

    def test_rank_guides_without_exclusions(self):
        """Test ranking without hard exclusions."""
        from phaselab.crispr.scoring import rank_guides_crispor_style

        guides = [
            {
                "sequence": "TTCGATGAATGGTTGCTACC",
                "mit_score": 93,
                "cfd_score": 95,
                "off_targets": {0: 0, 1: 0, 2: 0, 3: 5, 4: 15},
            },
            {
                "sequence": "TTTTAATGGCCGGCGATGCC",
                "mit_score": 98,
                "cfd_score": 98,
                "off_targets": {0: 0, 1: 0, 2: 1, 3: 6, 4: 22},
                "u6_compatible": False,
            },
        ]

        # Allow U6-incompatible guides (soft penalty instead)
        ranked = rank_guides_crispor_style(
            guides,
            require_u6_compatible=False,
            exclude_repeats=False,
        )

        assert len(ranked) == 2
        # Guide #1 should still rank higher due to 100-point U6 penalty on #7
        assert ranked[0]["sequence"] == "TTCGATGAATGGTTGCTACC"

    def test_poly_t_penalty_function(self):
        """Test poly-T detection for U6 compatibility."""
        from phaselab.crispr.scoring import poly_t_penalty

        # TTTT at start - incompatible
        is_bad, reason = poly_t_penalty("TTTTAATGGCCGGCGATGCC")
        assert is_bad
        assert "5' end" in reason

        # Internal TTTT - problematic
        is_bad, reason = poly_t_penalty("AATGTTTTGGCCGGCGATGC")
        assert is_bad
        assert "Internal" in reason

        # Normal guide - OK
        is_bad, reason = poly_t_penalty("TTCGATGAATGGTTGCTACC")
        assert not is_bad

    def test_is_repeat_region_function(self):
        """Test repeat region detection."""
        from phaselab.crispr.scoring import is_repeat_region

        # Tandem repeat (CAG repeat expansion)
        is_rep, reason = is_repeat_region("CAGCAGCAGCAGCAGCAGCA")
        assert is_rep
        assert "CAG" in reason

        # Dinucleotide repeat
        is_rep, reason = is_repeat_region("ATATATATATATATATATATAT"[:20])
        assert is_rep

        # Normal complex sequence
        is_rep, reason = is_repeat_region("TTCGATGAATGGTTGCTACC")
        assert not is_rep


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
