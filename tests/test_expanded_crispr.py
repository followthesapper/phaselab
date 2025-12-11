"""
Tests for expanded CRISPR modules: knockout, CRISPRi, prime editing, base editing.
"""

import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from phaselab.crispr.knockout import (
    cut_efficiency_score,
    frameshift_probability,
    repair_pathway_prediction,
    KnockoutConfig,
)

from phaselab.crispr.interference import (
    repression_efficiency_score,
    steric_hindrance_score,
    CRISPRiConfig,
)

from phaselab.crispr.prime_editing import (
    pbs_score,
    rt_template_score,
    reverse_complement,
    estimate_hairpin_dg,
    PrimeEditConfig,
)

from phaselab.crispr.base_editing import (
    editing_efficiency_at_position,
    sequence_context_score,
    find_bystanders,
    get_activity_window,
    BaseEditConfig,
)


class TestKnockoutModule:
    """Tests for CRISPR knockout module."""

    def test_cut_efficiency_basic(self):
        """Test basic cut efficiency calculation."""
        # Optimal guide should have high efficiency
        optimal = "GCGACTGCTACATAGCCAGG"
        eff = cut_efficiency_score(optimal)
        assert 0.5 <= eff <= 1.0

    def test_cut_efficiency_polyT_penalty(self):
        """Test that poly-T sequences are penalized."""
        with_polyT = "TTTTGACTGCTACATAGCCA"
        without_polyT = "GCGACTGCTACATAGCCAGG"
        eff_with = cut_efficiency_score(with_polyT)
        eff_without = cut_efficiency_score(without_polyT)
        assert eff_with < eff_without  # Poly-T should be penalized

    def test_cut_efficiency_G_at_20(self):
        """Test that G at position 20 is favored."""
        with_G = "GCGACTGCTACATAGCCAGG"
        with_A = "GCGACTGCTACATAGCCAGA"
        eff_G = cut_efficiency_score(with_G)
        eff_A = cut_efficiency_score(with_A)
        assert eff_G >= eff_A  # G at pos 20 should be equal or better

    def test_frameshift_probability_early_cut(self):
        """Test that early cuts have higher frameshift probability."""
        # frameshift_probability(guide_position, exon_length, cds_position)
        early = frameshift_probability(25, 300, 25)  # First 10% of CDS
        late = frameshift_probability(250, 300, 250)  # ~83% through CDS
        assert early > late

    def test_repair_pathway_nhej_dominant(self):
        """Test that NHEJ is dominant repair pathway."""
        guide = "GCGACTGCTACATAGCCAGG"
        repair = repair_pathway_prediction(guide)
        assert repair['NHEJ'] > repair['HDR']
        assert repair['preferred'] == 'NHEJ'


class TestCRISPRiModule:
    """Tests for CRISPRi module."""

    def test_repression_efficiency_optimal_position(self):
        """Test that TSS-proximal positions have high efficiency."""
        guide = "GCGACTGCTACATAGCCAGG"
        # Position 75bp from TSS on template strand - optimal
        eff_optimal = repression_efficiency_score(guide, 75, '+', 'KRAB')
        # Position 400bp - suboptimal
        eff_far = repression_efficiency_score(guide, 400, '+', 'KRAB')
        assert eff_optimal > eff_far

    def test_steric_hindrance_peak(self):
        """Test that steric hindrance peaks near TSS."""
        peak = steric_hindrance_score(50, '+')
        far = steric_hindrance_score(400, '+')
        assert peak > far

    def test_repressor_type_affects_efficiency(self):
        """Test that different repressors have different efficiencies."""
        guide = "GCGACTGCTACATAGCCAGG"
        krab = repression_efficiency_score(guide, 75, '+', 'KRAB')
        dcas9_only = repression_efficiency_score(guide, 75, '+', 'dCas9_only')
        assert krab > dcas9_only  # KRAB is more effective


class TestPrimeEditingModule:
    """Tests for prime editing module."""

    def test_reverse_complement(self):
        """Test reverse complement function."""
        seq = "ATGC"
        assert reverse_complement(seq) == "GCAT"

    def test_reverse_complement_longer(self):
        """Test reverse complement on longer sequence."""
        seq = "GCGACTGCTACATAGCCAGG"
        rc = reverse_complement(seq)
        assert len(rc) == len(seq)
        assert rc[0] == 'C'  # complement of G
        assert rc[-1] == 'C'  # complement of G

    def test_pbs_score_optimal_length(self):
        """Test that optimal PBS length scores higher."""
        pbs_14 = "CCTGGCTATGTAGC"  # 14bp - optimal
        pbs_8 = "CGATCGAT"  # 8bp - minimum
        score_14 = pbs_score(pbs_14)
        score_8 = pbs_score(pbs_8)
        assert score_14 > score_8

    def test_rt_template_score(self):
        """Test RT template scoring."""
        good_rt = "AGTCGATCGATCG"  # 13bp, moderate GC
        score = rt_template_score(good_rt)
        assert 0 < score <= 2.0

    def test_hairpin_detection(self):
        """Test secondary structure detection."""
        # Sequence with self-complementary regions
        hairpin_seq = "GCGCGCGCGCGCGC"  # Will form structure
        simple_seq = "AGTCAGTCAGTCAG"  # Less structure
        dg_hairpin = estimate_hairpin_dg(hairpin_seq)
        dg_simple = estimate_hairpin_dg(simple_seq)
        # More negative = more structure
        assert dg_hairpin <= dg_simple


class TestBaseEditingModule:
    """Tests for base editing module."""

    def test_activity_windows(self):
        """Test activity window retrieval."""
        abe_window = get_activity_window('ABE8e')
        cbe_window = get_activity_window('BE4')
        assert abe_window == (4, 8)
        assert cbe_window == (4, 8)

    def test_editing_efficiency_position_5(self):
        """Test that position 5 has peak efficiency."""
        eff_5 = editing_efficiency_at_position(5, 'ABE8e')
        eff_1 = editing_efficiency_at_position(1, 'ABE8e')
        eff_10 = editing_efficiency_at_position(10, 'ABE8e')
        assert eff_5 == 1.0
        assert eff_5 > eff_1
        assert eff_5 > eff_10

    def test_context_scoring_abe(self):
        """Test ABE context preferences."""
        # Position is 1-indexed, so position 4 means guide[3]
        # Context is guide[idx-1] + guide[idx] = guide[2] + guide[3]
        guide = "GCTACTGCTACATAGCCAGG"  # TA at positions 3-4 (0-indexed: 2-3)
        ctx = sequence_context_score(guide, 4, 'ABE8e')  # A at pos 4, preceded by T
        assert ctx == 1.1  # TA context preferred for ABE

    def test_bystander_detection(self):
        """Test bystander A detection in activity window."""
        guide = "AAAACTGCTACATAGCCAGG"  # Multiple As at positions 1-4
        bystanders = find_bystanders(guide, 5, 'ABE8e')  # Target at pos 5
        # Should find A at position 4 (in window)
        positions = [b['position'] for b in bystanders]
        assert 4 in positions


class TestConfigs:
    """Test configuration dataclasses."""

    def test_knockout_config_defaults(self):
        """Test KnockoutConfig default values."""
        config = KnockoutConfig()
        assert config.pam == "NGG"
        assert config.guide_length == 20
        assert config.min_cut_efficiency == 0.3

    def test_crispri_config_defaults(self):
        """Test CRISPRiConfig default values."""
        config = CRISPRiConfig()
        assert config.repressor == "KRAB"
        assert config.crispri_window == (-50, +300)

    def test_prime_edit_config_defaults(self):
        """Test PrimeEditConfig default values."""
        config = PrimeEditConfig()
        assert config.pbs_length_min == 8
        assert config.pbs_length_max == 17
        assert config.rt_length_min == 7

    def test_base_edit_config_defaults(self):
        """Test BaseEditConfig default values."""
        config = BaseEditConfig()
        assert config.editor == "ABE8e"
        assert config.target_base == "A"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
