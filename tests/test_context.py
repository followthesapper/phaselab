"""Tests for the biological context module."""

import pytest
import numpy as np

from phaselab.context import (
    ContextStack,
    CellType,
    ChromatinContext,
    MethylationContext,
    HistoneContext,
    AccessibilityScore,
    MethylationScore,
    HistoneMarks,
    BiologicalContext,
)
from phaselab.context.chromatin import AccessibilitySource
from phaselab.context.methylation import MethylationSource, CpGSite
from phaselab.context.histone import HistoneMark, ChromatinState


class TestChromatinContext:
    """Tests for ChromatinContext."""

    def test_init(self):
        """Test context initialization."""
        ctx = ChromatinContext("K562")
        assert ctx.cell_type == "K562"
        assert ctx.source == AccessibilitySource.ENCODE_ATAC

    def test_get_accessibility_no_data(self):
        """Test accessibility with no data loaded."""
        ctx = ChromatinContext("K562")
        score = ctx.get_accessibility("chr4", 55095264)

        # Should return neutral default
        assert isinstance(score, AccessibilityScore)
        assert score.score == 0.5
        assert score.confidence == 0.0  # Low confidence without data

    def test_accessibility_score_bounds(self):
        """Test AccessibilityScore clipping."""
        # Score should be clipped to [0, 1]
        score = AccessibilityScore(
            score=1.5,  # Out of bounds
            raw_signal=100.0,
            peak_overlap=True,
            source=AccessibilitySource.ENCODE_ATAC,
            cell_type="K562",
            confidence=2.0,  # Out of bounds
        )
        assert score.score == 1.0
        assert score.confidence == 1.0


class TestMethylationContext:
    """Tests for MethylationContext."""

    def test_init(self):
        """Test context initialization."""
        ctx = MethylationContext("K562")
        assert ctx.cell_type == "K562"
        assert ctx.source == MethylationSource.ENCODE_WGBS

    def test_get_methylation_no_data(self):
        """Test methylation with no data loaded."""
        ctx = MethylationContext("K562")
        meth = ctx.get_methylation("chr4", 55095264, 55095284)

        assert isinstance(meth, MethylationScore)
        assert meth.confidence == 0.0

    def test_predict_cpg_sites(self):
        """Test CpG site prediction from sequence."""
        ctx = MethylationContext("K562")

        # Sequence with known CpG sites
        sequence = "ATCGATCGATCGATCGATCG"  # CG at positions 2, 8, 14
        sites = ctx._predict_cpg_sites(sequence, 0)

        # Should find CpG dinucleotides
        assert len(sites) > 0
        assert all(isinstance(s, CpGSite) for s in sites)

    def test_methylation_with_sequence(self):
        """Test methylation context using sequence."""
        ctx = MethylationContext("K562")

        sequence = "CGCGCGCGCGATCGATCGATCG"  # CpG-rich region
        meth = ctx.get_methylation("chr4", 0, len(sequence), sequence=sequence)

        assert meth.cpg_count > 0
        assert meth.cpg_density > 0

    def test_cpg_effect_on_editing(self):
        """Test editing efficiency modifier."""
        ctx = MethylationContext("K562")

        # High methylation should reduce efficiency
        high_meth = MethylationScore(
            mean_methylation=0.9,
            cpg_count=5,
            cpg_density=5.0,
            is_cpg_island=False,
            cpg_sites=[],
            source=MethylationSource.ENCODE_WGBS,
            cell_type="K562",
        )
        modifier = ctx.get_cpg_effect_on_editing(high_meth)
        assert modifier < 1.0

        # CpG island with low methylation should be accessible
        low_meth = MethylationScore(
            mean_methylation=0.1,
            cpg_count=10,
            cpg_density=10.0,
            is_cpg_island=True,
            cpg_sites=[],
            source=MethylationSource.ENCODE_WGBS,
            cell_type="K562",
        )
        modifier = ctx.get_cpg_effect_on_editing(low_meth)
        assert modifier > 0.9


class TestHistoneContext:
    """Tests for HistoneContext."""

    def test_init(self):
        """Test context initialization."""
        ctx = HistoneContext("K562")
        assert ctx.cell_type == "K562"

    def test_get_histone_marks_no_data(self):
        """Test histone marks with no data loaded."""
        ctx = HistoneContext("K562")
        marks = ctx.get_histone_marks("chr4", 55095264, 55095284)

        assert isinstance(marks, HistoneMarks)
        assert marks.confidence < 1.0  # Low confidence without data

    def test_histone_marks_activity_score(self):
        """Test activity score calculation."""
        # Active marks present
        marks = HistoneMarks(
            marks={
                HistoneMark.H3K4ME3: 0.8,
                HistoneMark.H3K27AC: 0.7,
                HistoneMark.H3K27ME3: 0.1,
            },
            chromatin_state=ChromatinState.ACTIVE_TSS,
            is_active=True,
            is_repressed=False,
            cell_type="K562",
        )
        assert marks.activity_score > 0.5

        # Repressive marks present
        marks_rep = HistoneMarks(
            marks={
                HistoneMark.H3K4ME3: 0.1,
                HistoneMark.H3K27ME3: 0.8,
                HistoneMark.H3K9ME3: 0.7,
            },
            chromatin_state=ChromatinState.HETEROCHROMATIN,
            is_active=False,
            is_repressed=True,
            cell_type="K562",
        )
        assert marks_rep.activity_score < 0.5

    def test_editing_modifier(self):
        """Test editing efficiency modifier by chromatin state."""
        ctx = HistoneContext("K562")

        # Active TSS should boost editing
        active_marks = HistoneMarks(
            marks={},
            chromatin_state=ChromatinState.ACTIVE_TSS,
            is_active=True,
            is_repressed=False,
            cell_type="K562",
        )
        modifier = ctx.get_editing_modifier(active_marks)
        assert modifier > 1.0

        # Heterochromatin should reduce editing
        hetero_marks = HistoneMarks(
            marks={},
            chromatin_state=ChromatinState.HETEROCHROMATIN,
            is_active=False,
            is_repressed=True,
            cell_type="K562",
        )
        modifier = ctx.get_editing_modifier(hetero_marks)
        assert modifier < 0.5


class TestContextStack:
    """Tests for the unified ContextStack."""

    def test_init_with_enum(self):
        """Test initialization with CellType enum."""
        ctx = ContextStack(cell_type=CellType.K562)
        assert ctx.cell_type == CellType.K562

    def test_init_with_string(self):
        """Test initialization with string cell type."""
        ctx = ContextStack(cell_type="K562")
        assert ctx.cell_type == CellType.K562

    def test_init_with_custom_string(self):
        """Test initialization with custom cell type."""
        ctx = ContextStack(cell_type="CustomCellLine")
        assert ctx.cell_type == CellType.CUSTOM
        assert ctx._custom_cell_type == "CustomCellLine"

    def test_get_context(self):
        """Test getting complete biological context."""
        ctx = ContextStack(cell_type=CellType.K562)
        context = ctx.get_context("chr4", 55095264, 55095284)

        assert isinstance(context, BiologicalContext)
        assert context.chrom == "chr4"
        assert context.start == 55095264
        assert context.end == 55095284
        assert isinstance(context.accessibility, AccessibilityScore)
        assert isinstance(context.methylation, MethylationScore)
        assert isinstance(context.histones, HistoneMarks)

    def test_context_overall_score(self):
        """Test overall context score calculation."""
        ctx = ContextStack(cell_type=CellType.K562)
        context = ctx.get_context("chr4", 55095264, 55095284)

        # Score should be in [0, 1]
        assert 0.0 <= context.overall_score <= 1.0

    def test_context_editing_modifier(self):
        """Test editing modifier from context."""
        ctx = ContextStack(cell_type=CellType.K562)
        context = ctx.get_context("chr4", 55095264, 55095284)

        # Modifier should be in reasonable range
        assert 0.3 <= context.editing_modifier <= 1.3

    def test_context_to_dict(self):
        """Test context serialization."""
        ctx = ContextStack(cell_type=CellType.K562)
        context = ctx.get_context("chr4", 55095264, 55095284)
        d = context.to_dict()

        assert "chrom" in d
        assert "overall_score" in d
        assert "accessibility" in d
        assert "methylation" in d
        assert "histones" in d

    def test_batch_context(self):
        """Test batch context retrieval."""
        ctx = ContextStack(cell_type=CellType.K562)
        regions = [
            ("chr4", 55095264, 55095284),
            ("chr4", 55096000, 55096020),
        ]
        contexts = ctx.batch_context(regions)

        assert len(contexts) == 2
        assert all(isinstance(c, BiologicalContext) for c in contexts)

    def test_available_cell_types(self):
        """Test listing available cell types."""
        cell_types = ContextStack.available_cell_types()

        assert isinstance(cell_types, list)
        assert "K562" in cell_types
        assert "custom" not in cell_types

    def test_context_manager(self):
        """Test context manager protocol."""
        with ContextStack(cell_type=CellType.K562) as ctx:
            context = ctx.get_context("chr4", 55095264, 55095284)
            assert isinstance(context, BiologicalContext)
