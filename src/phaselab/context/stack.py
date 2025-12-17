"""
Unified biological context stack.

Combines chromatin accessibility, methylation, and histone marks
into a single interface for CRISPR guide scoring.
"""

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional, Dict, List, Any, Union
import logging

import numpy as np

from phaselab.context.chromatin import (
    ChromatinContext,
    AccessibilityScore,
    AccessibilitySource,
)
from phaselab.context.methylation import (
    MethylationContext,
    MethylationScore,
    MethylationSource,
)
from phaselab.context.histone import (
    HistoneContext,
    HistoneMarks,
    HistoneMark,
    ChromatinState,
)

logger = logging.getLogger(__name__)


class CellType(Enum):
    """Common cell types with ENCODE data."""
    K562 = "K562"                    # Chronic myelogenous leukemia
    HEPG2 = "HepG2"                  # Hepatocellular carcinoma
    GM12878 = "GM12878"              # B-lymphocyte
    HEK293 = "HEK293"                # Embryonic kidney
    IPSC = "iPSC"                    # Induced pluripotent stem cells
    H1_HESC = "H1-hESC"              # Human embryonic stem cells
    A549 = "A549"                    # Lung carcinoma
    MCFN7 = "MCF-7"                  # Breast cancer
    HELA = "HeLa"                    # Cervical cancer
    CUSTOM = "custom"                # User-provided data


@dataclass
class BiologicalContext:
    """
    Complete biological context for a genomic region.

    Combines accessibility, methylation, and histone information
    into a unified view of the chromatin landscape.

    Attributes
    ----------
    chrom : str
        Chromosome
    start : int
        Region start
    end : int
        Region end
    accessibility : AccessibilityScore
        Chromatin accessibility
    methylation : MethylationScore
        DNA methylation context
    histones : HistoneMarks
        Histone modification profile
    cell_type : CellType
        Cell type
    overall_score : float
        Combined context favorability score [0, 1]
    confidence : float
        Overall confidence in the context
    """
    chrom: str
    start: int
    end: int
    accessibility: AccessibilityScore
    methylation: MethylationScore
    histones: HistoneMarks
    cell_type: CellType
    overall_score: float = 0.5
    confidence: float = 0.0

    def __post_init__(self):
        # Calculate overall score
        self.overall_score = self._calculate_overall_score()
        # Calculate confidence
        self.confidence = self._calculate_confidence()

    def _calculate_overall_score(self) -> float:
        """
        Calculate overall context favorability.

        Higher scores indicate more favorable editing context.
        """
        # Accessibility is key - open chromatin is easier to edit
        accessibility_contrib = self.accessibility.score * 0.4

        # Low methylation is favorable
        methylation_contrib = (1.0 - self.methylation.mean_methylation) * 0.3

        # Active chromatin state is favorable
        histone_contrib = self.histones.activity_score * 0.3

        return float(accessibility_contrib + methylation_contrib + histone_contrib)

    def _calculate_confidence(self) -> float:
        """Calculate overall confidence from component confidences."""
        confidences = [
            self.accessibility.confidence,
            self.methylation.confidence,
            self.histones.confidence,
        ]
        return float(np.mean(confidences))

    @property
    def is_favorable(self) -> bool:
        """Whether context is favorable for editing."""
        return self.overall_score > 0.5

    @property
    def chromatin_state(self) -> ChromatinState:
        """Predicted chromatin state."""
        return self.histones.chromatin_state

    @property
    def editing_modifier(self) -> float:
        """
        Multiplier for editing efficiency based on context.

        Returns
        -------
        float
            Efficiency modifier in range [0.3, 1.3]
        """
        # Base from accessibility
        base = 0.5 + 0.5 * self.accessibility.score

        # Methylation penalty
        meth_penalty = 0.3 * self.methylation.mean_methylation

        # Chromatin state bonus/penalty
        state_effects = {
            ChromatinState.ACTIVE_TSS: 0.3,
            ChromatinState.FLANKING_TSS: 0.2,
            ChromatinState.STRONG_ENHANCER: 0.2,
            ChromatinState.WEAK_ENHANCER: 0.1,
            ChromatinState.TRANSCRIBED: 0.0,
            ChromatinState.REPRESSED_POLYCOMB: -0.3,
            ChromatinState.HETEROCHROMATIN: -0.4,
            ChromatinState.QUIESCENT: -0.1,
        }
        state_effect = state_effects.get(self.chromatin_state, 0.0)

        return float(np.clip(base - meth_penalty + state_effect, 0.3, 1.3))

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "cell_type": self.cell_type.value,
            "overall_score": self.overall_score,
            "confidence": self.confidence,
            "editing_modifier": self.editing_modifier,
            "accessibility": {
                "score": self.accessibility.score,
                "peak_overlap": self.accessibility.peak_overlap,
            },
            "methylation": {
                "mean": self.methylation.mean_methylation,
                "cpg_count": self.methylation.cpg_count,
                "is_cpg_island": self.methylation.is_cpg_island,
            },
            "histones": {
                "state": self.chromatin_state.value,
                "is_active": self.histones.is_active,
                "is_repressed": self.histones.is_repressed,
                "marks": {
                    k.value: v for k, v in self.histones.marks.items()
                },
            },
        }


class ContextStack:
    """
    Unified biological context provider.

    Combines chromatin accessibility, DNA methylation, and histone marks
    from ENCODE/Roadmap Epigenomics into a single interface.

    Parameters
    ----------
    cell_type : CellType or str
        Cell type for context
    cache_dir : Path, optional
        Cache directory for downloaded data
    lazy_load : bool
        Whether to defer data loading until first query

    Examples
    --------
    >>> from phaselab.context import ContextStack, CellType
    >>>
    >>> # Create context for K562 cells
    >>> ctx = ContextStack(cell_type=CellType.K562)
    >>>
    >>> # Get context for a guide target region
    >>> context = ctx.get_context("chr4", 55095264, 55095284)
    >>>
    >>> print(f"Overall score: {context.overall_score:.2f}")
    >>> print(f"Editing modifier: {context.editing_modifier:.2f}")
    >>> print(f"Chromatin state: {context.chromatin_state.value}")
    """

    def __init__(
        self,
        cell_type: Union[CellType, str] = CellType.K562,
        cache_dir: Optional[Path] = None,
        lazy_load: bool = True,
    ):
        # Initialize _custom_cell_type first
        self._custom_cell_type: Optional[str] = None

        if isinstance(cell_type, str):
            try:
                self.cell_type = CellType(cell_type)
            except ValueError:
                self.cell_type = CellType.CUSTOM
                self._custom_cell_type = cell_type
        else:
            self.cell_type = cell_type

        self.cache_dir = cache_dir or Path.home() / ".phaselab" / "cache"
        self.lazy_load = lazy_load

        # Cell type string for providers
        cell_str = (
            self._custom_cell_type
            if self._custom_cell_type
            else self.cell_type.value
        )

        # Initialize providers
        self._chromatin = ChromatinContext(
            cell_type=cell_str,
            cache_dir=self.cache_dir / "chromatin",
        )
        self._methylation = MethylationContext(
            cell_type=cell_str,
            cache_dir=self.cache_dir / "methylation",
        )
        self._histone = HistoneContext(
            cell_type=cell_str,
            cache_dir=self.cache_dir / "histone",
        )

        self._initialized = False

    def _ensure_initialized(self) -> None:
        """Initialize data providers if needed."""
        if self._initialized or self.lazy_load:
            return
        # Force load all data
        self._chromatin._ensure_data()
        self._methylation._ensure_data()
        self._histone._ensure_data()
        self._initialized = True

    def get_context(
        self,
        chrom: str,
        start: int,
        end: int,
        sequence: Optional[str] = None,
    ) -> BiologicalContext:
        """
        Get complete biological context for a genomic region.

        Parameters
        ----------
        chrom : str
            Chromosome (e.g., "chr4")
        start : int
            Region start position
        end : int
            Region end position
        sequence : str, optional
            Sequence for CpG prediction if no methylation data

        Returns
        -------
        BiologicalContext
            Complete context including accessibility, methylation, histones
        """
        # Get accessibility (use midpoint with window)
        midpoint = (start + end) // 2
        accessibility = self._chromatin.get_accessibility(
            chrom, midpoint, window=end - start
        )

        # Get methylation
        methylation = self._methylation.get_methylation(
            chrom, start, end, sequence=sequence
        )

        # Get histone marks
        histones = self._histone.get_histone_marks(chrom, start, end)

        return BiologicalContext(
            chrom=chrom,
            start=start,
            end=end,
            accessibility=accessibility,
            methylation=methylation,
            histones=histones,
            cell_type=self.cell_type,
        )

    def get_accessibility(
        self,
        chrom: str,
        position: int,
        window: int = 200,
    ) -> AccessibilityScore:
        """Get accessibility score for a position."""
        return self._chromatin.get_accessibility(chrom, position, window)

    def get_methylation(
        self,
        chrom: str,
        start: int,
        end: int,
        sequence: Optional[str] = None,
    ) -> MethylationScore:
        """Get methylation context for a region."""
        return self._methylation.get_methylation(chrom, start, end, sequence)

    def get_histone_marks(
        self,
        chrom: str,
        start: int,
        end: int,
    ) -> HistoneMarks:
        """Get histone marks for a region."""
        return self._histone.get_histone_marks(chrom, start, end)

    def batch_context(
        self,
        regions: List[tuple],
    ) -> List[BiologicalContext]:
        """
        Get context for multiple regions efficiently.

        Parameters
        ----------
        regions : List[tuple]
            List of (chrom, start, end) tuples

        Returns
        -------
        List[BiologicalContext]
            Context for each region
        """
        return [
            self.get_context(chrom, start, end)
            for chrom, start, end in regions
        ]

    def download_all_data(self, force: bool = False) -> Dict[str, bool]:
        """
        Download all context data from ENCODE.

        Parameters
        ----------
        force : bool
            Re-download even if cached

        Returns
        -------
        Dict[str, bool]
            Success status for each data type
        """
        results = {}
        results["chromatin"] = self._chromatin.download_data(force=force)
        results["methylation"] = self._methylation.download_data(force=force)
        results["histones"] = self._histone.download_data()
        return results

    @staticmethod
    def available_cell_types() -> List[str]:
        """List cell types with available ENCODE data."""
        return [ct.value for ct in CellType if ct != CellType.CUSTOM]

    def close(self) -> None:
        """Close all data handles."""
        self._chromatin.close()
        self._histone.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False
