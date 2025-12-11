"""
PhaseLab Chromatin Accessibility: Core accessibility scoring module.

Implements chromatin accessibility models for predicting CRISPR guide
efficacy based on local chromatin state.
"""

import numpy as np
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass
from enum import Enum


class ChromatinState(Enum):
    """Chromatin accessibility states."""
    OPEN = "open"
    PARTIALLY_OPEN = "partially_open"
    CLOSED = "closed"
    UNKNOWN = "unknown"


@dataclass
class ChromatinScore:
    """Container for chromatin accessibility scores."""

    # Overall score
    accessibility: float  # 0-1 scale
    state: ChromatinState

    # Component scores
    dnase_score: float
    atac_score: float
    histone_score: float

    # Position info
    chrom: str
    start: int
    end: int

    # Metadata
    tissue: str
    cell_type: Optional[str] = None
    data_source: str = "model"

    def to_dict(self) -> Dict:
        return {
            "accessibility": self.accessibility,
            "state": self.state.value,
            "dnase_score": self.dnase_score,
            "atac_score": self.atac_score,
            "histone_score": self.histone_score,
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "tissue": self.tissue,
            "cell_type": self.cell_type,
            "data_source": self.data_source,
        }


@dataclass
class ChromatinModel:
    """
    Chromatin accessibility prediction model.

    Combines multiple signals to predict accessibility:
    - DNase hypersensitivity
    - ATAC-seq peaks
    - Histone modifications (H3K4me3, H3K27ac)
    """

    name: str
    tissue: str
    cell_type: Optional[str] = None

    # Model weights
    dnase_weight: float = 0.4
    atac_weight: float = 0.3
    histone_weight: float = 0.3

    # Thresholds
    open_threshold: float = 0.7
    closed_threshold: float = 0.3

    # Optional precomputed data
    dnase_peaks: Optional[Dict] = None
    atac_peaks: Optional[Dict] = None
    histone_marks: Optional[Dict] = None

    def score_position(
        self,
        chrom: str,
        position: int,
        window: int = 500
    ) -> ChromatinScore:
        """
        Score chromatin accessibility at a genomic position.

        Args:
            chrom: Chromosome
            position: Genomic position
            window: Window size around position

        Returns:
            ChromatinScore for the position
        """
        start = position - window // 2
        end = position + window // 2

        # Score each component
        dnase = self._score_dnase(chrom, start, end)
        atac = self._score_atac(chrom, start, end)
        histone = self._score_histone(chrom, start, end)

        # Weighted combination
        accessibility = (
            self.dnase_weight * dnase +
            self.atac_weight * atac +
            self.histone_weight * histone
        )

        # Determine state
        if accessibility >= self.open_threshold:
            state = ChromatinState.OPEN
        elif accessibility <= self.closed_threshold:
            state = ChromatinState.CLOSED
        else:
            state = ChromatinState.PARTIALLY_OPEN

        return ChromatinScore(
            accessibility=accessibility,
            state=state,
            dnase_score=dnase,
            atac_score=atac,
            histone_score=histone,
            chrom=chrom,
            start=start,
            end=end,
            tissue=self.tissue,
            cell_type=self.cell_type,
            data_source=self.name
        )

    def _score_dnase(self, chrom: str, start: int, end: int) -> float:
        """Score DNase hypersensitivity."""
        if self.dnase_peaks is None:
            return self._estimate_dnase(chrom, start, end)

        # Check for peaks in region
        key = f"{chrom}:{start}-{end}"
        if key in self.dnase_peaks:
            return float(self.dnase_peaks[key])

        # Check for overlapping peaks
        for peak_key, score in self.dnase_peaks.items():
            if self._regions_overlap(key, peak_key):
                return float(score)

        return 0.3  # Default moderate accessibility

    def _score_atac(self, chrom: str, start: int, end: int) -> float:
        """Score ATAC-seq accessibility."""
        if self.atac_peaks is None:
            return self._estimate_atac(chrom, start, end)

        key = f"{chrom}:{start}-{end}"
        if key in self.atac_peaks:
            return float(self.atac_peaks[key])

        return 0.3

    def _score_histone(self, chrom: str, start: int, end: int) -> float:
        """Score histone modifications."""
        if self.histone_marks is None:
            return self._estimate_histone(chrom, start, end)

        key = f"{chrom}:{start}-{end}"
        if key in self.histone_marks:
            return float(self.histone_marks[key])

        return 0.3

    def _estimate_dnase(self, chrom: str, start: int, end: int) -> float:
        """
        Estimate DNase score from sequence features.

        Uses GC content and promoter proximity as proxies.
        """
        # Promoter regions tend to be more accessible
        # This is a simplified model - real implementation would use
        # actual DNase-seq data or machine learning predictions
        return 0.5  # Neutral estimate

    def _estimate_atac(self, chrom: str, start: int, end: int) -> float:
        """Estimate ATAC score from sequence features."""
        return 0.5

    def _estimate_histone(self, chrom: str, start: int, end: int) -> float:
        """Estimate histone modification score."""
        return 0.5

    def _regions_overlap(self, region1: str, region2: str) -> bool:
        """Check if two genomic regions overlap."""
        try:
            chrom1, coords1 = region1.split(":")
            start1, end1 = map(int, coords1.split("-"))

            chrom2, coords2 = region2.split(":")
            start2, end2 = map(int, coords2.split("-"))

            if chrom1 != chrom2:
                return False

            return not (end1 < start2 or end2 < start1)
        except:
            return False


def compute_accessibility(
    chrom: str,
    position: int,
    tissue: str = "generic",
    cell_type: Optional[str] = None,
    model: Optional[ChromatinModel] = None
) -> ChromatinScore:
    """
    Compute chromatin accessibility for a genomic position.

    Args:
        chrom: Chromosome (e.g., "chr2")
        position: Genomic coordinate
        tissue: Tissue type for model selection
        cell_type: Optional cell type specification
        model: Optional pre-loaded ChromatinModel

    Returns:
        ChromatinScore with accessibility prediction
    """
    if model is None:
        model = get_default_model(tissue, cell_type)

    return model.score_position(chrom, position)


def get_default_model(
    tissue: str = "generic",
    cell_type: Optional[str] = None
) -> ChromatinModel:
    """Get default chromatin model for a tissue."""
    return ChromatinModel(
        name=f"default_{tissue}",
        tissue=tissue,
        cell_type=cell_type
    )


def load_encode_model(
    cell_line: str,
    data_dir: Optional[str] = None
) -> ChromatinModel:
    """
    Load chromatin model from ENCODE data.

    Args:
        cell_line: ENCODE cell line (e.g., "K562", "HepG2")
        data_dir: Directory containing ENCODE data files

    Returns:
        ChromatinModel with ENCODE data

    Note:
        Requires downloaded ENCODE data files.
        See docs/CHROMATIN.md for data preparation.
    """
    # TODO: Implement ENCODE data loading
    # This would parse DNase-seq and histone ChIP-seq peaks

    return ChromatinModel(
        name=f"ENCODE_{cell_line}",
        tissue=cell_line,
        cell_type=cell_line,
        data_source="ENCODE"
    )


def load_psychencode_model(
    brain_region: str = "cortex",
    data_dir: Optional[str] = None
) -> ChromatinModel:
    """
    Load brain-specific chromatin model from PsychENCODE.

    Args:
        brain_region: Brain region (cortex, hippocampus, etc.)
        data_dir: Directory containing PsychENCODE data

    Returns:
        ChromatinModel with brain-specific data

    Note:
        PsychENCODE provides chromatin accessibility data
        for human brain tissues, critical for neurological
        disease targets like SCN2A and RAI1.
    """
    # TODO: Implement PsychENCODE data loading

    return ChromatinModel(
        name=f"PsychENCODE_{brain_region}",
        tissue="brain",
        cell_type=brain_region,
        # Brain-specific weights - DNase more important
        dnase_weight=0.5,
        atac_weight=0.25,
        histone_weight=0.25
    )
