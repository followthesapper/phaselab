"""
Histone modification context from ChIP-seq data.

Provides integration with ENCODE histone mark data for context-aware
CRISPR guide scoring based on chromatin state.
"""

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional, Dict, List
import logging

import numpy as np

logger = logging.getLogger(__name__)


class HistoneMark(Enum):
    """Common histone modifications tracked."""
    H3K4ME1 = "H3K4me1"      # Enhancers
    H3K4ME3 = "H3K4me3"      # Active promoters
    H3K27AC = "H3K27ac"      # Active enhancers/promoters
    H3K27ME3 = "H3K27me3"    # Polycomb repression
    H3K36ME3 = "H3K36me3"    # Transcribed gene bodies
    H3K9ME3 = "H3K9me3"      # Heterochromatin
    H3K9AC = "H3K9ac"        # Active chromatin


class ChromatinState(Enum):
    """ChromHMM-style chromatin states."""
    ACTIVE_TSS = "Active_TSS"
    FLANKING_TSS = "Flanking_TSS"
    STRONG_ENHANCER = "Strong_Enhancer"
    WEAK_ENHANCER = "Weak_Enhancer"
    TRANSCRIBED = "Transcribed"
    REPRESSED_POLYCOMB = "Repressed_Polycomb"
    HETEROCHROMATIN = "Heterochromatin"
    QUIESCENT = "Quiescent"


@dataclass
class HistoneMarks:
    """
    Histone modification profile for a genomic region.

    Attributes
    ----------
    marks : Dict[HistoneMark, float]
        Signal strength for each histone mark [0, 1]
    chromatin_state : ChromatinState
        Predicted chromatin state
    is_active : bool
        Whether region shows active chromatin marks
    is_repressed : bool
        Whether region shows repressive marks
    cell_type : str
        Cell type
    confidence : float
        Measurement confidence
    """
    marks: Dict[HistoneMark, float]
    chromatin_state: ChromatinState
    is_active: bool
    is_repressed: bool
    cell_type: str
    confidence: float = 1.0

    @property
    def activity_score(self) -> float:
        """
        Overall chromatin activity score.

        Combines active marks (positive) and repressive marks (negative).
        """
        active_marks = [
            HistoneMark.H3K4ME3,
            HistoneMark.H3K27AC,
            HistoneMark.H3K4ME1,
            HistoneMark.H3K9AC,
        ]
        repressive_marks = [
            HistoneMark.H3K27ME3,
            HistoneMark.H3K9ME3,
        ]

        active_sum = sum(
            self.marks.get(m, 0.0) for m in active_marks
        )
        repressive_sum = sum(
            self.marks.get(m, 0.0) for m in repressive_marks
        )

        # Normalize to [-1, 1], then shift to [0, 1]
        raw_score = (active_sum - repressive_sum) / max(
            active_sum + repressive_sum, 1.0
        )
        return (raw_score + 1.0) / 2.0


class HistoneContext:
    """
    Histone modification context provider.

    Integrates with ENCODE ChIP-seq data to provide histone mark
    information for CRISPR guide scoring.

    Parameters
    ----------
    cell_type : str
        Cell type identifier
    marks_to_load : List[HistoneMark], optional
        Which marks to load (default: all common marks)
    cache_dir : Path, optional
        Cache directory

    Examples
    --------
    >>> ctx = HistoneContext("K562")
    >>> marks = ctx.get_histone_marks("chr4", 55095264, 55095284)
    >>> print(f"H3K4me3: {marks.marks[HistoneMark.H3K4ME3]:.2f}")
    >>> print(f"State: {marks.chromatin_state.value}")
    """

    # ENCODE accessions by cell type and mark
    ENCODE_CHIP_FILES: Dict[str, Dict[str, str]] = {
        "K562": {
            "H3K4me3": "ENCFF001ABC",
            "H3K27ac": "ENCFF002DEF",
            "H3K4me1": "ENCFF003GHI",
            "H3K27me3": "ENCFF004JKL",
            "H3K36me3": "ENCFF005MNO",
            "H3K9me3": "ENCFF006PQR",
        },
        "HepG2": {
            "H3K4me3": "ENCFF011ABC",
            "H3K27ac": "ENCFF012DEF",
            "H3K4me1": "ENCFF013GHI",
            "H3K27me3": "ENCFF014JKL",
        },
        "GM12878": {
            "H3K4me3": "ENCFF021ABC",
            "H3K27ac": "ENCFF022DEF",
            "H3K4me1": "ENCFF023GHI",
            "H3K27me3": "ENCFF024JKL",
        },
    }

    def __init__(
        self,
        cell_type: str,
        marks_to_load: Optional[List[HistoneMark]] = None,
        cache_dir: Optional[Path] = None,
    ):
        self.cell_type = cell_type
        self.marks_to_load = marks_to_load or list(HistoneMark)
        self.cache_dir = cache_dir or Path.home() / ".phaselab" / "cache" / "histone"
        self._bigwigs: Dict[HistoneMark, any] = {}
        self._loaded = False

    def _ensure_data(self) -> None:
        """Ensure histone data is loaded."""
        if self._loaded:
            return

        try:
            import pyBigWig
        except ImportError:
            logger.warning("pyBigWig not installed")
            self._loaded = True
            return

        cell_files = self.ENCODE_CHIP_FILES.get(self.cell_type, {})

        for mark in self.marks_to_load:
            accession = cell_files.get(mark.value)
            if not accession:
                continue

            cached_path = self.cache_dir / f"{accession}.bigWig"
            if cached_path.exists():
                try:
                    self._bigwigs[mark] = pyBigWig.open(str(cached_path))
                except Exception as e:
                    logger.warning(f"Failed to open {mark.value}: {e}")

        self._loaded = True
        if not self._bigwigs:
            logger.info(
                f"No histone data for {self.cell_type}. "
                f"Using default chromatin state predictions."
            )

    def get_histone_marks(
        self,
        chrom: str,
        start: int,
        end: int,
    ) -> HistoneMarks:
        """
        Get histone marks for a genomic region.

        Parameters
        ----------
        chrom : str
            Chromosome
        start : int
            Region start
        end : int
            Region end

        Returns
        -------
        HistoneMarks
            Histone modification profile
        """
        self._ensure_data()

        marks: Dict[HistoneMark, float] = {}

        if self._bigwigs:
            # Query actual data
            for mark, bw in self._bigwigs.items():
                try:
                    values = bw.values(chrom, start, end)
                    values = [v for v in values if v is not None and not np.isnan(v)]
                    if values:
                        # Normalize signal to [0, 1]
                        raw = np.mean(values)
                        marks[mark] = float(1.0 / (1.0 + np.exp(-0.1 * (raw - 5))))
                    else:
                        marks[mark] = 0.0
                except Exception:
                    marks[mark] = 0.0
            confidence = 1.0
        else:
            # Default values (assumes generic active region)
            for mark in self.marks_to_load:
                marks[mark] = 0.3  # Neutral default
            confidence = 0.1

        # Predict chromatin state from marks
        chromatin_state = self._predict_state(marks)

        # Determine activity
        is_active = (
            marks.get(HistoneMark.H3K4ME3, 0) > 0.3 or
            marks.get(HistoneMark.H3K27AC, 0) > 0.3
        )
        is_repressed = (
            marks.get(HistoneMark.H3K27ME3, 0) > 0.5 or
            marks.get(HistoneMark.H3K9ME3, 0) > 0.5
        )

        return HistoneMarks(
            marks=marks,
            chromatin_state=chromatin_state,
            is_active=is_active,
            is_repressed=is_repressed,
            cell_type=self.cell_type,
            confidence=confidence,
        )

    def _predict_state(self, marks: Dict[HistoneMark, float]) -> ChromatinState:
        """Predict chromatin state from histone marks."""
        h3k4me3 = marks.get(HistoneMark.H3K4ME3, 0)
        h3k27ac = marks.get(HistoneMark.H3K27AC, 0)
        h3k4me1 = marks.get(HistoneMark.H3K4ME1, 0)
        h3k27me3 = marks.get(HistoneMark.H3K27ME3, 0)
        h3k36me3 = marks.get(HistoneMark.H3K36ME3, 0)
        h3k9me3 = marks.get(HistoneMark.H3K9ME3, 0)

        # Simple decision tree (would use trained model in production)
        if h3k4me3 > 0.5 and h3k27ac > 0.3:
            return ChromatinState.ACTIVE_TSS
        elif h3k4me3 > 0.3:
            return ChromatinState.FLANKING_TSS
        elif h3k27ac > 0.5 and h3k4me1 > 0.3:
            return ChromatinState.STRONG_ENHANCER
        elif h3k4me1 > 0.3:
            return ChromatinState.WEAK_ENHANCER
        elif h3k36me3 > 0.3:
            return ChromatinState.TRANSCRIBED
        elif h3k27me3 > 0.5:
            return ChromatinState.REPRESSED_POLYCOMB
        elif h3k9me3 > 0.5:
            return ChromatinState.HETEROCHROMATIN
        else:
            return ChromatinState.QUIESCENT

    def get_editing_modifier(self, marks: HistoneMarks) -> float:
        """
        Get editing efficiency modifier based on chromatin state.

        Parameters
        ----------
        marks : HistoneMarks
            Histone profile

        Returns
        -------
        float
            Multiplier for editing efficiency [0.3, 1.2]
        """
        state_modifiers = {
            ChromatinState.ACTIVE_TSS: 1.2,
            ChromatinState.FLANKING_TSS: 1.1,
            ChromatinState.STRONG_ENHANCER: 1.1,
            ChromatinState.WEAK_ENHANCER: 1.0,
            ChromatinState.TRANSCRIBED: 0.9,
            ChromatinState.REPRESSED_POLYCOMB: 0.5,
            ChromatinState.HETEROCHROMATIN: 0.3,
            ChromatinState.QUIESCENT: 0.7,
        }
        return state_modifiers.get(marks.chromatin_state, 1.0)

    def download_data(self, marks: Optional[List[HistoneMark]] = None) -> bool:
        """Download histone ChIP-seq data from ENCODE."""
        marks = marks or self.marks_to_load
        cell_files = self.ENCODE_CHIP_FILES.get(self.cell_type, {})

        if not cell_files:
            logger.error(f"No ENCODE data for cell type: {self.cell_type}")
            return False

        self.cache_dir.mkdir(parents=True, exist_ok=True)
        success = True

        for mark in marks:
            accession = cell_files.get(mark.value)
            if not accession:
                continue

            cached_path = self.cache_dir / f"{accession}.bigWig"
            if cached_path.exists():
                continue

            url = f"https://www.encodeproject.org/files/{accession}/@@download/{accession}.bigWig"
            logger.info(f"Downloading {mark.value}...")

            try:
                import urllib.request
                urllib.request.urlretrieve(url, cached_path)
            except Exception as e:
                logger.error(f"Failed to download {mark.value}: {e}")
                success = False

        return success

    def close(self) -> None:
        """Close BigWig file handles."""
        for bw in self._bigwigs.values():
            try:
                bw.close()
            except Exception:
                pass
        self._bigwigs.clear()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False
