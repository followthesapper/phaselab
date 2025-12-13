"""
Chromatin accessibility context from ATAC-seq and DNase-seq data.

Provides integration with ENCODE and Roadmap Epigenomics accessibility data
for context-aware CRISPR guide scoring.
"""

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional, Union, Dict, Any
import logging

import numpy as np

logger = logging.getLogger(__name__)


class AccessibilitySource(Enum):
    """Source of accessibility data."""
    ENCODE_DNASE = "encode_dnase"
    ENCODE_ATAC = "encode_atac"
    ROADMAP = "roadmap"
    CUSTOM = "custom"


@dataclass
class AccessibilityScore:
    """
    Chromatin accessibility score for a genomic region.

    Attributes
    ----------
    score : float
        Normalized accessibility score in [0, 1].
        0 = completely closed, 1 = maximally accessible
    raw_signal : float
        Raw signal value from source data
    peak_overlap : bool
        Whether region overlaps a called peak
    source : AccessibilitySource
        Data source
    cell_type : str
        Cell type for this measurement
    confidence : float
        Confidence in the measurement (based on coverage, replicates)
    """
    score: float
    raw_signal: float
    peak_overlap: bool
    source: AccessibilitySource
    cell_type: str
    confidence: float = 1.0

    def __post_init__(self):
        self.score = float(np.clip(self.score, 0.0, 1.0))
        self.confidence = float(np.clip(self.confidence, 0.0, 1.0))


class ChromatinContext:
    """
    Chromatin accessibility context provider.

    Integrates with ENCODE/Roadmap BigWig files to provide accessibility
    scores for genomic regions.

    Parameters
    ----------
    cell_type : str
        Cell type identifier (e.g., "K562", "HepG2", "GM12878")
    source : AccessibilitySource
        Data source preference
    cache_dir : Path, optional
        Directory for caching downloaded data
    bigwig_path : Path, optional
        Path to custom BigWig file (for CUSTOM source)

    Examples
    --------
    >>> ctx = ChromatinContext("K562", source=AccessibilitySource.ENCODE_ATAC)
    >>> score = ctx.get_accessibility("chr4", 55095264)
    >>> print(f"Accessibility: {score.score:.2f}")
    """

    # ENCODE file accessions for common cell types
    ENCODE_ATAC_FILES: Dict[str, str] = {
        "K562": "ENCFF123ABC",  # Placeholder accessions
        "HepG2": "ENCFF456DEF",
        "GM12878": "ENCFF789GHI",
        "HEK293": "ENCFF012JKL",
        "iPSC": "ENCFF345MNO",
    }

    ENCODE_DNASE_FILES: Dict[str, str] = {
        "K562": "ENCFF111AAA",
        "HepG2": "ENCFF222BBB",
        "GM12878": "ENCFF333CCC",
        "HEK293": "ENCFF444DDD",
        "iPSC": "ENCFF555EEE",
    }

    def __init__(
        self,
        cell_type: str,
        source: AccessibilitySource = AccessibilitySource.ENCODE_ATAC,
        cache_dir: Optional[Path] = None,
        bigwig_path: Optional[Path] = None,
    ):
        self.cell_type = cell_type
        self.source = source
        self.cache_dir = cache_dir or Path.home() / ".phaselab" / "cache" / "chromatin"
        self.bigwig_path = bigwig_path
        self._bigwig = None
        self._peaks = None

    def _ensure_data(self) -> None:
        """Ensure BigWig data is loaded."""
        if self._bigwig is not None:
            return

        try:
            import pyBigWig
        except ImportError:
            logger.warning(
                "pyBigWig not installed. Install with: pip install phaselab[chromatin]"
            )
            return

        if self.source == AccessibilitySource.CUSTOM:
            if self.bigwig_path and self.bigwig_path.exists():
                self._bigwig = pyBigWig.open(str(self.bigwig_path))
        else:
            # Try to load from cache or download
            cached_path = self._get_cached_path()
            if cached_path and cached_path.exists():
                self._bigwig = pyBigWig.open(str(cached_path))
            else:
                logger.info(
                    f"Accessibility data for {self.cell_type} not cached. "
                    f"Use download_data() to fetch from ENCODE."
                )

    def _get_cached_path(self) -> Optional[Path]:
        """Get path to cached BigWig file."""
        if self.source == AccessibilitySource.ENCODE_ATAC:
            accession = self.ENCODE_ATAC_FILES.get(self.cell_type)
        elif self.source == AccessibilitySource.ENCODE_DNASE:
            accession = self.ENCODE_DNASE_FILES.get(self.cell_type)
        else:
            return None

        if accession:
            return self.cache_dir / f"{accession}.bigWig"
        return None

    def download_data(self, force: bool = False) -> bool:
        """
        Download accessibility data from ENCODE.

        Parameters
        ----------
        force : bool
            Re-download even if cached

        Returns
        -------
        bool
            True if download successful
        """
        cached_path = self._get_cached_path()
        if not cached_path:
            logger.error(f"No ENCODE accession for cell type: {self.cell_type}")
            return False

        if cached_path.exists() and not force:
            logger.info(f"Data already cached at {cached_path}")
            return True

        # Create cache directory
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        # Get accession
        if self.source == AccessibilitySource.ENCODE_ATAC:
            accession = self.ENCODE_ATAC_FILES.get(self.cell_type)
        else:
            accession = self.ENCODE_DNASE_FILES.get(self.cell_type)

        if not accession:
            return False

        # Download from ENCODE
        url = f"https://www.encodeproject.org/files/{accession}/@@download/{accession}.bigWig"
        logger.info(f"Downloading {url}...")

        try:
            import urllib.request
            urllib.request.urlretrieve(url, cached_path)
            logger.info(f"Downloaded to {cached_path}")
            return True
        except Exception as e:
            logger.error(f"Download failed: {e}")
            return False

    def get_accessibility(
        self,
        chrom: str,
        position: int,
        window: int = 200,
    ) -> AccessibilityScore:
        """
        Get accessibility score for a genomic position.

        Parameters
        ----------
        chrom : str
            Chromosome (e.g., "chr4")
        position : int
            Genomic position (0-based)
        window : int
            Window size around position to average

        Returns
        -------
        AccessibilityScore
            Accessibility information
        """
        self._ensure_data()

        # Default for when data unavailable
        if self._bigwig is None:
            return AccessibilityScore(
                score=0.5,  # Neutral default
                raw_signal=0.0,
                peak_overlap=False,
                source=self.source,
                cell_type=self.cell_type,
                confidence=0.0,  # Low confidence when no data
            )

        # Query BigWig
        start = max(0, position - window // 2)
        end = position + window // 2

        try:
            values = self._bigwig.values(chrom, start, end)
            values = np.array([v for v in values if v is not None and not np.isnan(v)])

            if len(values) == 0:
                raw_signal = 0.0
            else:
                raw_signal = float(np.mean(values))

            # Normalize to [0, 1] using sigmoid-like transform
            # Typical ATAC-seq signals range 0-100
            score = 1.0 / (1.0 + np.exp(-0.1 * (raw_signal - 10)))

            # Check peak overlap (simplified - would use peak BED file)
            peak_overlap = raw_signal > 5.0

            return AccessibilityScore(
                score=score,
                raw_signal=raw_signal,
                peak_overlap=peak_overlap,
                source=self.source,
                cell_type=self.cell_type,
                confidence=1.0,
            )

        except Exception as e:
            logger.warning(f"Error querying BigWig: {e}")
            return AccessibilityScore(
                score=0.5,
                raw_signal=0.0,
                peak_overlap=False,
                source=self.source,
                cell_type=self.cell_type,
                confidence=0.0,
            )

    def get_region_accessibility(
        self,
        chrom: str,
        start: int,
        end: int,
    ) -> np.ndarray:
        """
        Get accessibility signal across a region.

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
        np.ndarray
            Signal values for each position
        """
        self._ensure_data()

        if self._bigwig is None:
            return np.full(end - start, 0.5)

        try:
            values = self._bigwig.values(chrom, start, end)
            values = np.array(values, dtype=float)
            values = np.nan_to_num(values, nan=0.0)
            return values
        except Exception as e:
            logger.warning(f"Error querying region: {e}")
            return np.full(end - start, 0.5)

    def close(self) -> None:
        """Close BigWig file handle."""
        if self._bigwig is not None:
            self._bigwig.close()
            self._bigwig = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False
