"""
PhaseLab ATAC-seq Integration: Real ATAC-seq data loading and processing.

Provides:
- BigWig file loading via pyBigWig
- BED file peak loading
- ENCODE/4DN data integration
- Tissue-specific ATAC-seq profiles
- Real-time accessibility scoring from experimental data

References:
- ENCODE ATAC-seq pipeline: https://www.encodeproject.org/atac-seq/
- pyBigWig: https://github.com/deeptools/pyBigWig
- Buenrostro et al. (2015) ATAC-seq method

Version: 0.5.0
"""

import numpy as np
from typing import Dict, List, Optional, Tuple, Union, Any
from dataclasses import dataclass, field
from pathlib import Path
import logging
import json

logger = logging.getLogger(__name__)


# Try to import pyBigWig - gracefully degrade if not available
try:
    import pyBigWig
    HAS_PYBIGWIG = True
except ImportError:
    HAS_PYBIGWIG = False
    logger.warning(
        "pyBigWig not installed. Install with: pip install pyBigWig\n"
        "BigWig functionality will use fallback estimation."
    )


@dataclass
class ATACPeak:
    """ATAC-seq peak from BED file."""
    chrom: str
    start: int
    end: int
    name: str = ""
    score: float = 0.0
    strand: str = "."
    signal_value: float = 0.0
    p_value: float = -1.0
    q_value: float = -1.0
    peak_offset: int = -1  # Offset of peak summit from start

    @property
    def center(self) -> int:
        """Peak center position."""
        if self.peak_offset >= 0:
            return self.start + self.peak_offset
        return (self.start + self.end) // 2

    @property
    def width(self) -> int:
        """Peak width."""
        return self.end - self.start


@dataclass
class ATACProfile:
    """
    ATAC-seq accessibility profile for a tissue/cell type.

    Can be loaded from:
    - BigWig files (continuous signal)
    - BED/narrowPeak files (called peaks)
    - Pre-computed accessibility scores
    """

    name: str
    tissue: str
    cell_type: Optional[str] = None

    # Data sources
    bigwig_path: Optional[str] = None
    peaks_path: Optional[str] = None

    # Loaded data
    _bigwig: Any = field(default=None, repr=False)
    _peaks: Dict[str, List[ATACPeak]] = field(default_factory=dict, repr=False)
    _peak_index: Dict[str, np.ndarray] = field(default_factory=dict, repr=False)

    # Quality metrics
    frip_score: float = 0.0  # Fraction of reads in peaks
    tss_enrichment: float = 0.0  # TSS enrichment score

    # Metadata
    assembly: str = "hg38"
    data_source: str = "unknown"

    def __post_init__(self):
        """Load data if paths provided."""
        if self.bigwig_path and HAS_PYBIGWIG:
            self.load_bigwig(self.bigwig_path)
        if self.peaks_path:
            self.load_peaks(self.peaks_path)

    def load_bigwig(self, path: str) -> bool:
        """
        Load BigWig file for continuous signal access.

        Args:
            path: Path to BigWig file (local or URL)

        Returns:
            True if successful
        """
        if not HAS_PYBIGWIG:
            logger.warning("pyBigWig not available, cannot load BigWig")
            return False

        try:
            self._bigwig = pyBigWig.open(path)
            self.bigwig_path = path
            logger.info(f"Loaded BigWig: {path}")
            return True
        except Exception as e:
            logger.error(f"Failed to load BigWig {path}: {e}")
            return False

    def load_peaks(self, path: str) -> int:
        """
        Load peaks from BED/narrowPeak file.

        Supports:
        - BED format (3-12 columns)
        - narrowPeak format (10 columns)
        - broadPeak format (9 columns)

        Args:
            path: Path to peaks file

        Returns:
            Number of peaks loaded
        """
        peaks_loaded = 0
        self._peaks = {}

        try:
            with open(path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#') or line.startswith('track'):
                        continue

                    fields = line.split('\t')
                    if len(fields) < 3:
                        continue

                    peak = self._parse_peak_line(fields)
                    if peak:
                        if peak.chrom not in self._peaks:
                            self._peaks[peak.chrom] = []
                        self._peaks[peak.chrom].append(peak)
                        peaks_loaded += 1

            # Build index for fast lookup
            self._build_peak_index()
            self.peaks_path = path
            logger.info(f"Loaded {peaks_loaded} peaks from {path}")

        except Exception as e:
            logger.error(f"Failed to load peaks {path}: {e}")

        return peaks_loaded

    def _parse_peak_line(self, fields: List[str]) -> Optional[ATACPeak]:
        """Parse a single peak line."""
        try:
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            peak = ATACPeak(chrom=chrom, start=start, end=end)

            if len(fields) >= 4:
                peak.name = fields[3]
            if len(fields) >= 5:
                try:
                    peak.score = float(fields[4])
                except ValueError:
                    pass
            if len(fields) >= 6:
                peak.strand = fields[5]

            # narrowPeak format
            if len(fields) >= 7:
                try:
                    peak.signal_value = float(fields[6])
                except ValueError:
                    pass
            if len(fields) >= 8:
                try:
                    peak.p_value = float(fields[7])
                except ValueError:
                    pass
            if len(fields) >= 9:
                try:
                    peak.q_value = float(fields[8])
                except ValueError:
                    pass
            if len(fields) >= 10:
                try:
                    peak.peak_offset = int(fields[9])
                except ValueError:
                    pass

            return peak

        except Exception:
            return None

    def _build_peak_index(self):
        """Build sorted index for fast peak lookup."""
        self._peak_index = {}
        for chrom, peaks in self._peaks.items():
            # Sort by start position
            peaks.sort(key=lambda p: p.start)
            # Create numpy arrays for binary search
            starts = np.array([p.start for p in peaks])
            ends = np.array([p.end for p in peaks])
            self._peak_index[chrom] = (starts, ends, peaks)

    def get_signal(
        self,
        chrom: str,
        start: int,
        end: int,
        stat: str = "mean"
    ) -> float:
        """
        Get ATAC-seq signal in a genomic region from BigWig.

        Args:
            chrom: Chromosome
            start: Start position (0-based)
            end: End position (exclusive)
            stat: Statistic to compute ("mean", "max", "min", "sum")

        Returns:
            Signal value (0 if no data)
        """
        if self._bigwig is None:
            return self._estimate_signal(chrom, start, end)

        try:
            if stat == "mean":
                val = self._bigwig.stats(chrom, start, end, type="mean")
            elif stat == "max":
                val = self._bigwig.stats(chrom, start, end, type="max")
            elif stat == "min":
                val = self._bigwig.stats(chrom, start, end, type="min")
            elif stat == "sum":
                val = self._bigwig.stats(chrom, start, end, type="sum")
            else:
                val = self._bigwig.stats(chrom, start, end, type="mean")

            if val and val[0] is not None:
                return float(val[0])
            return 0.0

        except Exception as e:
            logger.debug(f"BigWig query failed: {e}")
            return 0.0

    def get_signal_array(
        self,
        chrom: str,
        start: int,
        end: int
    ) -> np.ndarray:
        """
        Get per-base ATAC-seq signal array from BigWig.

        Args:
            chrom: Chromosome
            start: Start position
            end: End position

        Returns:
            Numpy array of signal values
        """
        if self._bigwig is None:
            return np.zeros(end - start)

        try:
            values = self._bigwig.values(chrom, start, end)
            arr = np.array(values, dtype=float)
            arr = np.nan_to_num(arr, nan=0.0)
            return arr
        except Exception:
            return np.zeros(end - start)

    def find_overlapping_peaks(
        self,
        chrom: str,
        start: int,
        end: int
    ) -> List[ATACPeak]:
        """
        Find peaks overlapping a genomic region.

        Uses binary search for efficiency.

        Args:
            chrom: Chromosome
            start: Region start
            end: Region end

        Returns:
            List of overlapping peaks
        """
        if chrom not in self._peak_index:
            return []

        starts, ends, peaks = self._peak_index[chrom]

        # Binary search for potential overlaps
        # Find peaks where peak.end > start
        left_idx = np.searchsorted(ends, start, side='right')
        # Find peaks where peak.start < end
        right_idx = np.searchsorted(starts, end, side='left')

        overlapping = []
        for i in range(left_idx, min(right_idx + 1, len(peaks))):
            if i < len(peaks):
                peak = peaks[i]
                # Check actual overlap
                if peak.start < end and peak.end > start:
                    overlapping.append(peak)

        return overlapping

    def get_peak_score(
        self,
        chrom: str,
        position: int,
        window: int = 500
    ) -> float:
        """
        Get peak-based accessibility score for a position.

        Args:
            chrom: Chromosome
            position: Genomic position
            window: Window around position

        Returns:
            Accessibility score 0-1
        """
        start = position - window // 2
        end = position + window // 2

        peaks = self.find_overlapping_peaks(chrom, start, end)

        if not peaks:
            return 0.0

        # Score based on distance to peak center and peak strength
        scores = []
        for peak in peaks:
            distance = abs(position - peak.center)
            half_width = peak.width / 2

            # Gaussian-like decay from peak center
            if half_width > 0:
                distance_score = np.exp(-0.5 * (distance / half_width) ** 2)
            else:
                distance_score = 1.0 if distance < 50 else 0.5

            # Combine with peak signal strength
            if peak.signal_value > 0:
                # Normalize signal (typical range 0-100)
                signal_score = min(1.0, peak.signal_value / 50.0)
            else:
                signal_score = 0.5

            scores.append(0.7 * distance_score + 0.3 * signal_score)

        return float(np.max(scores))

    def score_position(
        self,
        chrom: str,
        position: int,
        window: int = 500,
        use_signal: bool = True,
        use_peaks: bool = True
    ) -> Dict[str, float]:
        """
        Compute comprehensive accessibility score for a position.

        Combines:
        - BigWig signal (continuous)
        - Peak overlap (discrete)
        - Position relative to peak summit

        Args:
            chrom: Chromosome
            position: Genomic position
            window: Window around position
            use_signal: Include BigWig signal
            use_peaks: Include peak data

        Returns:
            Dictionary with scores
        """
        start = position - window // 2
        end = position + window // 2

        result = {
            'signal_score': 0.0,
            'peak_score': 0.0,
            'combined_score': 0.0,
            'in_peak': False,
            'peak_distance': -1,
            'peak_count': 0,
        }

        # BigWig signal score
        if use_signal and self._bigwig is not None:
            raw_signal = self.get_signal(chrom, start, end, stat="mean")
            # Normalize to 0-1 (typical ATAC signal ranges 0-10 for normalized data)
            result['signal_score'] = min(1.0, raw_signal / 5.0)

        # Peak-based score
        if use_peaks and self._peaks:
            peaks = self.find_overlapping_peaks(chrom, start, end)
            result['peak_count'] = len(peaks)

            if peaks:
                result['in_peak'] = any(
                    p.start <= position < p.end for p in peaks
                )

                # Distance to nearest peak center
                distances = [abs(position - p.center) for p in peaks]
                result['peak_distance'] = min(distances)

                result['peak_score'] = self.get_peak_score(chrom, position, window)

        # Combined score
        if use_signal and use_peaks:
            result['combined_score'] = (
                0.5 * result['signal_score'] +
                0.5 * result['peak_score']
            )
        elif use_signal:
            result['combined_score'] = result['signal_score']
        elif use_peaks:
            result['combined_score'] = result['peak_score']

        return result

    def _estimate_signal(
        self,
        chrom: str,
        start: int,
        end: int
    ) -> float:
        """Estimate signal when no BigWig available."""
        # Fall back to peak data if available
        if self._peaks:
            return self.get_peak_score(chrom, (start + end) // 2, end - start)
        return 0.3  # Default moderate accessibility

    def close(self):
        """Close BigWig file handle."""
        if self._bigwig is not None:
            try:
                self._bigwig.close()
            except Exception:
                pass
            self._bigwig = None


# Pre-defined ENCODE profiles (URLs for common cell types)
ENCODE_ATAC_URLS = {
    # Human cell lines
    "K562": {
        "bigwig": "https://www.encodeproject.org/files/ENCFF814FUP/@@download/ENCFF814FUP.bigWig",
        "peaks": "https://www.encodeproject.org/files/ENCFF804DKG/@@download/ENCFF804DKG.bed.gz",
        "description": "K562 chronic myelogenous leukemia cell line",
    },
    "HepG2": {
        "bigwig": "https://www.encodeproject.org/files/ENCFF821QEY/@@download/ENCFF821QEY.bigWig",
        "peaks": "https://www.encodeproject.org/files/ENCFF324ELW/@@download/ENCFF324ELW.bed.gz",
        "description": "HepG2 hepatocellular carcinoma cell line",
    },
    "GM12878": {
        "bigwig": "https://www.encodeproject.org/files/ENCFF757FRW/@@download/ENCFF757FRW.bigWig",
        "peaks": "https://www.encodeproject.org/files/ENCFF559YHR/@@download/ENCFF559YHR.bed.gz",
        "description": "GM12878 lymphoblastoid cell line",
    },
    # Brain tissue (from PsychENCODE/ENCODE)
    "brain_cortex": {
        "description": "Adult brain cortex ATAC-seq",
        "data_source": "PsychENCODE",
    },
    "brain_hippocampus": {
        "description": "Adult brain hippocampus ATAC-seq",
        "data_source": "PsychENCODE",
    },
    # iPSC-derived neurons
    "iPS_neuron": {
        "description": "iPSC-derived neurons",
        "data_source": "ENCODE",
    },
}


def load_encode_atac(
    cell_type: str,
    cache_dir: Optional[str] = None
) -> ATACProfile:
    """
    Load ENCODE ATAC-seq profile for a cell type.

    Args:
        cell_type: ENCODE cell type (K562, HepG2, GM12878, etc.)
        cache_dir: Directory to cache downloaded files

    Returns:
        ATACProfile with ENCODE data

    Example:
        >>> profile = load_encode_atac("K562")
        >>> score = profile.score_position("chr17", 17681458, window=1000)
        >>> print(f"Accessibility: {score['combined_score']:.3f}")
    """
    if cell_type not in ENCODE_ATAC_URLS:
        available = list(ENCODE_ATAC_URLS.keys())
        raise ValueError(
            f"Unknown cell type '{cell_type}'. "
            f"Available: {available}"
        )

    info = ENCODE_ATAC_URLS[cell_type]

    profile = ATACProfile(
        name=f"ENCODE_{cell_type}",
        tissue=cell_type,
        cell_type=cell_type,
        data_source="ENCODE",
    )

    # Load BigWig if URL available
    if "bigwig" in info and HAS_PYBIGWIG:
        try:
            # pyBigWig can open URLs directly
            profile.load_bigwig(info["bigwig"])
        except Exception as e:
            logger.warning(f"Could not load BigWig URL: {e}")

    return profile


def create_tissue_profile(
    tissue: str,
    atac_bigwig: Optional[str] = None,
    atac_peaks: Optional[str] = None,
) -> ATACProfile:
    """
    Create tissue-specific ATAC profile from user data.

    Args:
        tissue: Tissue name
        atac_bigwig: Path to ATAC-seq BigWig file
        atac_peaks: Path to ATAC-seq peaks file

    Returns:
        ATACProfile configured for the tissue

    Example:
        >>> profile = create_tissue_profile(
        ...     tissue="brain_pvn",
        ...     atac_bigwig="/path/to/pvn_atac.bw",
        ...     atac_peaks="/path/to/pvn_peaks.bed"
        ... )
    """
    profile = ATACProfile(
        name=f"custom_{tissue}",
        tissue=tissue,
        data_source="user",
    )

    if atac_bigwig:
        profile.load_bigwig(atac_bigwig)
    if atac_peaks:
        profile.load_peaks(atac_peaks)

    return profile


def score_guide_accessibility(
    guide_chrom: str,
    guide_position: int,
    profile: ATACProfile,
    window: int = 1000
) -> Dict[str, Any]:
    """
    Score CRISPR guide accessibility using ATAC-seq data.

    This function provides the interface between CRISPR guide design
    and real ATAC-seq accessibility data.

    Args:
        guide_chrom: Chromosome of guide target
        guide_position: Genomic position of guide
        profile: ATACProfile with tissue-specific data
        window: Window around guide to assess

    Returns:
        Accessibility assessment dictionary

    Example:
        >>> profile = load_encode_atac("K562")
        >>> result = score_guide_accessibility(
        ...     "chr17", 17681458, profile
        ... )
        >>> print(f"GO/NO-GO: {result['recommendation']}")
    """
    scores = profile.score_position(guide_chrom, guide_position, window)

    # Classify accessibility
    combined = scores['combined_score']

    if combined >= 0.7:
        classification = "HIGHLY_ACCESSIBLE"
        recommendation = "GO"
    elif combined >= 0.4:
        classification = "MODERATELY_ACCESSIBLE"
        recommendation = "GO"
    elif combined >= 0.2:
        classification = "LOW_ACCESSIBILITY"
        recommendation = "CAUTION"
    else:
        classification = "CLOSED"
        recommendation = "NO-GO"

    return {
        'accessibility_score': combined,
        'signal_score': scores['signal_score'],
        'peak_score': scores['peak_score'],
        'in_peak': scores['in_peak'],
        'peak_count': scores['peak_count'],
        'classification': classification,
        'recommendation': recommendation,
        'tissue': profile.tissue,
        'data_source': profile.data_source,
    }


def compare_tissue_accessibility(
    chrom: str,
    position: int,
    tissues: List[str],
    window: int = 1000
) -> Dict[str, Dict[str, float]]:
    """
    Compare accessibility across multiple tissues.

    Useful for identifying tissue-specific expression patterns.

    Args:
        chrom: Chromosome
        position: Genomic position
        tissues: List of tissue/cell type names
        window: Window around position

    Returns:
        Dictionary mapping tissue to accessibility scores
    """
    results = {}

    for tissue in tissues:
        if tissue in ENCODE_ATAC_URLS:
            try:
                profile = load_encode_atac(tissue)
                scores = profile.score_position(chrom, position, window)
                results[tissue] = scores
                profile.close()
            except Exception as e:
                logger.warning(f"Could not score {tissue}: {e}")
                results[tissue] = {'combined_score': 0.0, 'error': str(e)}
        else:
            results[tissue] = {'combined_score': 0.0, 'error': 'Profile not available'}

    return results
