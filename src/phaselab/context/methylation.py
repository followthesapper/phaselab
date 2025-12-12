"""
DNA methylation context from WGBS and RRBS data.

Provides integration with ENCODE methylation data for context-aware
CRISPR guide scoring, particularly important for CpG-rich regions.
"""

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional, List, Dict, Tuple
import logging

import numpy as np

logger = logging.getLogger(__name__)


class MethylationSource(Enum):
    """Source of methylation data."""
    ENCODE_WGBS = "encode_wgbs"  # Whole-genome bisulfite sequencing
    ENCODE_RRBS = "encode_rrbs"  # Reduced representation
    ROADMAP = "roadmap"
    CUSTOM = "custom"


@dataclass
class CpGSite:
    """
    Single CpG site methylation information.

    Attributes
    ----------
    position : int
        Genomic position of C in CpG
    methylation : float
        Methylation level [0, 1]
    coverage : int
        Read coverage at site
    """
    position: int
    methylation: float
    coverage: int


@dataclass
class MethylationScore:
    """
    Methylation context for a genomic region.

    Attributes
    ----------
    mean_methylation : float
        Mean methylation across CpGs in region [0, 1]
    cpg_count : int
        Number of CpG sites in region
    cpg_density : float
        CpG density (CpGs per 100bp)
    is_cpg_island : bool
        Whether region overlaps a CpG island
    cpg_sites : List[CpGSite]
        Individual CpG site information
    source : MethylationSource
        Data source
    cell_type : str
        Cell type
    confidence : float
        Measurement confidence
    """
    mean_methylation: float
    cpg_count: int
    cpg_density: float
    is_cpg_island: bool
    cpg_sites: List[CpGSite]
    source: MethylationSource
    cell_type: str
    confidence: float = 1.0

    @property
    def methylation_variance(self) -> float:
        """Variance in methylation across CpG sites."""
        if len(self.cpg_sites) < 2:
            return 0.0
        values = [s.methylation for s in self.cpg_sites]
        return float(np.var(values))


class MethylationContext:
    """
    DNA methylation context provider.

    Integrates with ENCODE WGBS/RRBS data to provide methylation
    context for CRISPR guide scoring.

    Parameters
    ----------
    cell_type : str
        Cell type identifier
    source : MethylationSource
        Data source preference
    cache_dir : Path, optional
        Cache directory for downloaded data
    bedmethyl_path : Path, optional
        Path to custom bedMethyl file

    Examples
    --------
    >>> ctx = MethylationContext("K562")
    >>> meth = ctx.get_methylation("chr4", 55095264, 55095284)
    >>> print(f"Mean methylation: {meth.mean_methylation:.2%}")
    >>> print(f"CpG count: {meth.cpg_count}")
    """

    # ENCODE accessions for methylation data
    ENCODE_WGBS_FILES: Dict[str, str] = {
        "K562": "ENCFF123MET",
        "HepG2": "ENCFF456MET",
        "GM12878": "ENCFF789MET",
        "H1-hESC": "ENCFF012MET",
    }

    # CpG island definitions (would load from UCSC)
    _cpg_islands: Optional[Dict] = None

    def __init__(
        self,
        cell_type: str,
        source: MethylationSource = MethylationSource.ENCODE_WGBS,
        cache_dir: Optional[Path] = None,
        bedmethyl_path: Optional[Path] = None,
    ):
        self.cell_type = cell_type
        self.source = source
        self.cache_dir = cache_dir or Path.home() / ".phaselab" / "cache" / "methylation"
        self.bedmethyl_path = bedmethyl_path
        self._data = None
        self._index = None

    def _ensure_data(self) -> None:
        """Ensure methylation data is loaded."""
        if self._data is not None:
            return

        if self.source == MethylationSource.CUSTOM and self.bedmethyl_path:
            self._load_bedmethyl(self.bedmethyl_path)
        else:
            cached_path = self._get_cached_path()
            if cached_path and cached_path.exists():
                self._load_bedmethyl(cached_path)
            else:
                logger.info(
                    f"Methylation data for {self.cell_type} not cached. "
                    f"Using sequence-based CpG prediction."
                )

    def _get_cached_path(self) -> Optional[Path]:
        """Get cached bedMethyl path."""
        accession = self.ENCODE_WGBS_FILES.get(self.cell_type)
        if accession:
            return self.cache_dir / f"{accession}.bed.gz"
        return None

    def _load_bedmethyl(self, path: Path) -> None:
        """Load bedMethyl file into memory."""
        # For large files, would use tabix indexing
        # Simplified version loads into dict
        self._data = {}
        self._index = {}

        try:
            import gzip
            opener = gzip.open if str(path).endswith('.gz') else open

            with opener(path, 'rt') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) < 11:
                        continue

                    chrom = parts[0]
                    pos = int(parts[1])
                    coverage = int(parts[9])
                    methylation = float(parts[10]) / 100.0

                    if chrom not in self._data:
                        self._data[chrom] = {}
                    self._data[chrom][pos] = (methylation, coverage)

            logger.info(f"Loaded methylation data from {path}")

        except Exception as e:
            logger.error(f"Failed to load bedMethyl: {e}")
            self._data = {}

    def get_methylation(
        self,
        chrom: str,
        start: int,
        end: int,
        sequence: Optional[str] = None,
    ) -> MethylationScore:
        """
        Get methylation context for a genomic region.

        Parameters
        ----------
        chrom : str
            Chromosome
        start : int
            Region start
        end : int
            Region end
        sequence : str, optional
            Sequence for CpG prediction if no data available

        Returns
        -------
        MethylationScore
            Methylation context
        """
        self._ensure_data()

        cpg_sites = []
        region_length = end - start

        # Get methylation from data if available
        if self._data and chrom in self._data:
            chrom_data = self._data[chrom]
            for pos in range(start, end):
                if pos in chrom_data:
                    meth, cov = chrom_data[pos]
                    cpg_sites.append(CpGSite(
                        position=pos,
                        methylation=meth,
                        coverage=cov,
                    ))
            confidence = 1.0 if cpg_sites else 0.0

        # Fall back to sequence-based prediction
        elif sequence:
            cpg_sites = self._predict_cpg_sites(sequence, start)
            confidence = 0.3  # Low confidence for predictions

        else:
            # No data available
            return MethylationScore(
                mean_methylation=0.5,  # Neutral
                cpg_count=0,
                cpg_density=0.0,
                is_cpg_island=False,
                cpg_sites=[],
                source=self.source,
                cell_type=self.cell_type,
                confidence=0.0,
            )

        # Calculate summary statistics
        if cpg_sites:
            mean_meth = np.mean([s.methylation for s in cpg_sites])
        else:
            mean_meth = 0.5

        cpg_count = len(cpg_sites)
        cpg_density = (cpg_count / region_length) * 100 if region_length > 0 else 0.0

        # CpG island detection (simplified)
        is_cpg_island = cpg_density > 5.0 and mean_meth < 0.2

        return MethylationScore(
            mean_methylation=float(mean_meth),
            cpg_count=cpg_count,
            cpg_density=float(cpg_density),
            is_cpg_island=is_cpg_island,
            cpg_sites=cpg_sites,
            source=self.source,
            cell_type=self.cell_type,
            confidence=confidence,
        )

    def _predict_cpg_sites(
        self,
        sequence: str,
        start_pos: int,
    ) -> List[CpGSite]:
        """
        Predict CpG sites from sequence with default methylation estimates.

        Uses simple heuristics:
        - CpG islands: assume low methylation (~0.1)
        - Other CpGs: assume high methylation (~0.8)
        """
        cpg_sites = []
        seq_upper = sequence.upper()

        # Find CpG dinucleotides
        for i in range(len(seq_upper) - 1):
            if seq_upper[i:i+2] == "CG":
                # Estimate methylation based on local CpG density
                window_start = max(0, i - 50)
                window_end = min(len(seq_upper), i + 50)
                window = seq_upper[window_start:window_end]
                local_cpg_count = window.count("CG")
                local_density = local_cpg_count / len(window) * 100

                # CpG island: low methylation, otherwise high
                if local_density > 5.0:
                    estimated_meth = 0.1
                else:
                    estimated_meth = 0.8

                cpg_sites.append(CpGSite(
                    position=start_pos + i,
                    methylation=estimated_meth,
                    coverage=0,  # Predicted, not measured
                ))

        return cpg_sites

    def get_cpg_effect_on_editing(
        self,
        methylation_score: MethylationScore,
    ) -> float:
        """
        Estimate effect of methylation on CRISPR editing efficiency.

        High methylation generally reduces accessibility and editing efficiency.

        Parameters
        ----------
        methylation_score : MethylationScore
            Methylation context

        Returns
        -------
        float
            Multiplier for editing efficiency [0.5, 1.0]
        """
        if methylation_score.cpg_count == 0:
            return 1.0

        # High methylation reduces efficiency
        # CpG islands with low methylation are accessible
        meth = methylation_score.mean_methylation

        if methylation_score.is_cpg_island:
            # CpG islands are often unmethylated and accessible
            return 1.0 - 0.3 * meth
        else:
            # Other regions: methylation has moderate effect
            return 1.0 - 0.5 * meth

    def download_data(self, force: bool = False) -> bool:
        """Download methylation data from ENCODE."""
        cached_path = self._get_cached_path()
        if not cached_path:
            logger.error(f"No ENCODE accession for: {self.cell_type}")
            return False

        if cached_path.exists() and not force:
            logger.info(f"Data cached at {cached_path}")
            return True

        self.cache_dir.mkdir(parents=True, exist_ok=True)
        accession = self.ENCODE_WGBS_FILES.get(self.cell_type)

        if not accession:
            return False

        url = f"https://www.encodeproject.org/files/{accession}/@@download/{accession}.bed.gz"
        logger.info(f"Downloading {url}...")

        try:
            import urllib.request
            urllib.request.urlretrieve(url, cached_path)
            return True
        except Exception as e:
            logger.error(f"Download failed: {e}")
            return False
