"""
PhaseLab SURF: Parser for CRISPR-SURF output files.

CRISPR-SURF outputs several file types:
- Beta profiles (deconvolved regulatory signal)
- Significant regions (peaks above threshold)
- QC metrics (deconvolution statistics)

This module parses these outputs into PhaseLab-compatible structures.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Union, Tuple
from pathlib import Path


@dataclass
class SURFRegion:
    """
    A significant region identified by CRISPR-SURF.

    Attributes:
        chrom: Chromosome
        start: Start position (genomic)
        end: End position (genomic)
        beta_mean: Mean beta value in region
        beta_max: Maximum beta value
        length: Region length in bp
        significance: Statistical significance (p-value or FDR)
        direction: "activating" or "repressing" based on beta sign
    """
    chrom: str
    start: int
    end: int
    beta_mean: float
    beta_max: float
    length: int
    significance: float
    direction: str
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def midpoint(self) -> int:
        """Region midpoint."""
        return (self.start + self.end) // 2

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'chrom': self.chrom,
            'start': self.start,
            'end': self.end,
            'beta_mean': self.beta_mean,
            'beta_max': self.beta_max,
            'length': self.length,
            'significance': self.significance,
            'direction': self.direction,
            **self.metadata,
        }


@dataclass
class SURFOutput:
    """
    Complete CRISPR-SURF output for a gene/locus.

    Attributes:
        gene_symbol: Gene name
        positions: Genomic positions (or TSS-relative)
        beta: Deconvolved beta values (regulatory signal)
        pvalues: Per-position p-values (if available)
        regions: List of significant regions
        qc_metrics: Quality control metrics from SURF
        raw_data: Original screen data (if retained)
        metadata: Additional metadata
    """
    gene_symbol: str
    positions: np.ndarray
    beta: np.ndarray
    pvalues: Optional[np.ndarray] = None
    regions: List[SURFRegion] = field(default_factory=list)
    qc_metrics: Dict[str, Any] = field(default_factory=dict)
    raw_data: Optional[np.ndarray] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_positions(self) -> int:
        """Number of positions."""
        return len(self.positions)

    @property
    def n_regions(self) -> int:
        """Number of significant regions."""
        return len(self.regions)

    @property
    def position_range(self) -> Tuple[int, int]:
        """Range of positions."""
        return int(self.positions.min()), int(self.positions.max())

    @property
    def activating_regions(self) -> List[SURFRegion]:
        """Regions with positive beta (activating)."""
        return [r for r in self.regions if r.direction == "activating"]

    @property
    def repressing_regions(self) -> List[SURFRegion]:
        """Regions with negative beta (repressing)."""
        return [r for r in self.regions if r.direction == "repressing"]

    def get_beta_at_position(self, pos: int) -> float:
        """Get beta value at nearest position."""
        idx = np.argmin(np.abs(self.positions - pos))
        return float(self.beta[idx])

    def get_region_at_position(self, pos: int) -> Optional[SURFRegion]:
        """Get region containing position (if any)."""
        for region in self.regions:
            if region.start <= pos <= region.end:
                return region
        return None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'gene_symbol': self.gene_symbol,
            'positions': self.positions.tolist(),
            'beta': self.beta.tolist(),
            'pvalues': self.pvalues.tolist() if self.pvalues is not None else None,
            'regions': [r.to_dict() for r in self.regions],
            'qc_metrics': self.qc_metrics,
            'metadata': self.metadata,
        }


def parse_surf_output(
    filepath: Union[str, Path],
    gene_symbol: str = "unknown",
    position_col: str = "position",
    beta_col: str = "beta",
    pvalue_col: Optional[str] = "pvalue",
    delimiter: str = "\t",
    tss_position: Optional[int] = None,
) -> SURFOutput:
    """
    Parse a CRISPR-SURF beta profile file.

    SURF typically outputs a TSV with columns:
    - position (or chr, start, end)
    - beta (deconvolved signal)
    - pvalue (optional)

    Args:
        filepath: Path to SURF output file.
        gene_symbol: Gene name for this analysis.
        position_col: Column name for position.
        beta_col: Column name for beta values.
        pvalue_col: Column name for p-values (optional).
        delimiter: Column delimiter.
        tss_position: TSS position for converting to relative coords.

    Returns:
        SURFOutput with parsed data.

    Example:
        >>> surf = parse_surf_output(
        ...     'cd69_surf_beta.tsv',
        ...     gene_symbol='CD69',
        ...     tss_position=9912000
        ... )
        >>> print(f"Beta range: {surf.beta.min():.3f} to {surf.beta.max():.3f}")
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"SURF output not found: {filepath}")

    try:
        import pandas as pd
        df = pd.read_csv(filepath, delimiter=delimiter)

        # Handle different position formats
        if position_col in df.columns:
            positions = df[position_col].values
        elif 'start' in df.columns:
            positions = df['start'].values
        elif 'pos' in df.columns:
            positions = df['pos'].values
        else:
            # Try to find any position-like column
            pos_cols = [c for c in df.columns if 'pos' in c.lower() or 'start' in c.lower()]
            if pos_cols:
                positions = df[pos_cols[0]].values
            else:
                raise ValueError(f"Could not find position column in {filepath}")

        # Get beta values
        if beta_col in df.columns:
            beta = df[beta_col].values
        elif 'signal' in df.columns:
            beta = df['signal'].values
        else:
            raise ValueError(f"Could not find beta column '{beta_col}' in {filepath}")

        # Get p-values if available
        pvalues = None
        if pvalue_col and pvalue_col in df.columns:
            pvalues = df[pvalue_col].values

        # Convert to TSS-relative if TSS provided
        if tss_position is not None:
            positions = positions - tss_position

        return SURFOutput(
            gene_symbol=gene_symbol,
            positions=positions.astype(float),
            beta=beta.astype(float),
            pvalues=pvalues,
            metadata={
                'source_file': str(filepath),
                'tss_position': tss_position,
            },
        )

    except ImportError:
        # Numpy fallback
        data = np.loadtxt(filepath, delimiter=delimiter, skiprows=1)
        positions = data[:, 0]
        beta = data[:, 1]
        pvalues = data[:, 2] if data.shape[1] > 2 else None

        if tss_position is not None:
            positions = positions - tss_position

        return SURFOutput(
            gene_symbol=gene_symbol,
            positions=positions,
            beta=beta,
            pvalues=pvalues,
            metadata={'source_file': str(filepath)},
        )


def parse_surf_regions(
    filepath: Union[str, Path],
    gene_symbol: str = "unknown",
    significance_threshold: float = 0.05,
) -> List[SURFRegion]:
    """
    Parse CRISPR-SURF significant regions file.

    SURF outputs a regions file with significant peaks above threshold.

    Args:
        filepath: Path to regions file.
        gene_symbol: Gene name.
        significance_threshold: FDR/p-value threshold for inclusion.

    Returns:
        List of SURFRegion objects.
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"SURF regions file not found: {filepath}")

    regions = []

    try:
        import pandas as pd
        df = pd.read_csv(filepath, delimiter='\t')

        for _, row in df.iterrows():
            # Determine significance column
            sig_col = None
            for col in ['FDR', 'fdr', 'pvalue', 'p_value', 'significance']:
                if col in df.columns:
                    sig_col = col
                    break

            significance = row[sig_col] if sig_col else 0.0

            if significance > significance_threshold:
                continue

            # Determine direction from beta
            beta_col = None
            for col in ['beta_mean', 'beta', 'signal', 'mean_beta']:
                if col in df.columns:
                    beta_col = col
                    break

            beta_mean = row[beta_col] if beta_col else 0.0
            direction = "activating" if beta_mean > 0 else "repressing"

            # Get coordinates
            chrom = row.get('chr', row.get('chrom', 'chr1'))
            start = int(row.get('start', row.get('chromStart', 0)))
            end = int(row.get('end', row.get('chromEnd', start + 100)))

            # Get max beta if available
            beta_max = row.get('beta_max', row.get('max_beta', abs(beta_mean)))

            region = SURFRegion(
                chrom=str(chrom),
                start=start,
                end=end,
                beta_mean=float(beta_mean),
                beta_max=float(beta_max),
                length=end - start,
                significance=float(significance),
                direction=direction,
                metadata={'gene': gene_symbol},
            )
            regions.append(region)

    except ImportError:
        # Basic parsing without pandas
        with open(filepath, 'r') as f:
            header = f.readline().strip().split('\t')
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 4:
                    region = SURFRegion(
                        chrom=fields[0],
                        start=int(fields[1]),
                        end=int(fields[2]),
                        beta_mean=float(fields[3]) if len(fields) > 3 else 0.0,
                        beta_max=float(fields[3]) if len(fields) > 3 else 0.0,
                        length=int(fields[2]) - int(fields[1]),
                        significance=float(fields[4]) if len(fields) > 4 else 0.0,
                        direction="activating" if float(fields[3]) > 0 else "repressing",
                    )
                    if region.significance <= significance_threshold:
                        regions.append(region)

    return regions


def load_surf_directory(
    directory: Union[str, Path],
    gene_symbol: str = "unknown",
    tss_position: Optional[int] = None,
) -> SURFOutput:
    """
    Load complete SURF output from a directory.

    CRISPR-SURF typically outputs multiple files:
    - *_beta.tsv: Beta profile
    - *_regions.tsv: Significant regions
    - *_qc.json: QC metrics

    This function loads all available files.

    Args:
        directory: Directory containing SURF output files.
        gene_symbol: Gene name.
        tss_position: TSS for coordinate conversion.

    Returns:
        SURFOutput with all available data.
    """
    directory = Path(directory)

    if not directory.exists():
        raise FileNotFoundError(f"SURF directory not found: {directory}")

    # Find beta profile file
    beta_files = list(directory.glob("*beta*.tsv")) + list(directory.glob("*beta*.txt"))
    if not beta_files:
        beta_files = list(directory.glob("*.tsv"))

    if not beta_files:
        raise FileNotFoundError(f"No SURF beta profile found in {directory}")

    # Parse main beta profile
    surf_output = parse_surf_output(
        beta_files[0],
        gene_symbol=gene_symbol,
        tss_position=tss_position,
    )

    # Try to load regions
    region_files = list(directory.glob("*region*.tsv")) + list(directory.glob("*peak*.tsv"))
    if region_files:
        regions = parse_surf_regions(region_files[0], gene_symbol=gene_symbol)
        surf_output.regions = regions

    # Try to load QC metrics
    qc_files = list(directory.glob("*qc*.json")) + list(directory.glob("*metrics*.json"))
    if qc_files:
        import json
        with open(qc_files[0], 'r') as f:
            surf_output.qc_metrics = json.load(f)

    surf_output.metadata['source_directory'] = str(directory)

    return surf_output


def convert_surf_to_tss_relative(
    surf: SURFOutput,
    tss_position: int,
) -> SURFOutput:
    """
    Convert genomic coordinates to TSS-relative.

    Args:
        surf: SURFOutput with genomic coordinates.
        tss_position: TSS position in same coordinate system.

    Returns:
        New SURFOutput with TSS-relative positions.
    """
    new_positions = surf.positions - tss_position

    new_regions = []
    for region in surf.regions:
        new_region = SURFRegion(
            chrom=region.chrom,
            start=region.start - tss_position,
            end=region.end - tss_position,
            beta_mean=region.beta_mean,
            beta_max=region.beta_max,
            length=region.length,
            significance=region.significance,
            direction=region.direction,
            metadata=region.metadata,
        )
        new_regions.append(new_region)

    return SURFOutput(
        gene_symbol=surf.gene_symbol,
        positions=new_positions,
        beta=surf.beta.copy(),
        pvalues=surf.pvalues.copy() if surf.pvalues is not None else None,
        regions=new_regions,
        qc_metrics=surf.qc_metrics,
        raw_data=surf.raw_data,
        metadata={
            **surf.metadata,
            'tss_position': tss_position,
            'coordinate_system': 'tss_relative',
        },
    )
