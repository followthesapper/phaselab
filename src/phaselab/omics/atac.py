"""
PhaseLab Omics: ATAC-seq coherence analysis.

ATAC-seq measures chromatin accessibility. PhaseLab identifies:
- Stable accessibility regions (reliable open chromatin)
- Variable regions (condition-dependent accessibility)
- Integration with CRISPR targeting

Key insight:
- ATAC signal across genomic positions = response landscape
- Spatial coherence identifies reliably accessible regions
- Coherence predicts where CRISPRa/i will work consistently
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Union, Tuple
from pathlib import Path

from ..landscapes.core import (
    ResponseLandscape,
    CoherenceProfile,
    StabilityClass,
)
from ..landscapes.coherence import compute_spatial_coherence
from ..landscapes.classification import classify_regions


@dataclass
class ATACLandscape:
    """
    ATAC-seq signal landscape.

    Attributes:
        positions: Genomic positions
        signal: ATAC-seq signal (normalized counts, RPKM, etc.)
        chromosome: Chromosome
        gene_symbol: Associated gene (if gene-centric)
        cell_type: Cell type/line
        condition: Experimental condition
        replicates: Optional replicate signals
        metadata: Additional metadata
    """
    positions: np.ndarray
    signal: np.ndarray
    chromosome: str = "unknown"
    gene_symbol: Optional[str] = None
    cell_type: str = "unknown"
    condition: str = "default"
    replicates: Optional[np.ndarray] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_positions(self) -> int:
        """Number of positions."""
        return len(self.positions)

    @property
    def position_range(self) -> Tuple[int, int]:
        """Genomic position range."""
        return int(self.positions.min()), int(self.positions.max())

    @property
    def mean_signal(self) -> float:
        """Mean ATAC signal."""
        return float(np.nanmean(self.signal))

    @property
    def max_signal(self) -> float:
        """Maximum ATAC signal."""
        return float(np.nanmax(self.signal))

    def to_response_landscape(self) -> ResponseLandscape:
        """Convert to generic ResponseLandscape."""
        return ResponseLandscape(
            coords=self.positions,
            responses=self.signal if self.replicates is None else self.replicates,
            metadata={
                'type': 'atac',
                'chromosome': self.chromosome,
                'gene': self.gene_symbol,
                'cell_type': self.cell_type,
                'condition': self.condition,
                **self.metadata,
            },
        )


@dataclass
class ATACResult:
    """
    Result of ATAC coherence analysis.

    Attributes:
        landscape: Original ATAC landscape
        profile: Coherence profile
        stable_regions: Stably accessible regions
        peak_regions: High-signal peak regions
        accessible_stable: Regions that are both accessible and stable
        validation: Validation statistics
    """
    landscape: ATACLandscape
    profile: CoherenceProfile
    stable_regions: List[Dict[str, Any]]
    peak_regions: List[Dict[str, Any]]
    accessible_stable: List[Dict[str, Any]]
    validation: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_validated(self) -> bool:
        """Whether coherence model is validated."""
        return self.profile.is_validated

    @property
    def n_stable_peaks(self) -> int:
        """Number of stable accessible regions."""
        return len(self.accessible_stable)

    def summary(self) -> str:
        """Generate summary string."""
        lines = [
            "=" * 60,
            f"ATAC-seq ANALYSIS: {self.landscape.gene_symbol or self.landscape.chromosome}",
            "=" * 60,
            "",
            f"Cell type: {self.landscape.cell_type}",
            f"Condition: {self.landscape.condition}",
            f"Positions: {self.landscape.n_positions}",
            f"Position range: {self.landscape.position_range}",
            "",
            "COHERENCE:",
            f"  Correlation: {self.profile.correlation:.3f}",
            f"  Validated: {'YES' if self.is_validated else 'NO'}",
            "",
            "REGIONS:",
            f"  Total peaks: {len(self.peak_regions)}",
            f"  Stable regions: {len(self.stable_regions)}",
            f"  Stable accessible: {self.n_stable_peaks}",
        ]

        if self.accessible_stable:
            lines.extend(["", "TOP STABLE ACCESSIBLE REGIONS:"])
            for region in self.accessible_stable[:5]:
                lines.append(
                    f"  [{region['start']}-{region['end']}] "
                    f"signal={region['mean_signal']:.2f} "
                    f"coh={region['coherence']:.3f}"
                )

        lines.append("=" * 60)
        return "\n".join(lines)


def load_atac_data(
    filepath: Union[str, Path],
    position_col: str = "position",
    signal_col: str = "signal",
    replicate_cols: Optional[List[str]] = None,
    delimiter: str = "\t",
    chromosome: str = "unknown",
    gene_symbol: Optional[str] = None,
    cell_type: str = "unknown",
    condition: str = "default",
) -> ATACLandscape:
    """
    Load ATAC-seq data from a file.

    Supports common ATAC output formats:
    - BigWig converted to TSV
    - BedGraph
    - Peak counts

    Args:
        filepath: Path to data file.
        position_col: Position column name.
        signal_col: Signal column name.
        replicate_cols: Replicate column names.
        delimiter: Column delimiter.
        chromosome: Chromosome identifier.
        gene_symbol: Associated gene.
        cell_type: Cell type.
        condition: Experimental condition.

    Returns:
        ATACLandscape for analysis.
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"ATAC file not found: {filepath}")

    try:
        import pandas as pd
        df = pd.read_csv(filepath, delimiter=delimiter)

        positions = df[position_col].values
        signal = df[signal_col].values

        replicates = None
        if replicate_cols:
            replicates = df[replicate_cols].values

        return ATACLandscape(
            positions=positions.astype(float),
            signal=signal.astype(float),
            chromosome=chromosome,
            gene_symbol=gene_symbol,
            cell_type=cell_type,
            condition=condition,
            replicates=replicates,
            metadata={'source_file': str(filepath)},
        )

    except ImportError:
        raise ImportError("pandas required for ATAC loading")


def analyze_atac_coherence(
    landscape: ATACLandscape,
    window: int = 50,
    stable_threshold: float = 0.7,
    peak_threshold: float = 2.0,
) -> ATACResult:
    """
    Analyze spatial coherence in ATAC-seq data.

    Identifies:
    - Stable accessibility regions (coherent signal)
    - Peak regions (high signal)
    - Stable accessible regions (both)

    Args:
        landscape: ATAC landscape to analyze.
        window: Window size for coherence.
        stable_threshold: Coherence threshold for stable regions.
        peak_threshold: Signal threshold for peaks (in std above mean).

    Returns:
        ATACResult with complete analysis.

    Example:
        >>> landscape = load_atac_data('cd69_atac.tsv', gene_symbol='CD69')
        >>> result = analyze_atac_coherence(landscape)
        >>> stable_peaks = result.accessible_stable
    """
    response_landscape = landscape.to_response_landscape()

    # Compute coherence
    profile = compute_spatial_coherence(
        response_landscape,
        window=window,
    )

    # Classify regions
    classification = classify_regions(
        response_landscape,
        profile=profile,
        stable_threshold=stable_threshold,
    )

    # Extract stable regions
    stable_regions = []
    for region in classification.regions:
        start, end, stability, score = region
        if stability == StabilityClass.STABLE:
            mask = (landscape.positions >= start) & (landscape.positions <= end)
            stable_regions.append({
                'start': int(start),
                'end': int(end),
                'coherence': float(score),
                'mean_signal': float(np.nanmean(landscape.signal[mask])),
                'n_positions': int(np.sum(mask)),
            })

    # Identify peak regions (high signal)
    mean_signal = landscape.mean_signal
    std_signal = float(np.nanstd(landscape.signal))
    signal_cutoff = mean_signal + peak_threshold * std_signal

    peak_regions = identify_stable_peaks(
        landscape, profile, signal_cutoff
    )

    # Find stable accessible regions (intersection)
    accessible_stable = []
    for region in stable_regions:
        mask = (landscape.positions >= region['start']) & (landscape.positions <= region['end'])
        region_signal = landscape.signal[mask]

        if np.nanmean(region_signal) > signal_cutoff:
            accessible_stable.append({
                **region,
                'status': 'STABLE_ACCESSIBLE',
            })

    # Sort by signal
    accessible_stable.sort(key=lambda x: x['mean_signal'], reverse=True)

    validation = {
        'coherence_variance_correlation': profile.correlation,
        'p_value': profile.p_value,
        'is_validated': profile.is_validated,
    }

    return ATACResult(
        landscape=landscape,
        profile=profile,
        stable_regions=stable_regions,
        peak_regions=peak_regions,
        accessible_stable=accessible_stable,
        validation=validation,
    )


def identify_stable_peaks(
    landscape: ATACLandscape,
    profile: CoherenceProfile,
    signal_cutoff: float,
    coherence_threshold: float = 0.7,
) -> List[Dict[str, Any]]:
    """
    Identify peaks that are both high-signal and stable.

    Args:
        landscape: ATAC landscape.
        profile: Coherence profile.
        signal_cutoff: Signal threshold for peaks.
        coherence_threshold: Minimum coherence.

    Returns:
        List of stable peak dictionaries.
    """
    peaks = []

    # Find contiguous high-signal regions
    high_signal = landscape.signal > signal_cutoff
    changes = np.diff(high_signal.astype(int))
    starts = np.where(changes == 1)[0] + 1
    ends = np.where(changes == -1)[0] + 1

    # Handle edge cases
    if high_signal[0]:
        starts = np.concatenate([[0], starts])
    if high_signal[-1]:
        ends = np.concatenate([ends, [len(high_signal)]])

    for start_idx, end_idx in zip(starts, ends):
        if end_idx <= start_idx:
            continue

        start_pos = landscape.positions[start_idx]
        end_pos = landscape.positions[end_idx - 1]

        # Get coherence in this region
        coh_mask = (profile.coords >= start_pos) & (profile.coords <= end_pos)
        if not np.any(coh_mask):
            continue

        region_coherence = float(np.nanmean(profile.coherence[coh_mask]))

        if region_coherence < coherence_threshold:
            continue

        mask = (landscape.positions >= start_pos) & (landscape.positions <= end_pos)
        peaks.append({
            'start': int(start_pos),
            'end': int(end_pos),
            'coherence': region_coherence,
            'mean_signal': float(np.nanmean(landscape.signal[mask])),
            'max_signal': float(np.nanmax(landscape.signal[mask])),
            'n_positions': int(end_idx - start_idx),
            'stable': region_coherence >= coherence_threshold,
        })

    return peaks


def compare_atac_conditions(
    landscape_a: ATACLandscape,
    landscape_b: ATACLandscape,
    window: int = 50,
) -> Dict[str, Any]:
    """
    Compare ATAC coherence between two conditions.

    Useful for:
    - Comparing treatment vs control
    - Comparing cell types
    - Comparing time points

    Args:
        landscape_a: First ATAC landscape.
        landscape_b: Second ATAC landscape.
        window: Coherence window size.

    Returns:
        Comparison results.
    """
    result_a = analyze_atac_coherence(landscape_a, window=window)
    result_b = analyze_atac_coherence(landscape_b, window=window)

    # Find common positions
    common_positions = np.intersect1d(landscape_a.positions, landscape_b.positions)

    mask_a = np.isin(landscape_a.positions, common_positions)
    mask_b = np.isin(landscape_b.positions, common_positions)

    if len(common_positions) > 2:
        signal_correlation = np.corrcoef(
            landscape_a.signal[mask_a],
            landscape_b.signal[mask_b],
        )[0, 1]
    else:
        signal_correlation = np.nan

    return {
        'condition_a': landscape_a.condition,
        'condition_b': landscape_b.condition,
        'comparison': {
            'n_common_positions': len(common_positions),
            'signal_correlation': float(signal_correlation),
            'coherence_a': result_a.profile.correlation,
            'coherence_b': result_b.profile.correlation,
            'stable_peaks_a': result_a.n_stable_peaks,
            'stable_peaks_b': result_b.n_stable_peaks,
        },
        'differential_accessibility': {
            'regions_more_accessible_a': [
                r for r in result_a.accessible_stable
                if not any(
                    r['start'] <= s['end'] and r['end'] >= s['start']
                    for s in result_b.accessible_stable
                )
            ],
            'regions_more_accessible_b': [
                r for r in result_b.accessible_stable
                if not any(
                    r['start'] <= s['end'] and r['end'] >= s['start']
                    for s in result_a.accessible_stable
                )
            ],
        },
    }
