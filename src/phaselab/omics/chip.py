"""
PhaseLab Omics: ChIP-seq coherence analysis.

ChIP-seq measures protein-DNA binding (TFs, histones). PhaseLab identifies:
- Stable binding sites (reliable ChIP signal)
- Variable binding (condition-dependent)
- High-confidence peak regions

Key insight:
- ChIP signal across positions = response landscape
- Spatial coherence identifies reproducible binding
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
class ChIPLandscape:
    """
    ChIP-seq signal landscape.

    Attributes:
        positions: Genomic positions
        signal: ChIP signal (fold enrichment, RPKM, etc.)
        chromosome: Chromosome
        target: ChIP target (TF name, histone mark)
        gene_symbol: Associated gene
        cell_type: Cell type
        condition: Experimental condition
        metadata: Additional metadata
    """
    positions: np.ndarray
    signal: np.ndarray
    chromosome: str = "unknown"
    target: str = "unknown"
    gene_symbol: Optional[str] = None
    cell_type: str = "unknown"
    condition: str = "default"
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_positions(self) -> int:
        return len(self.positions)

    @property
    def mean_signal(self) -> float:
        return float(np.nanmean(self.signal))

    def to_response_landscape(self) -> ResponseLandscape:
        return ResponseLandscape(
            coords=self.positions,
            responses=self.signal,
            metadata={
                'type': 'chip',
                'target': self.target,
                'gene': self.gene_symbol,
                'cell_type': self.cell_type,
                **self.metadata,
            },
        )


@dataclass
class ChIPResult:
    """
    Result of ChIP coherence analysis.

    Attributes:
        landscape: Original ChIP landscape
        profile: Coherence profile
        stable_binding: Stable binding regions
        peak_regions: High-signal peaks
        validation: Validation statistics
    """
    landscape: ChIPLandscape
    profile: CoherenceProfile
    stable_binding: List[Dict[str, Any]]
    peak_regions: List[Dict[str, Any]]
    validation: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_validated(self) -> bool:
        return self.profile.is_validated

    @property
    def n_stable_peaks(self) -> int:
        return len(self.stable_binding)

    def summary(self) -> str:
        lines = [
            "=" * 60,
            f"ChIP-seq ANALYSIS: {self.landscape.target}",
            "=" * 60,
            "",
            f"Gene: {self.landscape.gene_symbol or 'N/A'}",
            f"Cell type: {self.landscape.cell_type}",
            f"Positions: {self.landscape.n_positions}",
            "",
            "COHERENCE:",
            f"  Correlation: {self.profile.correlation:.3f}",
            f"  Validated: {'YES' if self.is_validated else 'NO'}",
            "",
            f"Stable binding sites: {self.n_stable_peaks}",
            f"Total peaks: {len(self.peak_regions)}",
            "=" * 60,
        ]
        return "\n".join(lines)


def load_chip_data(
    filepath: Union[str, Path],
    position_col: str = "position",
    signal_col: str = "signal",
    delimiter: str = "\t",
    target: str = "unknown",
    gene_symbol: Optional[str] = None,
    cell_type: str = "unknown",
    condition: str = "default",
) -> ChIPLandscape:
    """Load ChIP-seq data from a file."""
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"ChIP file not found: {filepath}")

    try:
        import pandas as pd
        df = pd.read_csv(filepath, delimiter=delimiter)

        positions = df[position_col].values
        signal = df[signal_col].values

        return ChIPLandscape(
            positions=positions.astype(float),
            signal=signal.astype(float),
            target=target,
            gene_symbol=gene_symbol,
            cell_type=cell_type,
            condition=condition,
            metadata={'source_file': str(filepath)},
        )

    except ImportError:
        raise ImportError("pandas required for ChIP loading")


def analyze_chip_coherence(
    landscape: ChIPLandscape,
    window: int = 50,
    stable_threshold: float = 0.7,
    peak_threshold: float = 2.0,
) -> ChIPResult:
    """
    Analyze spatial coherence in ChIP-seq data.

    Args:
        landscape: ChIP landscape to analyze.
        window: Window size for coherence.
        stable_threshold: Coherence threshold for stable binding.
        peak_threshold: Signal threshold for peaks (std above mean).

    Returns:
        ChIPResult with analysis.
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

    # Extract stable binding regions
    stable_binding = []
    for region in classification.regions:
        start, end, stability, score = region
        if stability == StabilityClass.STABLE:
            mask = (landscape.positions >= start) & (landscape.positions <= end)
            mean_signal = float(np.nanmean(landscape.signal[mask]))

            # Only include if signal is above background
            if mean_signal > landscape.mean_signal:
                stable_binding.append({
                    'start': int(start),
                    'end': int(end),
                    'coherence': float(score),
                    'mean_signal': mean_signal,
                    'n_positions': int(np.sum(mask)),
                })

    # Identify peak regions
    peak_threshold_value = landscape.mean_signal + peak_threshold * float(np.nanstd(landscape.signal))
    peak_regions = identify_stable_binding(landscape, profile, peak_threshold_value)

    validation = {
        'coherence_variance_correlation': profile.correlation,
        'p_value': profile.p_value,
        'is_validated': profile.is_validated,
    }

    return ChIPResult(
        landscape=landscape,
        profile=profile,
        stable_binding=stable_binding,
        peak_regions=peak_regions,
        validation=validation,
    )


def identify_stable_binding(
    landscape: ChIPLandscape,
    profile: CoherenceProfile,
    signal_cutoff: float,
    coherence_threshold: float = 0.7,
) -> List[Dict[str, Any]]:
    """Identify stable binding sites above signal threshold."""
    peaks = []

    high_signal = landscape.signal > signal_cutoff
    changes = np.diff(high_signal.astype(int))
    starts = np.where(changes == 1)[0] + 1
    ends = np.where(changes == -1)[0] + 1

    if high_signal[0]:
        starts = np.concatenate([[0], starts])
    if high_signal[-1]:
        ends = np.concatenate([ends, [len(high_signal)]])

    for start_idx, end_idx in zip(starts, ends):
        if end_idx <= start_idx:
            continue

        start_pos = landscape.positions[start_idx]
        end_pos = landscape.positions[end_idx - 1]

        coh_mask = (profile.coords >= start_pos) & (profile.coords <= end_pos)
        if not np.any(coh_mask):
            continue

        region_coherence = float(np.nanmean(profile.coherence[coh_mask]))

        mask = (landscape.positions >= start_pos) & (landscape.positions <= end_pos)
        peaks.append({
            'start': int(start_pos),
            'end': int(end_pos),
            'coherence': region_coherence,
            'mean_signal': float(np.nanmean(landscape.signal[mask])),
            'stable': region_coherence >= coherence_threshold,
        })

    return peaks
