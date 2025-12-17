"""
PhaseLab Chem: Binding landscape analysis.

Analyzes protein-ligand and protein-protein binding data to identify:
- Stable binding modes (reliable affinity predictions)
- Binding hot spots (residues with consistent effects)
- Coherent regions in binding site

The key insight:
- Mutations/modifications at different positions = perturbation coordinates
- Binding affinity (Kd, Ki, IC50) = response signal
- Spatial coherence identifies where binding is predictable

Applications:
- Drug design: Which modifications reliably improve affinity?
- Protein engineering: Which mutations predictably affect binding?
- Docking validation: Which poses show coherent scoring?
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
class BindingLandscape:
    """
    Binding affinity landscape.

    Represents how binding affinity varies across:
    - Mutation positions (alanine scanning, saturation mutagenesis)
    - Ligand modifications (SAR)
    - Docking poses

    Attributes:
        positions: Position coordinates (residue number, modification site, etc.)
        affinities: Binding affinity values (Kd, Ki, IC50, or scores)
        affinity_type: Type of affinity measurement
        target: Target protein/receptor name
        ligand: Ligand/partner name
        mutations: Optional list of mutation labels
        metadata: Additional metadata
    """
    positions: np.ndarray
    affinities: np.ndarray
    affinity_type: str = "Kd"
    target: str = "unknown"
    ligand: str = "unknown"
    mutations: Optional[List[str]] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_positions(self) -> int:
        """Number of positions."""
        return len(self.positions)

    @property
    def affinity_range(self) -> Tuple[float, float]:
        """Range of affinities."""
        return float(np.nanmin(self.affinities)), float(np.nanmax(self.affinities))

    @property
    def mean_affinity(self) -> float:
        """Mean affinity."""
        return float(np.nanmean(self.affinities))

    def to_response_landscape(self) -> ResponseLandscape:
        """Convert to generic ResponseLandscape."""
        return ResponseLandscape(
            coords=self.positions,
            responses=self.affinities,
            coord_labels=self.mutations,
            metadata={
                'type': 'binding',
                'affinity_type': self.affinity_type,
                'target': self.target,
                'ligand': self.ligand,
                **self.metadata,
            },
        )


@dataclass
class BindingResult:
    """
    Result of binding coherence analysis.

    Attributes:
        landscape: Original binding landscape
        profile: Coherence profile
        stable_regions: Stable binding regions
        hot_spots: Identified binding hot spots
        validation: Validation statistics
    """
    landscape: BindingLandscape
    profile: CoherenceProfile
    stable_regions: List[Dict[str, Any]]
    hot_spots: List[Dict[str, Any]]
    validation: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_validated(self) -> bool:
        """Whether coherence model is validated."""
        return self.profile.is_validated

    @property
    def n_hot_spots(self) -> int:
        """Number of hot spots identified."""
        return len(self.hot_spots)

    def summary(self) -> str:
        """Generate summary string."""
        lines = [
            "=" * 60,
            f"BINDING ANALYSIS: {self.landscape.target} + {self.landscape.ligand}",
            "=" * 60,
            "",
            f"Affinity type: {self.landscape.affinity_type}",
            f"Positions analyzed: {self.landscape.n_positions}",
            f"Affinity range: {self.landscape.affinity_range}",
            "",
            "COHERENCE:",
            f"  Correlation: {self.profile.correlation:.3f}",
            f"  Validated: {'YES' if self.is_validated else 'NO'}",
            f"  Variance reduction: {self.profile.variance_reduction_estimate:.1%}",
            "",
            f"Stable regions: {len(self.stable_regions)}",
            f"Hot spots: {self.n_hot_spots}",
        ]

        if self.hot_spots:
            lines.extend(["", "TOP HOT SPOTS:"])
            for hs in self.hot_spots[:5]:
                lines.append(
                    f"  Position {hs['position']}: "
                    f"effect={hs['effect']:.3f}, coherence={hs['coherence']:.3f}"
                )

        lines.append("=" * 60)
        return "\n".join(lines)


def load_binding_data(
    filepath: Union[str, Path],
    position_col: str = "position",
    affinity_col: str = "affinity",
    mutation_col: Optional[str] = "mutation",
    delimiter: str = "\t",
    target: str = "unknown",
    ligand: str = "unknown",
    affinity_type: str = "Kd",
    log_transform: bool = False,
) -> BindingLandscape:
    """
    Load binding data from a file.

    Supports common binding data formats:
    - Alanine scanning results
    - Saturation mutagenesis data
    - Structure-activity relationship (SAR) data

    Args:
        filepath: Path to data file.
        position_col: Position column name.
        affinity_col: Affinity column name.
        mutation_col: Mutation label column name.
        delimiter: Column delimiter.
        target: Target name.
        ligand: Ligand name.
        affinity_type: Type of affinity (Kd, Ki, IC50, score).
        log_transform: Whether to log-transform affinities.

    Returns:
        BindingLandscape for analysis.

    Example:
        >>> landscape = load_binding_data(
        ...     'alanine_scan.tsv',
        ...     target='BCL2',
        ...     ligand='BH3_peptide',
        ...     affinity_type='Kd',
        ... )
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"Binding data not found: {filepath}")

    try:
        import pandas as pd
        df = pd.read_csv(filepath, delimiter=delimiter)

        positions = df[position_col].values
        affinities = df[affinity_col].values

        if log_transform and np.all(affinities > 0):
            affinities = np.log10(affinities)

        mutations = None
        if mutation_col and mutation_col in df.columns:
            mutations = df[mutation_col].tolist()

        return BindingLandscape(
            positions=positions.astype(float),
            affinities=affinities.astype(float),
            affinity_type=affinity_type,
            target=target,
            ligand=ligand,
            mutations=mutations,
            metadata={
                'source_file': str(filepath),
                'log_transformed': log_transform,
            },
        )

    except ImportError:
        raise ImportError("pandas required for binding data loading")


def analyze_binding_coherence(
    landscape: BindingLandscape,
    window: int = 10,
    stable_threshold: float = 0.7,
) -> BindingResult:
    """
    Analyze spatial coherence in binding landscape.

    Identifies:
    - Stable binding regions (predictable affinity effects)
    - Binding hot spots (key residues/modifications)
    - Coherent structural elements

    Args:
        landscape: Binding landscape to analyze.
        window: Window size for coherence.
        stable_threshold: Threshold for stable regions.

    Returns:
        BindingResult with complete analysis.

    Example:
        >>> result = analyze_binding_coherence(landscape)
        >>> hot_spots = result.hot_spots
        >>> stable = result.stable_regions
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
                'mean_affinity': float(np.nanmean(landscape.affinities[mask])),
                'n_positions': int(np.sum(mask)),
            })

    # Identify hot spots
    hot_spots = hot_spot_analysis(landscape, profile)

    validation = {
        'coherence_variance_correlation': profile.correlation,
        'p_value': profile.p_value,
        'is_validated': profile.is_validated,
    }

    return BindingResult(
        landscape=landscape,
        profile=profile,
        stable_regions=stable_regions,
        hot_spots=hot_spots,
        validation=validation,
    )


def hot_spot_analysis(
    landscape: BindingLandscape,
    profile: CoherenceProfile,
    effect_threshold: float = 1.0,
    coherence_threshold: float = 0.7,
) -> List[Dict[str, Any]]:
    """
    Identify binding hot spots.

    Hot spots are positions that:
    1. Have large affinity effects (strong perturbation)
    2. Are in coherent regions (reliable measurement)

    Args:
        landscape: Binding landscape.
        profile: Coherence profile.
        effect_threshold: Minimum affinity effect magnitude.
        coherence_threshold: Minimum local coherence.

    Returns:
        List of hot spot dictionaries.
    """
    # Calculate effect relative to mean
    mean_affinity = landscape.mean_affinity
    effects = landscape.affinities - mean_affinity

    hot_spots = []

    for i, (pos, effect) in enumerate(zip(landscape.positions, effects)):
        # Find local coherence (handle empty profile edge case)
        if len(profile.coords) == 0:
            local_coherence = 0.0
        else:
            coh_idx = np.argmin(np.abs(profile.coords - pos))
            local_coherence = float(profile.coherence[coh_idx])

        if abs(effect) < effect_threshold:
            continue
        if local_coherence < coherence_threshold:
            continue

        hot_spot = {
            'position': int(pos),
            'effect': float(effect),
            'affinity': float(landscape.affinities[i]),
            'coherence': local_coherence,
            'is_stabilizing': effect < 0,  # Lower affinity = tighter binding
            'confidence': 'high' if local_coherence > 0.8 else 'medium',
        }

        if landscape.mutations and i < len(landscape.mutations):
            hot_spot['mutation'] = landscape.mutations[i]

        hot_spots.append(hot_spot)

    # Sort by effect magnitude
    hot_spots.sort(key=lambda x: abs(x['effect']), reverse=True)

    return hot_spots


def identify_stable_binding_modes(
    result: BindingResult,
    min_coherence: float = 0.7,
    min_size: int = 3,
) -> List[Dict[str, Any]]:
    """
    Identify stable binding modes from coherence analysis.

    Stable binding modes are contiguous regions where:
    1. Affinity effects are coherent (spatially consistent)
    2. Predictions are reliable (low variance)

    Args:
        result: BindingResult from analyze_binding_coherence.
        min_coherence: Minimum coherence score.
        min_size: Minimum number of positions.

    Returns:
        List of stable binding mode dictionaries.
    """
    modes = []

    for region in result.stable_regions:
        if region['coherence'] < min_coherence:
            continue
        if region['n_positions'] < min_size:
            continue

        # Count hot spots in this region
        hot_spots_in_region = [
            hs for hs in result.hot_spots
            if region['start'] <= hs['position'] <= region['end']
        ]

        mode = {
            **region,
            'n_hot_spots': len(hot_spots_in_region),
            'hot_spots': hot_spots_in_region,
            'recommendation': (
                "HIGH CONFIDENCE: Stable binding mode with predictable effects"
                if region['coherence'] > 0.8
                else "MEDIUM CONFIDENCE: Binding mode identified"
            ),
        }
        modes.append(mode)

    # Sort by coherence
    modes.sort(key=lambda x: x['coherence'], reverse=True)

    return modes
