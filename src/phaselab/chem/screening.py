"""
PhaseLab Chem: High-throughput screening analysis.

Analyzes HTS data to identify reliable hits:
- Compounds with coherent activity profiles
- Assay regions with stable readouts
- Reproducible vs artifactual hits

The key insight:
- Screening position (well, compound index) = perturbation coordinate
- Activity readout = response signal
- Spatial coherence identifies reliable signals vs noise

Applications:
- Drug discovery: Which hits are reproducible?
- Assay development: Which regions have reliable readouts?
- Quality control: Identifying systematic artifacts
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
class ScreeningLandscape:
    """
    High-throughput screening landscape.

    Represents screening data where:
    - Position can be well index, compound index, or concentration
    - Activity is the measured response

    Attributes:
        positions: Position indices or coordinates
        activities: Activity values (% inhibition, signal, etc.)
        compound_ids: Optional compound identifiers
        activity_type: Type of activity measurement
        assay_name: Assay identifier
        plate_id: Plate identifier
        metadata: Additional metadata
    """
    positions: np.ndarray
    activities: np.ndarray
    compound_ids: Optional[List[str]] = None
    activity_type: str = "inhibition"
    assay_name: str = "unknown"
    plate_id: str = "unknown"
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_compounds(self) -> int:
        """Number of compounds/wells."""
        return len(self.positions)

    @property
    def activity_range(self) -> Tuple[float, float]:
        """Activity range."""
        return float(np.nanmin(self.activities)), float(np.nanmax(self.activities))

    @property
    def mean_activity(self) -> float:
        """Mean activity."""
        return float(np.nanmean(self.activities))

    @property
    def activity_std(self) -> float:
        """Activity standard deviation."""
        return float(np.nanstd(self.activities))

    def to_response_landscape(self) -> ResponseLandscape:
        """Convert to generic ResponseLandscape."""
        return ResponseLandscape(
            coords=self.positions,
            responses=self.activities,
            coord_labels=self.compound_ids,
            metadata={
                'type': 'screening',
                'activity_type': self.activity_type,
                'assay': self.assay_name,
                'plate': self.plate_id,
                **self.metadata,
            },
        )


@dataclass
class ScreeningResult:
    """
    Result of screening coherence analysis.

    Attributes:
        landscape: Original screening landscape
        profile: Coherence profile
        reliable_regions: Regions with reliable readouts
        reliable_hits: Compounds flagged as reliable hits
        artifact_regions: Potential artifact regions
        validation: Validation statistics
    """
    landscape: ScreeningLandscape
    profile: CoherenceProfile
    reliable_regions: List[Dict[str, Any]]
    reliable_hits: List[Dict[str, Any]]
    artifact_regions: List[Dict[str, Any]]
    validation: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_validated(self) -> bool:
        """Whether coherence model is validated."""
        return self.profile.is_validated

    @property
    def n_reliable_hits(self) -> int:
        """Number of reliable hits."""
        return len(self.reliable_hits)

    @property
    def hit_rate(self) -> float:
        """Reliable hit rate."""
        return self.n_reliable_hits / self.landscape.n_compounds

    def summary(self) -> str:
        """Generate summary string."""
        lines = [
            "=" * 60,
            f"SCREENING ANALYSIS: {self.landscape.assay_name}",
            "=" * 60,
            "",
            f"Plate: {self.landscape.plate_id}",
            f"Activity type: {self.landscape.activity_type}",
            f"Compounds screened: {self.landscape.n_compounds}",
            "",
            "COHERENCE:",
            f"  Correlation: {self.profile.correlation:.3f}",
            f"  Validated: {'YES' if self.is_validated else 'NO'}",
            "",
            "RESULTS:",
            f"  Reliable hits: {self.n_reliable_hits}",
            f"  Hit rate: {self.hit_rate:.1%}",
            f"  Artifact regions: {len(self.artifact_regions)}",
        ]

        if self.reliable_hits:
            lines.extend(["", "TOP RELIABLE HITS:"])
            for hit in self.reliable_hits[:5]:
                compound_id = hit.get('compound_id', f"pos_{hit['position']}")
                lines.append(
                    f"  {compound_id} "
                    f"activity={hit['activity']:.3f} "
                    f"coherence={hit['coherence']:.3f}"
                )

        lines.append("=" * 60)
        return "\n".join(lines)


def load_screening_data(
    filepath: Union[str, Path],
    position_col: str = "well",
    activity_col: str = "activity",
    compound_col: Optional[str] = "compound_id",
    delimiter: str = "\t",
    assay_name: str = "unknown",
    plate_id: str = "unknown",
    activity_type: str = "inhibition",
) -> ScreeningLandscape:
    """
    Load HTS screening data from a file.

    Args:
        filepath: Path to data file.
        position_col: Position/well column name.
        activity_col: Activity column name.
        compound_col: Compound ID column name.
        delimiter: Column delimiter.
        assay_name: Assay identifier.
        plate_id: Plate identifier.
        activity_type: Type of activity.

    Returns:
        ScreeningLandscape for analysis.
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"Screening data not found: {filepath}")

    try:
        import pandas as pd
        df = pd.read_csv(filepath, delimiter=delimiter)

        # Handle well position (might be A1, B2 format or numeric)
        if position_col in df.columns:
            positions = df[position_col].values
            # Convert well names to numeric if needed
            if isinstance(positions[0], str):
                positions = np.arange(len(positions))
        else:
            positions = np.arange(len(df))

        activities = df[activity_col].values

        compound_ids = None
        if compound_col and compound_col in df.columns:
            compound_ids = df[compound_col].tolist()

        return ScreeningLandscape(
            positions=positions.astype(float),
            activities=activities.astype(float),
            compound_ids=compound_ids,
            activity_type=activity_type,
            assay_name=assay_name,
            plate_id=plate_id,
            metadata={'source_file': str(filepath)},
        )

    except ImportError:
        raise ImportError("pandas required for screening data loading")


def analyze_screening_coherence(
    landscape: ScreeningLandscape,
    window: int = 10,
    stable_threshold: float = 0.7,
    hit_threshold: float = 3.0,
) -> ScreeningResult:
    """
    Analyze spatial coherence in screening data.

    Identifies:
    - Reliable signal regions (coherent activity)
    - Artifact regions (incoherent, possibly edge effects)
    - Reliable hits (active + in coherent region)

    Args:
        landscape: Screening landscape to analyze.
        window: Window size for coherence.
        stable_threshold: Threshold for stable regions.
        hit_threshold: Activity threshold for hits (in std devs above mean).

    Returns:
        ScreeningResult with complete analysis.

    Example:
        >>> landscape = load_screening_data('plate_001.tsv')
        >>> result = analyze_screening_coherence(landscape)
        >>> reliable_hits = result.reliable_hits
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

    # Extract reliable and artifact regions
    reliable_regions = []
    artifact_regions = []

    for region in classification.regions:
        start, end, stability, score = region
        mask = (landscape.positions >= start) & (landscape.positions <= end)

        region_info = {
            'start': int(start),
            'end': int(end),
            'coherence': float(score),
            'mean_activity': float(np.nanmean(landscape.activities[mask])),
            'n_compounds': int(np.sum(mask)),
        }

        if stability == StabilityClass.STABLE:
            reliable_regions.append(region_info)
        elif stability == StabilityClass.AMPLIFYING:
            artifact_regions.append(region_info)

    # Identify reliable hits
    reliable_hits = identify_reliable_hits(
        landscape,
        profile,
        reliable_regions,
        hit_threshold=hit_threshold,
    )

    validation = {
        'coherence_variance_correlation': profile.correlation,
        'p_value': profile.p_value,
        'is_validated': profile.is_validated,
    }

    return ScreeningResult(
        landscape=landscape,
        profile=profile,
        reliable_regions=reliable_regions,
        reliable_hits=reliable_hits,
        artifact_regions=artifact_regions,
        validation=validation,
    )


def identify_reliable_hits(
    landscape: ScreeningLandscape,
    profile: CoherenceProfile,
    reliable_regions: List[Dict[str, Any]],
    hit_threshold: float = 3.0,
    coherence_threshold: float = 0.6,
) -> List[Dict[str, Any]]:
    """
    Identify reliable hits from screening data.

    A reliable hit must:
    1. Have activity above threshold (statistical hit)
    2. Be in a coherent region (reliable measurement)

    Args:
        landscape: Screening landscape.
        profile: Coherence profile.
        reliable_regions: List of reliable region dictionaries.
        hit_threshold: Activity threshold in std devs.
        coherence_threshold: Minimum local coherence.

    Returns:
        List of reliable hit dictionaries.
    """
    # Calculate hit threshold
    mean_act = landscape.mean_activity
    std_act = landscape.activity_std
    activity_cutoff = mean_act + hit_threshold * std_act

    hits = []

    for i, (pos, activity) in enumerate(zip(landscape.positions, landscape.activities)):
        # Check if activity is above threshold
        if activity < activity_cutoff:
            continue

        # Find local coherence
        coh_idx = np.argmin(np.abs(profile.coords - pos))
        local_coherence = float(profile.coherence[coh_idx])

        # Check if in reliable region
        in_reliable_region = any(
            region['start'] <= pos <= region['end']
            for region in reliable_regions
        )

        if local_coherence < coherence_threshold:
            continue

        hit = {
            'position': int(pos),
            'activity': float(activity),
            'z_score': float((activity - mean_act) / std_act),
            'coherence': local_coherence,
            'in_reliable_region': in_reliable_region,
            'confidence': (
                'high' if in_reliable_region and local_coherence > 0.8
                else 'medium' if in_reliable_region or local_coherence > 0.7
                else 'low'
            ),
        }

        if landscape.compound_ids and i < len(landscape.compound_ids):
            hit['compound_id'] = landscape.compound_ids[i]

        hits.append(hit)

    # Sort by activity
    hits.sort(key=lambda x: x['activity'], reverse=True)

    return hits
