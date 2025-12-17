"""
PhaseLab Landscapes: Core data structures.

This module defines the fundamental data structures for perturbation-response
coherence analysis. These structures are domain-agnostic and can represent
data from genomics, microbiology, chemistry, or any other perturbation science.

The key insight from E212-E213:
    The probe is not the structure - the response manifold is.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Union, Tuple
from enum import Enum
import json


class StabilityClass(Enum):
    """
    Classification of regions based on perturbation response stability.

    Based on E213 spatial coherence analysis:
    - STABLE: High coherence, low variance - safe to perturb, predictable outcomes
    - MIXED: Moderate coherence - context-dependent, use with caution
    - AMPLIFYING: High variance despite perturbation - MYC-like, avoid
    - IRRELEVANT: Low signal, no meaningful response detected

    The classification is determined by the coherence-variance relationship,
    not by arbitrary thresholds.
    """
    STABLE = "stable"
    MIXED = "mixed"
    AMPLIFYING = "amplifying"
    IRRELEVANT = "irrelevant"

    @property
    def is_safe(self) -> bool:
        """Whether this classification indicates safe perturbation."""
        return self == StabilityClass.STABLE

    @property
    def recommendation(self) -> str:
        """Human-readable recommendation for this class."""
        recommendations = {
            StabilityClass.STABLE: "PRIMARY - Safe to perturb, predictable outcomes expected",
            StabilityClass.MIXED: "SECONDARY - Context-dependent, validate before use",
            StabilityClass.AMPLIFYING: "AVOID - High variance amplification risk (MYC-like)",
            StabilityClass.IRRELEVANT: "SKIP - No meaningful response detected",
        }
        return recommendations[self]


@dataclass
class ResponseLandscape:
    """
    A structured representation of perturbation-response data.

    This is the core data structure for all PhaseLab landscape analysis.
    It represents a mapping from coordinates (positions, conditions, etc.)
    to responses (expression changes, fitness effects, binding energies, etc.).

    Attributes:
        coords: Array of coordinates where perturbations occurred.
                Shape: (n_positions,) or (n_positions, n_dims) for multi-dimensional.
        responses: Array of response values at each coordinate.
                   Shape: (n_positions,) for single measurement,
                          (n_positions, n_replicates) for replicated data,
                          (n_positions, n_conditions) for multi-condition.
        coord_labels: Optional labels for coordinates (e.g., gene names, positions).
        response_labels: Optional labels for response dimensions.
        metadata: Additional metadata about the landscape.

    Example:
        >>> # CRISPRa tiling data
        >>> landscape = ResponseLandscape(
        ...     coords=np.array([-500, -450, -400, ...]),  # TSS-relative positions
        ...     responses=expression_changes,  # logFC values
        ...     coord_labels=['pos_' + str(p) for p in positions],
        ...     metadata={'gene': 'RAI1', 'modality': 'CRISPRa'}
        ... )

        >>> # Drug binding landscape
        >>> landscape = ResponseLandscape(
        ...     coords=residue_indices,
        ...     responses=binding_energies,
        ...     metadata={'protein': 'PTEN', 'domain': 'pocket_coherence'}
        ... )
    """
    coords: np.ndarray
    responses: np.ndarray
    coord_labels: Optional[List[str]] = None
    response_labels: Optional[List[str]] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Validate and normalize inputs."""
        self.coords = np.asarray(self.coords)
        self.responses = np.asarray(self.responses)

        # Ensure coords is 1D or 2D
        if self.coords.ndim == 0:
            self.coords = self.coords.reshape(1)
        elif self.coords.ndim > 2:
            raise ValueError(f"coords must be 1D or 2D, got {self.coords.ndim}D")

        # Ensure responses matches coords length
        if len(self.responses) != len(self.coords):
            raise ValueError(
                f"coords length ({len(self.coords)}) must match "
                f"responses length ({len(self.responses)})"
            )

    @property
    def n_positions(self) -> int:
        """Number of positions/coordinates in the landscape."""
        return len(self.coords)

    @property
    def n_replicates(self) -> int:
        """Number of replicates (1 if single measurement)."""
        if self.responses.ndim == 1:
            return 1
        return self.responses.shape[1]

    @property
    def has_replicates(self) -> bool:
        """Whether the landscape has replicate measurements."""
        return self.n_replicates > 1

    @property
    def mean_response(self) -> np.ndarray:
        """Mean response across replicates (or raw response if no replicates)."""
        if self.responses.ndim == 1:
            return self.responses
        return np.nanmean(self.responses, axis=1)

    @property
    def response_variance(self) -> np.ndarray:
        """Variance of response across replicates."""
        if self.responses.ndim == 1:
            return np.zeros(len(self.responses))
        return np.nanvar(self.responses, axis=1)

    @property
    def response_std(self) -> np.ndarray:
        """Standard deviation of response across replicates."""
        return np.sqrt(self.response_variance)

    @property
    def coord_range(self) -> Tuple[float, float]:
        """Range of coordinates (min, max)."""
        if self.coords.ndim == 1:
            return float(np.min(self.coords)), float(np.max(self.coords))
        return float(np.min(self.coords[:, 0])), float(np.max(self.coords[:, 0]))

    @property
    def response_range(self) -> Tuple[float, float]:
        """Range of mean responses (min, max)."""
        mean = self.mean_response
        return float(np.nanmin(mean)), float(np.nanmax(mean))

    def subset(self, indices: np.ndarray) -> 'ResponseLandscape':
        """Create a subset landscape from given indices."""
        new_coords = self.coords[indices]
        new_responses = self.responses[indices]
        new_labels = None
        if self.coord_labels is not None:
            new_labels = [self.coord_labels[i] for i in indices]

        return ResponseLandscape(
            coords=new_coords,
            responses=new_responses,
            coord_labels=new_labels,
            response_labels=self.response_labels,
            metadata={**self.metadata, 'subset_of': 'parent'},
        )

    def window(self, start: float, end: float) -> 'ResponseLandscape':
        """Extract a window of the landscape by coordinate range."""
        if self.coords.ndim == 1:
            mask = (self.coords >= start) & (self.coords <= end)
        else:
            mask = (self.coords[:, 0] >= start) & (self.coords[:, 0] <= end)

        indices = np.where(mask)[0]
        return self.subset(indices)

    def sort_by_coords(self) -> 'ResponseLandscape':
        """Return a new landscape sorted by coordinates."""
        if self.coords.ndim == 1:
            order = np.argsort(self.coords)
        else:
            order = np.lexsort(self.coords.T[::-1])

        return self.subset(order)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            'coords': self.coords.tolist(),
            'responses': self.responses.tolist(),
            'coord_labels': self.coord_labels,
            'response_labels': self.response_labels,
            'metadata': self.metadata,
            'n_positions': self.n_positions,
            'n_replicates': self.n_replicates,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'ResponseLandscape':
        """Create from dictionary."""
        return cls(
            coords=np.array(data['coords']),
            responses=np.array(data['responses']),
            coord_labels=data.get('coord_labels'),
            response_labels=data.get('response_labels'),
            metadata=data.get('metadata', {}),
        )


@dataclass
class CoherenceProfile:
    """
    Spatial coherence profile across a landscape.

    This represents the output of sliding-window coherence analysis,
    following the E213 methodology.

    Attributes:
        coords: Center coordinates of each window.
        coherence: R̄ value at each position (spatial coherence).
        local_variance: Variance of responses within each window.
        window_size: Size of sliding window used.
        correlation: Pearson correlation between coherence and variance.
        p_value: Statistical significance of the correlation.
        variance_reduction_estimate: Expected variance reduction in high-coherence regions.
    """
    coords: np.ndarray
    coherence: np.ndarray
    local_variance: np.ndarray
    window_size: int
    correlation: float
    p_value: float
    variance_reduction_estimate: float

    def __post_init__(self):
        """Validate inputs."""
        self.coords = np.asarray(self.coords)
        self.coherence = np.asarray(self.coherence)
        self.local_variance = np.asarray(self.local_variance)

    @property
    def n_positions(self) -> int:
        """Number of positions in the profile."""
        return len(self.coords)

    @property
    def mean_coherence(self) -> float:
        """Mean coherence across all positions."""
        return float(np.nanmean(self.coherence))

    @property
    def mean_variance(self) -> float:
        """Mean local variance across all positions."""
        return float(np.nanmean(self.local_variance))

    @property
    def coherence_range(self) -> Tuple[float, float]:
        """Range of coherence values."""
        return float(np.nanmin(self.coherence)), float(np.nanmax(self.coherence))

    @property
    def is_validated(self) -> bool:
        """
        Whether the coherence-variance relationship is validated.

        E213 showed that spatial coherence should NEGATIVELY correlate
        with variance (r < 0). If correlation is positive, the methodology
        may not apply to this system (e.g., amplifying super-enhancers).
        """
        return self.correlation < 0 and self.p_value < 0.05

    @property
    def interpretation(self) -> str:
        """Human-readable interpretation of the profile."""
        if self.correlation < -0.2 and self.p_value < 0.05:
            return (
                f"VALIDATED: Negative correlation (r={self.correlation:.3f}, p={self.p_value:.2e}) "
                f"supports spatial coherence → stability relationship. "
                f"Expected variance reduction: {self.variance_reduction_estimate:.0%}"
            )
        elif self.correlation > 0.2 and self.p_value < 0.05:
            return (
                f"WARNING: Positive correlation (r={self.correlation:.3f}) indicates "
                f"AMPLIFYING behavior. This region may be like MYC - perturbations "
                f"amplify rather than stabilize. Exercise extreme caution."
            )
        else:
            return (
                f"INCONCLUSIVE: Weak or non-significant correlation (r={self.correlation:.3f}, "
                f"p={self.p_value:.2e}). Spatial coherence may not be predictive for this system."
            )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            'coords': self.coords.tolist(),
            'coherence': self.coherence.tolist(),
            'local_variance': self.local_variance.tolist(),
            'window_size': self.window_size,
            'correlation': self.correlation,
            'p_value': self.p_value,
            'variance_reduction_estimate': self.variance_reduction_estimate,
            'mean_coherence': self.mean_coherence,
            'is_validated': self.is_validated,
        }


@dataclass
class RegionClassification:
    """
    Classification of landscape regions by stability.

    This is the primary output of landscape analysis - regions classified
    as STABLE (safe to perturb), MIXED (context-dependent), AMPLIFYING
    (dangerous), or IRRELEVANT (no signal).

    Attributes:
        regions: List of (start, end, class, score) tuples.
        landscape: The source ResponseLandscape.
        profile: The CoherenceProfile used for classification.
        classification_params: Parameters used for classification.
    """
    regions: List[Tuple[float, float, StabilityClass, float]]
    landscape: ResponseLandscape
    profile: CoherenceProfile
    classification_params: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_regions(self) -> int:
        """Total number of classified regions."""
        return len(self.regions)

    @property
    def stable_regions(self) -> List[Tuple[float, float, float]]:
        """Regions classified as STABLE (start, end, score)."""
        return [
            (start, end, score)
            for start, end, cls, score in self.regions
            if cls == StabilityClass.STABLE
        ]

    @property
    def mixed_regions(self) -> List[Tuple[float, float, float]]:
        """Regions classified as MIXED (start, end, score)."""
        return [
            (start, end, score)
            for start, end, cls, score in self.regions
            if cls == StabilityClass.MIXED
        ]

    @property
    def amplifying_regions(self) -> List[Tuple[float, float, float]]:
        """Regions classified as AMPLIFYING (start, end, score)."""
        return [
            (start, end, score)
            for start, end, cls, score in self.regions
            if cls == StabilityClass.AMPLIFYING
        ]

    @property
    def irrelevant_regions(self) -> List[Tuple[float, float, float]]:
        """Regions classified as IRRELEVANT (start, end, score)."""
        return [
            (start, end, score)
            for start, end, cls, score in self.regions
            if cls == StabilityClass.IRRELEVANT
        ]

    @property
    def n_stable(self) -> int:
        """Number of stable regions."""
        return len(self.stable_regions)

    @property
    def n_amplifying(self) -> int:
        """Number of amplifying regions."""
        return len(self.amplifying_regions)

    @property
    def fraction_stable(self) -> float:
        """Fraction of total coordinate range that is stable."""
        if not self.stable_regions:
            return 0.0

        stable_length = sum(end - start for start, end, _ in self.stable_regions)
        total_length = self.landscape.coord_range[1] - self.landscape.coord_range[0]
        return stable_length / total_length if total_length > 0 else 0.0

    def get_class_at(self, coord: float) -> StabilityClass:
        """Get the stability class at a specific coordinate."""
        for start, end, cls, _ in self.regions:
            if start <= coord <= end:
                return cls
        return StabilityClass.IRRELEVANT

    def summary(self) -> str:
        """Generate a human-readable summary."""
        lines = [
            "=" * 60,
            "REGION CLASSIFICATION SUMMARY",
            "=" * 60,
            f"Total regions: {self.n_regions}",
            f"  STABLE: {self.n_stable} ({self.fraction_stable:.1%} of range)",
            f"  MIXED: {len(self.mixed_regions)}",
            f"  AMPLIFYING: {self.n_amplifying}",
            f"  IRRELEVANT: {len(self.irrelevant_regions)}",
            "",
            f"Coherence-variance correlation: r={self.profile.correlation:.3f}",
            f"Expected variance reduction: {self.profile.variance_reduction_estimate:.1%}",
            "",
        ]

        if self.stable_regions:
            lines.append("STABLE REGIONS (recommended for perturbation):")
            for start, end, score in self.stable_regions[:5]:
                lines.append(f"  [{start:.0f}, {end:.0f}] score={score:.3f}")
            if len(self.stable_regions) > 5:
                lines.append(f"  ... and {len(self.stable_regions) - 5} more")

        if self.amplifying_regions:
            lines.append("")
            lines.append("AMPLIFYING REGIONS (AVOID - MYC-like behavior):")
            for start, end, score in self.amplifying_regions[:3]:
                lines.append(f"  [{start:.0f}, {end:.0f}] score={score:.3f}")

        lines.append("=" * 60)
        return "\n".join(lines)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            'regions': [
                {'start': s, 'end': e, 'class': c.value, 'score': sc}
                for s, e, c, sc in self.regions
            ],
            'n_regions': self.n_regions,
            'n_stable': self.n_stable,
            'n_amplifying': self.n_amplifying,
            'fraction_stable': self.fraction_stable,
            'profile_summary': {
                'correlation': self.profile.correlation,
                'p_value': self.profile.p_value,
                'variance_reduction': self.profile.variance_reduction_estimate,
            },
            'classification_params': self.classification_params,
        }


@dataclass
class LandscapeMetrics:
    """
    Summary metrics for a response landscape.

    Provides quick statistics and quality indicators for a landscape
    before detailed analysis.
    """
    n_positions: int
    n_replicates: int
    coord_range: Tuple[float, float]
    response_range: Tuple[float, float]
    mean_response: float
    response_variance: float
    coverage: float  # Fraction of coordinate range with data
    signal_to_noise: float  # Mean / std ratio

    @classmethod
    def from_landscape(cls, landscape: ResponseLandscape) -> 'LandscapeMetrics':
        """Compute metrics from a landscape."""
        mean_resp = landscape.mean_response
        variance = np.nanvar(mean_resp)
        std = np.sqrt(variance)

        # Signal-to-noise ratio
        mean_val = np.nanmean(np.abs(mean_resp))
        snr = mean_val / std if std > 0 else float('inf')

        # Coverage (positions with non-NaN responses)
        n_valid = np.sum(~np.isnan(mean_resp))
        coverage = n_valid / landscape.n_positions

        return cls(
            n_positions=landscape.n_positions,
            n_replicates=landscape.n_replicates,
            coord_range=landscape.coord_range,
            response_range=landscape.response_range,
            mean_response=float(np.nanmean(mean_resp)),
            response_variance=float(variance),
            coverage=coverage,
            signal_to_noise=snr,
        )

    def is_suitable_for_analysis(self) -> Tuple[bool, str]:
        """Check if the landscape is suitable for coherence analysis."""
        issues = []

        if self.n_positions < 20:
            issues.append(f"Too few positions ({self.n_positions} < 20)")

        if self.coverage < 0.5:
            issues.append(f"Low coverage ({self.coverage:.1%} < 50%)")

        if self.signal_to_noise < 1.0:
            issues.append(f"Low signal-to-noise ({self.signal_to_noise:.2f} < 1.0)")

        if issues:
            return False, "; ".join(issues)
        return True, "Landscape suitable for analysis"

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'n_positions': self.n_positions,
            'n_replicates': self.n_replicates,
            'coord_range': self.coord_range,
            'response_range': self.response_range,
            'mean_response': self.mean_response,
            'response_variance': self.response_variance,
            'coverage': self.coverage,
            'signal_to_noise': self.signal_to_noise,
        }
