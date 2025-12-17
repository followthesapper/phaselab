"""
PhaseLab SURF: Spatial coherence analysis of CRISPR-SURF output.

This module applies the E213-validated spatial coherence methodology
to SURF-deconvolved regulatory signals.

Key insight:
- SURF deconvolves raw screen data → cleaner regulatory signal
- Spatial coherence identifies stable regions in that signal
- Combination provides highest-confidence targeting zones
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Tuple

from .parser import SURFOutput, SURFRegion
from ..spatial.regulatory import RegulatoryLandscape, RegulatoryRegion, classify_regulatory_regions
from ..landscapes.core import CoherenceProfile, StabilityClass
from ..landscapes.coherence import compute_spatial_coherence


@dataclass
class SURFCoherenceResult:
    """
    Result of coherence analysis on SURF output.

    Attributes:
        surf: Original SURF output
        landscape: Converted regulatory landscape
        profile: Computed coherence profile
        regions: Classified regulatory regions
        stable_surf_regions: SURF regions that overlap stable coherence zones
        validation: Validation statistics
    """
    surf: SURFOutput
    landscape: RegulatoryLandscape
    profile: CoherenceProfile
    regions: List[RegulatoryRegion]
    stable_surf_regions: List[SURFRegion] = field(default_factory=list)
    validation: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_stable_regions(self) -> int:
        """Number of stable regions."""
        return sum(1 for r in self.regions if r.is_safe)

    @property
    def correlation(self) -> float:
        """Coherence-variance correlation."""
        return self.profile.correlation

    @property
    def is_validated(self) -> bool:
        """Whether coherence model is validated for this locus."""
        return self.profile.is_validated

    def get_best_targeting_zones(self, n: int = 5) -> List[RegulatoryRegion]:
        """
        Get best regions for targeting, prioritizing:
        1. Stable coherence regions
        2. Overlap with SURF significant peaks
        3. Higher beta values (stronger regulatory signal)

        Args:
            n: Number of regions to return.

        Returns:
            Top n targeting zones.
        """
        stable = [r for r in self.regions if r.is_safe]

        # Score by coherence and variance reduction
        def score_region(region):
            base_score = region.coherence_score + region.variance_reduction

            # Bonus for overlap with SURF significant regions
            for surf_region in self.surf.regions:
                if self._regions_overlap(region, surf_region):
                    base_score += 0.5 * abs(surf_region.beta_mean)

            return base_score

        stable_scored = sorted(stable, key=score_region, reverse=True)
        return stable_scored[:n]

    def _regions_overlap(
        self,
        coherence_region: RegulatoryRegion,
        surf_region: SURFRegion,
    ) -> bool:
        """Check if two regions overlap."""
        return not (
            coherence_region.end < surf_region.start or
            coherence_region.start > surf_region.end
        )

    def summary(self) -> str:
        """Generate summary string."""
        lines = [
            "=" * 60,
            f"SURF COHERENCE ANALYSIS: {self.surf.gene_symbol}",
            "=" * 60,
            "",
            "SURF DECONVOLUTION:",
            f"  Positions analyzed: {self.surf.n_positions}",
            f"  Beta range: [{self.surf.beta.min():.3f}, {self.surf.beta.max():.3f}]",
            f"  SURF significant regions: {self.surf.n_regions}",
            "",
            "SPATIAL COHERENCE:",
            f"  Coherence-variance correlation: {self.correlation:.3f}",
            f"  Validated: {'YES' if self.is_validated else 'NO'}",
            f"  Variance reduction estimate: {self.profile.variance_reduction_estimate:.1%}",
            "",
            "REGION CLASSIFICATION:",
            f"  Stable regions: {self.n_stable_regions}",
            f"  Total regions: {len(self.regions)}",
        ]

        if self.stable_surf_regions:
            lines.extend([
                "",
                "STABLE SURF REGIONS (highest confidence):",
            ])
            for region in self.stable_surf_regions[:5]:
                lines.append(
                    f"  [{region.start}, {region.end}] beta={region.beta_mean:.3f}"
                )

        lines.append("=" * 60)
        return "\n".join(lines)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'gene_symbol': self.surf.gene_symbol,
            'surf': self.surf.to_dict(),
            'profile': self.profile.to_dict(),
            'regions': [r.to_dict() for r in self.regions],
            'stable_surf_regions': [r.to_dict() for r in self.stable_surf_regions],
            'validation': self.validation,
            'summary_stats': {
                'correlation': self.correlation,
                'is_validated': self.is_validated,
                'n_stable_regions': self.n_stable_regions,
            },
        }


def surf_to_regulatory_landscape(
    surf: SURFOutput,
    modality: str = "CRISPRa",
) -> RegulatoryLandscape:
    """
    Convert SURF output to a RegulatoryLandscape for coherence analysis.

    The SURF beta profile becomes the response signal for spatial
    coherence analysis.

    Args:
        surf: SURF output to convert.
        modality: Screen modality ("CRISPRa" or "CRISPRi").

    Returns:
        RegulatoryLandscape for coherence analysis.
    """
    return RegulatoryLandscape(
        positions=surf.positions,
        responses=surf.beta,
        gene_symbol=surf.gene_symbol,
        modality=modality,
        metadata={
            'source': 'CRISPR-SURF',
            'n_surf_regions': surf.n_regions,
            **surf.metadata,
        },
    )


def compute_surf_coherence(
    surf: SURFOutput,
    window: int = 50,
    stable_threshold: float = 0.7,
    modality: str = "CRISPRa",
) -> SURFCoherenceResult:
    """
    Compute spatial coherence on SURF-deconvolved data.

    This is the main function for combining SURF deconvolution with
    PhaseLab spatial coherence analysis.

    Args:
        surf: SURF output with beta profile.
        window: Window size for coherence computation.
        stable_threshold: Threshold for stable region classification.
        modality: Screen modality.

    Returns:
        SURFCoherenceResult with complete analysis.

    Example:
        >>> surf = parse_surf_output('cd69_surf.tsv', gene_symbol='CD69')
        >>> result = compute_surf_coherence(surf)
        >>> print(f"Correlation: {result.correlation:.3f}")
        >>> best_zones = result.get_best_targeting_zones()
    """
    # Convert to landscape
    landscape = surf_to_regulatory_landscape(surf, modality=modality)

    # Compute coherence profile using landscapes module
    from ..landscapes.core import ResponseLandscape

    response_landscape = ResponseLandscape(
        coords=landscape.positions,
        responses=landscape.responses,
        metadata=landscape.metadata,
    )

    profile = compute_spatial_coherence(
        response_landscape,
        window=window,
    )

    # Classify regions
    regions = classify_regulatory_regions(landscape, window=window)

    # Find SURF regions that overlap stable coherence zones
    stable_surf_regions = []
    stable_coherence_regions = [r for r in regions if r.is_safe]

    for surf_region in surf.regions:
        for coherence_region in stable_coherence_regions:
            # Check overlap
            if not (coherence_region.end < surf_region.start or
                    coherence_region.start > surf_region.end):
                stable_surf_regions.append(surf_region)
                break

    # Validation statistics
    validation = {
        'coherence_variance_correlation': profile.correlation,
        'p_value': profile.p_value,
        'is_validated': profile.is_validated,
        'variance_reduction_estimate': profile.variance_reduction_estimate,
        'n_surf_regions_in_stable_zones': len(stable_surf_regions),
        'fraction_surf_regions_stable': (
            len(stable_surf_regions) / len(surf.regions)
            if surf.regions else 0.0
        ),
    }

    return SURFCoherenceResult(
        surf=surf,
        landscape=landscape,
        profile=profile,
        regions=regions,
        stable_surf_regions=stable_surf_regions,
        validation=validation,
    )


def compare_raw_vs_surf(
    raw_positions: np.ndarray,
    raw_responses: np.ndarray,
    surf: SURFOutput,
    window: int = 50,
    gene_symbol: str = "unknown",
) -> Dict[str, Any]:
    """
    Compare coherence analysis on raw screen data vs SURF-deconvolved data.

    This validates that SURF deconvolution improves the coherence signal,
    which is expected because SURF removes technical noise.

    Args:
        raw_positions: Raw screen positions.
        raw_responses: Raw screen responses (log fold change).
        surf: SURF output for same screen.
        window: Coherence window size.
        gene_symbol: Gene name.

    Returns:
        Dictionary comparing raw vs SURF coherence metrics.

    Example:
        >>> # Load raw data and SURF output
        >>> comparison = compare_raw_vs_surf(
        ...     raw_positions=screen_data['position'],
        ...     raw_responses=screen_data['logFC'],
        ...     surf=surf_output,
        ... )
        >>> print(f"SURF improvement: {comparison['correlation_improvement']:.3f}")
    """
    from ..landscapes.core import ResponseLandscape

    # Analyze raw data
    raw_landscape = ResponseLandscape(
        coords=raw_positions,
        responses=raw_responses,
        metadata={'source': 'raw', 'gene': gene_symbol},
    )
    raw_profile = compute_spatial_coherence(raw_landscape, window=window)

    # Analyze SURF data
    surf_landscape = ResponseLandscape(
        coords=surf.positions,
        responses=surf.beta,
        metadata={'source': 'SURF', 'gene': gene_symbol},
    )
    surf_profile = compute_spatial_coherence(surf_landscape, window=window)

    return {
        'gene_symbol': gene_symbol,
        'raw': {
            'correlation': raw_profile.correlation,
            'p_value': raw_profile.p_value,
            'is_validated': raw_profile.is_validated,
            'variance_reduction': raw_profile.variance_reduction_estimate,
            'mean_coherence': float(np.nanmean(raw_profile.coherence)),
        },
        'surf': {
            'correlation': surf_profile.correlation,
            'p_value': surf_profile.p_value,
            'is_validated': surf_profile.is_validated,
            'variance_reduction': surf_profile.variance_reduction_estimate,
            'mean_coherence': float(np.nanmean(surf_profile.coherence)),
        },
        'comparison': {
            'correlation_improvement': (
                abs(surf_profile.correlation) - abs(raw_profile.correlation)
            ),
            'variance_reduction_improvement': (
                surf_profile.variance_reduction_estimate -
                raw_profile.variance_reduction_estimate
            ),
            'surf_better': abs(surf_profile.correlation) > abs(raw_profile.correlation),
        },
    }


def identify_high_confidence_targets(
    result: SURFCoherenceResult,
    min_beta: float = 0.5,
    min_coherence: float = 0.7,
) -> List[Dict[str, Any]]:
    """
    Identify highest-confidence targeting zones.

    High-confidence targets have:
    1. Stable coherence (low variance)
    2. Strong SURF signal (significant beta)
    3. Overlap between coherence and SURF regions

    Args:
        result: SURFCoherenceResult from compute_surf_coherence.
        min_beta: Minimum absolute beta value.
        min_coherence: Minimum coherence score.

    Returns:
        List of high-confidence target regions.
    """
    targets = []

    stable_regions = [r for r in result.regions if r.is_safe]

    for region in stable_regions:
        if region.coherence_score < min_coherence:
            continue

        # Find overlapping SURF regions
        overlapping_surf = []
        for surf_region in result.surf.regions:
            if not (region.end < surf_region.start or region.start > surf_region.end):
                if abs(surf_region.beta_mean) >= min_beta:
                    overlapping_surf.append(surf_region)

        if overlapping_surf:
            # Compute combined confidence
            best_surf = max(overlapping_surf, key=lambda r: abs(r.beta_mean))

            target = {
                'start': region.start,
                'end': region.end,
                'coherence_score': region.coherence_score,
                'variance_reduction': region.variance_reduction,
                'surf_beta': best_surf.beta_mean,
                'surf_significance': best_surf.significance,
                'direction': best_surf.direction,
                'confidence_score': (
                    region.coherence_score *
                    region.variance_reduction *
                    abs(best_surf.beta_mean)
                ),
            }
            targets.append(target)

    # Sort by confidence
    targets.sort(key=lambda x: x['confidence_score'], reverse=True)

    return targets
