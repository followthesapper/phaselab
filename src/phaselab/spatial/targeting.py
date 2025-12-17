"""
PhaseLab Spatial: Guide targeting within classified regions.

This module connects region classification to guide selection.
Guides are only ranked WITHIN stable regions - this is the key insight
from E212-E216.

The workflow:
1. Classify regions (spatial module)
2. Get guide candidates (CRISPOR or enumeration)
3. Filter guides to stable regions only
4. Rank guides by CRISPOR metrics within each region
5. Output prioritized guide list with region context
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Tuple
from enum import Enum

from .regulatory import (
    RegulatoryLandscape,
    RegulatoryRegion,
    classify_regulatory_regions,
)
from ..landscapes.core import StabilityClass


class TargetPriority(Enum):
    """Priority level for targeting recommendations."""
    PRIMARY = "primary"      # Stable region, high confidence
    SECONDARY = "secondary"  # Mixed region or lower confidence
    EXPLORATORY = "exploratory"  # Not validated but potentially useful
    AVOID = "avoid"          # Amplifying region


@dataclass
class TargetRecommendation:
    """
    A specific targeting recommendation.

    Combines region stability with guide-level information.
    """
    position: int  # TSS-relative position
    region_start: int
    region_end: int
    region_stability: StabilityClass
    coherence_score: float
    variance_reduction: float
    priority: TargetPriority
    guide_sequence: Optional[str] = None
    guide_score: Optional[float] = None  # CRISPOR or other score
    off_targets: Optional[Dict[str, int]] = None
    recommendation: str = ""

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'position': self.position,
            'region_start': self.region_start,
            'region_end': self.region_end,
            'region_stability': self.region_stability.value,
            'coherence_score': self.coherence_score,
            'variance_reduction': self.variance_reduction,
            'priority': self.priority.value,
            'guide_sequence': self.guide_sequence,
            'guide_score': self.guide_score,
            'off_targets': self.off_targets,
            'recommendation': self.recommendation,
        }


class CRISPRaTargetFinder:
    """
    Find optimal CRISPRa targets within classified regulatory regions.

    This class implements the complete E215 methodology:
    1. Classify regions by spatial coherence
    2. Filter to stable/mixed regions
    3. Rank positions within regions
    4. Generate prioritized recommendations

    Example:
        >>> finder = CRISPRaTargetFinder(landscape)
        >>> targets = finder.recommend_targets(n=10)
        >>> for t in targets:
        ...     print(f"{t.position}: {t.priority.value} - {t.recommendation}")
    """

    def __init__(
        self,
        landscape: RegulatoryLandscape,
        regions: Optional[List[RegulatoryRegion]] = None,
        window: int = 50,
    ):
        """
        Initialize the target finder.

        Args:
            landscape: RegulatoryLandscape with perturbation response data.
            regions: Pre-computed regions (computed if not provided).
            window: Window size for coherence computation.
        """
        self.landscape = landscape
        self.window = window

        if regions is None:
            self.regions = classify_regulatory_regions(landscape, window=window)
        else:
            self.regions = regions

        # Index regions for fast lookup
        self._build_region_index()

    def _build_region_index(self):
        """Build index for fast position → region lookup."""
        self._region_index = {}
        for region in self.regions:
            for pos in range(region.start, region.end + 1):
                self._region_index[pos] = region

    def get_region_at(self, position: int) -> Optional[RegulatoryRegion]:
        """Get the region containing a specific position."""
        return self._region_index.get(position)

    @property
    def stable_regions(self) -> List[RegulatoryRegion]:
        """Get all stable regions."""
        return [r for r in self.regions if r.is_safe]

    @property
    def mixed_regions(self) -> List[RegulatoryRegion]:
        """Get all mixed regions."""
        return [r for r in self.regions if r.stability == StabilityClass.MIXED]

    @property
    def amplifying_regions(self) -> List[RegulatoryRegion]:
        """Get all amplifying regions (to avoid)."""
        return [r for r in self.regions if r.is_dangerous]

    def recommend_targets(
        self,
        n: int = 10,
        include_mixed: bool = False,
        avoid_boundaries: int = 20,
        min_coherence: float = 0.5,
    ) -> List[TargetRecommendation]:
        """
        Generate prioritized target recommendations.

        Args:
            n: Maximum number of recommendations.
            include_mixed: Whether to include mixed regions as SECONDARY.
            avoid_boundaries: Distance from region boundaries to avoid.
            min_coherence: Minimum coherence score for recommendation.

        Returns:
            List of TargetRecommendation objects, sorted by priority.
        """
        recommendations = []

        # Get candidate regions
        candidate_regions = self.stable_regions.copy()
        if include_mixed:
            candidate_regions.extend(self.mixed_regions)

        for region in candidate_regions:
            if region.coherence_score < min_coherence:
                continue

            # Calculate safe targeting range within region
            safe_start = region.start + avoid_boundaries
            safe_end = region.end - avoid_boundaries

            if safe_end <= safe_start:
                continue

            # Best position is center of region (most stable)
            center = (safe_start + safe_end) // 2

            # Determine priority
            if region.stability == StabilityClass.STABLE:
                priority = TargetPriority.PRIMARY
            else:
                priority = TargetPriority.SECONDARY

            # Generate recommendation text
            rec_text = self._generate_recommendation(region, priority)

            recommendations.append(TargetRecommendation(
                position=center,
                region_start=region.start,
                region_end=region.end,
                region_stability=region.stability,
                coherence_score=region.coherence_score,
                variance_reduction=region.variance_reduction,
                priority=priority,
                recommendation=rec_text,
            ))

        # Sort by priority and coherence score
        recommendations.sort(
            key=lambda x: (
                0 if x.priority == TargetPriority.PRIMARY else 1,
                -x.coherence_score,
            )
        )

        return recommendations[:n]

    def _generate_recommendation(
        self,
        region: RegulatoryRegion,
        priority: TargetPriority,
    ) -> str:
        """Generate recommendation text for a region."""
        if priority == TargetPriority.PRIMARY:
            return (
                f"Target region [{region.start}, {region.end}] relative to TSS. "
                f"Spatial coherence {region.coherence_score:.2f} predicts "
                f"{region.variance_reduction:.0%} variance reduction. "
                f"RECOMMENDED for {self.landscape.modality}."
            )
        elif priority == TargetPriority.SECONDARY:
            return (
                f"Region [{region.start}, {region.end}] shows moderate stability. "
                f"Consider as backup if primary targets fail. "
                f"Validate experimentally before large-scale use."
            )
        else:
            return f"Region [{region.start}, {region.end}] - exercise caution."

    def filter_guides_to_safe_regions(
        self,
        guides: List[Dict[str, Any]],
        position_key: str = 'position',
    ) -> List[Dict[str, Any]]:
        """
        Filter a list of guides to only those in stable regions.

        Args:
            guides: List of guide dictionaries with position information.
            position_key: Key for position in guide dictionaries.

        Returns:
            Filtered list of guides with region information added.
        """
        filtered = []

        for guide in guides:
            pos = guide.get(position_key)
            if pos is None:
                continue

            region = self.get_region_at(pos)
            if region is None:
                continue

            if region.is_safe or region.stability == StabilityClass.MIXED:
                guide_with_region = guide.copy()
                guide_with_region['region_start'] = region.start
                guide_with_region['region_end'] = region.end
                guide_with_region['region_stability'] = region.stability.value
                guide_with_region['coherence_score'] = region.coherence_score
                guide_with_region['variance_reduction'] = region.variance_reduction
                guide_with_region['is_in_stable_region'] = region.is_safe
                filtered.append(guide_with_region)

        return filtered

    def annotate_guides_with_regions(
        self,
        guides: List[Dict[str, Any]],
        position_key: str = 'position',
    ) -> List[Dict[str, Any]]:
        """
        Annotate all guides with region information (without filtering).

        Args:
            guides: List of guide dictionaries.
            position_key: Key for position in guide dictionaries.

        Returns:
            Annotated guides with region information.
        """
        annotated = []

        for guide in guides:
            pos = guide.get(position_key)
            guide_annotated = guide.copy()

            if pos is not None:
                region = self.get_region_at(pos)
                if region is not None:
                    guide_annotated['region_start'] = region.start
                    guide_annotated['region_end'] = region.end
                    guide_annotated['region_stability'] = region.stability.value
                    guide_annotated['coherence_score'] = region.coherence_score
                    guide_annotated['is_in_stable_region'] = region.is_safe
                    guide_annotated['is_in_amplifying_region'] = region.is_dangerous
                else:
                    guide_annotated['region_stability'] = 'unknown'
                    guide_annotated['is_in_stable_region'] = False
                    guide_annotated['is_in_amplifying_region'] = False

            annotated.append(guide_annotated)

        return annotated

    def summary(self) -> str:
        """Generate a summary of the targeting analysis."""
        lines = [
            "=" * 60,
            f"CRISPRa TARGET FINDER: {self.landscape.gene_symbol}",
            "=" * 60,
            f"Total regions: {len(self.regions)}",
            f"  STABLE (primary): {len(self.stable_regions)}",
            f"  MIXED (secondary): {len(self.mixed_regions)}",
            f"  AMPLIFYING (avoid): {len(self.amplifying_regions)}",
            "",
        ]

        if self.stable_regions:
            lines.append("PRIMARY TARGET REGIONS:")
            for r in self.stable_regions[:3]:
                lines.append(
                    f"  [{r.start}, {r.end}] - coherence={r.coherence_score:.2f}, "
                    f"var_reduction={r.variance_reduction:.0%}"
                )

        if self.amplifying_regions:
            lines.append("")
            lines.append("REGIONS TO AVOID (amplifying):")
            for r in self.amplifying_regions[:3]:
                lines.append(f"  [{r.start}, {r.end}]")

        lines.append("=" * 60)
        return "\n".join(lines)


class CRISPRiTargetFinder(CRISPRaTargetFinder):
    """
    Find optimal CRISPRi targets within classified regulatory regions.

    CRISPRi has similar region-based targeting logic to CRISPRa,
    but with different optimal positioning (closer to TSS typically).
    """

    def __init__(
        self,
        landscape: RegulatoryLandscape,
        regions: Optional[List[RegulatoryRegion]] = None,
        window: int = 50,
    ):
        """Initialize CRISPRi target finder."""
        # Ensure modality is CRISPRi
        if landscape.modality != "CRISPRi":
            landscape = RegulatoryLandscape(
                positions=landscape.positions,
                responses=landscape.responses,
                gene_symbol=landscape.gene_symbol,
                tss_position=landscape.tss_position,
                chromosome=landscape.chromosome,
                genome_build=landscape.genome_build,
                modality="CRISPRi",
                replicates=landscape.replicates,
                metadata=landscape.metadata,
            )

        super().__init__(landscape, regions, window)

    def recommend_targets(
        self,
        n: int = 10,
        include_mixed: bool = False,
        avoid_boundaries: int = 15,  # Smaller buffer for CRISPRi
        min_coherence: float = 0.5,
        prefer_proximal: bool = True,  # CRISPRi often works better near TSS
    ) -> List[TargetRecommendation]:
        """
        Generate CRISPRi target recommendations.

        CRISPRi typically works best in the region -50 to +300 from TSS.
        """
        recommendations = super().recommend_targets(
            n=n * 2,  # Get more, then filter
            include_mixed=include_mixed,
            avoid_boundaries=avoid_boundaries,
            min_coherence=min_coherence,
        )

        if prefer_proximal:
            # Prefer positions closer to TSS
            recommendations.sort(
                key=lambda x: (
                    0 if x.priority == TargetPriority.PRIMARY else 1,
                    abs(x.position),  # Prefer closer to TSS
                    -x.coherence_score,
                )
            )

        return recommendations[:n]


def rank_guides_in_regions(
    guides: List[Dict[str, Any]],
    regions: List[RegulatoryRegion],
    position_key: str = 'position',
    score_key: str = 'mit_specificity',
    higher_is_better: bool = True,
) -> List[Dict[str, Any]]:
    """
    Rank guides within their respective regions.

    Guides in STABLE regions are ranked first, then MIXED, then others.
    Within each region, guides are ranked by the specified score.

    Args:
        guides: List of guide dictionaries.
        regions: List of classified RegulatoryRegion objects.
        position_key: Key for position in guide dictionaries.
        score_key: Key for score to rank by (e.g., 'mit_specificity').
        higher_is_better: Whether higher scores are better.

    Returns:
        Ranked list of guides with region annotations.
    """
    # Build region index
    region_index = {}
    for region in regions:
        for pos in range(region.start, region.end + 1):
            region_index[pos] = region

    # Annotate guides with regions
    annotated = []
    for guide in guides:
        pos = guide.get(position_key)
        guide_copy = guide.copy()

        region = region_index.get(pos) if pos is not None else None
        if region is not None:
            guide_copy['_region'] = region
            guide_copy['_stability_rank'] = {
                StabilityClass.STABLE: 0,
                StabilityClass.MIXED: 1,
                StabilityClass.IRRELEVANT: 2,
                StabilityClass.AMPLIFYING: 3,
            }.get(region.stability, 4)
        else:
            guide_copy['_region'] = None
            guide_copy['_stability_rank'] = 5

        annotated.append(guide_copy)

    # Sort by stability rank, then by score
    def sort_key(g):
        stability_rank = g.get('_stability_rank', 5)
        score = g.get(score_key, 0)
        if not higher_is_better:
            score = -score
        return (stability_rank, -score)

    annotated.sort(key=sort_key)

    # Clean up and add final annotations
    result = []
    for i, guide in enumerate(annotated):
        region = guide.pop('_region', None)
        guide.pop('_stability_rank', None)

        guide['rank'] = i + 1
        if region is not None:
            guide['region_stability'] = region.stability.value
            guide['in_stable_region'] = region.is_safe
            guide['coherence_score'] = region.coherence_score
        else:
            guide['region_stability'] = 'unknown'
            guide['in_stable_region'] = False

        result.append(guide)

    return result
