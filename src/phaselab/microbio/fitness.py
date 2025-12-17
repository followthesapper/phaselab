"""
PhaseLab Microbio: General fitness landscape analysis.

This module provides generic fitness landscape tools that work across:
- TnSeq data
- CRISPRi screens
- Growth rate measurements
- Competitive fitness assays

These are the building blocks used by more specific modules.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Union
from pathlib import Path

from ..landscapes.core import (
    ResponseLandscape,
    CoherenceProfile,
    StabilityClass,
    RegionClassification,
)
from ..landscapes.coherence import compute_spatial_coherence
from ..landscapes.classification import classify_regions


@dataclass
class FitnessLandscape:
    """
    Generic fitness landscape.

    Can represent any perturbation-fitness relationship:
    - Position → fitness (TnSeq, CRISPRi)
    - Mutation → fitness (deep mutational scanning)
    - Condition → fitness (environment screens)

    Attributes:
        coordinates: Perturbation coordinates
        fitness: Fitness values
        labels: Optional labels for coordinates
        experiment_type: Type of experiment
        organism: Organism identifier
        condition: Experimental condition
        metadata: Additional metadata
    """
    coordinates: np.ndarray
    fitness: np.ndarray
    labels: Optional[List[str]] = None
    experiment_type: str = "generic"
    organism: str = "unknown"
    condition: str = "default"
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_points(self) -> int:
        """Number of data points."""
        return len(self.coordinates)

    @property
    def mean_fitness(self) -> float:
        """Mean fitness."""
        return float(np.nanmean(self.fitness))

    @property
    def fitness_variance(self) -> float:
        """Fitness variance."""
        return float(np.nanvar(self.fitness))

    def to_response_landscape(self) -> ResponseLandscape:
        """Convert to generic ResponseLandscape."""
        return ResponseLandscape(
            coords=self.coordinates,
            responses=self.fitness,
            coord_labels=self.labels,
            metadata={
                'type': self.experiment_type,
                'organism': self.organism,
                'condition': self.condition,
                **self.metadata,
            },
        )


@dataclass
class FunctionalDomain:
    """
    A functional domain identified by coherence analysis.

    Attributes:
        start: Start coordinate
        end: End coordinate
        stability: Stability classification
        coherence: Coherence score
        mean_fitness: Mean fitness in domain
        fitness_effect: Classification of fitness effect
        n_points: Number of data points
    """
    start: float
    end: float
    stability: StabilityClass
    coherence: float
    mean_fitness: float
    fitness_effect: str  # "essential", "neutral", "beneficial"
    n_points: int

    @property
    def is_essential(self) -> bool:
        """Whether domain is essential (strong fitness defect)."""
        return self.fitness_effect == "essential"

    @property
    def is_stable(self) -> bool:
        """Whether domain is stable."""
        return self.stability == StabilityClass.STABLE

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'start': self.start,
            'end': self.end,
            'stability': self.stability.value,
            'coherence': self.coherence,
            'mean_fitness': self.mean_fitness,
            'fitness_effect': self.fitness_effect,
            'n_points': self.n_points,
            'is_essential': self.is_essential,
            'is_stable': self.is_stable,
        }


def compute_fitness_coherence(
    landscape: FitnessLandscape,
    window: int = 20,
) -> CoherenceProfile:
    """
    Compute spatial coherence for a fitness landscape.

    Args:
        landscape: Fitness landscape to analyze.
        window: Window size for coherence computation.

    Returns:
        CoherenceProfile with spatial coherence metrics.

    Example:
        >>> profile = compute_fitness_coherence(landscape)
        >>> print(f"Correlation: {profile.correlation:.3f}")
        >>> print(f"Variance reduction: {profile.variance_reduction_estimate:.1%}")
    """
    response_landscape = landscape.to_response_landscape()
    return compute_spatial_coherence(response_landscape, window=window)


def identify_functional_domains(
    landscape: FitnessLandscape,
    window: int = 20,
    stable_threshold: float = 0.7,
    essential_threshold: float = -1.0,
    beneficial_threshold: float = 0.5,
) -> List[FunctionalDomain]:
    """
    Identify functional domains in a fitness landscape.

    Domains are classified by:
    1. Stability (spatial coherence)
    2. Fitness effect (essential, neutral, beneficial)

    Args:
        landscape: Fitness landscape to analyze.
        window: Window size for coherence.
        stable_threshold: Coherence threshold for stable regions.
        essential_threshold: Fitness threshold for essential (below = essential).
        beneficial_threshold: Fitness threshold for beneficial (above = beneficial).

    Returns:
        List of FunctionalDomain objects.

    Example:
        >>> domains = identify_functional_domains(landscape)
        >>> essential = [d for d in domains if d.is_essential and d.is_stable]
        >>> print(f"High-confidence essential domains: {len(essential)}")
    """
    response_landscape = landscape.to_response_landscape()

    # Compute coherence and classify regions
    profile = compute_spatial_coherence(response_landscape, window=window)
    classification = classify_regions(
        response_landscape,
        profile=profile,
        stable_threshold=stable_threshold,
    )

    domains = []

    for region in classification.regions:
        start, end, stability, score = region

        # Get fitness values in this region
        mask = (landscape.coordinates >= start) & (landscape.coordinates <= end)
        region_fitness = landscape.fitness[mask]

        if len(region_fitness) == 0:
            continue

        mean_fitness = float(np.nanmean(region_fitness))

        # Classify fitness effect
        if mean_fitness < essential_threshold:
            fitness_effect = "essential"
        elif mean_fitness > beneficial_threshold:
            fitness_effect = "beneficial"
        else:
            fitness_effect = "neutral"

        domain = FunctionalDomain(
            start=float(start),
            end=float(end),
            stability=stability,
            coherence=float(score),
            mean_fitness=mean_fitness,
            fitness_effect=fitness_effect,
            n_points=int(np.sum(mask)),
        )
        domains.append(domain)

    return domains


def essential_gene_analysis(
    landscape: FitnessLandscape,
    window: int = 20,
    stable_threshold: float = 0.7,
    essential_threshold: float = -1.0,
) -> Dict[str, Any]:
    """
    Comprehensive essential gene analysis.

    Determines:
    1. Is the gene essential (overall)?
    2. Which domains are essential?
    3. How confident is the essentiality call?

    Args:
        landscape: Fitness landscape for a gene.
        window: Coherence window size.
        stable_threshold: Threshold for stable regions.
        essential_threshold: Fitness threshold for essential.

    Returns:
        Dictionary with essentiality assessment.

    Example:
        >>> result = essential_gene_analysis(landscape)
        >>> if result['is_essential']:
        ...     print(f"Gene is essential with confidence: {result['confidence']}")
        ...     print(f"Essential domains: {result['n_essential_domains']}")
    """
    domains = identify_functional_domains(
        landscape,
        window=window,
        stable_threshold=stable_threshold,
        essential_threshold=essential_threshold,
    )

    # Find essential domains
    essential_domains = [d for d in domains if d.is_essential]
    stable_essential = [d for d in essential_domains if d.is_stable]

    # Overall essentiality
    total_essential_points = sum(d.n_points for d in essential_domains)
    total_points = landscape.n_points

    # Confidence based on:
    # - Fraction of gene that's essential
    # - Whether essential regions are stable (reproducible)
    essential_fraction = total_essential_points / total_points if total_points > 0 else 0

    if stable_essential:
        stable_essential_fraction = sum(d.n_points for d in stable_essential) / total_points
        confidence = "high" if stable_essential_fraction > 0.3 else "medium"
    elif essential_domains:
        confidence = "low"  # Essential but not stable
    else:
        confidence = "none"  # Not essential

    is_essential = essential_fraction > 0.2 or (stable_essential and essential_fraction > 0.1)

    return {
        'is_essential': is_essential,
        'confidence': confidence,
        'essential_fraction': essential_fraction,
        'n_essential_domains': len(essential_domains),
        'n_stable_essential_domains': len(stable_essential),
        'mean_fitness': landscape.mean_fitness,
        'essential_domains': [d.to_dict() for d in stable_essential[:5]],
        'recommendation': (
            f"ESSENTIAL (confidence: {confidence})"
            if is_essential
            else "NOT ESSENTIAL or inconclusive"
        ),
    }


def compare_fitness_landscapes(
    landscape_a: FitnessLandscape,
    landscape_b: FitnessLandscape,
    window: int = 20,
) -> Dict[str, Any]:
    """
    Compare two fitness landscapes.

    Useful for:
    - Comparing conditions (e.g., +/- drug)
    - Comparing strains
    - Comparing replicates

    Args:
        landscape_a: First fitness landscape.
        landscape_b: Second fitness landscape.
        window: Coherence window size.

    Returns:
        Comparison results.
    """
    profile_a = compute_fitness_coherence(landscape_a, window=window)
    profile_b = compute_fitness_coherence(landscape_b, window=window)

    domains_a = identify_functional_domains(landscape_a, window=window)
    domains_b = identify_functional_domains(landscape_b, window=window)

    # Find common coordinates
    common = np.intersect1d(landscape_a.coordinates, landscape_b.coordinates)

    mask_a = np.isin(landscape_a.coordinates, common)
    mask_b = np.isin(landscape_b.coordinates, common)

    if len(common) > 2:
        fitness_correlation = np.corrcoef(
            landscape_a.fitness[mask_a],
            landscape_b.fitness[mask_b],
        )[0, 1]
    else:
        fitness_correlation = np.nan

    return {
        'landscape_a': {
            'condition': landscape_a.condition,
            'coherence': profile_a.correlation,
            'n_essential_domains': sum(1 for d in domains_a if d.is_essential),
        },
        'landscape_b': {
            'condition': landscape_b.condition,
            'coherence': profile_b.correlation,
            'n_essential_domains': sum(1 for d in domains_b if d.is_essential),
        },
        'comparison': {
            'n_common_points': len(common),
            'fitness_correlation': float(fitness_correlation),
            'coherence_difference': abs(profile_a.correlation - profile_b.correlation),
        },
    }
