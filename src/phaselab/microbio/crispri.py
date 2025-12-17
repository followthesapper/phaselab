"""
PhaseLab Microbio: CRISPRi screen analysis for bacteria.

CRISPRi (CRISPR interference) screens in bacteria measure:
- Gene essentiality under selection
- Guide position effects on knockdown
- Fitness consequences of expression reduction

PhaseLab applies spatial coherence to:
- Identify stable knockdown positions
- Rank guides by reliability
- Compare conditions

This extends the E213 methodology to bacterial gene knockdown.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Union
from pathlib import Path

from ..landscapes.core import (
    ResponseLandscape,
    CoherenceProfile,
    StabilityClass,
)
from ..landscapes.coherence import compute_spatial_coherence
from ..landscapes.classification import classify_regions


@dataclass
class CRISPRiLandscape:
    """
    CRISPRi screen landscape for a gene.

    Attributes:
        positions: Guide positions (TSS-relative, negative = upstream)
        fitness: Fitness scores (log2 fold change)
        gene_id: Gene identifier
        gene_name: Gene name
        guides: Guide sequences (optional)
        strain: Bacterial strain
        condition: Growth condition
        metadata: Additional metadata
    """
    positions: np.ndarray
    fitness: np.ndarray
    gene_id: str
    gene_name: Optional[str] = None
    guides: Optional[List[str]] = None
    strain: str = "unknown"
    condition: str = "default"
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_guides(self) -> int:
        """Number of guides."""
        return len(self.positions)

    @property
    def mean_fitness(self) -> float:
        """Mean fitness effect."""
        return float(np.nanmean(self.fitness))

    @property
    def is_essential(self) -> bool:
        """Whether gene appears essential (mean fitness < -1)."""
        return self.mean_fitness < -1.0

    def to_response_landscape(self) -> ResponseLandscape:
        """Convert to generic ResponseLandscape."""
        return ResponseLandscape(
            coords=self.positions,
            responses=self.fitness,
            metadata={
                'type': 'crispri',
                'gene_id': self.gene_id,
                'gene_name': self.gene_name,
                'strain': self.strain,
                'condition': self.condition,
                **self.metadata,
            },
        )


@dataclass
class CRISPRiResult:
    """
    Result of CRISPRi coherence analysis.

    Attributes:
        landscape: Original CRISPRi landscape
        profile: Coherence profile
        stable_regions: Stable knockdown regions
        ranked_guides: Guides ranked by reliability
        essentiality_score: Overall essentiality assessment
        validation: Validation statistics
    """
    landscape: CRISPRiLandscape
    profile: CoherenceProfile
    stable_regions: List[Dict[str, Any]]
    ranked_guides: List[Dict[str, Any]]
    essentiality_score: float
    validation: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_validated(self) -> bool:
        """Whether coherence model is validated."""
        return self.profile.is_validated

    @property
    def top_guides(self) -> List[Dict[str, Any]]:
        """Top 5 most reliable guides."""
        return self.ranked_guides[:5]

    def summary(self) -> str:
        """Generate summary string."""
        lines = [
            "=" * 60,
            f"CRISPRi ANALYSIS: {self.landscape.gene_name or self.landscape.gene_id}",
            "=" * 60,
            "",
            f"Strain: {self.landscape.strain}",
            f"Condition: {self.landscape.condition}",
            f"Total guides: {self.landscape.n_guides}",
            "",
            "ESSENTIALITY:",
            f"  Mean fitness: {self.landscape.mean_fitness:.3f}",
            f"  Essential: {'YES' if self.landscape.is_essential else 'NO'}",
            f"  Essentiality score: {self.essentiality_score:.3f}",
            "",
            "COHERENCE:",
            f"  Correlation: {self.profile.correlation:.3f}",
            f"  Validated: {'YES' if self.is_validated else 'NO'}",
            "",
            f"Stable regions: {len(self.stable_regions)}",
            f"Reliable guides: {len([g for g in self.ranked_guides if g.get('reliable')])}",
        ]

        if self.top_guides:
            lines.extend(["", "TOP GUIDES:"])
            for guide in self.top_guides:
                lines.append(
                    f"  pos={guide['position']:+d} "
                    f"fitness={guide['fitness']:.3f} "
                    f"reliability={guide.get('reliability', 'N/A')}"
                )

        lines.append("=" * 60)
        return "\n".join(lines)


def load_crispri_screen(
    filepath: Union[str, Path],
    gene_id: str,
    position_col: str = "position",
    fitness_col: str = "fitness",
    guide_col: Optional[str] = "guide",
    delimiter: str = "\t",
    gene_name: Optional[str] = None,
    strain: str = "unknown",
    condition: str = "default",
) -> CRISPRiLandscape:
    """
    Load CRISPRi screen data for a gene.

    Args:
        filepath: Path to screen data file.
        gene_id: Gene identifier.
        position_col: Position column name.
        fitness_col: Fitness column name.
        guide_col: Guide sequence column name.
        delimiter: Column delimiter.
        gene_name: Optional gene name.
        strain: Bacterial strain.
        condition: Growth condition.

    Returns:
        CRISPRiLandscape for analysis.
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"CRISPRi file not found: {filepath}")

    try:
        import pandas as pd
        df = pd.read_csv(filepath, delimiter=delimiter)

        # Filter by gene
        if 'gene' in df.columns:
            df = df[df['gene'] == gene_id]
        elif 'gene_id' in df.columns:
            df = df[df['gene_id'] == gene_id]
        elif 'target' in df.columns:
            df = df[df['target'] == gene_id]

        if len(df) == 0:
            raise ValueError(f"No data found for gene {gene_id}")

        positions = df[position_col].values
        fitness = df[fitness_col].values

        guides = None
        if guide_col and guide_col in df.columns:
            guides = df[guide_col].tolist()

        return CRISPRiLandscape(
            positions=positions.astype(float),
            fitness=fitness.astype(float),
            gene_id=gene_id,
            gene_name=gene_name,
            guides=guides,
            strain=strain,
            condition=condition,
            metadata={'source_file': str(filepath)},
        )

    except ImportError:
        raise ImportError("pandas required for CRISPRi loading")


def analyze_crispri_coherence(
    landscape: CRISPRiLandscape,
    window: int = 10,
    stable_threshold: float = 0.7,
) -> CRISPRiResult:
    """
    Analyze spatial coherence in CRISPRi screen.

    Identifies:
    - Stable knockdown regions
    - Reliable guide positions
    - Overall gene essentiality with confidence

    Args:
        landscape: CRISPRi landscape to analyze.
        window: Window size for coherence.
        stable_threshold: Threshold for stable regions.

    Returns:
        CRISPRiResult with complete analysis.

    Example:
        >>> landscape = load_crispri_screen('screen.tsv', 'dnaA')
        >>> result = analyze_crispri_coherence(landscape)
        >>> reliable_guides = result.top_guides
    """
    # Convert to response landscape
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
                'mean_fitness': float(np.nanmean(landscape.fitness[mask])),
                'n_guides': int(np.sum(mask)),
            })

    # Rank guides
    ranked_guides = rank_crispri_guides(landscape, profile, stable_regions)

    # Compute essentiality score
    # Combines fitness magnitude with reliability (coherence)
    if stable_regions:
        stable_fitness = np.mean([r['mean_fitness'] for r in stable_regions])
        stable_coherence = np.mean([r['coherence'] for r in stable_regions])
        essentiality_score = -stable_fitness * stable_coherence
    else:
        essentiality_score = -landscape.mean_fitness * 0.5  # Discount without stability

    validation = {
        'coherence_variance_correlation': profile.correlation,
        'p_value': profile.p_value,
        'is_validated': profile.is_validated,
    }

    return CRISPRiResult(
        landscape=landscape,
        profile=profile,
        stable_regions=stable_regions,
        ranked_guides=ranked_guides,
        essentiality_score=essentiality_score,
        validation=validation,
    )


def rank_crispri_guides(
    landscape: CRISPRiLandscape,
    profile: CoherenceProfile,
    stable_regions: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """
    Rank CRISPRi guides by reliability.

    Guides in stable regions are prioritized.
    Score combines fitness effect and local coherence.

    Args:
        landscape: CRISPRi landscape.
        profile: Coherence profile.
        stable_regions: List of stable region dictionaries.

    Returns:
        List of ranked guide dictionaries.
    """
    ranked = []

    for i, (pos, fitness) in enumerate(zip(landscape.positions, landscape.fitness)):
        # Find local coherence (handle empty profile edge case)
        if len(profile.coords) == 0:
            local_coherence = 0.0
            local_variance = 1.0
        else:
            coh_idx = np.argmin(np.abs(profile.coords - pos))
            local_coherence = float(profile.coherence[coh_idx])
            local_variance = float(profile.local_variance[coh_idx])

        # Check if in stable region
        in_stable = False
        for region in stable_regions:
            if region['start'] <= pos <= region['end']:
                in_stable = True
                break

        guide_info = {
            'position': int(pos),
            'fitness': float(fitness),
            'coherence': local_coherence,
            'local_variance': local_variance,
            'in_stable_region': in_stable,
            'reliable': in_stable and local_coherence > 0.7,
        }

        if landscape.guides and i < len(landscape.guides):
            guide_info['sequence'] = landscape.guides[i]

        # Compute reliability score
        # Strong fitness + high coherence + in stable region
        reliability = abs(fitness) * local_coherence
        if in_stable:
            reliability *= 1.5
        guide_info['reliability_score'] = reliability

        ranked.append(guide_info)

    # Sort by reliability
    ranked.sort(key=lambda x: x['reliability_score'], reverse=True)

    return ranked
