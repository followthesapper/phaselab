"""
PhaseLab Microbio: TnSeq (transposon sequencing) analysis.

TnSeq experiments measure gene fitness by:
1. Creating a library of transposon insertions
2. Growing under selection
3. Counting insertion abundance changes

PhaseLab analyzes TnSeq data to identify:
- Which gene regions show stable fitness effects
- Which insertions are reliable reporters
- Domain-level functional organization

This is the microbial analog of CRISPRa tiling screens.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Union, Tuple
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
class TnSeqLandscape:
    """
    TnSeq fitness landscape for a gene or genomic region.

    Attributes:
        positions: Insertion positions (within gene or genomic)
        fitness: Fitness scores (log2 fold change or similar)
        gene_id: Gene identifier
        gene_name: Gene name (if known)
        strand: Gene strand (+/-)
        gene_start: Gene start position
        gene_end: Gene end position
        condition: Experimental condition
        replicates: Optional replicate fitness values
        metadata: Additional metadata
    """
    positions: np.ndarray
    fitness: np.ndarray
    gene_id: str
    gene_name: Optional[str] = None
    strand: str = "+"
    gene_start: Optional[int] = None
    gene_end: Optional[int] = None
    condition: str = "default"
    replicates: Optional[np.ndarray] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_insertions(self) -> int:
        """Number of insertion sites."""
        return len(self.positions)

    @property
    def gene_length(self) -> Optional[int]:
        """Gene length in bp."""
        if self.gene_start is not None and self.gene_end is not None:
            return abs(self.gene_end - self.gene_start)
        return None

    @property
    def insertion_density(self) -> Optional[float]:
        """Insertions per bp."""
        if self.gene_length:
            return self.n_insertions / self.gene_length
        return None

    @property
    def mean_fitness(self) -> float:
        """Mean fitness score."""
        return float(np.nanmean(self.fitness))

    @property
    def fitness_variance(self) -> float:
        """Fitness variance."""
        return float(np.nanvar(self.fitness))

    def to_response_landscape(self) -> ResponseLandscape:
        """Convert to generic ResponseLandscape."""
        return ResponseLandscape(
            coords=self.positions,
            responses=self.fitness if self.replicates is None else self.replicates,
            metadata={
                'type': 'tnseq',
                'gene_id': self.gene_id,
                'gene_name': self.gene_name,
                'condition': self.condition,
                **self.metadata,
            },
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'gene_id': self.gene_id,
            'gene_name': self.gene_name,
            'n_insertions': self.n_insertions,
            'gene_length': self.gene_length,
            'mean_fitness': self.mean_fitness,
            'fitness_variance': self.fitness_variance,
            'condition': self.condition,
            'positions': self.positions.tolist(),
            'fitness': self.fitness.tolist(),
        }


@dataclass
class TnSeqResult:
    """
    Result of TnSeq coherence analysis.

    Attributes:
        landscape: Original TnSeq landscape
        profile: Coherence profile
        classification: Region classification
        essential_domains: Regions with stable fitness effects
        dispensable_domains: Regions with neutral fitness
        validation: Validation statistics
    """
    landscape: TnSeqLandscape
    profile: CoherenceProfile
    classification: RegionClassification
    essential_domains: List[Dict[str, Any]] = field(default_factory=list)
    dispensable_domains: List[Dict[str, Any]] = field(default_factory=list)
    validation: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_validated(self) -> bool:
        """Whether coherence model is validated."""
        return self.profile.is_validated

    @property
    def correlation(self) -> float:
        """Coherence-variance correlation."""
        return self.profile.correlation

    def summary(self) -> str:
        """Generate summary string."""
        lines = [
            "=" * 60,
            f"TnSeq ANALYSIS: {self.landscape.gene_name or self.landscape.gene_id}",
            "=" * 60,
            "",
            "INSERTIONS:",
            f"  Total insertions: {self.landscape.n_insertions}",
            f"  Gene length: {self.landscape.gene_length or 'unknown'} bp",
            f"  Mean fitness: {self.landscape.mean_fitness:.3f}",
            "",
            "COHERENCE:",
            f"  Correlation: {self.correlation:.3f}",
            f"  Validated: {'YES' if self.is_validated else 'NO'}",
            f"  Variance reduction: {self.profile.variance_reduction_estimate:.1%}",
            "",
            "DOMAINS:",
            f"  Essential domains: {len(self.essential_domains)}",
            f"  Dispensable domains: {len(self.dispensable_domains)}",
        ]

        if self.essential_domains:
            lines.append("")
            lines.append("TOP ESSENTIAL DOMAINS:")
            for domain in self.essential_domains[:3]:
                lines.append(
                    f"  [{domain['start']}-{domain['end']}] "
                    f"fitness={domain['mean_fitness']:.3f}"
                )

        lines.append("=" * 60)
        return "\n".join(lines)


def load_tnseq_data(
    filepath: Union[str, Path],
    gene_id: str,
    position_col: str = "position",
    fitness_col: str = "fitness",
    replicate_cols: Optional[List[str]] = None,
    delimiter: str = "\t",
    gene_name: Optional[str] = None,
    condition: str = "default",
) -> TnSeqLandscape:
    """
    Load TnSeq data from a file.

    Supports common TnSeq output formats:
    - TRANSIT output
    - Tn-seq analysis pipelines
    - Generic TSV with position/fitness columns

    Args:
        filepath: Path to TnSeq data file.
        gene_id: Gene identifier to filter for.
        position_col: Column name for insertion position.
        fitness_col: Column name for fitness score.
        replicate_cols: Column names for replicate fitness values.
        delimiter: Column delimiter.
        gene_name: Optional gene name.
        condition: Experimental condition label.

    Returns:
        TnSeqLandscape for the specified gene.

    Example:
        >>> landscape = load_tnseq_data(
        ...     'tnseq_results.tsv',
        ...     gene_id='b0001',
        ...     gene_name='thrL',
        ... )
        >>> print(f"Insertions: {landscape.n_insertions}")
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"TnSeq file not found: {filepath}")

    try:
        import pandas as pd
        df = pd.read_csv(filepath, delimiter=delimiter)

        # Filter by gene if gene column exists
        if 'gene' in df.columns:
            df = df[df['gene'] == gene_id]
        elif 'gene_id' in df.columns:
            df = df[df['gene_id'] == gene_id]
        elif 'locus_tag' in df.columns:
            df = df[df['locus_tag'] == gene_id]

        if len(df) == 0:
            raise ValueError(f"No data found for gene {gene_id}")

        positions = df[position_col].values
        fitness = df[fitness_col].values

        replicates = None
        if replicate_cols:
            replicates = df[replicate_cols].values

        return TnSeqLandscape(
            positions=positions.astype(float),
            fitness=fitness.astype(float),
            gene_id=gene_id,
            gene_name=gene_name,
            condition=condition,
            replicates=replicates,
            metadata={
                'source_file': str(filepath),
            },
        )

    except ImportError:
        raise ImportError("pandas required for TnSeq loading")


def analyze_tnseq_coherence(
    landscape: TnSeqLandscape,
    window: int = 20,
    stable_threshold: float = 0.7,
    fitness_threshold: float = -1.0,
) -> TnSeqResult:
    """
    Analyze spatial coherence in TnSeq fitness landscape.

    Identifies:
    - Regions with stable (reproducible) fitness effects
    - Essential domains (stable + strong fitness defect)
    - Dispensable domains (stable + neutral fitness)

    Args:
        landscape: TnSeq landscape to analyze.
        window: Window size for coherence (in insertion sites).
        stable_threshold: Coherence threshold for stable regions.
        fitness_threshold: Fitness threshold for "essential" (negative = defect).

    Returns:
        TnSeqResult with complete analysis.

    Example:
        >>> result = analyze_tnseq_coherence(landscape)
        >>> print(f"Essential domains: {len(result.essential_domains)}")
        >>> print(result.summary())
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
        min_region_size=max(5, window // 4),
    )

    # Identify essential and dispensable domains
    essential_domains = []
    dispensable_domains = []

    for region in classification.regions:
        start, end, stability, score = region

        if stability != StabilityClass.STABLE:
            continue

        # Get fitness values in this region
        mask = (landscape.positions >= start) & (landscape.positions <= end)
        region_fitness = landscape.fitness[mask]

        if len(region_fitness) == 0:
            continue

        mean_fitness = float(np.nanmean(region_fitness))

        domain_info = {
            'start': int(start),
            'end': int(end),
            'coherence': float(score),
            'mean_fitness': mean_fitness,
            'n_insertions': int(np.sum(mask)),
        }

        if mean_fitness < fitness_threshold:
            domain_info['classification'] = 'essential'
            essential_domains.append(domain_info)
        else:
            domain_info['classification'] = 'dispensable'
            dispensable_domains.append(domain_info)

    # Sort by fitness effect
    essential_domains.sort(key=lambda x: x['mean_fitness'])
    dispensable_domains.sort(key=lambda x: -x['coherence'])

    # Validation
    validation = {
        'coherence_variance_correlation': profile.correlation,
        'p_value': profile.p_value,
        'is_validated': profile.is_validated,
        'variance_reduction': profile.variance_reduction_estimate,
    }

    return TnSeqResult(
        landscape=landscape,
        profile=profile,
        classification=classification,
        essential_domains=essential_domains,
        dispensable_domains=dispensable_domains,
        validation=validation,
    )


def identify_essential_domains(
    result: TnSeqResult,
    min_fitness_defect: float = -1.0,
    min_coherence: float = 0.7,
    min_size: int = 5,
) -> List[Dict[str, Any]]:
    """
    Extract high-confidence essential domains.

    Essential domains have:
    1. High spatial coherence (reproducible effect)
    2. Strong fitness defect (insertions are deleterious)
    3. Sufficient size (not single-insertion noise)

    Args:
        result: TnSeqResult from analyze_tnseq_coherence.
        min_fitness_defect: Maximum fitness for "essential" (more negative = more essential).
        min_coherence: Minimum coherence score.
        min_size: Minimum number of insertions in domain.

    Returns:
        List of essential domain dictionaries.
    """
    domains = []

    for domain in result.essential_domains:
        if domain['mean_fitness'] > min_fitness_defect:
            continue
        if domain['coherence'] < min_coherence:
            continue
        if domain['n_insertions'] < min_size:
            continue

        domains.append({
            **domain,
            'confidence': 'high' if domain['coherence'] > 0.8 else 'medium',
        })

    return domains


def compare_conditions(
    landscape_a: TnSeqLandscape,
    landscape_b: TnSeqLandscape,
    window: int = 20,
) -> Dict[str, Any]:
    """
    Compare TnSeq coherence between two conditions.

    Useful for identifying:
    - Condition-specific essential genes
    - Regions that become stable/unstable under selection
    - Fitness landscape changes

    Args:
        landscape_a: First condition landscape.
        landscape_b: Second condition landscape.
        window: Coherence window size.

    Returns:
        Comparison results dictionary.
    """
    result_a = analyze_tnseq_coherence(landscape_a, window=window)
    result_b = analyze_tnseq_coherence(landscape_b, window=window)

    # Find common positions
    common_positions = np.intersect1d(landscape_a.positions, landscape_b.positions)

    # Compare fitness
    mask_a = np.isin(landscape_a.positions, common_positions)
    mask_b = np.isin(landscape_b.positions, common_positions)

    fitness_correlation = np.corrcoef(
        landscape_a.fitness[mask_a],
        landscape_b.fitness[mask_b],
    )[0, 1]

    return {
        'condition_a': landscape_a.condition,
        'condition_b': landscape_b.condition,
        'gene_id': landscape_a.gene_id,
        'comparison': {
            'n_common_insertions': len(common_positions),
            'fitness_correlation': float(fitness_correlation),
            'coherence_a': result_a.correlation,
            'coherence_b': result_b.correlation,
            'essential_domains_a': len(result_a.essential_domains),
            'essential_domains_b': len(result_b.essential_domains),
        },
        'condition_specific': {
            'a_essential_only': [
                d for d in result_a.essential_domains
                if not any(
                    d['start'] <= e['end'] and d['end'] >= e['start']
                    for e in result_b.essential_domains
                )
            ],
            'b_essential_only': [
                d for d in result_b.essential_domains
                if not any(
                    d['start'] <= e['end'] and d['end'] >= e['start']
                    for e in result_a.essential_domains
                )
            ],
        },
    }
