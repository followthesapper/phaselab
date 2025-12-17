"""
PhaseLab Spatial: Regulatory landscape analysis.

Implements the E213 methodology specifically for gene regulatory regions.
This is the genomics-specific implementation of the general landscapes framework.

The key concept:
    A regulatory landscape maps positions (relative to TSS) to expression responses
    under perturbation (CRISPRa, CRISPRi, etc.). Spatial coherence of these responses
    predicts whether perturbations at those positions will have reproducible effects.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Tuple
from enum import Enum

from ..landscapes.core import ResponseLandscape, CoherenceProfile, StabilityClass
from ..landscapes.coherence import compute_spatial_coherence
from ..landscapes.classification import classify_regions, RegionClassification


class RegionType(Enum):
    """Types of regulatory regions."""
    PROMOTER = "promoter"
    ENHANCER = "enhancer"
    SILENCER = "silencer"
    INSULATOR = "insulator"
    UTR_5 = "5'UTR"
    UTR_3 = "3'UTR"
    INTRON = "intron"
    INTERGENIC = "intergenic"
    UNKNOWN = "unknown"


@dataclass
class RegulatoryRegion:
    """
    A classified regulatory region.

    Attributes:
        start: Start position (TSS-relative).
        end: End position (TSS-relative).
        stability: Stability classification (STABLE/MIXED/AMPLIFYING/IRRELEVANT).
        coherence_score: Mean coherence in this region.
        variance_reduction: Expected variance reduction vs random targeting.
        region_type: Type of regulatory element (if known).
        contains_tss: Whether this region contains the TSS.
        gene_symbol: Associated gene.
    """
    start: int
    end: int
    stability: StabilityClass
    coherence_score: float
    variance_reduction: float
    region_type: RegionType = RegionType.UNKNOWN
    contains_tss: bool = False
    gene_symbol: Optional[str] = None

    @property
    def size(self) -> int:
        """Size of the region in bp."""
        return abs(self.end - self.start)

    @property
    def is_safe(self) -> bool:
        """Whether this region is safe for targeting."""
        return self.stability == StabilityClass.STABLE

    @property
    def is_dangerous(self) -> bool:
        """Whether this region shows amplifying behavior."""
        return self.stability == StabilityClass.AMPLIFYING

    @property
    def recommendation(self) -> str:
        """Human-readable recommendation for this region."""
        if self.is_safe:
            return f"PRIMARY: Target region [{self.start}, {self.end}] - {self.variance_reduction:.0%} expected variance reduction"
        elif self.stability == StabilityClass.MIXED:
            return f"SECONDARY: Region [{self.start}, {self.end}] - context-dependent, validate before use"
        elif self.is_dangerous:
            return f"AVOID: Region [{self.start}, {self.end}] - amplifying behavior detected"
        else:
            return f"SKIP: Region [{self.start}, {self.end}] - no significant response"

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'start': self.start,
            'end': self.end,
            'stability': self.stability.value,
            'coherence_score': self.coherence_score,
            'variance_reduction': self.variance_reduction,
            'region_type': self.region_type.value,
            'contains_tss': self.contains_tss,
            'gene_symbol': self.gene_symbol,
            'size': self.size,
            'is_safe': self.is_safe,
        }


@dataclass
class RegulatoryLandscape:
    """
    A regulatory landscape for a gene.

    This extends ResponseLandscape with gene-specific metadata and
    TSS-relative coordinate handling.

    Attributes:
        positions: Array of positions relative to TSS.
        responses: Array of expression responses (logFC, etc.).
        gene_symbol: Gene name.
        tss_position: TSS position in genomic coordinates (if known).
        chromosome: Chromosome (if known).
        genome_build: Genome build (hg38, mm10, etc.).
        modality: Perturbation modality (CRISPRa, CRISPRi, etc.).
        replicates: Optional replicate data.
        metadata: Additional metadata.
    """
    positions: np.ndarray  # TSS-relative positions
    responses: np.ndarray  # Expression responses
    gene_symbol: str
    tss_position: Optional[int] = None
    chromosome: Optional[str] = None
    genome_build: str = "hg38"
    modality: str = "CRISPRa"
    replicates: Optional[np.ndarray] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Validate and normalize inputs."""
        self.positions = np.asarray(self.positions)
        self.responses = np.asarray(self.responses)

        if len(self.positions) != len(self.responses):
            raise ValueError("positions and responses must have same length")

        # Sort by position
        order = np.argsort(self.positions)
        self.positions = self.positions[order]
        self.responses = self.responses[order]
        if self.replicates is not None:
            self.replicates = np.asarray(self.replicates)[order]

    @property
    def n_positions(self) -> int:
        """Number of positions."""
        return len(self.positions)

    @property
    def position_range(self) -> Tuple[int, int]:
        """Range of positions (min, max)."""
        return int(self.positions.min()), int(self.positions.max())

    @property
    def mean_response(self) -> np.ndarray:
        """Mean response (across replicates if present)."""
        if self.replicates is not None:
            return np.nanmean(self.replicates, axis=1)
        return self.responses

    @property
    def response_variance(self) -> np.ndarray:
        """Variance of responses (requires replicates)."""
        if self.replicates is not None:
            return np.nanvar(self.replicates, axis=1)
        return np.zeros(len(self.responses))

    @property
    def coverage(self) -> float:
        """Average positions per bp in the screened region."""
        pos_range = self.position_range[1] - self.position_range[0]
        if pos_range > 0:
            return self.n_positions / pos_range
        return 0.0

    def to_response_landscape(self) -> ResponseLandscape:
        """Convert to generic ResponseLandscape."""
        responses = self.replicates if self.replicates is not None else self.responses

        return ResponseLandscape(
            coords=self.positions,
            responses=responses,
            metadata={
                'gene_symbol': self.gene_symbol,
                'tss_position': self.tss_position,
                'chromosome': self.chromosome,
                'genome_build': self.genome_build,
                'modality': self.modality,
                **self.metadata,
            },
        )

    def window(self, start: int, end: int) -> 'RegulatoryLandscape':
        """Extract a window of the landscape."""
        mask = (self.positions >= start) & (self.positions <= end)

        new_replicates = None
        if self.replicates is not None:
            new_replicates = self.replicates[mask]

        return RegulatoryLandscape(
            positions=self.positions[mask],
            responses=self.responses[mask],
            gene_symbol=self.gene_symbol,
            tss_position=self.tss_position,
            chromosome=self.chromosome,
            genome_build=self.genome_build,
            modality=self.modality,
            replicates=new_replicates,
            metadata={**self.metadata, 'window': (start, end)},
        )

    @classmethod
    def from_tiling_screen(
        cls,
        positions: np.ndarray,
        responses: np.ndarray,
        gene: str,
        tss: int = 0,
        replicates: Optional[np.ndarray] = None,
        modality: str = "CRISPRa",
        **kwargs,
    ) -> 'RegulatoryLandscape':
        """
        Create RegulatoryLandscape from tiling screen data.

        Args:
            positions: Genomic positions of guides/perturbations.
            responses: Expression responses (logFC, etc.).
            gene: Gene symbol.
            tss: TSS position (to compute relative coordinates).
            replicates: Optional replicate data (n_positions x n_replicates).
            modality: Perturbation type.
            **kwargs: Additional metadata.

        Returns:
            RegulatoryLandscape instance.
        """
        # Convert to TSS-relative positions
        relative_positions = np.asarray(positions) - tss

        return cls(
            positions=relative_positions,
            responses=responses,
            gene_symbol=gene,
            tss_position=tss,
            modality=modality,
            replicates=replicates,
            metadata=kwargs,
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        data = {
            'positions': self.positions.tolist(),
            'responses': self.responses.tolist(),
            'gene_symbol': self.gene_symbol,
            'tss_position': self.tss_position,
            'chromosome': self.chromosome,
            'genome_build': self.genome_build,
            'modality': self.modality,
            'n_positions': self.n_positions,
            'position_range': self.position_range,
            'coverage': self.coverage,
            'metadata': self.metadata,
        }
        if self.replicates is not None:
            data['replicates'] = self.replicates.tolist()
        return data


def compute_regulatory_coherence(
    landscape: RegulatoryLandscape,
    window: int = 50,
    step: int = 1,
) -> CoherenceProfile:
    """
    Compute spatial coherence profile for a regulatory landscape.

    This is the E213 algorithm applied to gene regulatory regions.

    Args:
        landscape: RegulatoryLandscape to analyze.
        window: Window size in bp for coherence computation.
        step: Step size for sliding window.

    Returns:
        CoherenceProfile with spatial coherence at each position.

    Example:
        >>> profile = compute_regulatory_coherence(landscape, window=50)
        >>> print(f"Coherence-variance correlation: {profile.correlation:.3f}")
        >>> if profile.is_validated:
        ...     print("Spatial coherence predicts stability!")
    """
    response_landscape = landscape.to_response_landscape()
    return compute_spatial_coherence(response_landscape, window=window, step=step)


def classify_regulatory_regions(
    landscape: RegulatoryLandscape,
    window: int = 50,
    stable_threshold: float = 0.7,
    mixed_threshold: float = 0.4,
    min_region_size: int = 20,
) -> List[RegulatoryRegion]:
    """
    Classify regions of a regulatory landscape by stability.

    This is the main analysis function that identifies safe targeting regions.

    Args:
        landscape: RegulatoryLandscape to analyze.
        window: Window size for coherence computation.
        stable_threshold: Coherence threshold for STABLE classification.
        mixed_threshold: Coherence threshold for MIXED.
        min_region_size: Minimum region size in bp.

    Returns:
        List of RegulatoryRegion objects with classifications.

    Example:
        >>> regions = classify_regulatory_regions(landscape)
        >>> safe_regions = [r for r in regions if r.is_safe]
        >>> print(f"Found {len(safe_regions)} safe targeting regions")
    """
    response_landscape = landscape.to_response_landscape()
    profile = compute_spatial_coherence(response_landscape, window=window)

    classification = classify_regions(
        response_landscape,
        profile=profile,
        window=window,
        stable_threshold=stable_threshold,
        mixed_threshold=mixed_threshold,
        min_region_size=min_region_size,
    )

    # Convert to RegulatoryRegion objects
    regulatory_regions = []

    for start, end, stability_class, score in classification.regions:
        # Determine region type based on position
        region_type = _infer_region_type(start, end)

        # Check if contains TSS
        contains_tss = start <= 0 <= end

        # Estimate variance reduction
        if stability_class == StabilityClass.STABLE:
            variance_reduction = profile.variance_reduction_estimate
        else:
            variance_reduction = 0.0

        regulatory_regions.append(RegulatoryRegion(
            start=int(start),
            end=int(end),
            stability=stability_class,
            coherence_score=score,
            variance_reduction=variance_reduction,
            region_type=region_type,
            contains_tss=contains_tss,
            gene_symbol=landscape.gene_symbol,
        ))

    return regulatory_regions


def _infer_region_type(start: int, end: int) -> RegionType:
    """Infer regulatory region type from position relative to TSS."""
    center = (start + end) / 2

    # Typical CRISPRa window: -400 to -50
    if -500 <= center <= 0:
        return RegionType.PROMOTER
    elif center < -500:
        return RegionType.ENHANCER  # Could be distal enhancer
    elif 0 < center <= 200:
        return RegionType.UTR_5
    else:
        return RegionType.UNKNOWN


def optimal_crispra_window(
    landscape: RegulatoryLandscape,
    candidate_windows: Optional[List[int]] = None,
) -> Dict[str, Any]:
    """
    Find the optimal CRISPRa targeting window.

    The optimal window maximizes the negative correlation between
    spatial coherence and response variance.

    Args:
        landscape: RegulatoryLandscape to analyze.
        candidate_windows: Windows to test (default: standard CRISPRa windows).

    Returns:
        Dictionary with optimal window and analysis results.
    """
    if candidate_windows is None:
        # Standard CRISPRa windows relative to TSS
        candidate_windows = [
            (-400, -50),   # Standard CRISPRa
            (-300, -50),   # Narrower
            (-500, -100),  # Wider
            (-200, 0),     # Proximal
        ]

    results = []

    for start, end in candidate_windows:
        try:
            window_landscape = landscape.window(start, end)
            if window_landscape.n_positions < 20:
                continue

            profile = compute_regulatory_coherence(window_landscape, window=30)

            results.append({
                'window': (start, end),
                'correlation': profile.correlation,
                'p_value': profile.p_value,
                'variance_reduction': profile.variance_reduction_estimate,
                'n_positions': window_landscape.n_positions,
                'is_validated': profile.is_validated,
            })
        except Exception:
            continue

    if not results:
        return {'error': 'No valid windows found'}

    # Find best window (most negative correlation with significance)
    validated = [r for r in results if r['is_validated']]
    if validated:
        best = min(validated, key=lambda x: x['correlation'])
    else:
        best = min(results, key=lambda x: x['correlation'])

    return {
        'optimal_window': best['window'],
        'all_results': results,
        'best_correlation': best['correlation'],
        'best_variance_reduction': best['variance_reduction'],
    }


def summarize_landscape(landscape: RegulatoryLandscape) -> str:
    """
    Generate a human-readable summary of a regulatory landscape.

    Args:
        landscape: RegulatoryLandscape to summarize.

    Returns:
        Formatted summary string.
    """
    regions = classify_regulatory_regions(landscape)
    profile = compute_regulatory_coherence(landscape)

    safe_regions = [r for r in regions if r.is_safe]
    dangerous_regions = [r for r in regions if r.is_dangerous]

    lines = [
        "=" * 60,
        f"REGULATORY LANDSCAPE SUMMARY: {landscape.gene_symbol}",
        "=" * 60,
        f"Modality: {landscape.modality}",
        f"Genome: {landscape.genome_build}",
        f"Positions: {landscape.n_positions} ({landscape.position_range[0]} to {landscape.position_range[1]})",
        f"Coverage: {landscape.coverage:.2f} positions/bp",
        "",
        "COHERENCE ANALYSIS:",
        f"  Coherence-variance correlation: r = {profile.correlation:.3f} (p = {profile.p_value:.2e})",
        f"  Expected variance reduction: {profile.variance_reduction_estimate:.1%}",
        f"  Validation: {'PASSED' if profile.is_validated else 'NOT VALIDATED'}",
        "",
        "REGION CLASSIFICATION:",
        f"  Total regions: {len(regions)}",
        f"  STABLE (safe): {len(safe_regions)}",
        f"  AMPLIFYING (avoid): {len(dangerous_regions)}",
    ]

    if safe_regions:
        lines.append("")
        lines.append("RECOMMENDED TARGET REGIONS:")
        for r in safe_regions[:3]:
            lines.append(f"  [{r.start}, {r.end}] - score={r.coherence_score:.3f}, {r.variance_reduction:.0%} var reduction")

    if dangerous_regions:
        lines.append("")
        lines.append("WARNING - AMPLIFYING REGIONS (avoid):")
        for r in dangerous_regions[:3]:
            lines.append(f"  [{r.start}, {r.end}] - MYC-like behavior detected")

    lines.append("=" * 60)
    return "\n".join(lines)
