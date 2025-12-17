"""
PhaseLab Spatial: Regulatory landscape analysis for CRISPR targeting.

This module implements the E213/E215/E216 validated methodology for
identifying stable regulatory regions for CRISPRa/CRISPRi targeting.

Key insight from E212:
    The guide is the probe, not the structure.
    IR should characterize REGIONS, not score guides.

Validated results (E213, E216):
    - 4/4 datasets passed cross-validation (CD69, IL2RA, GATA1, MYC)
    - 115,251 sgRNAs analyzed
    - 32-49% variance reduction in high-coherence regions
    - MYC correctly identified as AMPLIFYING (super-enhancer)

Usage:
    >>> from phaselab.spatial import (
    ...     RegulatoryLandscape,
    ...     compute_regulatory_coherence,
    ...     classify_regulatory_regions,
    ...     CRISPRaTargetFinder,
    ... )
    >>>
    >>> # Load tiling screen data
    >>> landscape = RegulatoryLandscape.from_tiling_screen(
    ...     positions, responses, gene='RAI1', tss=0
    ... )
    >>>
    >>> # Identify stable regions
    >>> regions = classify_regulatory_regions(landscape)
    >>> print(regions.stable_regions)  # Safe for CRISPRa targeting
    >>>
    >>> # Get guide recommendations within stable regions
    >>> finder = CRISPRaTargetFinder(landscape, regions)
    >>> targets = finder.recommend_targets(n=10)

Version: 1.0.0
"""

from .regulatory import (
    RegulatoryLandscape,
    RegulatoryRegion,
    RegionType,
    compute_regulatory_coherence,
    classify_regulatory_regions,
)

from .targeting import (
    CRISPRaTargetFinder,
    CRISPRiTargetFinder,
    TargetRecommendation,
    rank_guides_in_regions,
)

from .validation import (
    validate_coherence_model,
    cross_validate_regions,
    falsification_tests,
    ValidationResult,
)

from .integration import (
    integrate_with_crispor,
    integrate_with_chromatin,
    integrate_with_biogrid,
    IntegratedAnalysis,
)

__all__ = [
    # Core regulatory analysis
    "RegulatoryLandscape",
    "RegulatoryRegion",
    "RegionType",
    "compute_regulatory_coherence",
    "classify_regulatory_regions",

    # Targeting
    "CRISPRaTargetFinder",
    "CRISPRiTargetFinder",
    "TargetRecommendation",
    "rank_guides_in_regions",

    # Validation
    "validate_coherence_model",
    "cross_validate_regions",
    "falsification_tests",
    "ValidationResult",

    # Integration
    "integrate_with_crispor",
    "integrate_with_chromatin",
    "integrate_with_biogrid",
    "IntegratedAnalysis",
]

__version__ = "1.0.0"
