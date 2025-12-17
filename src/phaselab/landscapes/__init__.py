"""
PhaseLab Landscapes: General perturbation-response coherence analysis.

This module implements the core insight from E212-E213:
    IR measures coherence of the SYSTEM'S RESPONSE, not the probe.

The landscapes module provides generic primitives for analyzing ANY
perturbation-response dataset, whether from:
- CRISPRa/CRISPRi tiling screens
- TnSeq/RB-TnSeq fitness data
- Deep mutational scanning
- Drug binding landscapes
- Reaction parameter sweeps
- Chromatin accessibility profiles

Core Concepts:
    ResponseLandscape: A structured representation of perturbation → response data
    SpatialCoherence: Sliding-window coherence analysis
    StabilityClass: Classification of regions as STABLE/MIXED/AMPLIFYING
    AmplificationDetector: Detection of MYC-like amplifying regimes

The key equation remains R̄ = exp(-V_φ/2), but applied to RESPONSE variance,
not probe properties.

Example:
    >>> from phaselab.landscapes import ResponseLandscape, compute_spatial_coherence
    >>>
    >>> # Load your perturbation-response data
    >>> landscape = ResponseLandscape(
    ...     coords=positions,           # Where perturbations occurred
    ...     responses=expression_changes,  # What happened
    ...     replicates=replicate_data,  # Optional: multiple measurements
    ... )
    >>>
    >>> # Compute spatial coherence
    >>> profile = compute_spatial_coherence(landscape, window=50)
    >>>
    >>> # Classify regions
    >>> regions = classify_regions(profile)
    >>> print(regions.stable_regions)  # Safe to perturb
    >>> print(regions.amplifying_regions)  # Avoid - MYC-like behavior

Version: 1.0.0
"""

from .core import (
    ResponseLandscape,
    CoherenceProfile,
    StabilityClass,
    RegionClassification,
    LandscapeMetrics,
)

from .coherence import (
    compute_spatial_coherence,
    compute_local_variance,
    coherence_variance_correlation,
    estimate_variance_reduction,
)

from .classification import (
    classify_regions,
    detect_amplifying_regions,
    stability_boundaries,
    region_summary,
)

from .amplification import (
    AmplificationDetector,
    amplification_score,
    is_amplifying,
)

from .io import (
    load_tiling_data,
    load_response_matrix,
    export_landscape,
    export_regions,
)

__all__ = [
    # Core data structures
    "ResponseLandscape",
    "CoherenceProfile",
    "StabilityClass",
    "RegionClassification",
    "LandscapeMetrics",

    # Coherence analysis
    "compute_spatial_coherence",
    "compute_local_variance",
    "coherence_variance_correlation",
    "estimate_variance_reduction",

    # Classification
    "classify_regions",
    "detect_amplifying_regions",
    "stability_boundaries",
    "region_summary",

    # Amplification detection
    "AmplificationDetector",
    "amplification_score",
    "is_amplifying",

    # I/O
    "load_tiling_data",
    "load_response_matrix",
    "export_landscape",
    "export_regions",
]

__version__ = "1.0.0"
