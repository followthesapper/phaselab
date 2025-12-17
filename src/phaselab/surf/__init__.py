"""
PhaseLab SURF: CRISPR-SURF integration for regulatory signal deconvolution.

CRISPR-SURF (Systematic Unbiased Regulatory Factor discovery) performs
deconvolution on CRISPRa/CRISPRi tiling screen data to extract the
underlying regulatory signal.

This module bridges CRISPR-SURF output with PhaseLab spatial coherence:
1. Parse SURF beta profiles
2. Apply spatial coherence analysis
3. Identify stable regulatory elements
4. Integrate with guide design

Key insight from E213:
- SURF deconvolves the noisy screen data into a cleaner regulatory signal
- Spatial coherence then identifies which regions of that signal are stable
- The combination provides higher-confidence target regions

Version: 0.1.0
"""

from .parser import (
    SURFOutput,
    parse_surf_output,
    parse_surf_regions,
    load_surf_directory,
)

from .coherence import (
    SURFCoherenceResult,
    compute_surf_coherence,
    surf_to_regulatory_landscape,
    compare_raw_vs_surf,
)

from .pipeline import (
    SURFPipeline,
    SURFPipelineConfig,
    run_surf_analysis,
    integrated_surf_pipeline,
)

from .visualization import (
    plot_surf_coherence,
    plot_raw_vs_deconvolved,
    plot_surf_regions,
)

__all__ = [
    # Parser
    "SURFOutput",
    "parse_surf_output",
    "parse_surf_regions",
    "load_surf_directory",
    # Coherence
    "SURFCoherenceResult",
    "compute_surf_coherence",
    "surf_to_regulatory_landscape",
    "compare_raw_vs_surf",
    # Pipeline
    "SURFPipeline",
    "SURFPipelineConfig",
    "run_surf_analysis",
    "integrated_surf_pipeline",
    # Visualization
    "plot_surf_coherence",
    "plot_raw_vs_deconvolved",
    "plot_surf_regions",
]
