"""
PhaseLab SURF: End-to-end pipeline for SURF + coherence analysis.

This module provides a complete pipeline that:
1. Runs CRISPR-SURF deconvolution (via subprocess or API)
2. Applies spatial coherence analysis
3. Integrates with guide design
4. Generates targeting recommendations
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Union
from pathlib import Path
import subprocess
import tempfile
import json

from .parser import SURFOutput, parse_surf_output, parse_surf_regions
from .coherence import (
    SURFCoherenceResult,
    compute_surf_coherence,
    identify_high_confidence_targets,
)


@dataclass
class SURFPipelineConfig:
    """
    Configuration for SURF + coherence pipeline.

    Attributes:
        surf_path: Path to CRISPR-SURF executable (if running locally).
        lambda_param: SURF regularization parameter.
        scale: SURF scale parameter (guide density).
        window: Coherence window size.
        stable_threshold: Threshold for stable regions.
        min_beta: Minimum beta for high-confidence targets.
        modality: Screen modality ("CRISPRa" or "CRISPRi").
        output_dir: Directory for output files.
    """
    surf_path: Optional[str] = None
    lambda_param: float = 1.0
    scale: int = 10
    window: int = 50
    stable_threshold: float = 0.7
    min_beta: float = 0.5
    modality: str = "CRISPRa"
    output_dir: Optional[str] = None


@dataclass
class SURFPipelineResult:
    """
    Complete pipeline result.

    Attributes:
        gene_symbol: Gene analyzed.
        surf_output: SURF deconvolution results.
        coherence_result: Spatial coherence analysis.
        high_confidence_targets: Top targeting zones.
        recommendations: Text recommendations.
        files_generated: List of output files.
    """
    gene_symbol: str
    surf_output: SURFOutput
    coherence_result: SURFCoherenceResult
    high_confidence_targets: List[Dict[str, Any]]
    recommendations: List[str] = field(default_factory=list)
    files_generated: List[str] = field(default_factory=list)

    def summary(self) -> str:
        """Generate summary string."""
        lines = [
            "=" * 70,
            f"SURF PIPELINE RESULTS: {self.gene_symbol}",
            "=" * 70,
            "",
            self.coherence_result.summary(),
            "",
            "HIGH CONFIDENCE TARGETS:",
        ]

        for i, target in enumerate(self.high_confidence_targets[:5], 1):
            lines.append(
                f"  {i}. [{target['start']}, {target['end']}] "
                f"beta={target['surf_beta']:.3f} "
                f"coh={target['coherence_score']:.3f} "
                f"({target['direction']})"
            )

        if self.recommendations:
            lines.extend(["", "RECOMMENDATIONS:"])
            for rec in self.recommendations:
                lines.append(f"  - {rec}")

        lines.append("=" * 70)
        return "\n".join(lines)


class SURFPipeline:
    """
    End-to-end SURF + coherence pipeline.

    This class manages the complete workflow:
    1. (Optional) Run CRISPR-SURF deconvolution
    2. Parse SURF output
    3. Compute spatial coherence
    4. Identify high-confidence targets
    5. Generate recommendations

    Example:
        >>> config = SURFPipelineConfig(
        ...     window=50,
        ...     modality="CRISPRa",
        ... )
        >>> pipeline = SURFPipeline(config)
        >>> result = pipeline.run_from_surf_output('cd69_surf.tsv', gene_symbol='CD69')
        >>> print(result.summary())
    """

    def __init__(self, config: Optional[SURFPipelineConfig] = None):
        """Initialize pipeline with configuration."""
        self.config = config or SURFPipelineConfig()

    def run_from_surf_output(
        self,
        surf_file: Union[str, Path],
        gene_symbol: str,
        tss_position: Optional[int] = None,
        regions_file: Optional[Union[str, Path]] = None,
    ) -> SURFPipelineResult:
        """
        Run pipeline starting from existing SURF output.

        Args:
            surf_file: Path to SURF beta profile.
            gene_symbol: Gene name.
            tss_position: TSS for coordinate conversion.
            regions_file: Optional path to SURF regions file.

        Returns:
            SURFPipelineResult with complete analysis.
        """
        # Parse SURF output
        surf_output = parse_surf_output(
            surf_file,
            gene_symbol=gene_symbol,
            tss_position=tss_position,
        )

        # Parse regions if provided
        if regions_file:
            regions = parse_surf_regions(regions_file, gene_symbol=gene_symbol)
            surf_output.regions = regions

        return self._run_coherence_analysis(surf_output)

    def run_from_raw_data(
        self,
        positions: np.ndarray,
        responses: np.ndarray,
        gene_symbol: str,
        run_surf: bool = True,
    ) -> SURFPipelineResult:
        """
        Run pipeline starting from raw screen data.

        If run_surf=True, will attempt to run CRISPR-SURF deconvolution.
        Otherwise, uses raw data directly.

        Args:
            positions: Screen positions.
            responses: Screen responses (log fold change).
            gene_symbol: Gene name.
            run_surf: Whether to run SURF deconvolution.

        Returns:
            SURFPipelineResult with complete analysis.
        """
        if run_surf and self.config.surf_path:
            # Run SURF deconvolution
            surf_output = self._run_surf(positions, responses, gene_symbol)
        else:
            # Use raw data as "pseudo-SURF" output
            surf_output = SURFOutput(
                gene_symbol=gene_symbol,
                positions=positions,
                beta=responses,  # Raw responses as beta
                metadata={'source': 'raw_data', 'deconvolved': False},
            )

        return self._run_coherence_analysis(surf_output)

    def _run_coherence_analysis(
        self,
        surf_output: SURFOutput,
    ) -> SURFPipelineResult:
        """Run coherence analysis on SURF output."""

        # Compute spatial coherence
        coherence_result = compute_surf_coherence(
            surf_output,
            window=self.config.window,
            stable_threshold=self.config.stable_threshold,
            modality=self.config.modality,
        )

        # Identify high-confidence targets
        high_confidence = identify_high_confidence_targets(
            coherence_result,
            min_beta=self.config.min_beta,
            min_coherence=self.config.stable_threshold,
        )

        # Generate recommendations
        recommendations = self._generate_recommendations(
            coherence_result,
            high_confidence,
        )

        # Save outputs if output_dir specified
        files_generated = []
        if self.config.output_dir:
            files_generated = self._save_outputs(
                surf_output.gene_symbol,
                coherence_result,
                high_confidence,
            )

        return SURFPipelineResult(
            gene_symbol=surf_output.gene_symbol,
            surf_output=surf_output,
            coherence_result=coherence_result,
            high_confidence_targets=high_confidence,
            recommendations=recommendations,
            files_generated=files_generated,
        )

    def _run_surf(
        self,
        positions: np.ndarray,
        responses: np.ndarray,
        gene_symbol: str,
    ) -> SURFOutput:
        """
        Run CRISPR-SURF deconvolution.

        Requires CRISPR-SURF to be installed and surf_path configured.
        """
        if not self.config.surf_path:
            raise ValueError("SURF path not configured")

        # Create temporary input file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("position\tlogFC\n")
            for pos, resp in zip(positions, responses):
                f.write(f"{pos}\t{resp}\n")
            input_file = f.name

        # Create temporary output directory
        with tempfile.TemporaryDirectory() as output_dir:
            # Run SURF
            cmd = [
                self.config.surf_path,
                '-i', input_file,
                '-o', output_dir,
                '-l', str(self.config.lambda_param),
                '-s', str(self.config.scale),
            ]

            try:
                subprocess.run(cmd, check=True, capture_output=True)
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"SURF failed: {e.stderr.decode()}")

            # Parse output
            output_file = Path(output_dir) / "beta.tsv"
            if not output_file.exists():
                # Try to find any output file
                output_files = list(Path(output_dir).glob("*.tsv"))
                if output_files:
                    output_file = output_files[0]
                else:
                    raise FileNotFoundError("SURF did not produce output")

            return parse_surf_output(output_file, gene_symbol=gene_symbol)

    def _generate_recommendations(
        self,
        coherence_result: SURFCoherenceResult,
        high_confidence: List[Dict[str, Any]],
    ) -> List[str]:
        """Generate text recommendations from analysis."""
        recommendations = []

        # Validation status
        if coherence_result.is_validated:
            recommendations.append(
                f"VALIDATED: Spatial coherence model confirmed "
                f"(r = {coherence_result.correlation:.3f})"
            )
        else:
            recommendations.append(
                "CAUTION: Spatial coherence not strongly validated for this locus"
            )

        # High-confidence targets
        if high_confidence:
            best = high_confidence[0]
            recommendations.append(
                f"TOP TARGET: [{best['start']}, {best['end']}] "
                f"({best['direction']}, beta={best['surf_beta']:.3f})"
            )

            # Count by direction
            activating = sum(1 for t in high_confidence if t['direction'] == 'activating')
            repressing = sum(1 for t in high_confidence if t['direction'] == 'repressing')

            if self.config.modality == "CRISPRa" and activating > 0:
                recommendations.append(
                    f"CRISPRa TARGETS: {activating} activating enhancer(s) identified"
                )
            if self.config.modality == "CRISPRi" and repressing > 0:
                recommendations.append(
                    f"CRISPRi TARGETS: {repressing} repressing region(s) identified"
                )
        else:
            recommendations.append(
                "WARNING: No high-confidence targets found - consider relaxing thresholds"
            )

        # SURF region overlap
        n_stable_surf = len(coherence_result.stable_surf_regions)
        if n_stable_surf > 0:
            recommendations.append(
                f"SURF VALIDATION: {n_stable_surf} significant region(s) in stable zones"
            )

        return recommendations

    def _save_outputs(
        self,
        gene_symbol: str,
        coherence_result: SURFCoherenceResult,
        high_confidence: List[Dict[str, Any]],
    ) -> List[str]:
        """Save analysis outputs to files."""
        output_dir = Path(self.config.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        files = []

        # Save coherence profile
        profile_file = output_dir / f"{gene_symbol}_coherence_profile.tsv"
        with open(profile_file, 'w') as f:
            f.write("position\tcoherence\tlocal_variance\n")
            for pos, coh, var in zip(
                coherence_result.profile.coords,
                coherence_result.profile.coherence,
                coherence_result.profile.local_variance,
            ):
                f.write(f"{pos}\t{coh}\t{var}\n")
        files.append(str(profile_file))

        # Save regions
        regions_file = output_dir / f"{gene_symbol}_regions.tsv"
        with open(regions_file, 'w') as f:
            f.write("start\tend\tstability\tcoherence_score\tvariance_reduction\n")
            for region in coherence_result.regions:
                f.write(
                    f"{region.start}\t{region.end}\t{region.stability.value}\t"
                    f"{region.coherence_score}\t{region.variance_reduction}\n"
                )
        files.append(str(regions_file))

        # Save high-confidence targets
        targets_file = output_dir / f"{gene_symbol}_targets.json"
        with open(targets_file, 'w') as f:
            json.dump(high_confidence, f, indent=2)
        files.append(str(targets_file))

        # Save complete result
        result_file = output_dir / f"{gene_symbol}_complete.json"
        with open(result_file, 'w') as f:
            json.dump(coherence_result.to_dict(), f, indent=2, default=str)
        files.append(str(result_file))

        return files


def run_surf_analysis(
    surf_file: Union[str, Path],
    gene_symbol: str,
    window: int = 50,
    modality: str = "CRISPRa",
    **kwargs,
) -> SURFPipelineResult:
    """
    Convenience function to run SURF + coherence analysis.

    Args:
        surf_file: Path to SURF output file.
        gene_symbol: Gene name.
        window: Coherence window size.
        modality: Screen modality.
        **kwargs: Additional config parameters.

    Returns:
        SURFPipelineResult with analysis.

    Example:
        >>> result = run_surf_analysis(
        ...     'cd69_surf_beta.tsv',
        ...     gene_symbol='CD69',
        ...     modality='CRISPRa',
        ... )
        >>> print(result.summary())
    """
    config = SURFPipelineConfig(
        window=window,
        modality=modality,
        **kwargs,
    )

    pipeline = SURFPipeline(config)
    return pipeline.run_from_surf_output(surf_file, gene_symbol)


def integrated_surf_pipeline(
    screen_data: Dict[str, Any],
    gene_symbol: str,
    crispor_guides: Optional[List[Dict[str, Any]]] = None,
    modality: str = "CRISPRa",
    window: int = 50,
) -> Dict[str, Any]:
    """
    Run integrated pipeline combining SURF, coherence, and guide design.

    This is the full integration that:
    1. Takes raw screen data or SURF output
    2. Computes spatial coherence
    3. Matches with CRISPOR guides
    4. Returns comprehensive targeting recommendations

    Args:
        screen_data: Dictionary with 'positions' and 'responses' or 'surf_file'.
        gene_symbol: Gene name.
        crispor_guides: Optional CRISPOR guide data.
        modality: Screen modality.
        window: Coherence window size.

    Returns:
        Dictionary with complete integrated analysis.
    """
    config = SURFPipelineConfig(
        window=window,
        modality=modality,
    )
    pipeline = SURFPipeline(config)

    # Run SURF analysis
    if 'surf_file' in screen_data:
        result = pipeline.run_from_surf_output(
            screen_data['surf_file'],
            gene_symbol,
            tss_position=screen_data.get('tss_position'),
        )
    else:
        result = pipeline.run_from_raw_data(
            screen_data['positions'],
            screen_data['responses'],
            gene_symbol,
            run_surf=screen_data.get('run_surf', False),
        )

    # Integrate with CRISPOR guides if provided
    integrated_guides = []
    if crispor_guides:
        from ..spatial.integration import integrate_with_crispor
        from ..spatial.regulatory import RegulatoryLandscape

        landscape = RegulatoryLandscape(
            positions=result.surf_output.positions,
            responses=result.surf_output.beta,
            gene_symbol=gene_symbol,
            modality=modality,
        )

        integrated_guides = integrate_with_crispor(
            landscape,
            crispor_guides,
            regions=result.coherence_result.regions,
        )

    return {
        'gene_symbol': gene_symbol,
        'surf_result': result,
        'coherence': {
            'correlation': result.coherence_result.correlation,
            'is_validated': result.coherence_result.is_validated,
            'n_stable_regions': result.coherence_result.n_stable_regions,
        },
        'high_confidence_targets': result.high_confidence_targets,
        'integrated_guides': integrated_guides,
        'recommendations': result.recommendations,
    }
