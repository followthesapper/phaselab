"""
PhaseLab Microbio: Drug response landscape analysis.

This module analyzes dose-response and drug exposure data to identify:
- Stable dosing windows (reproducible pharmacodynamics)
- Unstable concentration ranges (high variability)
- Therapeutic windows with coherent response

The key insight:
- Drug concentration is a perturbation coordinate
- Cell/organism response is the measured signal
- Spatial coherence identifies stable pharmacology regions

This applies the E213 methodology to pharmacology.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Union, Tuple
from pathlib import Path

from ..landscapes.core import (
    ResponseLandscape,
    CoherenceProfile,
    StabilityClass,
)
from ..landscapes.coherence import compute_spatial_coherence
from ..landscapes.classification import classify_regions


@dataclass
class DrugResponseLandscape:
    """
    Drug dose-response landscape.

    Attributes:
        concentrations: Drug concentrations (log-scale recommended)
        responses: Response values (viability, activity, etc.)
        drug_name: Drug identifier
        cell_line: Cell line or organism
        response_type: Type of response measured
        replicates: Optional replicate measurements
        metadata: Additional metadata
    """
    concentrations: np.ndarray
    responses: np.ndarray
    drug_name: str
    cell_line: str = "unknown"
    response_type: str = "viability"
    replicates: Optional[np.ndarray] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_concentrations(self) -> int:
        """Number of concentration points."""
        return len(self.concentrations)

    @property
    def concentration_range(self) -> Tuple[float, float]:
        """Min and max concentration."""
        return float(self.concentrations.min()), float(self.concentrations.max())

    @property
    def response_range(self) -> Tuple[float, float]:
        """Min and max response."""
        return float(np.nanmin(self.responses)), float(np.nanmax(self.responses))

    def to_response_landscape(self) -> ResponseLandscape:
        """Convert to generic ResponseLandscape."""
        return ResponseLandscape(
            coords=self.concentrations,
            responses=self.responses if self.replicates is None else self.replicates,
            metadata={
                'type': 'drug_response',
                'drug': self.drug_name,
                'cell_line': self.cell_line,
                'response_type': self.response_type,
                **self.metadata,
            },
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'drug_name': self.drug_name,
            'cell_line': self.cell_line,
            'response_type': self.response_type,
            'n_concentrations': self.n_concentrations,
            'concentration_range': self.concentration_range,
            'response_range': self.response_range,
        }


@dataclass
class DrugResponseResult:
    """
    Result of drug response coherence analysis.

    Attributes:
        landscape: Original drug response landscape
        profile: Coherence profile
        stable_windows: Stable dosing windows
        therapeutic_window: Identified therapeutic window
        ic50_estimate: IC50 estimate with confidence
        validation: Validation statistics
    """
    landscape: DrugResponseLandscape
    profile: CoherenceProfile
    stable_windows: List[Dict[str, Any]]
    therapeutic_window: Optional[Dict[str, Any]] = None
    ic50_estimate: Optional[Dict[str, Any]] = None
    validation: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_validated(self) -> bool:
        """Whether coherence model is validated."""
        return self.profile.is_validated

    @property
    def has_therapeutic_window(self) -> bool:
        """Whether a therapeutic window was identified."""
        return self.therapeutic_window is not None

    def summary(self) -> str:
        """Generate summary string."""
        lines = [
            "=" * 60,
            f"DRUG RESPONSE ANALYSIS: {self.landscape.drug_name}",
            "=" * 60,
            "",
            f"Cell line: {self.landscape.cell_line}",
            f"Response type: {self.landscape.response_type}",
            f"Concentrations: {self.landscape.n_concentrations}",
            "",
            "COHERENCE:",
            f"  Correlation: {self.profile.correlation:.3f}",
            f"  Validated: {'YES' if self.is_validated else 'NO'}",
            "",
            f"Stable windows: {len(self.stable_windows)}",
        ]

        if self.therapeutic_window:
            tw = self.therapeutic_window
            lines.extend([
                "",
                "THERAPEUTIC WINDOW:",
                f"  Range: [{tw['start']:.2e}, {tw['end']:.2e}]",
                f"  Coherence: {tw['coherence']:.3f}",
                f"  Confidence: {tw['confidence']}",
            ])

        if self.ic50_estimate:
            ic50 = self.ic50_estimate
            lines.extend([
                "",
                "IC50 ESTIMATE:",
                f"  Value: {ic50['value']:.2e}",
                f"  In stable window: {'YES' if ic50.get('in_stable_window') else 'NO'}",
                f"  Confidence: {ic50.get('confidence', 'N/A')}",
            ])

        lines.append("=" * 60)
        return "\n".join(lines)


def load_dose_response(
    filepath: Union[str, Path],
    drug_name: str,
    concentration_col: str = "concentration",
    response_col: str = "response",
    replicate_cols: Optional[List[str]] = None,
    delimiter: str = "\t",
    log_transform: bool = True,
    cell_line: str = "unknown",
    response_type: str = "viability",
) -> DrugResponseLandscape:
    """
    Load dose-response data from a file.

    Args:
        filepath: Path to data file.
        drug_name: Drug identifier.
        concentration_col: Concentration column name.
        response_col: Response column name.
        replicate_cols: Replicate column names.
        delimiter: Column delimiter.
        log_transform: Whether to log-transform concentrations.
        cell_line: Cell line identifier.
        response_type: Response type (viability, activity, etc.).

    Returns:
        DrugResponseLandscape for analysis.
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"Dose-response file not found: {filepath}")

    try:
        import pandas as pd
        df = pd.read_csv(filepath, delimiter=delimiter)

        # Filter by drug if multiple drugs in file
        if 'drug' in df.columns:
            df = df[df['drug'] == drug_name]
        elif 'compound' in df.columns:
            df = df[df['compound'] == drug_name]

        concentrations = df[concentration_col].values
        responses = df[response_col].values

        if log_transform and np.all(concentrations > 0):
            concentrations = np.log10(concentrations)

        replicates = None
        if replicate_cols:
            replicates = df[replicate_cols].values

        return DrugResponseLandscape(
            concentrations=concentrations.astype(float),
            responses=responses.astype(float),
            drug_name=drug_name,
            cell_line=cell_line,
            response_type=response_type,
            replicates=replicates,
            metadata={
                'source_file': str(filepath),
                'log_transformed': log_transform,
            },
        )

    except ImportError:
        raise ImportError("pandas required for dose-response loading")


def analyze_drug_coherence(
    landscape: DrugResponseLandscape,
    window: int = 5,
    stable_threshold: float = 0.7,
) -> DrugResponseResult:
    """
    Analyze spatial coherence in dose-response curve.

    Identifies:
    - Stable dosing windows (reproducible pharmacodynamics)
    - Therapeutic window with highest confidence
    - IC50 estimate with reliability assessment

    Args:
        landscape: Drug response landscape.
        window: Window size for coherence (in concentration points).
        stable_threshold: Threshold for stable windows.

    Returns:
        DrugResponseResult with complete analysis.

    Example:
        >>> landscape = load_dose_response('drug_data.tsv', 'CompoundA')
        >>> result = analyze_drug_coherence(landscape)
        >>> if result.has_therapeutic_window:
        ...     print(f"Therapeutic window: {result.therapeutic_window}")
    """
    # Convert to response landscape
    response_landscape = landscape.to_response_landscape()

    # Compute coherence
    profile = compute_spatial_coherence(
        response_landscape,
        window=window,
    )

    # Classify regions (concentration windows)
    classification = classify_regions(
        response_landscape,
        profile=profile,
        stable_threshold=stable_threshold,
    )

    # Extract stable windows
    stable_windows = []
    for region in classification.regions:
        start, end, stability, score = region
        if stability == StabilityClass.STABLE:
            mask = (landscape.concentrations >= start) & (landscape.concentrations <= end)
            stable_windows.append({
                'start': float(start),
                'end': float(end),
                'coherence': float(score),
                'mean_response': float(np.nanmean(landscape.responses[mask])),
                'response_variance': float(np.nanvar(landscape.responses[mask])),
                'n_points': int(np.sum(mask)),
            })

    # Identify therapeutic window
    therapeutic_window = identify_stable_dosing_window(
        landscape, stable_windows
    )

    # Estimate IC50 with confidence
    ic50_estimate = _estimate_ic50_with_confidence(landscape, stable_windows)

    validation = {
        'coherence_variance_correlation': profile.correlation,
        'p_value': profile.p_value,
        'is_validated': profile.is_validated,
    }

    return DrugResponseResult(
        landscape=landscape,
        profile=profile,
        stable_windows=stable_windows,
        therapeutic_window=therapeutic_window,
        ic50_estimate=ic50_estimate,
        validation=validation,
    )


def identify_stable_dosing_window(
    landscape: DrugResponseLandscape,
    stable_windows: List[Dict[str, Any]],
) -> Optional[Dict[str, Any]]:
    """
    Identify the best therapeutic window.

    Criteria:
    1. High coherence (stable)
    2. Intermediate response (not maxed out)
    3. Largest span

    Args:
        landscape: Drug response landscape.
        stable_windows: List of stable windows from coherence analysis.

    Returns:
        Best therapeutic window or None.
    """
    if not stable_windows:
        return None

    # Score windows by:
    # - Coherence
    # - Response in therapeutic range (20-80% of max)
    # - Window size
    max_response = np.nanmax(landscape.responses)
    min_response = np.nanmin(landscape.responses)

    best_window = None
    best_score = 0

    for window in stable_windows:
        # Prefer intermediate response ranges
        mean_resp = window['mean_response']
        resp_range = max_response - min_response
        if resp_range > 0:
            normalized_resp = (mean_resp - min_response) / resp_range
            # Score highest for 30-70% range
            range_score = 1 - 2 * abs(normalized_resp - 0.5)
        else:
            range_score = 0.5

        # Size score (log scale)
        size_score = np.log10(max(1, window['end'] - window['start'] + 1))

        # Combined score
        score = window['coherence'] * (1 + range_score) * size_score

        if score > best_score:
            best_score = score
            best_window = window

    if best_window:
        return {
            **best_window,
            'confidence': 'high' if best_window['coherence'] > 0.8 else 'medium',
        }

    return None


def _estimate_ic50_with_confidence(
    landscape: DrugResponseLandscape,
    stable_windows: List[Dict[str, Any]],
) -> Optional[Dict[str, Any]]:
    """Estimate IC50 with coherence-based confidence."""
    responses = landscape.responses
    concentrations = landscape.concentrations

    # Simple IC50 estimation (linear interpolation at 50% response)
    max_resp = np.nanmax(responses)
    min_resp = np.nanmin(responses)
    mid_resp = (max_resp + min_resp) / 2

    # Find where response crosses 50%
    sorted_idx = np.argsort(concentrations)
    sorted_conc = concentrations[sorted_idx]
    sorted_resp = responses[sorted_idx]

    ic50 = None
    for i in range(len(sorted_resp) - 1):
        if (sorted_resp[i] >= mid_resp >= sorted_resp[i + 1] or
            sorted_resp[i] <= mid_resp <= sorted_resp[i + 1]):
            # Linear interpolation
            if sorted_resp[i + 1] != sorted_resp[i]:
                frac = (mid_resp - sorted_resp[i]) / (sorted_resp[i + 1] - sorted_resp[i])
                ic50 = sorted_conc[i] + frac * (sorted_conc[i + 1] - sorted_conc[i])
            break

    if ic50 is None:
        return None

    # Check if IC50 is in a stable window
    in_stable = False
    for window in stable_windows:
        if window['start'] <= ic50 <= window['end']:
            in_stable = True
            break

    return {
        'value': float(ic50),
        'in_stable_window': in_stable,
        'confidence': 'high' if in_stable else 'low',
        'mid_response': float(mid_resp),
    }


def pharmacodynamic_stability(
    landscape: DrugResponseLandscape,
    target_response: float,
    tolerance: float = 0.1,
) -> Dict[str, Any]:
    """
    Assess pharmacodynamic stability around a target response.

    Useful for determining:
    - Is the target response achievable in a stable region?
    - What concentration range maintains target response?
    - How sensitive is response to concentration changes?

    Args:
        landscape: Drug response landscape.
        target_response: Desired response level.
        tolerance: Acceptable deviation from target.

    Returns:
        Stability assessment dictionary.
    """
    result = analyze_drug_coherence(landscape)

    # Find concentrations that achieve target response
    target_mask = np.abs(landscape.responses - target_response) <= tolerance
    target_concentrations = landscape.concentrations[target_mask]

    if len(target_concentrations) == 0:
        return {
            'achievable': False,
            'target_response': target_response,
            'message': f"Target response {target_response} not achievable",
        }

    # Check if any target concentrations are in stable windows
    stable_target_concs = []
    for conc in target_concentrations:
        for window in result.stable_windows:
            if window['start'] <= conc <= window['end']:
                stable_target_concs.append({
                    'concentration': float(conc),
                    'window_coherence': window['coherence'],
                })
                break

    return {
        'achievable': True,
        'target_response': target_response,
        'n_concentrations_achieving_target': len(target_concentrations),
        'concentration_range': [
            float(target_concentrations.min()),
            float(target_concentrations.max()),
        ],
        'stable_concentrations': stable_target_concs,
        'in_stable_window': len(stable_target_concs) > 0,
        'recommendation': (
            f"Target achievable in stable window"
            if stable_target_concs
            else "Target achievable but NOT in stable window - expect variability"
        ),
    }
