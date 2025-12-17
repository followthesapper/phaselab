"""
PhaseLab Chem: Reaction landscape analysis.

Analyzes reaction optimization data to identify:
- Stable reaction conditions (reproducible yields/selectivity)
- Optimal parameter ranges
- Coherent response regions

The key insight:
- Reaction parameters (temperature, concentration, catalyst) = perturbation coordinates
- Yield, selectivity, rate = response signals
- Spatial coherence identifies where optimization is reliable

Applications:
- Process chemistry: Which conditions give reproducible results?
- Enzyme engineering: Which mutations predictably affect activity?
- Catalyst screening: Which catalysts behave consistently?
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
class ReactionLandscape:
    """
    Reaction parameter-response landscape.

    Represents how reaction outcome varies across:
    - Temperature
    - Concentration
    - Catalyst loading
    - Time
    - Other parameters

    Attributes:
        parameters: Parameter values (single dimension for now)
        responses: Response values (yield, selectivity, rate)
        parameter_name: Name of the parameter
        response_name: Name of the response
        reaction_name: Reaction identifier
        conditions: Fixed conditions dictionary
        metadata: Additional metadata
    """
    parameters: np.ndarray
    responses: np.ndarray
    parameter_name: str = "parameter"
    response_name: str = "yield"
    reaction_name: str = "unknown"
    conditions: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_points(self) -> int:
        """Number of data points."""
        return len(self.parameters)

    @property
    def parameter_range(self) -> Tuple[float, float]:
        """Parameter range."""
        return float(self.parameters.min()), float(self.parameters.max())

    @property
    def response_range(self) -> Tuple[float, float]:
        """Response range."""
        return float(np.nanmin(self.responses)), float(np.nanmax(self.responses))

    @property
    def best_response(self) -> Tuple[float, float]:
        """Best response and its parameter value."""
        idx = np.nanargmax(self.responses)
        return float(self.parameters[idx]), float(self.responses[idx])

    def to_response_landscape(self) -> ResponseLandscape:
        """Convert to generic ResponseLandscape."""
        return ResponseLandscape(
            coords=self.parameters,
            responses=self.responses,
            metadata={
                'type': 'reaction',
                'parameter': self.parameter_name,
                'response': self.response_name,
                'reaction': self.reaction_name,
                **self.metadata,
            },
        )


@dataclass
class ReactionResult:
    """
    Result of reaction coherence analysis.

    Attributes:
        landscape: Original reaction landscape
        profile: Coherence profile
        stable_windows: Stable parameter windows
        optimal_window: Best stable window for optimization
        recommendations: Optimization recommendations
        validation: Validation statistics
    """
    landscape: ReactionLandscape
    profile: CoherenceProfile
    stable_windows: List[Dict[str, Any]]
    optimal_window: Optional[Dict[str, Any]] = None
    recommendations: List[str] = field(default_factory=list)
    validation: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_validated(self) -> bool:
        """Whether coherence model is validated."""
        return self.profile.is_validated

    @property
    def has_stable_optimum(self) -> bool:
        """Whether optimal conditions are in a stable window."""
        if not self.optimal_window:
            return False
        return self.optimal_window.get('contains_optimum', False)

    def summary(self) -> str:
        """Generate summary string."""
        lines = [
            "=" * 60,
            f"REACTION ANALYSIS: {self.landscape.reaction_name}",
            "=" * 60,
            "",
            f"Parameter: {self.landscape.parameter_name}",
            f"Response: {self.landscape.response_name}",
            f"Data points: {self.landscape.n_points}",
            "",
            f"Best response: {self.landscape.best_response[1]:.3f}",
            f"  at {self.landscape.parameter_name} = {self.landscape.best_response[0]:.3f}",
            "",
            "COHERENCE:",
            f"  Correlation: {self.profile.correlation:.3f}",
            f"  Validated: {'YES' if self.is_validated else 'NO'}",
            "",
            f"Stable windows: {len(self.stable_windows)}",
            f"Optimum in stable window: {'YES' if self.has_stable_optimum else 'NO'}",
        ]

        if self.optimal_window:
            ow = self.optimal_window
            lines.extend([
                "",
                "RECOMMENDED WINDOW:",
                f"  Range: [{ow['start']:.3f}, {ow['end']:.3f}]",
                f"  Mean response: {ow['mean_response']:.3f}",
                f"  Coherence: {ow['coherence']:.3f}",
            ])

        if self.recommendations:
            lines.extend(["", "RECOMMENDATIONS:"])
            for rec in self.recommendations:
                lines.append(f"  - {rec}")

        lines.append("=" * 60)
        return "\n".join(lines)


def load_reaction_data(
    filepath: Union[str, Path],
    parameter_col: str = "parameter",
    response_col: str = "response",
    delimiter: str = "\t",
    parameter_name: str = "parameter",
    response_name: str = "yield",
    reaction_name: str = "unknown",
) -> ReactionLandscape:
    """
    Load reaction data from a file.

    Args:
        filepath: Path to data file.
        parameter_col: Parameter column name.
        response_col: Response column name.
        delimiter: Column delimiter.
        parameter_name: Name of the parameter.
        response_name: Name of the response.
        reaction_name: Reaction identifier.

    Returns:
        ReactionLandscape for analysis.
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"Reaction data not found: {filepath}")

    try:
        import pandas as pd
        df = pd.read_csv(filepath, delimiter=delimiter)

        parameters = df[parameter_col].values
        responses = df[response_col].values

        return ReactionLandscape(
            parameters=parameters.astype(float),
            responses=responses.astype(float),
            parameter_name=parameter_name,
            response_name=response_name,
            reaction_name=reaction_name,
            metadata={'source_file': str(filepath)},
        )

    except ImportError:
        raise ImportError("pandas required for reaction data loading")


def analyze_reaction_coherence(
    landscape: ReactionLandscape,
    window: int = 5,
    stable_threshold: float = 0.7,
) -> ReactionResult:
    """
    Analyze spatial coherence in reaction landscape.

    Identifies:
    - Stable parameter windows (reproducible outcomes)
    - Optimal window containing best results
    - Recommendations for robust optimization

    Args:
        landscape: Reaction landscape to analyze.
        window: Window size for coherence.
        stable_threshold: Threshold for stable windows.

    Returns:
        ReactionResult with complete analysis.

    Example:
        >>> landscape = load_reaction_data('temperature_screen.tsv')
        >>> result = analyze_reaction_coherence(landscape)
        >>> if result.has_stable_optimum:
        ...     print(f"Optimize in: {result.optimal_window}")
    """
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

    # Extract stable windows
    stable_windows = []
    for region in classification.regions:
        start, end, stability, score = region
        if stability == StabilityClass.STABLE:
            mask = (landscape.parameters >= start) & (landscape.parameters <= end)
            region_responses = landscape.responses[mask]

            stable_windows.append({
                'start': float(start),
                'end': float(end),
                'coherence': float(score),
                'mean_response': float(np.nanmean(region_responses)),
                'max_response': float(np.nanmax(region_responses)),
                'variance': float(np.nanvar(region_responses)),
                'n_points': int(np.sum(mask)),
            })

    # Find optimal window
    optimal_window = identify_stable_conditions(landscape, stable_windows)

    # Generate recommendations
    recommendations = _generate_recommendations(landscape, stable_windows, optimal_window)

    validation = {
        'coherence_variance_correlation': profile.correlation,
        'p_value': profile.p_value,
        'is_validated': profile.is_validated,
    }

    return ReactionResult(
        landscape=landscape,
        profile=profile,
        stable_windows=stable_windows,
        optimal_window=optimal_window,
        recommendations=recommendations,
        validation=validation,
    )


def identify_stable_conditions(
    landscape: ReactionLandscape,
    stable_windows: List[Dict[str, Any]],
) -> Optional[Dict[str, Any]]:
    """
    Identify the best stable window for optimization.

    Prioritizes windows that:
    1. Contain high response values
    2. Have high coherence
    3. Are large enough for practical optimization

    Args:
        landscape: Reaction landscape.
        stable_windows: List of stable windows.

    Returns:
        Best stable window or None.
    """
    if not stable_windows:
        return None

    # Find global optimum
    best_param, best_response = landscape.best_response

    # Score windows
    best_window = None
    best_score = -np.inf

    for window in stable_windows:
        # Check if window contains optimum
        contains_optimum = window['start'] <= best_param <= window['end']

        # Score by:
        # - Mean response
        # - Coherence
        # - Contains optimum (bonus)
        score = (
            window['mean_response'] *
            window['coherence'] *
            (2.0 if contains_optimum else 1.0)
        )

        if score > best_score:
            best_score = score
            best_window = {
                **window,
                'contains_optimum': contains_optimum,
                'optimum_param': best_param if contains_optimum else None,
                'optimum_response': best_response if contains_optimum else None,
            }

    return best_window


def _generate_recommendations(
    landscape: ReactionLandscape,
    stable_windows: List[Dict[str, Any]],
    optimal_window: Optional[Dict[str, Any]],
) -> List[str]:
    """Generate optimization recommendations."""
    recommendations = []

    best_param, best_response = landscape.best_response

    if not stable_windows:
        recommendations.append(
            f"WARNING: No stable windows identified. "
            f"Results at {landscape.parameter_name}={best_param:.3f} may not be reproducible."
        )
        return recommendations

    if optimal_window and optimal_window.get('contains_optimum'):
        recommendations.append(
            f"OPTIMAL: Best {landscape.response_name} ({best_response:.3f}) "
            f"is in a stable window at {landscape.parameter_name}={best_param:.3f}"
        )
        recommendations.append(
            f"RECOMMENDED RANGE: {landscape.parameter_name} = "
            f"[{optimal_window['start']:.3f}, {optimal_window['end']:.3f}]"
        )
    else:
        # Optimum not in stable window
        recommendations.append(
            f"CAUTION: Global optimum at {landscape.parameter_name}={best_param:.3f} "
            f"is NOT in a stable window. Results may be irreproducible."
        )

        if stable_windows:
            best_stable = max(stable_windows, key=lambda w: w['mean_response'])
            recommendations.append(
                f"ALTERNATIVE: Consider {landscape.parameter_name} = "
                f"[{best_stable['start']:.3f}, {best_stable['end']:.3f}] "
                f"(mean {landscape.response_name}={best_stable['mean_response']:.3f}, "
                f"coherence={best_stable['coherence']:.3f})"
            )

    return recommendations


def optimize_with_stability(
    landscape: ReactionLandscape,
    target_response: Optional[float] = None,
    min_coherence: float = 0.7,
) -> Dict[str, Any]:
    """
    Find optimal conditions that are also stable.

    Balances:
    - Response magnitude (higher is better)
    - Stability (coherence)
    - Reproducibility (low variance)

    Args:
        landscape: Reaction landscape.
        target_response: Target response value (if None, maximize).
        min_coherence: Minimum coherence for stable windows.

    Returns:
        Optimization recommendation dictionary.
    """
    result = analyze_reaction_coherence(landscape)

    # Filter to sufficiently stable windows
    stable_windows = [w for w in result.stable_windows if w['coherence'] >= min_coherence]

    if not stable_windows:
        return {
            'success': False,
            'message': "No sufficiently stable windows found",
            'recommendation': None,
        }

    if target_response is not None:
        # Find window closest to target
        best_window = min(
            stable_windows,
            key=lambda w: abs(w['mean_response'] - target_response)
        )
    else:
        # Maximize response
        best_window = max(stable_windows, key=lambda w: w['mean_response'])

    return {
        'success': True,
        'window': best_window,
        'parameter_range': [best_window['start'], best_window['end']],
        'expected_response': best_window['mean_response'],
        'coherence': best_window['coherence'],
        'confidence': 'high' if best_window['coherence'] > 0.8 else 'medium',
        'message': (
            f"Recommend {landscape.parameter_name} in "
            f"[{best_window['start']:.3f}, {best_window['end']:.3f}] "
            f"for stable {landscape.response_name}={best_window['mean_response']:.3f}"
        ),
    }
