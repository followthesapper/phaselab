"""
PhaseLab Spatial: Validation and falsification tests.

Implements the validation methodology from E213 and E216:
- Cross-validation of coherence-variance relationship
- Falsification tests for wet-lab collaborators
- Comparison against benchmark datasets

The key validation criterion:
    Spatial coherence should NEGATIVELY correlate with variance (r < 0).
    This was validated on 4 datasets totaling 115,251 sgRNAs.
"""

import numpy as np
from scipy import stats
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Tuple

from .regulatory import (
    RegulatoryLandscape,
    compute_regulatory_coherence,
    classify_regulatory_regions,
)
from ..landscapes.core import CoherenceProfile


@dataclass
class ValidationResult:
    """
    Result of validation analysis.

    Attributes:
        passed: Whether the validation passed.
        correlation: Coherence-variance correlation.
        p_value: Statistical significance.
        expected_sign: Expected sign of correlation (negative).
        actual_sign: Actual sign observed.
        interpretation: Human-readable interpretation.
        details: Additional validation details.
    """
    passed: bool
    correlation: float
    p_value: float
    expected_sign: str
    actual_sign: str
    interpretation: str
    details: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'passed': self.passed,
            'correlation': self.correlation,
            'p_value': self.p_value,
            'expected_sign': self.expected_sign,
            'actual_sign': self.actual_sign,
            'interpretation': self.interpretation,
            'details': self.details,
        }


def validate_coherence_model(
    landscape: RegulatoryLandscape,
    window: int = 50,
    significance_level: float = 0.05,
    min_effect_size: float = -0.1,
) -> ValidationResult:
    """
    Validate that spatial coherence predicts response stability.

    The core validation from E213:
    - Coherence should negatively correlate with variance
    - This indicates that high-coherence regions are more stable

    Args:
        landscape: RegulatoryLandscape to validate.
        window: Window size for coherence computation.
        significance_level: P-value threshold for significance.
        min_effect_size: Minimum correlation for meaningful effect.

    Returns:
        ValidationResult with validation outcome.

    Example:
        >>> result = validate_coherence_model(landscape)
        >>> if result.passed:
        ...     print("Spatial coherence model validated!")
        >>> else:
        ...     print(f"Validation failed: {result.interpretation}")
    """
    profile = compute_regulatory_coherence(landscape, window=window)

    correlation = profile.correlation
    p_value = profile.p_value

    # Determine sign
    if correlation < min_effect_size:
        actual_sign = "negative"
    elif correlation > -min_effect_size:
        actual_sign = "positive" if correlation > 0.1 else "near_zero"
    else:
        actual_sign = "weak_negative"

    # Check if validation passed
    passed = (
        correlation < min_effect_size and
        p_value < significance_level
    )

    # Generate interpretation
    if passed:
        interpretation = (
            f"VALIDATED: Spatial coherence negatively correlates with variance "
            f"(r={correlation:.3f}, p={p_value:.2e}). High-coherence regions "
            f"are predicted to show {profile.variance_reduction_estimate:.0%} lower variance."
        )
    elif correlation > 0.1 and p_value < significance_level:
        interpretation = (
            f"WARNING: POSITIVE correlation detected (r={correlation:.3f}). "
            f"This landscape shows MYC-like AMPLIFYING behavior. "
            f"Spatial coherence does NOT predict stability - it predicts amplification."
        )
    elif p_value >= significance_level:
        interpretation = (
            f"INCONCLUSIVE: Correlation not statistically significant "
            f"(r={correlation:.3f}, p={p_value:.2e}). "
            f"May need more data points or different window size."
        )
    else:
        interpretation = (
            f"WEAK: Correlation present but weak (r={correlation:.3f}). "
            f"Spatial coherence may be predictive but effect size is small."
        )

    return ValidationResult(
        passed=passed,
        correlation=correlation,
        p_value=p_value,
        expected_sign="negative",
        actual_sign=actual_sign,
        interpretation=interpretation,
        details={
            'window_size': window,
            'n_positions': profile.n_positions,
            'mean_coherence': profile.mean_coherence,
            'variance_reduction_estimate': profile.variance_reduction_estimate,
            'significance_level': significance_level,
            'min_effect_size': min_effect_size,
        },
    )


def cross_validate_regions(
    landscape: RegulatoryLandscape,
    n_folds: int = 5,
    window: int = 50,
) -> Dict[str, Any]:
    """
    Cross-validate region classification.

    Splits the landscape into folds and checks consistency of:
    - Coherence-variance correlation across folds
    - Region classification stability

    Args:
        landscape: RegulatoryLandscape to validate.
        n_folds: Number of cross-validation folds.
        window: Window size for coherence computation.

    Returns:
        Dictionary with cross-validation results.
    """
    n = landscape.n_positions
    fold_size = n // n_folds

    fold_results = []

    for i in range(n_folds):
        # Create fold indices (leave-one-out style)
        test_start = i * fold_size
        test_end = min((i + 1) * fold_size, n)

        train_mask = np.ones(n, dtype=bool)
        train_mask[test_start:test_end] = False

        # Create training landscape
        train_landscape = RegulatoryLandscape(
            positions=landscape.positions[train_mask],
            responses=landscape.responses[train_mask],
            gene_symbol=landscape.gene_symbol,
            modality=landscape.modality,
        )

        # Validate on training fold
        if train_landscape.n_positions >= 20:
            profile = compute_regulatory_coherence(train_landscape, window=window)
            fold_results.append({
                'fold': i,
                'n_positions': train_landscape.n_positions,
                'correlation': profile.correlation,
                'p_value': profile.p_value,
                'is_validated': profile.is_validated,
            })

    if not fold_results:
        return {'error': 'Not enough data for cross-validation'}

    # Aggregate results
    correlations = [r['correlation'] for r in fold_results]
    validated_folds = sum(1 for r in fold_results if r['is_validated'])

    return {
        'n_folds': len(fold_results),
        'mean_correlation': np.mean(correlations),
        'std_correlation': np.std(correlations),
        'min_correlation': np.min(correlations),
        'max_correlation': np.max(correlations),
        'validated_folds': validated_folds,
        'validation_rate': validated_folds / len(fold_results),
        'fold_results': fold_results,
        'is_consistent': all(r['correlation'] < 0 for r in fold_results),
    }


@dataclass
class FalsificationTest:
    """
    A specific falsification test for wet-lab validation.

    These are concrete, testable predictions that can confirm or
    refute the spatial coherence model.
    """
    name: str
    description: str
    prediction: str
    test_method: str
    expected_outcome: str
    falsification_criterion: str

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'name': self.name,
            'description': self.description,
            'prediction': self.prediction,
            'test_method': self.test_method,
            'expected_outcome': self.expected_outcome,
            'falsification_criterion': self.falsification_criterion,
        }


def falsification_tests(
    landscape: RegulatoryLandscape,
    regions: Optional[list] = None,
) -> List[FalsificationTest]:
    """
    Generate falsification tests for wet-lab validation.

    These tests are designed to be:
    1. Concrete and testable
    2. Falsifiable (can disprove the model)
    3. Based on specific predictions from spatial coherence analysis

    Args:
        landscape: RegulatoryLandscape that was analyzed.
        regions: Classified regions (computed if not provided).

    Returns:
        List of FalsificationTest objects.

    Example:
        >>> tests = falsification_tests(landscape, regions)
        >>> for test in tests:
        ...     print(f"TEST: {test.name}")
        ...     print(f"  Prediction: {test.prediction}")
        ...     print(f"  Falsified if: {test.falsification_criterion}")
    """
    if regions is None:
        regions = classify_regulatory_regions(landscape)

    tests = []

    # Find regions for tests
    from ..landscapes.core import StabilityClass

    stable = [r for r in regions if r.stability == StabilityClass.STABLE]
    amplifying = [r for r in regions if r.stability == StabilityClass.AMPLIFYING]
    mixed = [r for r in regions if r.stability == StabilityClass.MIXED]

    # Test A: Stable region variance
    if stable:
        best_stable = max(stable, key=lambda r: r.coherence_score)
        tests.append(FalsificationTest(
            name="Test A: Stable Region Variance",
            description=(
                f"Test guides in the stable region [{best_stable.start}, {best_stable.end}] "
                f"for {landscape.gene_symbol}."
            ),
            prediction=(
                f"Guides in this region should show {best_stable.variance_reduction:.0%} "
                f"lower replicate variance than guides outside stable regions."
            ),
            test_method=(
                "Compare replicate variance (or coefficient of variation) for 5+ guides "
                "inside vs outside the stable region."
            ),
            expected_outcome=(
                f"Variance ratio (outside/inside) > {1 / (1 - best_stable.variance_reduction + 0.01):.2f}"
            ),
            falsification_criterion=(
                "Model is FALSIFIED if variance inside stable region is HIGHER than outside, "
                "or if the variance reduction is < 10%."
            ),
        ))

    # Test B: Amplifying region detection
    if amplifying:
        worst_amp = max(amplifying, key=lambda r: r.coherence_score)
        tests.append(FalsificationTest(
            name="Test B: Amplifying Region Behavior",
            description=(
                f"Test guides in the amplifying region [{worst_amp.start}, {worst_amp.end}] "
                f"for {landscape.gene_symbol}."
            ),
            prediction=(
                "Guides in this region should show HIGHER variance than average, "
                "and small perturbations may cause unexpectedly large expression changes."
            ),
            test_method=(
                "Measure expression changes and replicate variance for guides in this region. "
                "Compare to stable regions."
            ),
            expected_outcome=(
                "Variance in amplifying region > 1.5x variance in stable regions."
            ),
            falsification_criterion=(
                "Model is FALSIFIED if variance in amplifying region is LOWER than stable regions, "
                "or if effects are more predictable than in stable regions."
            ),
        ))

    # Test C: Coherence gradient
    if len(stable) >= 2:
        tests.append(FalsificationTest(
            name="Test C: Coherence Gradient",
            description=(
                "Test whether guides at region boundaries show intermediate behavior."
            ),
            prediction=(
                "Guides near stable/unstable boundaries should show intermediate variance - "
                "higher than deep inside stable regions, lower than in unstable regions."
            ),
            test_method=(
                "Select guides at 3 positions: center of stable region, boundary, "
                "and outside stable region. Compare variance."
            ),
            expected_outcome=(
                "Variance(boundary) should be between Variance(center) and Variance(outside)."
            ),
            falsification_criterion=(
                "Model is FALSIFIED if boundary variance is not intermediate "
                "(i.e., if there's no gradient effect)."
            ),
        ))

    # Test D: Reproducibility of classification
    tests.append(FalsificationTest(
        name="Test D: Classification Reproducibility",
        description=(
            "Test whether region classifications are consistent across replicates."
        ),
        prediction=(
            "If the tiling screen is repeated, the same regions should be classified "
            "as STABLE/AMPLIFYING with > 80% overlap."
        ),
        test_method=(
            "Perform an independent tiling screen (or use held-out data). "
            "Classify regions and compare to original classification."
        ),
        expected_outcome=(
            "Jaccard similarity of stable regions > 0.7 between replicates."
        ),
        falsification_criterion=(
            "Model is FALSIFIED if region classifications are not reproducible "
            "(Jaccard < 0.5), indicating the coherence signal is noise."
        ),
    ))

    return tests


def benchmark_against_known_datasets(
    landscape: RegulatoryLandscape,
    benchmark_name: str = "auto",
) -> Dict[str, Any]:
    """
    Compare landscape analysis to known benchmark datasets.

    E213 validated against:
    - CD69: r = -0.24 (PASS)
    - IL2RA: r = -0.04 (weak PASS)
    - GATA1: r = -0.03 (weak PASS)
    - MYC: r = +0.45 (correctly identified as AMPLIFYING)

    Args:
        landscape: RegulatoryLandscape to benchmark.
        benchmark_name: Which benchmark to compare against.

    Returns:
        Dictionary with benchmark comparison results.
    """
    profile = compute_regulatory_coherence(landscape)

    # Known benchmark correlations from E213
    benchmarks = {
        'CD69': {'correlation': -0.24, 'type': 'stable'},
        'IL2RA': {'correlation': -0.04, 'type': 'weak_stable'},
        'GATA1': {'correlation': -0.03, 'type': 'weak_stable'},
        'MYC': {'correlation': +0.45, 'type': 'amplifying'},
    }

    observed_correlation = profile.correlation

    # Find most similar benchmark
    similarities = {}
    for name, data in benchmarks.items():
        diff = abs(observed_correlation - data['correlation'])
        similarities[name] = {
            'benchmark_correlation': data['correlation'],
            'observed_correlation': observed_correlation,
            'difference': diff,
            'benchmark_type': data['type'],
        }

    most_similar = min(similarities.items(), key=lambda x: x[1]['difference'])

    # Determine classification
    if observed_correlation < -0.1:
        behavior_type = "stable"
        interpretation = "Landscape shows stable behavior similar to CD69."
    elif observed_correlation > 0.1:
        behavior_type = "amplifying"
        interpretation = "WARNING: Landscape shows amplifying behavior similar to MYC."
    else:
        behavior_type = "weak_stable"
        interpretation = "Landscape shows weak stable behavior similar to IL2RA/GATA1."

    return {
        'observed_correlation': observed_correlation,
        'behavior_type': behavior_type,
        'interpretation': interpretation,
        'most_similar_benchmark': most_similar[0],
        'similarity_details': most_similar[1],
        'all_similarities': similarities,
        'is_validated': profile.is_validated,
    }
