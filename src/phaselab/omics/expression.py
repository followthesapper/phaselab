"""
PhaseLab Omics: Expression coherence analysis.

Analyzes gene expression data to identify:
- Reliable expression changes
- Co-expression patterns
- Coherent transcriptional responses

Key insight:
- Gene/sample index = perturbation coordinate
- Expression level = response
- Coherence identifies reproducible expression patterns
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
class ExpressionLandscape:
    """
    Expression landscape across genes or samples.

    Attributes:
        indices: Gene or sample indices
        expression: Expression values (log2, TPM, etc.)
        labels: Gene names or sample IDs
        expression_type: Type of expression data
        condition_a: First condition
        condition_b: Second condition (for differential)
        metadata: Additional metadata
    """
    indices: np.ndarray
    expression: np.ndarray
    labels: Optional[List[str]] = None
    expression_type: str = "logFC"
    condition_a: str = "control"
    condition_b: str = "treatment"
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_genes(self) -> int:
        return len(self.indices)

    @property
    def mean_expression(self) -> float:
        return float(np.nanmean(self.expression))

    @property
    def expression_variance(self) -> float:
        return float(np.nanvar(self.expression))

    def to_response_landscape(self) -> ResponseLandscape:
        return ResponseLandscape(
            coords=self.indices,
            responses=self.expression,
            coord_labels=self.labels,
            metadata={
                'type': 'expression',
                'expression_type': self.expression_type,
                'condition_a': self.condition_a,
                'condition_b': self.condition_b,
                **self.metadata,
            },
        )


@dataclass
class ExpressionResult:
    """
    Result of expression coherence analysis.

    Attributes:
        landscape: Original expression landscape
        profile: Coherence profile
        reliable_changes: Genes with reliable expression changes
        unreliable_changes: Genes with variable expression
        validation: Validation statistics
    """
    landscape: ExpressionLandscape
    profile: CoherenceProfile
    reliable_changes: List[Dict[str, Any]]
    unreliable_changes: List[Dict[str, Any]]
    validation: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_validated(self) -> bool:
        return self.profile.is_validated

    @property
    def n_reliable(self) -> int:
        return len(self.reliable_changes)

    def summary(self) -> str:
        lines = [
            "=" * 60,
            f"EXPRESSION ANALYSIS: {self.landscape.condition_a} vs {self.landscape.condition_b}",
            "=" * 60,
            "",
            f"Total genes: {self.landscape.n_genes}",
            f"Expression type: {self.landscape.expression_type}",
            "",
            "COHERENCE:",
            f"  Correlation: {self.profile.correlation:.3f}",
            f"  Validated: {'YES' if self.is_validated else 'NO'}",
            "",
            f"Reliable changes: {self.n_reliable}",
            f"Unreliable changes: {len(self.unreliable_changes)}",
            "=" * 60,
        ]
        return "\n".join(lines)


def load_expression_data(
    filepath: Union[str, Path],
    gene_col: str = "gene",
    expression_col: str = "logFC",
    delimiter: str = "\t",
    condition_a: str = "control",
    condition_b: str = "treatment",
) -> ExpressionLandscape:
    """
    Load expression data from a file.

    Args:
        filepath: Path to expression file.
        gene_col: Gene name column.
        expression_col: Expression value column.
        delimiter: Column delimiter.
        condition_a: First condition name.
        condition_b: Second condition name.

    Returns:
        ExpressionLandscape for analysis.
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"Expression file not found: {filepath}")

    try:
        import pandas as pd
        df = pd.read_csv(filepath, delimiter=delimiter)

        indices = np.arange(len(df))
        expression = df[expression_col].values

        labels = None
        if gene_col in df.columns:
            labels = df[gene_col].tolist()

        return ExpressionLandscape(
            indices=indices.astype(float),
            expression=expression.astype(float),
            labels=labels,
            expression_type=expression_col,
            condition_a=condition_a,
            condition_b=condition_b,
            metadata={'source_file': str(filepath)},
        )

    except ImportError:
        raise ImportError("pandas required for expression loading")


def analyze_expression_coherence(
    landscape: ExpressionLandscape,
    window: int = 20,
    stable_threshold: float = 0.7,
    change_threshold: float = 1.0,
) -> ExpressionResult:
    """
    Analyze coherence in expression changes.

    For sorted gene lists (by expression or significance),
    identifies regions of coherent expression behavior.

    Args:
        landscape: Expression landscape.
        window: Window size for coherence.
        stable_threshold: Coherence threshold for reliable.
        change_threshold: Expression change threshold.

    Returns:
        ExpressionResult with analysis.
    """
    response_landscape = landscape.to_response_landscape()

    # Compute coherence
    profile = compute_spatial_coherence(
        response_landscape,
        window=window,
    )

    # Identify reliable vs unreliable changes
    reliable_changes, unreliable_changes = identify_reliable_changes(
        landscape,
        profile,
        stable_threshold=stable_threshold,
        change_threshold=change_threshold,
    )

    validation = {
        'coherence_variance_correlation': profile.correlation,
        'p_value': profile.p_value,
        'is_validated': profile.is_validated,
    }

    return ExpressionResult(
        landscape=landscape,
        profile=profile,
        reliable_changes=reliable_changes,
        unreliable_changes=unreliable_changes,
        validation=validation,
    )


def identify_reliable_changes(
    landscape: ExpressionLandscape,
    profile: CoherenceProfile,
    stable_threshold: float = 0.7,
    change_threshold: float = 1.0,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    Identify reliable vs unreliable expression changes.

    Reliable changes have:
    1. Significant expression change (above threshold)
    2. High local coherence (reproducible)

    Args:
        landscape: Expression landscape.
        profile: Coherence profile.
        stable_threshold: Minimum coherence for reliable.
        change_threshold: Minimum expression change.

    Returns:
        Tuple of (reliable_changes, unreliable_changes) lists.
    """
    reliable = []
    unreliable = []

    for i, (idx, expr) in enumerate(zip(landscape.indices, landscape.expression)):
        # Skip if below change threshold
        if abs(expr) < change_threshold:
            continue

        # Find local coherence (handle empty profile edge case)
        if len(profile.coords) == 0:
            local_coherence = 0.0
            local_variance = 1.0
        else:
            coh_idx = np.argmin(np.abs(profile.coords - idx))
            local_coherence = float(profile.coherence[coh_idx])
            local_variance = float(profile.local_variance[coh_idx])

        gene_info = {
            'index': int(idx),
            'expression': float(expr),
            'coherence': local_coherence,
            'local_variance': local_variance,
            'direction': 'up' if expr > 0 else 'down',
        }

        if landscape.labels and i < len(landscape.labels):
            gene_info['gene'] = landscape.labels[i]

        if local_coherence >= stable_threshold:
            gene_info['reliable'] = True
            reliable.append(gene_info)
        else:
            gene_info['reliable'] = False
            unreliable.append(gene_info)

    # Sort by expression magnitude
    reliable.sort(key=lambda x: abs(x['expression']), reverse=True)
    unreliable.sort(key=lambda x: abs(x['expression']), reverse=True)

    return reliable, unreliable
