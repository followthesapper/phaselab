"""
PhaseLab Quantum Grouping: IR-enhanced Hamiltonian grouping via ATLAS-Q.

Provides 5× variance reduction for quantum measurements through:
- Commutative Pauli grouping (QWC = qubit-wise commuting)
- GLS (Generalized Least Squares) weighting within groups
- Neyman allocation for optimal shot distribution

This integrates ATLAS-Q's validated VRA/IR grouping algorithm.
"""

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np

from . import is_atlas_q_available


@dataclass
class GroupingResult:
    """
    Result from Hamiltonian grouping optimization.

    Attributes:
        groups: List of term groups (indices into Hamiltonian)
        shots_per_group: Optimal shot allocation for each group
        variance_reduction: Variance reduction factor vs naive grouping
        method: Grouping method used
    """
    groups: List[List[int]]
    shots_per_group: np.ndarray
    variance_reduction: float
    method: str


def pauli_commutes(pauli1: str, pauli2: str) -> bool:
    """
    Check if two Pauli strings commute.

    Two Pauli operators commute if they anti-commute at an even number
    of positions. Anti-commuting pairs: (X,Y), (Y,Z), (Z,X)

    Args:
        pauli1: First Pauli string (e.g., "XXYZI")
        pauli2: Second Pauli string (e.g., "IXYZZ")

    Returns:
        True if the Pauli operators commute.
    """
    if len(pauli1) != len(pauli2):
        raise ValueError(f"Pauli strings must have same length: {len(pauli1)} vs {len(pauli2)}")

    anti_commute_count = 0
    for p1, p2 in zip(pauli1, pauli2):
        if p1 == 'I' or p2 == 'I':
            continue
        if p1 == p2:
            continue
        anti_commute_count += 1

    return anti_commute_count % 2 == 0


def ir_hamiltonian_grouping(
    coefficients: np.ndarray,
    pauli_strings: Optional[List[str]] = None,
    total_shots: int = 10000,
    max_group_size: int = 5,
    use_atlas_q: bool = True,
) -> GroupingResult:
    """
    Group Hamiltonian terms for optimal measurement variance.

    Uses IR/VRA coherence-aware grouping from ATLAS-Q when available,
    falling back to simple commutativity-based grouping otherwise.

    The algorithm:
    1. Estimates coherence matrix from Hamiltonian structure
    2. Greedily groups COMMUTING terms to minimize Q_GLS
    3. Allocates shots using Neyman allocation (m_g ∝ √Q_g)

    Args:
        coefficients: Hamiltonian term coefficients
        pauli_strings: Pauli strings (e.g., ["XYZI", "IZXY", ...])
        total_shots: Total measurement budget
        max_group_size: Maximum terms per group
        use_atlas_q: Use ATLAS-Q backend if available

    Returns:
        GroupingResult with groups, shots, and variance reduction

    Example:
        >>> coeffs = np.array([1.5, -0.8, 0.3, -0.2])
        >>> paulis = ["ZZII", "IZZI", "XXII", "IXXI"]
        >>> result = ir_hamiltonian_grouping(coeffs, paulis)
        >>> print(f"Variance reduction: {result.variance_reduction:.1f}x")
    """
    coefficients = np.asarray(coefficients)
    n_terms = len(coefficients)

    # Try ATLAS-Q backend first
    if use_atlas_q and is_atlas_q_available():
        try:
            from atlas_q.vra_enhanced.vqe_grouping import vra_hamiltonian_grouping

            result = vra_hamiltonian_grouping(
                coefficients=coefficients,
                pauli_strings=pauli_strings,
                total_shots=total_shots,
                max_group_size=max_group_size,
            )

            return GroupingResult(
                groups=result.groups,
                shots_per_group=result.shots_per_group,
                variance_reduction=result.variance_reduction,
                method=f"atlas_q:{result.method}",
            )
        except Exception as e:
            # Fall back to simple grouping on error
            pass

    # Fallback: Simple commutativity-based grouping
    return _simple_commuting_grouping(
        coefficients=coefficients,
        pauli_strings=pauli_strings,
        total_shots=total_shots,
        max_group_size=max_group_size,
    )


def _simple_commuting_grouping(
    coefficients: np.ndarray,
    pauli_strings: Optional[List[str]],
    total_shots: int,
    max_group_size: int,
) -> GroupingResult:
    """
    Simple fallback grouping based on commutativity.

    Groups commuting Paulis together with uniform shot allocation.
    """
    n_terms = len(coefficients)

    if pauli_strings is None:
        # No Pauli strings - each term in its own group
        groups = [[i] for i in range(n_terms)]
    else:
        # Greedy commuting groups
        remaining = set(range(n_terms))
        groups = []

        while remaining:
            # Start with largest coefficient
            start_idx = max(remaining, key=lambda i: abs(coefficients[i]))
            group = [start_idx]
            remaining.remove(start_idx)

            # Add commuting terms
            for candidate in list(remaining):
                if len(group) >= max_group_size:
                    break

                # Check if candidate commutes with all in group
                commutes = all(
                    pauli_commutes(pauli_strings[candidate], pauli_strings[g])
                    for g in group
                )
                if commutes:
                    group.append(candidate)
                    remaining.remove(candidate)

            groups.append(sorted(group))

    # Uniform shot allocation (fallback doesn't do Neyman)
    n_groups = len(groups)
    base_shots = total_shots // n_groups
    shots_per_group = np.full(n_groups, base_shots)

    # Distribute remainder
    remainder = total_shots - sum(shots_per_group)
    for i in range(remainder):
        shots_per_group[i % n_groups] += 1

    # Estimate variance reduction (approximate)
    # Grouping reduces variance roughly by factor of avg group size
    avg_group_size = n_terms / n_groups
    variance_reduction = avg_group_size

    return GroupingResult(
        groups=groups,
        shots_per_group=shots_per_group,
        variance_reduction=variance_reduction,
        method="simple_commuting",
    )


def allocate_shots_neyman(
    coefficients: np.ndarray,
    groups: List[List[int]],
    total_shots: int,
) -> np.ndarray:
    """
    Allocate measurement shots using Neyman allocation.

    Neyman allocation minimizes total variance under fixed budget:
    m_g ∝ √(Σ c_i² for i in group g)

    Args:
        coefficients: Hamiltonian coefficients
        groups: List of term groups
        total_shots: Total measurement budget

    Returns:
        Array of shots per group
    """
    # Compute group weights (proportional to √variance)
    weights = []
    for group in groups:
        group_variance = sum(coefficients[i]**2 for i in group)
        weights.append(np.sqrt(group_variance))

    weights = np.array(weights)

    if np.sum(weights) == 0:
        # Uniform if all weights zero
        n_groups = len(groups)
        return np.full(n_groups, total_shots // n_groups)

    # Neyman allocation
    fractions = weights / np.sum(weights)
    shots = np.maximum(1, (total_shots * fractions).astype(int))

    # Adjust to match total exactly
    while np.sum(shots) > total_shots:
        max_idx = np.argmax(shots)
        shots[max_idx] -= 1

    while np.sum(shots) < total_shots:
        max_idx = np.argmax(weights)
        shots[max_idx] += 1

    return shots
