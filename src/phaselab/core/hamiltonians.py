"""
PhaseLab Hamiltonians: Pauli Hamiltonian builders for various systems.

Supports:
- DNA/RNA base-pairing Hamiltonians (for gRNA binding)
- Generic spin Hamiltonians
- VQE-style molecular Hamiltonians
"""

import numpy as np
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass


@dataclass
class PauliTerm:
    """A single term in a Pauli Hamiltonian."""
    coefficient: float
    pauli_string: str  # e.g., "ZZII", "XXII"

    def __repr__(self):
        return f"{self.coefficient:+.4f} * {self.pauli_string}"


class PauliHamiltonian:
    """
    A Hamiltonian expressed as a sum of Pauli terms.

    H = Σ_i c_i * P_i where P_i is a Pauli string (e.g., ZZII, XIXI).
    """

    def __init__(self, n_qubits: int):
        self.n_qubits = n_qubits
        self.terms: List[PauliTerm] = []

    def add_term(self, coefficient: float, pauli_string: str):
        """Add a Pauli term to the Hamiltonian."""
        if len(pauli_string) != self.n_qubits:
            raise ValueError(
                f"Pauli string length {len(pauli_string)} != n_qubits {self.n_qubits}"
            )
        self.terms.append(PauliTerm(coefficient, pauli_string))

    def get_terms(self) -> List[Tuple[float, str]]:
        """Return list of (coefficient, pauli_string) tuples."""
        return [(t.coefficient, t.pauli_string) for t in self.terms]

    def __repr__(self):
        lines = [f"PauliHamiltonian (n_qubits={self.n_qubits}):"]
        for term in self.terms[:10]:  # Show first 10
            lines.append(f"  {term}")
        if len(self.terms) > 10:
            lines.append(f"  ... and {len(self.terms) - 10} more terms")
        return "\n".join(lines)


# Watson-Crick base pairing energies (simplified)
WC_PAIRING = {
    ('A', 'T'): -1.0, ('T', 'A'): -1.0,
    ('A', 'U'): -1.0, ('U', 'A'): -1.0,
    ('G', 'C'): -1.5, ('C', 'G'): -1.5,
    ('A', 'A'): +0.5, ('T', 'T'): +0.5, ('U', 'U'): +0.5,
    ('G', 'G'): +0.3, ('C', 'C'): +0.3,
    ('A', 'G'): +0.4, ('G', 'A'): +0.4,
    ('A', 'C'): +0.6, ('C', 'A'): +0.6,
    ('T', 'G'): +0.4, ('G', 'T'): +0.4,
    ('T', 'C'): +0.5, ('C', 'T'): +0.5,
    ('U', 'G'): +0.3, ('G', 'U'): +0.3,  # G-U wobble
    ('U', 'C'): +0.5, ('C', 'U'): +0.5,
}

# SantaLucia nearest-neighbor stacking energies (kcal/mol, simplified)
NN_STACKING = {
    'AA': -1.0, 'TT': -1.0, 'UU': -1.0,
    'AT': -0.9, 'TA': -0.6,
    'CA': -1.9, 'TG': -1.9, 'UG': -1.9,
    'GT': -1.3, 'AC': -1.3,
    'CT': -1.6, 'AG': -1.6,
    'GA': -1.5, 'TC': -1.5, 'UC': -1.5,
    'CG': -3.6,
    'GC': -3.1,
    'GG': -3.1, 'CC': -3.1,
}


def build_pauli_hamiltonian(
    n_qubits: int,
    zz_terms: Optional[List[Tuple[int, int, float]]] = None,
    xx_terms: Optional[List[Tuple[int, int, float]]] = None,
    z_terms: Optional[List[Tuple[int, float]]] = None,
) -> PauliHamiltonian:
    """
    Build a generic Pauli Hamiltonian.

    Args:
        n_qubits: Number of qubits.
        zz_terms: List of (i, j, coeff) for ZZ interactions.
        xx_terms: List of (i, j, coeff) for XX interactions.
        z_terms: List of (i, coeff) for single-qubit Z terms.

    Returns:
        PauliHamiltonian object.
    """
    H = PauliHamiltonian(n_qubits)

    if zz_terms:
        for i, j, coeff in zz_terms:
            pauli = ['I'] * n_qubits
            pauli[i] = 'Z'
            pauli[j] = 'Z'
            H.add_term(coeff, ''.join(pauli))

    if xx_terms:
        for i, j, coeff in xx_terms:
            pauli = ['I'] * n_qubits
            pauli[i] = 'X'
            pauli[j] = 'X'
            H.add_term(coeff, ''.join(pauli))

    if z_terms:
        for i, coeff in z_terms:
            pauli = ['I'] * n_qubits
            pauli[i] = 'Z'
            H.add_term(coeff, ''.join(pauli))

    return H


def build_grna_hamiltonian(
    guide_seq: str,
    target_seq: Optional[str] = None,
    include_stacking: bool = True,
) -> PauliHamiltonian:
    """
    Build a Hamiltonian encoding gRNA-DNA binding interactions.

    The Hamiltonian encodes:
    - Watson-Crick base pairing (ZZ terms)
    - Stacking interactions (XX/YY terms)
    - GC content bonus (Z terms)

    Args:
        guide_seq: 20bp gRNA sequence (RNA, uses U).
        target_seq: Optional DNA target. If None, assumes perfect WC complement.
        include_stacking: Whether to include nearest-neighbor stacking.

    Returns:
        PauliHamiltonian for VQE/simulation.
    """
    guide_seq = guide_seq.upper().replace('T', 'U')
    n = len(guide_seq)

    if target_seq is None:
        # Generate perfect Watson-Crick complement (DNA)
        complement = {'A': 'T', 'U': 'A', 'G': 'C', 'C': 'G'}
        target_seq = ''.join(complement.get(b, 'N') for b in guide_seq)
    else:
        target_seq = target_seq.upper()

    H = PauliHamiltonian(n)

    # Base pairing terms (ZZ)
    for i in range(n):
        g_base = guide_seq[i]
        t_base = target_seq[i]
        pair = (g_base, t_base)
        energy = WC_PAIRING.get(pair, 0.0)
        if abs(energy) > 1e-6:
            pauli = ['I'] * n
            pauli[i] = 'Z'
            # Single Z term for binding energy at position i
            H.add_term(energy, ''.join(pauli))

    # Stacking terms (XX between neighbors)
    if include_stacking:
        for i in range(n - 1):
            dinuc = guide_seq[i:i+2].replace('U', 'T')
            stack_energy = NN_STACKING.get(dinuc, -0.5)
            pauli = ['I'] * n
            pauli[i] = 'X'
            pauli[i + 1] = 'X'
            H.add_term(stack_energy * 0.1, ''.join(pauli))  # Scale down stacking

    # GC content bonus (identity term effect via ZZ on all GC positions)
    gc_count = sum(1 for b in guide_seq if b in 'GC')
    gc_bonus = -0.1 * gc_count  # Favorable for higher GC
    H.add_term(gc_bonus, 'I' * n)

    return H


def hamiltonian_to_qiskit(H: PauliHamiltonian):
    """
    Convert PauliHamiltonian to Qiskit SparsePauliOp.

    Requires qiskit to be installed.

    Returns:
        qiskit.quantum_info.SparsePauliOp
    """
    try:
        from qiskit.quantum_info import SparsePauliOp
    except ImportError:
        raise ImportError("Qiskit is required for this function. Install with: pip install qiskit")

    paulis = []
    coeffs = []
    for term in H.terms:
        paulis.append(term.pauli_string)
        coeffs.append(term.coefficient)

    return SparsePauliOp.from_list(list(zip(paulis, coeffs)))
