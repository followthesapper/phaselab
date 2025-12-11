"""
PhaseLab Protein Folding: IR coherence for structure prediction.

Applies the R̄ = exp(-V_φ/2) framework to assess protein folding simulation
reliability, using structural metrics as phase-like observables.

Key insight: Protein backbone dihedral angles (φ, ψ) are literal phases that
can be analyzed using Kuramoto-style order parameters.
"""

import numpy as np
from typing import List, Tuple, Optional, Dict, Union
from dataclasses import dataclass

from ..core.constants import E_MINUS_2
from ..core.coherence import coherence_score, go_no_go


@dataclass
class FoldingCoherence:
    """Container for protein folding coherence metrics."""

    # Overall coherence
    R_bar: float
    go_no_go: str

    # Component scores
    ramachandran_R: float
    contact_map_R: float
    secondary_structure_R: float

    # Metadata
    n_residues: int
    method: str  # "alphafold", "rosetta", "md", etc.
    confidence: Optional[float] = None  # e.g., pLDDT for AlphaFold

    def to_dict(self) -> Dict:
        return {
            "R_bar": self.R_bar,
            "go_no_go": self.go_no_go,
            "ramachandran_R": self.ramachandran_R,
            "contact_map_R": self.contact_map_R,
            "secondary_structure_R": self.secondary_structure_R,
            "n_residues": self.n_residues,
            "method": self.method,
            "confidence": self.confidence,
        }


def ramachandran_coherence(
    phi_angles: np.ndarray,
    psi_angles: np.ndarray,
    weights: Optional[np.ndarray] = None
) -> float:
    """
    Compute coherence from Ramachandran angles (φ, ψ).

    Uses the backbone dihedral angles as literal phase observables.
    Well-folded proteins cluster in specific Ramachandran regions,
    giving high coherence. Disordered/unfolded regions scatter.

    Args:
        phi_angles: Backbone φ angles in radians (N-Cα bond rotation)
        psi_angles: Backbone ψ angles in radians (Cα-C bond rotation)
        weights: Optional per-residue weights (e.g., from pLDDT)

    Returns:
        R̄ coherence score for Ramachandran distribution
    """
    if len(phi_angles) == 0 or len(psi_angles) == 0:
        return 0.0

    if len(phi_angles) != len(psi_angles):
        raise ValueError("phi and psi arrays must have same length")

    # Compute order parameters for each angle type
    z_phi = np.exp(1j * phi_angles)
    z_psi = np.exp(1j * psi_angles)

    if weights is not None:
        weights = np.array(weights)
        weights = weights / np.sum(weights)
        R_phi = np.abs(np.sum(weights * z_phi))
        R_psi = np.abs(np.sum(weights * z_psi))
    else:
        R_phi = np.abs(np.mean(z_phi))
        R_psi = np.abs(np.mean(z_psi))

    # Combined coherence (geometric mean)
    R_rama = np.sqrt(R_phi * R_psi)

    return float(R_rama)


def contact_map_coherence(
    contact_map: np.ndarray,
    reference_map: Optional[np.ndarray] = None,
    distance_threshold: float = 8.0
) -> float:
    """
    Compute coherence from residue contact maps.

    Contact maps encode the 3D structure as a 2D matrix of distances.
    Coherence measures how consistently contacts are predicted/observed.

    Args:
        contact_map: N×N matrix of Cα-Cα distances or contact probabilities
        reference_map: Optional reference structure for comparison
        distance_threshold: Threshold for defining contacts (Å)

    Returns:
        R̄ coherence score for contact map
    """
    if contact_map.size == 0:
        return 0.0

    # Binarize if distances
    if contact_map.max() > 1.0:
        contacts = (contact_map < distance_threshold).astype(float)
    else:
        contacts = contact_map

    # Mask diagonal (self-contacts)
    n = contacts.shape[0]
    mask = ~np.eye(n, dtype=bool)

    if reference_map is not None:
        # Compare to reference
        if reference_map.max() > 1.0:
            ref_contacts = (reference_map < distance_threshold).astype(float)
        else:
            ref_contacts = reference_map

        # Agreement as coherence
        agreement = np.mean(contacts[mask] == ref_contacts[mask])
        return float(agreement)
    else:
        # Self-consistency: variance in contact probabilities
        contact_values = contacts[mask]

        # High variance = uncertain contacts = low coherence
        # Convert to phase-like metric
        variance = np.var(contact_values)
        R_contact = np.exp(-variance * 2)  # Scale factor for contact variance

        return float(R_contact)


def secondary_structure_coherence(
    ss_sequence: str,
    confidence: Optional[np.ndarray] = None
) -> float:
    """
    Compute coherence from secondary structure assignment.

    Args:
        ss_sequence: String of H (helix), E (sheet), C (coil)
        confidence: Optional per-residue confidence scores

    Returns:
        R̄ coherence for secondary structure
    """
    if not ss_sequence:
        return 0.0

    # Map SS to phases (arbitrary but consistent)
    ss_to_phase = {'H': 0.0, 'E': 2*np.pi/3, 'C': 4*np.pi/3}

    phases = np.array([ss_to_phase.get(s, 4*np.pi/3) for s in ss_sequence])

    if confidence is not None:
        # Weighted coherence
        weights = np.array(confidence)
        z = np.sum(weights * np.exp(1j * phases)) / np.sum(weights)
    else:
        z = np.mean(np.exp(1j * phases))

    return float(np.abs(z))


def compute_structure_coherence(
    phi_angles: Optional[np.ndarray] = None,
    psi_angles: Optional[np.ndarray] = None,
    contact_map: Optional[np.ndarray] = None,
    ss_sequence: Optional[str] = None,
    plddt: Optional[np.ndarray] = None,
    method: str = "unknown"
) -> FoldingCoherence:
    """
    Compute comprehensive structure coherence from available data.

    Combines multiple structural metrics into overall coherence score.

    Args:
        phi_angles: Backbone φ angles (radians)
        psi_angles: Backbone ψ angles (radians)
        contact_map: Distance or contact probability matrix
        ss_sequence: Secondary structure string
        plddt: AlphaFold pLDDT or similar confidence
        method: Prediction method name

    Returns:
        FoldingCoherence dataclass with all metrics
    """
    scores = []

    # Ramachandran coherence
    if phi_angles is not None and psi_angles is not None:
        rama_R = ramachandran_coherence(phi_angles, psi_angles, weights=plddt)
        scores.append(rama_R)
    else:
        rama_R = 0.0

    # Contact map coherence
    if contact_map is not None:
        contact_R = contact_map_coherence(contact_map)
        scores.append(contact_R)
    else:
        contact_R = 0.0

    # Secondary structure coherence
    if ss_sequence is not None:
        ss_R = secondary_structure_coherence(ss_sequence, confidence=plddt)
        scores.append(ss_R)
    else:
        ss_R = 0.0

    # Overall coherence (geometric mean of available scores)
    if scores:
        R_bar = float(np.exp(np.mean(np.log(np.clip(scores, 1e-10, 1.0)))))
    else:
        R_bar = 0.0

    # Determine n_residues
    if phi_angles is not None:
        n_residues = len(phi_angles)
    elif ss_sequence is not None:
        n_residues = len(ss_sequence)
    elif contact_map is not None:
        n_residues = contact_map.shape[0]
    else:
        n_residues = 0

    # Mean pLDDT as confidence
    mean_plddt = float(np.mean(plddt)) if plddt is not None else None

    return FoldingCoherence(
        R_bar=R_bar,
        go_no_go=go_no_go(R_bar),
        ramachandran_R=rama_R,
        contact_map_R=contact_R,
        secondary_structure_R=ss_R,
        n_residues=n_residues,
        method=method,
        confidence=mean_plddt
    )


def go_no_go_structure(
    coherence: Union[float, FoldingCoherence],
    threshold: float = E_MINUS_2,
    min_plddt: float = 70.0
) -> str:
    """
    GO/NO-GO decision for structure prediction.

    Uses both IR coherence and method-specific confidence.

    Args:
        coherence: R̄ value or FoldingCoherence object
        threshold: IR coherence threshold
        min_plddt: Minimum pLDDT for GO (AlphaFold)

    Returns:
        "GO" or "NO-GO"
    """
    if isinstance(coherence, FoldingCoherence):
        R_bar = coherence.R_bar
        plddt = coherence.confidence
    else:
        R_bar = coherence
        plddt = None

    # Check IR coherence
    if R_bar < threshold:
        return "NO-GO"

    # Check method-specific confidence if available
    if plddt is not None and plddt < min_plddt:
        return "NO-GO"

    return "GO"


def alphafold_coherence(
    pdb_path: str,
    confidence_path: Optional[str] = None
) -> FoldingCoherence:
    """
    Compute coherence from AlphaFold prediction.

    Extracts backbone angles and pLDDT from AlphaFold output.

    Args:
        pdb_path: Path to AlphaFold PDB output
        confidence_path: Path to pLDDT JSON (optional)

    Returns:
        FoldingCoherence for AlphaFold prediction

    Note:
        Requires BioPython for PDB parsing.
    """
    try:
        from Bio.PDB import PDBParser
        from Bio.PDB.Polypeptide import PPBuilder
    except ImportError:
        raise ImportError("BioPython required: pip install biopython")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    ppb = PPBuilder()
    phi_list = []
    psi_list = []

    for pp in ppb.build_peptides(structure):
        for phi, psi in pp.get_phi_psi_list():
            if phi is not None:
                phi_list.append(phi)
            if psi is not None:
                psi_list.append(psi)

    phi_angles = np.array(phi_list) if phi_list else None
    psi_angles = np.array(psi_list) if psi_list else None

    # Extract pLDDT from B-factors (AlphaFold convention)
    plddt = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.name == "CA":
                        plddt.append(atom.bfactor)
                        break

    plddt = np.array(plddt) if plddt else None

    return compute_structure_coherence(
        phi_angles=phi_angles,
        psi_angles=psi_angles,
        plddt=plddt,
        method="alphafold"
    )
