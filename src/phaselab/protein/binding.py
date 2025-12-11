"""
PhaseLab Protein Binding: IR coherence for docking and affinity prediction.

Applies coherence metrics to:
- Molecular docking pose reliability
- Binding affinity prediction confidence
- Protein-protein interaction assessment
"""

import numpy as np
from typing import List, Optional, Dict, Tuple
from dataclasses import dataclass

from ..core.constants import E_MINUS_2
from ..core.coherence import coherence_score, go_no_go


@dataclass
class BindingCoherence:
    """Container for binding prediction coherence metrics."""

    R_bar: float
    go_no_go: str

    # Component scores
    pose_coherence_R: float
    affinity_coherence_R: float
    ensemble_coherence_R: float

    # Predictions
    predicted_affinity_kcal: Optional[float] = None
    affinity_std: Optional[float] = None
    n_poses: int = 0

    def to_dict(self) -> Dict:
        return {
            "R_bar": self.R_bar,
            "go_no_go": self.go_no_go,
            "pose_coherence_R": self.pose_coherence_R,
            "affinity_coherence_R": self.affinity_coherence_R,
            "ensemble_coherence_R": self.ensemble_coherence_R,
            "predicted_affinity_kcal": self.predicted_affinity_kcal,
            "affinity_std": self.affinity_std,
            "n_poses": self.n_poses,
        }


def docking_coherence(
    poses: List[np.ndarray],
    scores: Optional[List[float]] = None,
    rmsd_threshold: float = 2.0
) -> Tuple[float, int]:
    """
    Compute coherence from docking pose ensemble.

    Multiple docking runs should converge to similar poses if the
    binding mode is reliable. High pose diversity = low coherence.

    Args:
        poses: List of ligand coordinate arrays (n_atoms, 3)
        scores: Optional docking scores for each pose
        rmsd_threshold: RMSD threshold for clustering (Å)

    Returns:
        (R̄ coherence, n_clusters)
    """
    if len(poses) < 2:
        return 1.0, 1

    # Compute pairwise RMSD
    n_poses = len(poses)
    rmsd_matrix = np.zeros((n_poses, n_poses))

    for i in range(n_poses):
        for j in range(i + 1, n_poses):
            # Simple RMSD (assuming aligned)
            rmsd = np.sqrt(np.mean((poses[i] - poses[j]) ** 2))
            rmsd_matrix[i, j] = rmsd
            rmsd_matrix[j, i] = rmsd

    # Cluster by RMSD
    # Simple single-linkage clustering
    clusters = list(range(n_poses))
    for i in range(n_poses):
        for j in range(i + 1, n_poses):
            if rmsd_matrix[i, j] < rmsd_threshold:
                # Merge clusters
                old_cluster = clusters[j]
                new_cluster = clusters[i]
                clusters = [new_cluster if c == old_cluster else c for c in clusters]

    n_clusters = len(set(clusters))

    # Coherence based on clustering
    # All poses in one cluster = high coherence
    R_pose = 1.0 / n_clusters

    # Weight by scores if available
    if scores is not None:
        scores = np.array(scores)
        # Lower scores are better in docking
        # Variance in top scores indicates uncertainty
        top_scores = np.sort(scores)[:max(3, len(scores) // 4)]
        score_variance = np.var(top_scores)

        # Combine pose and score coherence
        R_score = np.exp(-score_variance / 2.0)
        R_pose = (R_pose + R_score) / 2

    return float(R_pose), n_clusters


def affinity_confidence(
    predictions: np.ndarray,
    method: str = "ensemble"
) -> Tuple[float, float, float]:
    """
    Compute confidence in binding affinity predictions.

    Args:
        predictions: Array of predicted affinities (kcal/mol or pKd)
        method: "ensemble" (multiple models) or "bootstrap"

    Returns:
        (mean_affinity, std, R̄ coherence)
    """
    if len(predictions) < 2:
        return float(predictions[0]) if len(predictions) == 1 else 0.0, 0.0, 1.0

    mean_aff = float(np.mean(predictions))
    std_aff = float(np.std(predictions))

    # Coherence from prediction variance
    # Typical affinity uncertainty ~ 1-2 kcal/mol
    # Lower variance = higher coherence
    R_aff = np.exp(-std_aff / 1.5)

    return mean_aff, std_aff, float(np.clip(R_aff, 0.0, 1.0))


def interaction_coherence(
    contacts: np.ndarray,
    reference_contacts: Optional[np.ndarray] = None
) -> float:
    """
    Compute coherence from protein-ligand interactions.

    Args:
        contacts: Binary contact matrix or interaction fingerprint
        reference_contacts: Optional reference for comparison

    Returns:
        R̄ coherence for interaction pattern
    """
    if contacts.size == 0:
        return 0.0

    if reference_contacts is not None:
        # Tanimoto similarity as coherence
        intersection = np.sum(contacts & reference_contacts)
        union = np.sum(contacts | reference_contacts)
        if union == 0:
            return 1.0
        return float(intersection / union)
    else:
        # Self-consistency of contacts
        # Check symmetry and completeness
        contact_density = np.mean(contacts)
        # Medium density is most reliable
        # Very sparse or very dense contacts are uncertain
        R_int = 1.0 - abs(contact_density - 0.3) / 0.3

        return float(np.clip(R_int, 0.0, 1.0))


def compute_binding_coherence(
    poses: Optional[List[np.ndarray]] = None,
    docking_scores: Optional[List[float]] = None,
    affinity_predictions: Optional[np.ndarray] = None,
    contacts: Optional[np.ndarray] = None
) -> BindingCoherence:
    """
    Compute comprehensive binding prediction coherence.

    Args:
        poses: Docking pose coordinates
        docking_scores: Scores for each pose
        affinity_predictions: Ensemble affinity predictions
        contacts: Interaction fingerprint

    Returns:
        BindingCoherence dataclass
    """
    scores = []

    # Pose coherence
    if poses is not None:
        pose_R, n_clusters = docking_coherence(poses, docking_scores)
        scores.append(pose_R)
        n_poses = len(poses)
    else:
        pose_R = 0.0
        n_poses = 0

    # Affinity coherence
    if affinity_predictions is not None:
        mean_aff, std_aff, aff_R = affinity_confidence(affinity_predictions)
        scores.append(aff_R)
    else:
        aff_R = 0.0
        mean_aff = None
        std_aff = None

    # Interaction coherence
    if contacts is not None:
        int_R = interaction_coherence(contacts)
        scores.append(int_R)
    else:
        int_R = 0.0

    # Overall coherence
    if scores:
        R_bar = float(np.mean(scores))
    else:
        R_bar = 0.0

    return BindingCoherence(
        R_bar=R_bar,
        go_no_go=go_no_go(R_bar),
        pose_coherence_R=pose_R,
        affinity_coherence_R=aff_R,
        ensemble_coherence_R=int_R,
        predicted_affinity_kcal=mean_aff,
        affinity_std=std_aff,
        n_poses=n_poses
    )
