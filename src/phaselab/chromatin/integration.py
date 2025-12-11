"""
PhaseLab Chromatin Integration: Combine chromatin with CRISPR scoring.

Provides functions to integrate chromatin accessibility into the
CRISPR guide design pipeline, improving predictions for tissue-specific
efficacy.
"""

import numpy as np
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

from .accessibility import ChromatinModel, ChromatinScore, ChromatinState
from .tissue_models import get_tissue_model
from ..core.constants import E_MINUS_2
from ..core.coherence import coherence_score


@dataclass
class IntegratedGuideScore:
    """Combined CRISPR and chromatin score."""

    # Guide info
    sequence: str
    position: int

    # CRISPR scores
    mit_score: float
    cfd_score: float
    coherence_R: float
    go_no_go: str

    # Chromatin scores
    chromatin_accessibility: float
    chromatin_state: ChromatinState

    # Integrated score
    integrated_score: float
    tissue_adjusted_R: float

    def to_dict(self) -> Dict:
        return {
            "sequence": self.sequence,
            "position": self.position,
            "mit_score": self.mit_score,
            "cfd_score": self.cfd_score,
            "coherence_R": self.coherence_R,
            "go_no_go": self.go_no_go,
            "chromatin_accessibility": self.chromatin_accessibility,
            "chromatin_state": self.chromatin_state.value,
            "integrated_score": self.integrated_score,
            "tissue_adjusted_R": self.tissue_adjusted_R,
        }


def combine_with_crispr(
    guide_sequence: str,
    guide_position: int,
    chrom: str,
    tss_position: int,
    mit_score: float,
    cfd_score: float,
    coherence_R: float,
    tissue: str = "generic",
    chromatin_weight: float = 0.2
) -> IntegratedGuideScore:
    """
    Combine CRISPR scores with chromatin accessibility.

    The integrated score accounts for both guide quality and
    local chromatin state, improving tissue-specific predictions.

    Args:
        guide_sequence: 20bp guide sequence
        guide_position: Position relative to TSS
        chrom: Chromosome
        tss_position: TSS genomic coordinate
        mit_score: MIT specificity score (0-100)
        cfd_score: CFD score (0-100)
        coherence_R: IR coherence score (0-1)
        tissue: Target tissue for chromatin model
        chromatin_weight: Weight for chromatin in integrated score

    Returns:
        IntegratedGuideScore with combined metrics
    """
    # Get chromatin accessibility
    model = get_tissue_model(tissue)
    genomic_pos = tss_position + guide_position
    chromatin = model.score_position(chrom, genomic_pos)

    # Normalize CRISPR scores to 0-1
    mit_norm = mit_score / 100.0
    cfd_norm = cfd_score / 100.0

    # Compute integrated score
    crispr_score = (mit_norm + cfd_norm + coherence_R) / 3.0
    integrated = (
        (1 - chromatin_weight) * crispr_score +
        chromatin_weight * chromatin.accessibility
    )

    # Tissue-adjusted coherence
    # Closed chromatin reduces effective coherence
    if chromatin.state == ChromatinState.CLOSED:
        tissue_adjusted_R = coherence_R * 0.5
    elif chromatin.state == ChromatinState.PARTIALLY_OPEN:
        tissue_adjusted_R = coherence_R * 0.8
    else:
        tissue_adjusted_R = coherence_R

    # GO/NO-GO based on tissue-adjusted coherence
    if tissue_adjusted_R > E_MINUS_2:
        go_status = "GO"
    else:
        go_status = "NO-GO"

    return IntegratedGuideScore(
        sequence=guide_sequence,
        position=guide_position,
        mit_score=mit_score,
        cfd_score=cfd_score,
        coherence_R=coherence_R,
        go_no_go=go_status,
        chromatin_accessibility=chromatin.accessibility,
        chromatin_state=chromatin.state,
        integrated_score=integrated,
        tissue_adjusted_R=tissue_adjusted_R
    )


def accessibility_weighted_coherence(
    coherence_R: float,
    chromatin_score: ChromatinScore,
    method: str = "multiplicative"
) -> float:
    """
    Adjust coherence by chromatin accessibility.

    Args:
        coherence_R: Base IR coherence
        chromatin_score: ChromatinScore for the position
        method: "multiplicative" or "penalty"

    Returns:
        Adjusted coherence score
    """
    accessibility = chromatin_score.accessibility

    if method == "multiplicative":
        # Scale coherence by accessibility
        # Open chromatin = full coherence
        # Closed chromatin = reduced coherence
        adjusted = coherence_R * (0.5 + 0.5 * accessibility)

    elif method == "penalty":
        # Apply penalty for low accessibility
        if accessibility < 0.3:
            penalty = 0.3  # Severe penalty for closed chromatin
        elif accessibility < 0.6:
            penalty = 0.15  # Moderate penalty
        else:
            penalty = 0.0  # No penalty for open chromatin

        adjusted = coherence_R * (1 - penalty)

    else:
        adjusted = coherence_R

    return float(np.clip(adjusted, 0.0, 1.0))


def rank_guides_by_tissue(
    guides: List[Dict],
    chrom: str,
    tss_position: int,
    tissue: str = "brain"
) -> List[IntegratedGuideScore]:
    """
    Re-rank guide candidates by tissue-specific scores.

    Args:
        guides: List of guide dicts with sequence, position, scores
        chrom: Chromosome
        tss_position: TSS genomic coordinate
        tissue: Target tissue

    Returns:
        List of IntegratedGuideScore, sorted by integrated_score
    """
    integrated_guides = []

    for guide in guides:
        integrated = combine_with_crispr(
            guide_sequence=guide["sequence"],
            guide_position=guide["position"],
            chrom=chrom,
            tss_position=tss_position,
            mit_score=guide.get("mit_score", 50),
            cfd_score=guide.get("cfd_score", 50),
            coherence_R=guide.get("coherence_R", 0.5),
            tissue=tissue
        )
        integrated_guides.append(integrated)

    # Sort by integrated score (descending)
    integrated_guides.sort(key=lambda x: x.integrated_score, reverse=True)

    return integrated_guides


def filter_by_chromatin(
    guides: List[Dict],
    chrom: str,
    tss_position: int,
    tissue: str = "brain",
    min_accessibility: float = 0.5
) -> List[Dict]:
    """
    Filter guides by chromatin accessibility.

    Args:
        guides: List of guide dicts
        chrom: Chromosome
        tss_position: TSS coordinate
        tissue: Target tissue
        min_accessibility: Minimum accessibility threshold

    Returns:
        Filtered list of guides in open/accessible chromatin
    """
    model = get_tissue_model(tissue)
    filtered = []

    for guide in guides:
        genomic_pos = tss_position + guide["position"]
        chromatin = model.score_position(chrom, genomic_pos)

        if chromatin.accessibility >= min_accessibility:
            # Add chromatin info to guide
            guide_copy = guide.copy()
            guide_copy["chromatin_accessibility"] = chromatin.accessibility
            guide_copy["chromatin_state"] = chromatin.state.value
            filtered.append(guide_copy)

    return filtered
