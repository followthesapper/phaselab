"""
PhaseLab Multi-Evidence Scorer: Unified CRISPRa guide ranking.

Combines all three breakthrough paths:
- Path A: Binding Energy Landscape (quantum chemistry)
- Path B: Transcriptional Phase Alignment (IR dynamics)
- Path C: Off-Target Landscape Geometry (coherence contrast)

Evidence Fusion Strategy:
Each path provides a score on [0, 1] scale. The fusion strategy uses
weighted geometric mean to ensure all paths must be favorable:

    S_combined = (S_A^w_A * S_B^w_B * S_C^w_C)^(1/(w_A + w_B + w_C))

This means a guide scoring 0 on ANY path will have S_combined = 0,
preventing "one good score masks bad ones" failure mode.

Default weights (tunable):
- Path A (binding): 0.3 (quantum chemistry is resource-intensive)
- Path B (phase): 0.4 (IR dynamics is fast and validated)
- Path C (geometry): 0.3 (depends on off-target data quality)

Author: PhaseLab
Date: December 2025
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any, Tuple
import logging

from .binding_landscape import compute_binding_energy, BindingEnergyResult
from .transcriptional_phase import compute_phase_alignment, PhaseAlignmentResult
from .offtarget_geometry import compute_offtarget_geometry, OffTargetGeometryResult

logger = logging.getLogger(__name__)

# Default evidence weights
DEFAULT_WEIGHTS = {
    'binding': 0.3,    # Path A
    'phase': 0.4,      # Path B
    'geometry': 0.3,   # Path C
}


@dataclass
class MultiEvidenceResult:
    """
    Result from multi-evidence scoring.

    Attributes:
        combined_score: Fused score from all paths (0-1)
        binding_result: Result from Path A
        phase_result: Result from Path B
        geometry_result: Result from Path C
        individual_scores: Normalized scores per path
        evidence_levels: Evidence level per path
        is_go: Whether combined analysis passes GO threshold
        recommendation: Human-readable recommendation
    """
    combined_score: float
    binding_result: Optional[BindingEnergyResult]
    phase_result: Optional[PhaseAlignmentResult]
    geometry_result: Optional[OffTargetGeometryResult]
    individual_scores: Dict[str, float]
    evidence_levels: Dict[str, str]
    is_go: bool
    recommendation: str
    details: Dict[str, Any] = field(default_factory=dict)

    def __repr__(self):
        status = "GO" if self.is_go else "NO-GO"
        return (
            f"MultiEvidenceResult("
            f"combined={self.combined_score:.3f} [{status}], "
            f"binding={self.individual_scores.get('binding', 0):.2f}, "
            f"phase={self.individual_scores.get('phase', 0):.2f}, "
            f"geometry={self.individual_scores.get('geometry', 0):.2f})"
        )


def _normalize_binding_score(result: BindingEnergyResult) -> float:
    """
    Normalize binding energy to [0, 1] score.

    More negative ΔE = better binding = higher score.
    """
    # Typical ΔE range: -0.1 (excellent) to +0.1 (poor)
    # Map to 0-1 where lower (more negative) is better
    delta_E = result.delta_E

    # Sigmoid-like normalization
    # ΔE = -0.1 → score ≈ 0.9
    # ΔE = 0 → score = 0.5
    # ΔE = +0.1 → score ≈ 0.1
    score = 1 / (1 + np.exp(20 * delta_E))

    # Coherence bonus
    if result.coherence > 0.135:  # GO threshold
        score = min(1.0, score * 1.1)

    return float(score)


def _normalize_phase_score(result: PhaseAlignmentResult) -> float:
    """
    Normalize phase alignment to [0, 1] score.

    Higher coherence and enhancement = better.
    """
    # Phase coherence already in [0, 1]
    base_score = result.phase_coherence

    # Enhancement factor bonus
    # Enhancement > 3× is good
    enhancement_bonus = min(0.2, (result.enhancement_factor - 1) / 10)

    # Critical window bonus
    window_bonus = 0.1 if result.critical_window else 0.0

    score = base_score + enhancement_bonus + window_bonus
    return float(min(1.0, max(0.0, score)))


def _normalize_geometry_score(result: OffTargetGeometryResult) -> float:
    """
    Normalize off-target geometry to [0, 1] score.

    Higher ΔR̄ (coherence contrast) = better specificity = higher score.
    """
    # ΔR̄ typically in [-0.5, 0.8]
    # Map to 0-1 where higher is better
    delta_R = result.delta_R

    # Linear mapping with clipping
    # ΔR̄ = 0.5 → score = 1.0
    # ΔR̄ = 0 → score = 0.5
    # ΔR̄ = -0.5 → score = 0.0
    score = (delta_R + 0.5) / 1.0
    score = max(0.0, min(1.0, score))

    # On-target coherence bonus
    if result.R_bar_on > 0.5:
        score = min(1.0, score * 1.1)

    return float(score)


def _fuse_scores(
    scores: Dict[str, float],
    weights: Dict[str, float],
) -> float:
    """
    Fuse individual scores using weighted geometric mean.

    Geometric mean ensures all scores must be good - one zero kills everything.
    """
    # Filter to available scores
    available = {k: v for k, v in scores.items() if k in weights and v is not None}

    if not available:
        return 0.0

    # Weighted geometric mean
    total_weight = sum(weights[k] for k in available)
    if total_weight == 0:
        return 0.0

    log_sum = sum(
        weights[k] * np.log(max(v, 1e-10))  # Avoid log(0)
        for k, v in available.items()
    )

    combined = np.exp(log_sum / total_weight)
    return float(combined)


def compute_multi_evidence_score(
    guide_sequence: str,
    promoter_sequence: Optional[str] = None,
    tss_position: Optional[int] = None,
    guide_position: Optional[int] = None,
    offtargets: Optional[List[Dict]] = None,
    weights: Optional[Dict[str, float]] = None,
    run_binding: bool = True,
    run_phase: bool = True,
    run_geometry: bool = True,
    use_quantum: bool = True,
) -> MultiEvidenceResult:
    """
    Compute multi-evidence score for a CRISPRa guide.

    This is the main entry point for unified guide evaluation.

    Args:
        guide_sequence: 20bp guide RNA sequence
        promoter_sequence: Full promoter DNA (required for Path B)
        tss_position: TSS position (required for Path B)
        guide_position: Guide binding position (required for Path B)
        offtargets: Off-target list for Path C (optional)
        weights: Custom weights for evidence fusion
        run_binding: Run Path A (quantum binding)
        run_phase: Run Path B (phase alignment)
        run_geometry: Run Path C (geometry)
        use_quantum: Use ATLAS-Q where available

    Returns:
        MultiEvidenceResult with fused score and individual results

    Example:
        >>> result = compute_multi_evidence_score(
        ...     guide_sequence="ATCGATCGATCGATCGATCG",
        ...     promoter_sequence=promoter,
        ...     tss_position=500,
        ...     guide_position=200,
        ... )
        >>> print(f"Combined: {result.combined_score:.3f}")
        >>> if result.is_go:
        ...     print(f"RECOMMENDATION: {result.recommendation}")
    """
    if weights is None:
        weights = DEFAULT_WEIGHTS.copy()

    individual_scores = {}
    evidence_levels = {}

    binding_result = None
    phase_result = None
    geometry_result = None

    # Path A: Binding Energy Landscape
    if run_binding:
        try:
            binding_result = compute_binding_energy(
                guide_sequence,
                use_quantum=use_quantum,
            )
            individual_scores['binding'] = _normalize_binding_score(binding_result)
            evidence_levels['binding'] = binding_result.evidence
        except Exception as e:
            logger.warning(f"Path A (binding) failed: {e}")
            individual_scores['binding'] = 0.5  # Neutral fallback
            evidence_levels['binding'] = "FAILED"

    # Path B: Transcriptional Phase Alignment
    if run_phase and promoter_sequence and tss_position is not None and guide_position is not None:
        try:
            phase_result = compute_phase_alignment(
                guide_sequence=guide_sequence,
                promoter_sequence=promoter_sequence,
                tss_position=tss_position,
                guide_position=guide_position,
            )
            individual_scores['phase'] = _normalize_phase_score(phase_result)
            evidence_levels['phase'] = phase_result.evidence
        except Exception as e:
            logger.warning(f"Path B (phase) failed: {e}")
            individual_scores['phase'] = 0.5
            evidence_levels['phase'] = "FAILED"
    elif run_phase:
        logger.info("Path B skipped: promoter_sequence, tss_position, or guide_position missing")
        evidence_levels['phase'] = "SKIPPED"

    # Path C: Off-Target Landscape Geometry
    if run_geometry:
        try:
            geometry_result = compute_offtarget_geometry(
                guide_sequence,
                offtargets=offtargets,
                use_quantum=use_quantum,
            )
            individual_scores['geometry'] = _normalize_geometry_score(geometry_result)
            evidence_levels['geometry'] = geometry_result.evidence
        except Exception as e:
            logger.warning(f"Path C (geometry) failed: {e}")
            individual_scores['geometry'] = 0.5
            evidence_levels['geometry'] = "FAILED"

    # Fuse scores
    combined_score = _fuse_scores(individual_scores, weights)

    # Determine GO/NO-GO
    # GO requires: combined > 0.5 AND all individual scores > 0.3
    min_individual = min(individual_scores.values()) if individual_scores else 0
    is_go = combined_score > 0.5 and min_individual > 0.3

    # Generate recommendation
    if is_go:
        if combined_score > 0.8:
            recommendation = "STRONG CANDIDATE - Proceed with validation"
        elif combined_score > 0.6:
            recommendation = "GOOD CANDIDATE - Consider for wet lab testing"
        else:
            recommendation = "MARGINAL - May work but not optimal"
    else:
        weak_paths = [k for k, v in individual_scores.items() if v < 0.4]
        if weak_paths:
            recommendation = f"NOT RECOMMENDED - Weak on: {', '.join(weak_paths)}"
        else:
            recommendation = "NOT RECOMMENDED - Combined evidence insufficient"

    return MultiEvidenceResult(
        combined_score=combined_score,
        binding_result=binding_result,
        phase_result=phase_result,
        geometry_result=geometry_result,
        individual_scores=individual_scores,
        evidence_levels=evidence_levels,
        is_go=is_go,
        recommendation=recommendation,
        details={
            'weights': weights,
            'guide_sequence': guide_sequence,
        }
    )


def rank_guides_multi_evidence(
    guides: List[Dict[str, Any]],
    promoter_sequence: Optional[str] = None,
    tss_position: Optional[int] = None,
    offtargets_per_guide: Optional[Dict[str, List[Dict]]] = None,
    weights: Optional[Dict[str, float]] = None,
    use_quantum: bool = True,
) -> List[Dict[str, Any]]:
    """
    Rank multiple guides using multi-evidence scoring.

    Args:
        guides: List of guide dicts with 'sequence' and optionally 'position'
        promoter_sequence: Full promoter DNA
        tss_position: TSS position
        offtargets_per_guide: Dict mapping guide sequence to off-targets
        weights: Custom weights
        use_quantum: Use ATLAS-Q

    Returns:
        Guides sorted by combined score (highest first)

    Example:
        >>> from phaselab.crispr import design_crispra_guides
        >>> result = design_crispra_guides(...)
        >>> ranked = rank_guides_multi_evidence(
        ...     result.ranked_guides,
        ...     promoter_sequence=promoter,
        ...     tss_position=500,
        ... )
        >>> for g in ranked[:5]:
        ...     print(f"{g['sequence']}: {g['multi_evidence']['combined_score']:.3f}")
    """
    for guide in guides:
        seq = guide.get('sequence', '')
        if not seq:
            continue

        # Get position (try different keys)
        pos = guide.get('position')
        if pos is None and 'tss_relative_position' in guide and tss_position is not None:
            pos = tss_position + guide['tss_relative_position']

        # Get off-targets for this guide
        offtargets = None
        if offtargets_per_guide and seq in offtargets_per_guide:
            offtargets = offtargets_per_guide[seq]

        try:
            result = compute_multi_evidence_score(
                guide_sequence=seq,
                promoter_sequence=promoter_sequence,
                tss_position=tss_position,
                guide_position=pos,
                offtargets=offtargets,
                weights=weights,
                use_quantum=use_quantum,
            )

            guide['multi_evidence'] = {
                'combined_score': result.combined_score,
                'individual_scores': result.individual_scores,
                'evidence_levels': result.evidence_levels,
                'is_go': result.is_go,
                'recommendation': result.recommendation,
            }

        except Exception as e:
            logger.warning(f"Multi-evidence scoring failed for {seq[:10]}...: {e}")
            guide['multi_evidence'] = {
                'combined_score': 0.0,
                'individual_scores': {},
                'evidence_levels': {'error': str(e)},
                'is_go': False,
                'recommendation': f"FAILED: {e}",
            }

    # Sort by combined score (descending)
    guides.sort(
        key=lambda g: g.get('multi_evidence', {}).get('combined_score', 0),
        reverse=True
    )

    return guides


def print_multi_evidence_report(
    guides: List[Dict[str, Any]],
    max_guides: int = 15,
) -> None:
    """
    Print formatted multi-evidence ranking report.

    Args:
        guides: Guides with 'multi_evidence' field
        max_guides: Maximum guides to display
    """
    print("\n" + "=" * 100)
    print("MULTI-EVIDENCE GUIDE RANKING (v0.9.4)")
    print("=" * 100)
    print(f"{'Rank':<5} {'Sequence':<22} {'Combined':<10} {'Binding':<10} "
          f"{'Phase':<10} {'Geometry':<10} {'GO?':<5} {'Recommendation'}")
    print("-" * 100)

    for i, guide in enumerate(guides[:max_guides], 1):
        me = guide.get('multi_evidence', {})
        scores = me.get('individual_scores', {})

        seq = guide.get('sequence', 'N/A')[:20]
        combined = me.get('combined_score', 0)
        binding = scores.get('binding', 0)
        phase = scores.get('phase', 0)
        geometry = scores.get('geometry', 0)
        is_go = "GO" if me.get('is_go') else "NO"
        rec = me.get('recommendation', '')[:30]

        print(f"{i:<5} {seq:<22} {combined:<10.3f} {binding:<10.2f} "
              f"{phase:<10.2f} {geometry:<10.2f} {is_go:<5} {rec}")

    print("-" * 100)
    print("Scoring: Binding (Path A), Phase (Path B), Geometry (Path C)")
    print("Combined = Weighted geometric mean of individual scores")
    print("GO requires combined > 0.5 AND all individual > 0.3")
    print()


__all__ = [
    'MultiEvidenceResult',
    'compute_multi_evidence_score',
    'rank_guides_multi_evidence',
    'print_multi_evidence_report',
    'DEFAULT_WEIGHTS',
]
