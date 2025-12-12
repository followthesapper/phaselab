"""
CRISPOR + PhaseLab Integrated Pipeline.

Combines IR phase coherence scoring with CRISPOR off-target analysis
for multi-objective guide design.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field, asdict
from enum import Enum
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple

import numpy as np

from .client import CrisporClient, CrisporConfig, CrisporOutput
from .parsers import CrisporGuideRow, parse_guides_tsv, index_guides_by_sequence

logger = logging.getLogger(__name__)


class EvidenceLevel(Enum):
    """
    Evidence level for guide validation.

    A = CRISPOR + experimental off-target + wet-lab outcome (gold standard)
    B = CRISPOR metrics + IR add-ons (entropy, ΔR̄) (validated)
    C = Missing CRISPOR, or relying only on heuristic coherence (suggestions only)
    """
    A = "A"  # Fully validated (CRISPOR + experimental)
    B = "B"  # CRISPOR validated + IR analysis
    C = "C"  # Unvalidated / suggestions only

    def __str__(self) -> str:
        descriptions = {
            "A": "VALIDATED (CRISPOR + experimental)",
            "B": "CRISPOR + IR analysis",
            "C": "UNVALIDATED (suggestions only)"
        }
        return descriptions.get(self.value, self.value)


class CoherenceSource(Enum):
    """
    Source of coherence calculation.

    Makes it clear where R̄ came from and how much to trust it.
    """
    HEURISTIC = "heuristic"  # Hamiltonian coefficient variance proxy
    ATLAS_Q_EXPECTATIONS = "atlas_q_expectations"  # ATLAS-Q from expectation values
    ATLAS_Q_VQE = "atlas_q_vqe"  # ATLAS-Q from VQE simulation
    HARDWARE = "hardware"  # Real quantum hardware measurement
    UNKNOWN = "unknown"

    @property
    def is_quantum_measured(self) -> bool:
        """True if coherence came from actual quantum measurement/simulation."""
        return self in (CoherenceSource.ATLAS_Q_VQE, CoherenceSource.HARDWARE)


@dataclass
class GuideCandidate:
    """
    Combined guide candidate with both PhaseLab IR and CRISPOR metrics.

    This unifies quantum phase coherence with classical off-target analysis.
    Now includes IR-enhanced off-target metrics (v0.6.0+).

    IMPORTANT: Check evidence_level before trusting scores:
    - Level A: Fully validated (can trust score)
    - Level B: CRISPOR validated + IR (can compare within category)
    - Level C: Unvalidated suggestions only (do NOT compare to A/B)
    """
    # Identification
    sequence: str  # 20bp guide
    pam: str  # NGG typically
    position: int  # Relative to TSS
    strand: str  # + or -

    # IR Phase Coherence (from PhaseLab ATLAS-Q)
    r_bar: float = 0.0  # Phase coherence R̄ = |⟨e^(iφ)⟩|
    r_bar_std: float = 0.0  # Uncertainty in R̄
    r_bar_zscore: float = 0.0  # Z-score within locus (relative ranking)
    phase_mean: float = 0.0  # Mean phase angle
    phase_concentration: float = 0.0  # von Mises κ
    coherence_source: CoherenceSource = CoherenceSource.UNKNOWN  # Where R̄ came from

    # CRISPOR Specificity
    mit_specificity: Optional[float] = None  # 0-100, higher = fewer off-targets
    cfd_specificity: Optional[float] = None

    # CRISPOR Activity
    doench_2016: Optional[float] = None  # On-target activity score
    moreno_mateos: Optional[float] = None  # T7 activity
    out_of_frame: Optional[float] = None

    # Off-target counts
    ot_0mm: int = 0
    ot_1mm: int = 0
    ot_2mm: int = 0
    ot_3mm: int = 0
    ot_4mm: int = 0

    # IR-Enhanced Off-Target Metrics (NEW in v0.6.0)
    ot_entropy: float = 0.0  # Shannon entropy of risk distribution (higher = safer)
    ot_entropy_normalized: float = 0.0  # 0-1 scale
    delta_r_bar: float = 0.0  # Coherence contrast: R̄_on - max(R̄_off)
    ot_r_bar_max: float = 0.0  # Highest R̄ among off-targets
    ot_energy_tail: float = 0.0  # 95th percentile binding energy
    risk_concentrated: bool = False  # Low entropy warning
    ot_competitive: bool = False  # Off-targets too similar to on-target
    has_exonic_close_ot: bool = False  # Critical safety flag

    # Risk Mass Metrics (NEW in v0.6.1)
    risk_mass_close: float = 0.0  # Sum(CFD) for 0-2mm off-targets
    risk_mass_distant: float = 0.0  # Sum(CFD) for 3-4mm off-targets
    risk_mass_exonic: float = 0.0  # Sum(CFD) for exonic off-targets
    tail_risk_score: float = 0.0  # max(CFD) in exonic regions
    gini_coefficient: float = 0.0  # Concentration measure (0=equal, 1=concentrated)
    herfindahl_index: float = 0.0  # HHI concentration measure

    # Computed metrics
    final_score: float = 0.0
    stage1_pass: bool = True  # Passed safety gate?
    rank: int = 0

    # Evidence & Metadata
    evidence_level: EvidenceLevel = EvidenceLevel.C  # Default to unvalidated
    crispor_validated: bool = False
    ir_analysis_done: bool = False  # Whether IR-enhanced analysis was performed
    has_experimental_data: bool = False  # Has wet-lab validation
    warnings: List[str] = field(default_factory=list)
    score_breakdown: Dict[str, float] = field(default_factory=dict)  # Score decomposition

    @property
    def total_close_offtargets(self) -> int:
        """Off-targets with 0-2 mismatches (most dangerous)."""
        return self.ot_0mm + self.ot_1mm + self.ot_2mm

    @property
    def total_offtargets(self) -> int:
        """All off-targets."""
        return self.ot_0mm + self.ot_1mm + self.ot_2mm + self.ot_3mm + self.ot_4mm

    @property
    def full_sequence(self) -> str:
        """Guide + PAM."""
        return self.sequence + self.pam

    def compute_evidence_level(self) -> EvidenceLevel:
        """Compute and return the evidence level based on available data."""
        if self.has_experimental_data and self.crispor_validated:
            return EvidenceLevel.A
        elif self.crispor_validated:
            return EvidenceLevel.B
        else:
            return EvidenceLevel.C

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        d = asdict(self)
        # Handle numpy types and enums
        for k, v in d.items():
            if isinstance(v, (np.floating, np.integer)):
                d[k] = float(v) if isinstance(v, np.floating) else int(v)
            elif isinstance(v, np.bool_):
                d[k] = bool(v)
            elif isinstance(v, Enum):
                d[k] = v.value
        return d


@dataclass
class ScoringWeights:
    """
    Weights for multi-objective guide scoring.

    Default weights balance safety (specificity) with efficacy (activity, coherence).
    Now includes IR-enhanced off-target metrics (v0.6.0+).

    TWO-STAGE SCORING MODEL:
    - Stage 1: Safety gate (hard filters) - rejects unsafe guides
    - Stage 2: Ranking (soft scores) - ranks safe guides

    This prevents unsafe guides from ever appearing in top rankings.

    COHERENCE MODE (v0.6.1):
    - Quantum mode (VQE): w_coherence = 0.30 (full weight)
    - Heuristic mode: w_coherence = 0.05 (tie-breaker only)

    Heuristic coherence is demoted because it's a structural proxy,
    not true quantum phase coherence. Use it to break ties among
    guides with similar CRISPOR scores, not as primary ranking signal.
    """
    # Phase coherence weight
    # NOTE: This is adjusted based on coherence_source in compute_final_score
    # - Quantum coherence: full weight (0.30)
    # - Heuristic coherence: demoted to tie-breaker (0.05)
    w_coherence: float = 0.30
    w_coherence_heuristic: float = 0.05  # Reduced weight for heuristic mode

    # Specificity weights (safety)
    w_mit: float = 0.25
    w_cfd: float = 0.10

    # Activity weights (efficacy)
    w_doench: float = 0.20

    # Standard off-target penalties
    penalty_per_close_ot: float = 0.02  # Per off-target with 0-2mm
    penalty_per_distant_ot: float = 0.005  # Per off-target with 3-4mm
    penalty_perfect_match_ot: float = 0.50  # Critical: perfect off-target match

    # IR-Enhanced Off-Target Penalties/Bonuses (NEW in v0.6.0)
    bonus_diffuse_risk: float = 0.05  # High entropy = diffuse risk = safer
    penalty_concentrated_risk: float = 0.10  # Low entropy = few bad sites
    penalty_competitive_ot: float = 0.08  # Off-targets with similar R̄
    penalty_exonic_close_ot: float = 0.15  # Critical: exonic off-targets
    penalty_high_affinity_ot: float = 0.05  # Off-targets with favorable binding

    # Risk mass penalties (NEW in v0.6.1)
    penalty_per_risk_mass_close: float = 0.01  # Per unit of sum(CFD) close
    penalty_per_risk_mass_exonic: float = 0.02  # Per unit of sum(CFD) exonic
    penalty_tail_risk: float = 0.10  # If max exonic CFD > 0.5

    # --- STAGE 1: SAFETY GATE (hard filters) ---
    # If ANY of these fail, guide is rejected (stage1_pass = False)
    min_mit_score: float = 30.0  # HARD GATE: Reject below this
    max_perfect_offtargets: int = 1  # HARD GATE: Only on-target allowed
    max_exonic_close_ot: int = 5  # HARD GATE: Exonic 0-2mm off-targets
    max_tail_risk: float = 0.7  # HARD GATE: max CFD in exon

    # --- STAGE 2: SOFT THRESHOLDS (for scoring) ---
    min_coherence: float = 0.80  # Soft: Minimum R̄ threshold
    max_close_offtargets: int = 50  # Soft: Flag if exceeded
    entropy_low_threshold: float = 0.3  # Below = concentrated risk
    entropy_high_threshold: float = 0.7  # Above = diffuse risk (safe)
    delta_r_bar_threshold: float = 0.05  # Below = competitive off-targets

    # --- SCORE CAPPING ---
    # Prevents unvalidated guides from competing with validated ones
    max_score_unvalidated: float = 0.20  # Cap for Level C guides


def check_safety_gate(
    guide: GuideCandidate,
    weights: ScoringWeights,
) -> Tuple[bool, List[str]]:
    """
    Stage 1: Safety gate - hard filters that reject unsafe guides.

    Returns:
        (passed, rejection_reasons)
    """
    rejections = []

    # HARD GATE 1: Minimum MIT specificity
    if guide.mit_specificity is not None and guide.mit_specificity < weights.min_mit_score:
        rejections.append(f"REJECTED: MIT={guide.mit_specificity:.0f} < {weights.min_mit_score}")

    # HARD GATE 2: Perfect off-targets (0mm beyond on-target)
    if guide.ot_0mm > weights.max_perfect_offtargets:
        rejections.append(f"REJECTED: {guide.ot_0mm - 1} perfect off-target(s)")

    # HARD GATE 3: Exonic close off-targets (if known)
    if guide.ir_analysis_done and guide.has_exonic_close_ot:
        # Count exonic close OTs if we have detailed data
        n_exonic_close = getattr(guide, '_n_exonic_close_ot', 0)
        if n_exonic_close > weights.max_exonic_close_ot:
            rejections.append(f"REJECTED: {n_exonic_close} exonic close off-targets")

    # HARD GATE 4: Tail risk (max CFD in exon)
    if guide.tail_risk_score > weights.max_tail_risk:
        rejections.append(f"REJECTED: tail risk {guide.tail_risk_score:.2f} > {weights.max_tail_risk}")

    passed = len(rejections) == 0
    return passed, rejections


def compute_final_score(
    guide: GuideCandidate,
    weights: ScoringWeights = None,
    coherence_contrast: float = 0.0,  # ΔR̄ = R̄_on - max(R̄_off)
) -> Tuple[float, List[str]]:
    """
    Two-stage scoring for guide candidates.

    STAGE 1: Safety gate (hard filters)
    - If failed, guide.stage1_pass = False and score is capped at 0.0

    STAGE 2: Ranking (soft scores)
    - Combines IR coherence, CRISPOR metrics, off-target penalties
    - Only applied to guides that pass Stage 1

    Returns:
        Tuple of (score, warnings)
    """
    if weights is None:
        weights = ScoringWeights()

    warnings = []
    breakdown = {}
    score = 0.0

    # =========================================================================
    # STAGE 1: SAFETY GATE (Hard Filters)
    # =========================================================================
    stage1_pass, rejections = check_safety_gate(guide, weights)
    guide.stage1_pass = stage1_pass

    if not stage1_pass:
        warnings.extend(rejections)
        guide.score_breakdown = {"stage1": "FAILED", "rejections": rejections}
        return 0.0, warnings

    # =========================================================================
    # STAGE 2: SOFT SCORING (Ranking)
    # =========================================================================

    # 1. Phase Coherence Component (0-1 scale)
    # IMPORTANT: Weight is adjusted based on coherence source (v0.6.1)
    # - Quantum coherence (VQE, hardware): full weight
    # - Heuristic coherence: demoted to tie-breaker weight
    if guide.coherence_source == CoherenceSource.HEURISTIC:
        effective_w_coherence = weights.w_coherence_heuristic
        breakdown['coherence_mode'] = 'heuristic (tie-breaker only)'
    else:
        effective_w_coherence = weights.w_coherence
        breakdown['coherence_mode'] = f'{guide.coherence_source.value} (full weight)'

    coherence_score = guide.r_bar * effective_w_coherence
    score += coherence_score
    breakdown['coherence'] = coherence_score
    breakdown['coherence_weight'] = effective_w_coherence

    # Bonus for coherence contrast (localized phase coherence)
    # Only apply bonus if we have quantum-grade coherence
    if coherence_contrast > 0 and guide.coherence_source != CoherenceSource.HEURISTIC:
        delta_bonus = 0.05 * coherence_contrast
        score += delta_bonus
        breakdown['delta_r_bar_bonus'] = delta_bonus

    # 2. MIT Specificity Component (0-100 scale -> 0-1)
    if guide.mit_specificity is not None:
        mit_normalized = guide.mit_specificity / 100.0
        mit_score = mit_normalized * weights.w_mit
        score += mit_score
        breakdown['mit'] = mit_score
    else:
        # No CRISPOR data - can't validate
        warnings.append("No CRISPOR validation available")
        breakdown['mit'] = 0.0

    # 3. CFD Specificity Component
    if guide.cfd_specificity is not None:
        cfd_normalized = guide.cfd_specificity / 100.0
        cfd_score = cfd_normalized * weights.w_cfd
        score += cfd_score
        breakdown['cfd'] = cfd_score

    # 4. Doench Activity Component (0-100 scale -> 0-1)
    if guide.doench_2016 is not None:
        doench_normalized = guide.doench_2016 / 100.0
        doench_score = doench_normalized * weights.w_doench
        score += doench_score
        breakdown['doench'] = doench_score

    # 5. Off-target Penalties (count-based)
    penalties = 0.0

    # Perfect match off-targets are critical
    if guide.ot_0mm > 1:  # More than the on-target site
        penalty = weights.penalty_perfect_match_ot * (guide.ot_0mm - 1)
        penalties += penalty
        breakdown['penalty_perfect_ot'] = -penalty
        warnings.append(f"CRITICAL: {guide.ot_0mm - 1} perfect off-target matches")

    # Close off-targets (1-2mm)
    close_ots = guide.ot_1mm + guide.ot_2mm
    close_penalty = close_ots * weights.penalty_per_close_ot
    penalties += close_penalty
    breakdown['penalty_close_ot'] = -close_penalty

    if close_ots > weights.max_close_offtargets:
        warnings.append(f"High close off-targets: {close_ots}")

    # Distant off-targets (3-4mm) - less penalty
    distant_ots = guide.ot_3mm + guide.ot_4mm
    distant_penalty = distant_ots * weights.penalty_per_distant_ot
    penalties += distant_penalty
    breakdown['penalty_distant_ot'] = -distant_penalty

    score -= penalties

    # 6. Risk Mass Penalties (sum(CFD) based - NEW in v0.6.1)
    if guide.risk_mass_close > 0:
        risk_mass_penalty = guide.risk_mass_close * weights.penalty_per_risk_mass_close
        score -= risk_mass_penalty
        breakdown['penalty_risk_mass_close'] = -risk_mass_penalty

    if guide.risk_mass_exonic > 0:
        risk_exonic_penalty = guide.risk_mass_exonic * weights.penalty_per_risk_mass_exonic
        score -= risk_exonic_penalty
        breakdown['penalty_risk_mass_exonic'] = -risk_exonic_penalty

    if guide.tail_risk_score > 0.5:
        score -= weights.penalty_tail_risk
        breakdown['penalty_tail_risk'] = -weights.penalty_tail_risk
        warnings.append(f"High tail risk: max exonic CFD = {guide.tail_risk_score:.2f}")

    # 7. IR-Enhanced Off-Target Analysis (v0.6.0)
    if guide.ir_analysis_done:
        # Entropy-based risk distribution
        if guide.ot_entropy_normalized > weights.entropy_high_threshold:
            # Diffuse risk = bonus (safer)
            score += weights.bonus_diffuse_risk
            breakdown['bonus_diffuse_risk'] = weights.bonus_diffuse_risk
        elif guide.ot_entropy_normalized < weights.entropy_low_threshold:
            # Concentrated risk = penalty (dangerous)
            score -= weights.penalty_concentrated_risk
            breakdown['penalty_concentrated_risk'] = -weights.penalty_concentrated_risk
            warnings.append(f"Concentrated off-target risk (entropy={guide.ot_entropy_normalized:.2f})")

        # Coherence contrast
        if guide.ot_competitive:
            score -= weights.penalty_competitive_ot
            breakdown['penalty_competitive_ot'] = -weights.penalty_competitive_ot
            warnings.append(f"Competitive off-targets (ΔR̄={guide.delta_r_bar:.3f})")

        # Exonic close off-targets (critical safety)
        if guide.has_exonic_close_ot:
            score -= weights.penalty_exonic_close_ot
            breakdown['penalty_exonic_close'] = -weights.penalty_exonic_close_ot
            warnings.append("CRITICAL: Exonic off-targets with 0-2mm")

        # High-affinity off-targets (binding energy)
        if guide.ot_energy_tail < -15:
            score -= weights.penalty_high_affinity_ot
            breakdown['penalty_high_affinity'] = -weights.penalty_high_affinity_ot
            warnings.append(f"High-affinity off-targets (E={guide.ot_energy_tail:.1f})")

    # 8. Ensure score is bounded
    score = max(0.0, min(1.0, score))

    # 9. Score capping for unvalidated guides
    # Evidence Level C guides cannot compete with validated ones
    guide.evidence_level = guide.compute_evidence_level()
    if guide.evidence_level == EvidenceLevel.C:
        original_score = score
        score = min(score, weights.max_score_unvalidated)
        if score < original_score:
            warnings.append(f"Score capped at {weights.max_score_unvalidated} (unvalidated)")
            breakdown['cap_applied'] = True

    breakdown['final'] = score
    guide.score_breakdown = breakdown

    return score, warnings


def merge_crispor_data(
    guides: List[GuideCandidate],
    crispor_output: CrisporOutput,
) -> List[GuideCandidate]:
    """
    Merge CRISPOR results into guide candidates.

    Matches guides by sequence and adds CRISPOR metrics.
    """
    if not crispor_output.guides_file or not crispor_output.guides_file.exists():
        logger.warning("No CRISPOR guides file available")
        return guides

    # Parse CRISPOR output
    crispor_guides = parse_guides_tsv(crispor_output.guides_file)
    crispor_index = index_guides_by_sequence(crispor_guides)

    matched = 0
    for guide in guides:
        seq_upper = guide.sequence.upper()

        if seq_upper in crispor_index:
            cg = crispor_index[seq_upper]

            # Transfer CRISPOR metrics
            guide.mit_specificity = cg.mit_specificity
            guide.cfd_specificity = cg.cfd_specificity
            guide.doench_2016 = cg.doench_2016
            guide.moreno_mateos = cg.moreno_mateos
            guide.out_of_frame = cg.out_of_frame

            guide.ot_0mm = cg.ot_0mm
            guide.ot_1mm = cg.ot_1mm
            guide.ot_2mm = cg.ot_2mm
            guide.ot_3mm = cg.ot_3mm
            guide.ot_4mm = cg.ot_4mm

            guide.crispor_validated = True
            matched += 1

    logger.info(f"Matched {matched}/{len(guides)} guides with CRISPOR data")
    return guides


def rank_guides(
    guides: List[GuideCandidate],
    weights: ScoringWeights = None,
) -> List[GuideCandidate]:
    """
    Score and rank all guide candidates.

    Returns guides sorted by final_score descending.
    """
    if weights is None:
        weights = ScoringWeights()

    for guide in guides:
        score, warnings = compute_final_score(guide, weights)
        guide.final_score = score
        guide.warnings = warnings

    # Sort by final score descending
    guides.sort(key=lambda g: g.final_score, reverse=True)

    # Assign ranks
    for i, guide in enumerate(guides, 1):
        guide.rank = i

    return guides


def design_guides_with_crispor(
    sequence: str,
    crispor_client: CrisporClient,
    name: str = "target",
    tss_position: int = 500,  # Position of TSS in sequence
    window_upstream: int = 400,
    window_downstream: int = 100,
    weights: ScoringWeights = None,
    ir_scores: Optional[Dict[str, float]] = None,  # Pre-computed R̄ values
) -> List[GuideCandidate]:
    """
    Full pipeline: find guides, run CRISPOR, score and rank.

    Args:
        sequence: Target DNA sequence
        crispor_client: Configured CRISPOR client
        name: Target name
        tss_position: Position of TSS in the sequence (for relative coords)
        window_upstream: bp upstream of TSS to search
        window_downstream: bp downstream of TSS to search
        weights: Scoring weights
        ir_scores: Dict mapping guide sequence -> R̄ value (from PhaseLab)

    Returns:
        Ranked list of GuideCandidate objects
    """
    if weights is None:
        weights = ScoringWeights()

    if ir_scores is None:
        ir_scores = {}

    # Define search window
    window_start = max(0, tss_position - window_upstream)
    window_end = min(len(sequence), tss_position + window_downstream)

    # Run CRISPOR on the full sequence
    logger.info(f"Running CRISPOR on {len(sequence)}bp sequence...")
    crispor_output = crispor_client.score_sequence(
        sequence,
        name=name,
        include_offtargets=True
    )

    if not crispor_output.success:
        logger.error(f"CRISPOR failed: {crispor_output.error}")
        return []

    # Parse CRISPOR results
    if not crispor_output.guides_file:
        logger.error("No guides file from CRISPOR")
        return []

    crispor_guides = parse_guides_tsv(crispor_output.guides_file)
    logger.info(f"CRISPOR found {len(crispor_guides)} potential guides")

    # Convert to GuideCandidate objects
    candidates = []
    for cg in crispor_guides:
        # Calculate position relative to TSS
        if cg.start is not None:
            rel_pos = cg.start - tss_position
        else:
            rel_pos = 0

        # Skip guides outside our window of interest
        if cg.start is not None:
            if cg.start < window_start or cg.start > window_end:
                continue

        candidate = GuideCandidate(
            sequence=cg.sequence,
            pam=cg.pam,
            position=rel_pos,
            strand=cg.strand or "+",

            # CRISPOR metrics
            mit_specificity=cg.mit_specificity,
            cfd_specificity=cg.cfd_specificity,
            doench_2016=cg.doench_2016,
            moreno_mateos=cg.moreno_mateos,
            out_of_frame=cg.out_of_frame,

            ot_0mm=cg.ot_0mm,
            ot_1mm=cg.ot_1mm,
            ot_2mm=cg.ot_2mm,
            ot_3mm=cg.ot_3mm,
            ot_4mm=cg.ot_4mm,

            crispor_validated=True,
        )

        # Add IR coherence if available
        seq_upper = cg.sequence.upper()
        if seq_upper in ir_scores:
            candidate.r_bar = ir_scores[seq_upper]

        candidates.append(candidate)

    logger.info(f"Found {len(candidates)} guides in target window")

    # Score and rank
    candidates = rank_guides(candidates, weights)

    return candidates


def generate_report(
    candidates: List[GuideCandidate],
    output_path: Path,
    top_n: int = 20,
) -> None:
    """
    Generate a detailed report of top guide candidates.
    """
    lines = []
    lines.append("=" * 80)
    lines.append("INTEGRATED GUIDE DESIGN REPORT")
    lines.append("PhaseLab IR Coherence + CRISPOR Off-Target Analysis")
    lines.append("=" * 80)
    lines.append("")

    lines.append(f"Total candidates analyzed: {len(candidates)}")
    lines.append(f"Showing top {min(top_n, len(candidates))} guides")
    lines.append("")

    # Summary table
    lines.append("-" * 80)
    lines.append(f"{'Rank':<5} {'Sequence':<22} {'Pos':<7} {'Score':<7} {'R̄':<6} {'MIT':<6} {'Doench':<7} {'OT(0-2)':<8}")
    lines.append("-" * 80)

    for c in candidates[:top_n]:
        mit_str = f"{c.mit_specificity:.0f}" if c.mit_specificity else "N/A"
        doench_str = f"{c.doench_2016:.0f}" if c.doench_2016 else "N/A"

        lines.append(
            f"{c.rank:<5} {c.sequence:<22} {c.position:<7} "
            f"{c.final_score:.4f} {c.r_bar:.3f} {mit_str:<6} {doench_str:<7} "
            f"{c.total_close_offtargets:<8}"
        )

        if c.warnings:
            for w in c.warnings:
                lines.append(f"       ⚠ {w}")

    lines.append("-" * 80)
    lines.append("")

    # Detailed analysis of top 5
    lines.append("DETAILED ANALYSIS - TOP 5 CANDIDATES")
    lines.append("=" * 80)

    for c in candidates[:5]:
        lines.append("")
        lines.append(f"Rank #{c.rank}: {c.sequence}")
        lines.append(f"  PAM: {c.pam} | Strand: {c.strand} | Position: {c.position}")
        lines.append(f"  ")
        lines.append(f"  IR Coherence:")
        lines.append(f"    R̄ = {c.r_bar:.4f}")
        lines.append(f"  ")
        lines.append(f"  CRISPOR Scores:")
        lines.append(f"    MIT Specificity: {c.mit_specificity or 'N/A'}")
        lines.append(f"    CFD Specificity: {c.cfd_specificity or 'N/A'}")
        lines.append(f"    Doench 2016: {c.doench_2016 or 'N/A'}")
        lines.append(f"  ")
        lines.append(f"  Off-targets:")
        lines.append(f"    0mm: {c.ot_0mm} | 1mm: {c.ot_1mm} | 2mm: {c.ot_2mm} | 3mm: {c.ot_3mm} | 4mm: {c.ot_4mm}")
        lines.append(f"    Total close (0-2mm): {c.total_close_offtargets}")
        lines.append(f"  ")
        lines.append(f"  Final Score: {c.final_score:.4f}")

        if c.warnings:
            lines.append(f"  ")
            lines.append(f"  Warnings:")
            for w in c.warnings:
                lines.append(f"    ⚠ {w}")

    # Write report
    output_path.write_text("\n".join(lines))
    logger.info(f"Report written to {output_path}")


def save_results_json(
    candidates: List[GuideCandidate],
    output_path: Path,
    metadata: Optional[Dict[str, Any]] = None,
) -> None:
    """Save results to JSON for programmatic access."""
    data = {
        "metadata": metadata or {},
        "total_candidates": len(candidates),
        "candidates": [c.to_dict() for c in candidates],
    }

    with output_path.open("w") as f:
        json.dump(data, f, indent=2)

    logger.info(f"Results saved to {output_path}")
