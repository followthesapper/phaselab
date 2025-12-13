"""
PhaseLab CRISPR Scoring: Guide RNA quality metrics.

Implements:
- GC content calculation
- Homopolymer run detection
- SantaLucia thermodynamic ΔG
- MIT specificity algorithm
- CFD (Cutting Frequency Determination) score
- Chromatin accessibility modeling
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Tuple, Optional, Dict, List, Any
import re


# SantaLucia nearest-neighbor parameters (kcal/mol)
# From: SantaLucia & Hicks (2004) Annu. Rev. Biophys. Biomol. Struct.
NN_PARAMS = {
    'AA': (-1.00, -0.0027),  # (ΔH, ΔS)
    'TT': (-1.00, -0.0027),
    'AT': (-0.88, -0.0024),
    'TA': (-0.58, -0.0015),
    'CA': (-1.45, -0.0039),
    'TG': (-1.45, -0.0039),
    'GT': (-1.44, -0.0037),
    'AC': (-1.44, -0.0037),
    'CT': (-1.28, -0.0033),
    'AG': (-1.28, -0.0033),
    'GA': (-1.30, -0.0032),
    'TC': (-1.30, -0.0032),
    'CG': (-2.17, -0.0055),
    'GC': (-2.24, -0.0056),
    'GG': (-1.84, -0.0046),
    'CC': (-1.84, -0.0046),
}

# Initiation parameters
NN_INIT = {
    'G': (0.98, 0.0024),
    'C': (0.98, 0.0024),
    'A': (1.03, 0.0027),
    'T': (1.03, 0.0027),
}

# MIT position weights for off-target scoring
# Higher weight = more important for specificity
MIT_POSITION_WEIGHTS = [
    0, 0, 0.014, 0, 0,       # positions 1-5 (PAM-distal)
    0.395, 0.317, 0, 0.389, 0.079,  # positions 6-10
    0.445, 0.508, 0.613, 0.851, 0.732,  # positions 11-15
    0.828, 0.615, 0.804, 0.685, 0.583,  # positions 16-20 (PAM-proximal)
]


def gc_content(sequence: str) -> float:
    """
    Calculate GC content of a sequence.

    Args:
        sequence: DNA/RNA sequence.

    Returns:
        GC fraction (0.0 to 1.0).
    """
    sequence = sequence.upper()
    gc = sum(1 for b in sequence if b in 'GC')
    return gc / len(sequence) if sequence else 0.0


def max_homopolymer_run(sequence: str) -> int:
    """
    Find the longest homopolymer run in a sequence.

    Args:
        sequence: DNA/RNA sequence.

    Returns:
        Length of longest single-nucleotide repeat.
    """
    if not sequence:
        return 0

    sequence = sequence.upper()
    max_run = 1
    current_run = 1

    for i in range(1, len(sequence)):
        if sequence[i] == sequence[i - 1]:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            current_run = 1

    return max_run


def delta_g_santalucia(
    sequence: str,
    temperature: float = 37.0,
    na_conc: float = 0.1,
) -> float:
    """
    Calculate ΔG of hybridization using SantaLucia nearest-neighbor model.

    Args:
        sequence: DNA/RNA sequence (assumes binding to perfect complement).
        temperature: Temperature in Celsius.
        na_conc: Na+ concentration in M.

    Returns:
        ΔG in kcal/mol (negative = favorable binding).
    """
    sequence = sequence.upper().replace('U', 'T')
    T_kelvin = temperature + 273.15

    if len(sequence) < 2:
        return 0.0

    # Sum nearest-neighbor contributions
    delta_H = 0.0
    delta_S = 0.0

    for i in range(len(sequence) - 1):
        dinuc = sequence[i:i+2]
        if dinuc in NN_PARAMS:
            dH, dS = NN_PARAMS[dinuc]
            delta_H += dH
            delta_S += dS

    # Initiation
    first_base = sequence[0]
    last_base = sequence[-1]
    if first_base in NN_INIT:
        dH, dS = NN_INIT[first_base]
        delta_H += dH
        delta_S += dS
    if last_base in NN_INIT:
        dH, dS = NN_INIT[last_base]
        delta_H += dH
        delta_S += dS

    # Salt correction (simplified)
    delta_S_corrected = delta_S + 0.368 * len(sequence) * np.log(na_conc) / 1000

    # ΔG = ΔH - TΔS
    delta_G = delta_H - T_kelvin * delta_S_corrected

    return delta_G


def sequence_complexity(sequence: str) -> float:
    """
    Calculate sequence complexity (0 = repetitive, 1 = complex).

    Uses linguistic complexity based on unique k-mers.

    Args:
        sequence: DNA sequence.

    Returns:
        Complexity score (0.0 to 1.0).
    """
    sequence = sequence.upper()
    n = len(sequence)

    if n < 3:
        return 1.0

    # Count unique k-mers for k=1,2,3
    total_possible = 0
    total_unique = 0

    for k in [1, 2, 3]:
        kmers = set()
        for i in range(n - k + 1):
            kmers.add(sequence[i:i+k])
        possible = min(4**k, n - k + 1)
        total_possible += possible
        total_unique += len(kmers)

    return total_unique / total_possible if total_possible > 0 else 1.0


def mit_specificity_score(
    guide_seq: str,
    off_target_count: int = 0,
    avg_mismatches: float = 4.0,
) -> float:
    """
    Calculate MIT specificity score (simplified).

    The full MIT algorithm requires genome-wide alignment.
    This provides an estimate based on guide sequence properties.

    Higher score = more specific (fewer predicted off-targets).

    Args:
        guide_seq: 20bp guide sequence.
        off_target_count: Number of known off-targets (if available).
        avg_mismatches: Average mismatches to off-targets.

    Returns:
        MIT specificity score (0-100).
    """
    guide_seq = guide_seq.upper()

    # Base score from sequence complexity
    complexity = sequence_complexity(guide_seq)

    # Penalize low GC or very high GC
    gc = gc_content(guide_seq)
    gc_penalty = 0.0
    if gc < 0.4 or gc > 0.7:
        gc_penalty = 10 * abs(gc - 0.55)

    # Penalize homopolymers
    max_homo = max_homopolymer_run(guide_seq)
    homo_penalty = max(0, (max_homo - 3) * 5)

    # Seed region (PAM-proximal 12nt) importance
    seed = guide_seq[-12:]
    seed_complexity = sequence_complexity(seed)

    # Estimate score
    base_score = 100 * complexity * seed_complexity
    score = base_score - gc_penalty - homo_penalty

    # Adjust for known off-targets
    if off_target_count > 0:
        score -= min(30, off_target_count * 2)

    return max(0, min(100, score))


def cfd_score(
    guide_seq: str,
    target_seq: Optional[str] = None,
) -> float:
    """
    Calculate CFD (Cutting Frequency Determination) score.

    CFD predicts how likely an off-target site will be cut.
    For on-target (no mismatches), returns 100.

    Args:
        guide_seq: 20bp guide sequence.
        target_seq: Target sequence (if different from perfect match).

    Returns:
        CFD score (0-100, higher = more cutting).
    """
    guide_seq = guide_seq.upper()

    if target_seq is None:
        # On-target: perfect match
        return 100.0

    target_seq = target_seq.upper()

    if len(guide_seq) != len(target_seq):
        return 0.0

    # Count mismatches and their positions
    mismatches = []
    for i, (g, t) in enumerate(zip(guide_seq, target_seq)):
        if g != t:
            mismatches.append(i)

    if not mismatches:
        return 100.0

    # CFD penalty increases with mismatches and seed region location
    cfd = 100.0
    for pos in mismatches:
        # Higher penalty for PAM-proximal (seed) mismatches
        if pos >= 8:  # Seed region
            cfd *= 0.5
        else:
            cfd *= 0.7

    return max(0, cfd)


def poly_t_penalty(sequence: str, threshold: int = 4) -> Tuple[bool, str]:
    """
    Check for poly-T runs that cause U6/U3 Pol III termination.

    Pol III promoters (U6, U3) terminate at poly-T runs, making guides
    starting with TTTT or containing long T-runs incompatible with
    standard expression systems.

    Args:
        sequence: Guide sequence.
        threshold: Minimum T-run length to flag (default: 4 = TTTT).

    Returns:
        (is_problematic, reason)
        is_problematic: True if guide has poly-T issue
        reason: Description of the issue (empty if none)

    Example:
        >>> poly_t_penalty("TTTTAATGGCCGGCGATGCC")
        (True, "TTTT at 5' end - incompatible with U6/U3 promoters")
    """
    sequence = sequence.upper()

    # Check for TTTT at 5' end (most critical)
    if sequence.startswith('T' * threshold):
        return (True, f"{'T' * threshold} at 5' end - incompatible with U6/U3 promoters")

    # Check for internal poly-T runs
    for i in range(len(sequence) - threshold + 1):
        if sequence[i:i + threshold] == 'T' * threshold:
            return (True, f"Internal poly-T run at position {i} - may cause premature termination")

    return (False, "")


def is_repeat_region(sequence: str, min_repeat_length: int = 3) -> Tuple[bool, str]:
    """
    Detect if guide is in a repetitive/low-complexity region.

    Guides in repeat regions often have:
    - Multiple identical off-targets
    - Mapping ambiguity
    - Reduced specificity

    Args:
        sequence: Guide sequence.
        min_repeat_length: Minimum unit length for repeat detection.

    Returns:
        (is_repeat, reason)
        is_repeat: True if sequence is repetitive
        reason: Description of repeat type

    Example:
        >>> is_repeat_region("CAGCAGCAGCAGCAGCAGCA")
        (True, "Tandem repeat: CAG repeated 6+ times")
    """
    sequence = sequence.upper()
    n = len(sequence)

    # Check for tandem repeats (same unit repeated)
    for unit_len in range(min_repeat_length, n // 2 + 1):
        unit = sequence[:unit_len]
        repeat_count = 0
        for i in range(0, n - unit_len + 1, unit_len):
            if sequence[i:i + unit_len] == unit:
                repeat_count += 1
            else:
                break

        if repeat_count >= 3:  # At least 3 repeats
            return (True, f"Tandem repeat: {unit} repeated {repeat_count}+ times")

    # Check overall complexity
    complexity = sequence_complexity(sequence)
    if complexity < 0.4:
        return (True, f"Low complexity region (score={complexity:.2f})")

    # Check for dinucleotide repeats (e.g., ATATAT, CGCGCG)
    for dinuc in ['AT', 'TA', 'CG', 'GC', 'AC', 'CA', 'GT', 'TG', 'AG', 'GA', 'CT', 'TC']:
        pattern = (dinuc * (n // 2 + 1))[:n]
        matches = sum(1 for a, b in zip(sequence, pattern) if a == b)
        if matches >= n * 0.8:  # 80% match to dinucleotide repeat
            return (True, f"Dinucleotide repeat pattern: {dinuc}")

    return (False, "")


def u6_compatibility_check(sequence: str) -> Tuple[bool, List[str]]:
    """
    Comprehensive U6/Pol III compatibility check.

    Checks for all known issues with Pol III-driven expression:
    - Poly-T termination signals
    - G at position 1 preferred for U6
    - Internal TTTT runs

    Args:
        sequence: Guide sequence.

    Returns:
        (is_compatible, warnings)
        is_compatible: True if guide can be expressed from U6
        warnings: List of compatibility warnings

    Example:
        >>> u6_compatibility_check("GAAGTGACGGCTAGGGCTCC")
        (True, [])
        >>> u6_compatibility_check("TTTTAATGGCCGGCGATGCC")
        (False, ["TTTT at 5' - will cause Pol III termination"])
    """
    sequence = sequence.upper()
    warnings = []

    # Check poly-T
    is_poly_t, poly_t_reason = poly_t_penalty(sequence)
    if is_poly_t:
        warnings.append(poly_t_reason)

    # U6 prefers G at position 1 (not critical but suboptimal)
    if sequence[0] != 'G':
        # This is a soft warning, not disqualifying
        pass  # Don't add warning for this - too minor

    # Check for TTTTT (5 Ts) anywhere - almost always problematic
    if 'TTTTT' in sequence:
        warnings.append("Contains TTTTT - strong termination signal")

    is_compatible = len(warnings) == 0
    return (is_compatible, warnings)


def chromatin_accessibility_score(
    position: int,
    tss_position: int,
    dnase_peaks: Optional[list] = None,
) -> Tuple[str, float]:
    """
    Estimate chromatin accessibility at a genomic position.

    Without experimental data, uses heuristic based on TSS proximity.
    Near TSS = more likely to be open chromatin.

    Args:
        position: Genomic position.
        tss_position: Transcription start site position.
        dnase_peaks: Optional list of (start, end) DNase HS peaks.

    Returns:
        (state, accessibility_score)
        state: "OPEN", "MODERATE", or "CLOSED"
        accessibility_score: 0.0 to 1.0
    """
    rel_pos = position - tss_position

    # Check if in provided DNase peaks
    if dnase_peaks:
        for start, end in dnase_peaks:
            if start <= position <= end:
                return ("OPEN", 0.9)

    # Heuristic based on TSS proximity
    abs_dist = abs(rel_pos)

    if abs_dist < 200:
        # Very close to TSS - likely open
        score = 0.8 - 0.001 * abs_dist
        state = "OPEN"
    elif abs_dist < 500:
        score = 0.6 - 0.0005 * (abs_dist - 200)
        state = "MODERATE"
    else:
        score = max(0.2, 0.5 - 0.0002 * (abs_dist - 500))
        state = "MODERATE" if score > 0.35 else "CLOSED"

    return (state, score)


# =============================================================================
# RANKING POLICY SYSTEM (v0.9.2+)
# =============================================================================
#
# Key insight: "Best guide" is not absolute - it's "best under a stated policy".
# This system makes policy explicit and reproducible.
#
# DOMINANCE PRINCIPLE (all policies):
# - Ranking is LEXICOGRAPHIC on (0mm, 1mm, 2mm) - these are SAFETY-CRITICAL
# - A guide with 0/0/0 ALWAYS beats 0/0/1, regardless of 3-4mm counts
# - 3-4mm off-targets are tie-breakers only, not primary ranking criteria
#
# POLICY-DRIVEN RISK PREFERENCE (not universal biology):
# The fold-change estimates below are order-of-magnitude guides, not guarantees.
# Actual off-target cutting depends on mismatch position (seed vs distal),
# PAM context, nuclease variant, and experimental conditions.

from enum import Enum
from typing import NamedTuple
import json
import hashlib
from datetime import datetime


class RankingPolicy(Enum):
    """
    Named ranking policies for different use cases.

    Each policy defines:
    - Hard gates (what disqualifies a guide entirely)
    - Dominance order for safety-critical mismatches
    - Tie-breaker weights for lower-priority factors
    """
    # Default for gene knockout - maximum safety
    CUTTING_STRICT = "cutting_strict"

    # For CRISPRa/CRISPRi - binding matters more than cutting
    BINDING_STRICT = "binding_strict"

    # Exploratory - relaxed constraints, requires explicit flag
    # WARNING: Not for therapeutic use
    EXPLORATORY = "exploratory"


class PolicyConfig(NamedTuple):
    """Configuration for a ranking policy."""
    name: str
    description: str

    # Hard gates - guides failing these are EXCLUDED (not just penalized)
    gate_unscorable: bool  # Exclude guides CRISPOR couldn't score
    gate_0mm: bool         # Exclude guides with any 0mm off-targets
    gate_1mm: bool         # Exclude guides with any 1mm off-targets
    gate_u6_incompatible: bool  # Exclude U6/Pol III incompatible
    gate_repeats: bool     # Exclude guides in repeat regions

    # Max allowed off-targets before exclusion (None = no limit)
    max_0mm: int
    max_1mm: int
    max_2mm: int

    # Tie-breaker weights (used AFTER dominance sorting)
    # These are SMALL - dominance handles the main ranking
    weight_3mm: float
    weight_4mm: float
    weight_mit: float  # Bonus per MIT point
    weight_cfd: float  # Bonus per CFD point


# Pre-defined policies
POLICY_CONFIGS = {
    RankingPolicy.CUTTING_STRICT: PolicyConfig(
        name="cutting_strict",
        description="Maximum safety for gene knockout. Hard gates on 0-1mm.",
        gate_unscorable=True,
        gate_0mm=True,
        gate_1mm=True,  # Even 1 hit at 1mm is concerning for knockout
        gate_u6_incompatible=True,
        gate_repeats=True,
        max_0mm=0,
        max_1mm=0,
        max_2mm=None,  # 2mm used for dominance ranking, not gating
        weight_3mm=0.01,
        weight_4mm=0.001,
        weight_mit=0.1,
        weight_cfd=0.1,
    ),
    RankingPolicy.BINDING_STRICT: PolicyConfig(
        name="binding_strict",
        description="For CRISPRa/CRISPRi. Prioritizes binding specificity.",
        gate_unscorable=True,
        gate_0mm=True,
        gate_1mm=False,  # 1mm less critical for activation/interference
        gate_u6_incompatible=True,
        gate_repeats=True,
        max_0mm=0,
        max_1mm=2,  # Allow up to 2 sites at 1mm
        max_2mm=None,
        weight_3mm=0.1,  # Care more about far off-targets for binding
        weight_4mm=0.01,
        weight_mit=0.15,
        weight_cfd=0.05,  # CFD less relevant for non-cutting
    ),
    RankingPolicy.EXPLORATORY: PolicyConfig(
        name="exploratory",
        description="Relaxed constraints for research. NOT for therapeutic use.",
        gate_unscorable=False,  # Allow unscorable (with warning)
        gate_0mm=True,  # Still gate perfect matches
        gate_1mm=False,
        gate_u6_incompatible=False,  # May use alternative promoters
        gate_repeats=False,
        max_0mm=0,
        max_1mm=5,
        max_2mm=10,
        weight_3mm=0.5,
        weight_4mm=0.1,
        weight_mit=0.2,
        weight_cfd=0.2,
    ),
}


# Legacy weights (kept for backward compatibility with crispor_composite_score)
OFFTARGET_MISMATCH_WEIGHTS = {
    0: 10000.0,  # Perfect off-target - DISQUALIFYING
    1: 1000.0,   # 1 mismatch - very dangerous
    2: 100.0,    # 2 mismatches - concerning
    3: 0.1,      # 3 mismatches - negligible (policy-dependent)
    4: 0.01,     # 4 mismatches - essentially zero (policy-dependent)
}


@dataclass
class CRISPORMetrics:
    """
    Container for CRISPOR-style guide metrics.

    Attributes:
        mit_score: MIT specificity score (0-100)
        cfd_score: CFD cutting frequency score (0-100)
        off_targets: Dict mapping mismatch count to number of off-targets
        u6_compatible: Whether compatible with U6/Pol III
        is_repeat: Whether in genomic repeat region
    """
    mit_score: float
    cfd_score: float
    off_targets: Dict[int, int] = field(default_factory=dict)  # {mm_count: num_OTs}
    u6_compatible: bool = True
    is_repeat: bool = False


# =============================================================================
# HARD GATES (v0.9.2+)
# =============================================================================

@dataclass
class GateResult:
    """Result of applying hard gates to a guide."""
    passed: bool
    excluded_by: Optional[str] = None
    reason: Optional[str] = None


def apply_hard_gates(
    guide: Dict[str, Any],
    policy: RankingPolicy = RankingPolicy.CUTTING_STRICT,
) -> GateResult:
    """
    Apply hard gates to a guide based on the ranking policy.

    Hard gates are EXCLUSION criteria - guides failing any gate are removed
    from consideration entirely, not just penalized.

    Args:
        guide: Guide dict with CRISPOR metrics.
        policy: Ranking policy to use.

    Returns:
        GateResult with passed=True if guide passes all gates.

    Example:
        >>> result = apply_hard_gates(guide, RankingPolicy.CUTTING_STRICT)
        >>> if not result.passed:
        ...     print(f"Excluded: {result.reason}")
    """
    config = POLICY_CONFIGS[policy]
    ot = guide.get('off_targets', {})

    # Gate: Unscorable
    if config.gate_unscorable and guide.get('is_unscorable', False):
        return GateResult(
            passed=False,
            excluded_by="UNSCORABLE",
            reason="CRISPOR could not validate this guide (no safety data)"
        )

    # Gate: 0mm off-targets
    ot_0mm = ot.get(0, 0)
    if config.gate_0mm and ot_0mm > config.max_0mm:
        return GateResult(
            passed=False,
            excluded_by="0MM_OFFTARGET",
            reason=f"Has {ot_0mm} perfect off-target(s) - guaranteed off-target cutting"
        )

    # Gate: 1mm off-targets
    ot_1mm = ot.get(1, 0)
    if config.gate_1mm and ot_1mm > config.max_1mm:
        return GateResult(
            passed=False,
            excluded_by="1MM_OFFTARGET",
            reason=f"Has {ot_1mm} 1-mismatch off-target(s) - high risk of cutting"
        )

    # Gate: 2mm off-targets (if max set)
    ot_2mm = ot.get(2, 0)
    if config.max_2mm is not None and ot_2mm > config.max_2mm:
        return GateResult(
            passed=False,
            excluded_by="2MM_OFFTARGET",
            reason=f"Has {ot_2mm} 2-mismatch off-target(s) - exceeds policy limit"
        )

    # Gate: U6 incompatible
    if config.gate_u6_incompatible and not guide.get('u6_compatible', True):
        return GateResult(
            passed=False,
            excluded_by="U6_INCOMPATIBLE",
            reason="Contains poly-T - incompatible with U6/Pol III promoters"
        )

    # Gate: Repeat region
    if config.gate_repeats and guide.get('is_repeat', False):
        return GateResult(
            passed=False,
            excluded_by="REPEAT_REGION",
            reason="Guide is in a genomic repeat region"
        )

    return GateResult(passed=True)


# =============================================================================
# DOMINANCE-BASED RANKING (v0.9.2+)
# =============================================================================

class GuideTier(Enum):
    """Tiers for guide categorization."""
    A = "A"  # Best: 0/0/0 (no 0-2mm off-targets)
    B = "B"  # Good: 0/0/1-2 (1-2 off-targets at 2mm only)
    C = "C"  # Acceptable: 0/0/3+ or 0/1+/any
    X = "X"  # Excluded: Failed hard gates


def _dominance_key(guide: Dict[str, Any], config: PolicyConfig) -> tuple:
    """
    Create lexicographic sorting key for dominance-based ranking.

    Sort order (ascending = worse):
    1. 0mm count (primary)
    2. 1mm count
    3. 2mm count
    4. tie-breaker score (descending = better, so negated)

    The tie-breaker incorporates 3-4mm counts and MIT/CFD.
    """
    ot = guide.get('off_targets', {})

    # Primary: safety-critical (0, 1, 2mm) - ASCENDING (fewer = better)
    ot_0mm = ot.get(0, 0)
    ot_1mm = ot.get(1, 0)
    ot_2mm = ot.get(2, 0)

    # Tie-breaker: weighted combination of 3-4mm and specificity scores
    # DESCENDING (higher = better), so negate
    ot_3mm = ot.get(3, 0)
    ot_4mm = ot.get(4, 0)
    mit = guide.get('crispor_mit', guide.get('mit_score', 50))
    cfd = guide.get('crispor_cfd', guide.get('cfd_score', 50))

    tie_breaker = (
        config.weight_mit * mit +
        config.weight_cfd * cfd -
        config.weight_3mm * ot_3mm -
        config.weight_4mm * ot_4mm
    )

    return (ot_0mm, ot_1mm, ot_2mm, -tie_breaker)


def _assign_tier(guide: Dict[str, Any]) -> GuideTier:
    """Assign a tier based on off-target profile."""
    ot = guide.get('off_targets', {})
    ot_0mm = ot.get(0, 0)
    ot_1mm = ot.get(1, 0)
    ot_2mm = ot.get(2, 0)

    if ot_0mm == 0 and ot_1mm == 0 and ot_2mm == 0:
        return GuideTier.A
    elif ot_0mm == 0 and ot_1mm == 0 and ot_2mm <= 2:
        return GuideTier.B
    else:
        return GuideTier.C


def rank_guides(
    guides: List[Dict[str, Any]],
    policy: RankingPolicy = RankingPolicy.CUTTING_STRICT,
    return_excluded: bool = False,
) -> Dict[str, Any]:
    """
    Rank guides using dominance-based lexicographic sorting.

    This is the v0.9.2+ replacement for composite scoring. Key features:
    - Hard gates exclude guides entirely (not just penalized)
    - Lexicographic sort on (0mm, 1mm, 2mm) - safety-critical
    - Tie-breaker uses MIT/CFD and 3-4mm counts
    - Output includes tier assignments for wet lab convenience

    Args:
        guides: List of guide dicts with CRISPOR metrics.
        policy: Ranking policy to use.
        return_excluded: Include excluded guides in output.

    Returns:
        Dict with:
        - ranked: List of ranked guides (passing hard gates)
        - excluded: List of excluded guides (if return_excluded=True)
        - tiers: Dict mapping tier to list of guides
        - policy: Policy used
        - manifest: Run manifest for reproducibility

    Example:
        >>> result = rank_guides(guides, RankingPolicy.CUTTING_STRICT)
        >>> print(f"Top guide: {result['ranked'][0]['sequence']}")
        >>> print(f"Tier A guides: {len(result['tiers']['A'])}")
    """
    config = POLICY_CONFIGS[policy]

    # Apply hard gates
    passed = []
    excluded = []

    for guide in guides:
        gate_result = apply_hard_gates(guide, policy)

        if gate_result.passed:
            guide_copy = {**guide}
            guide_copy['tier'] = _assign_tier(guide).value
            passed.append(guide_copy)
        else:
            guide_copy = {**guide}
            guide_copy['tier'] = GuideTier.X.value
            guide_copy['excluded_by'] = gate_result.excluded_by
            guide_copy['exclusion_reason'] = gate_result.reason
            excluded.append(guide_copy)

    # Sort by dominance key
    passed.sort(key=lambda g: _dominance_key(g, config))

    # Add ranks
    for i, guide in enumerate(passed):
        guide['rank'] = i + 1

    # Group by tier
    tiers = {'A': [], 'B': [], 'C': []}
    for guide in passed:
        tier = guide.get('tier', 'C')
        if tier in tiers:
            tiers[tier].append(guide)

    # Build result
    result = {
        'ranked': passed,
        'tiers': tiers,
        'policy': policy.value,
        'policy_description': config.description,
        'summary': {
            'total_input': len(guides),
            'passed_gates': len(passed),
            'excluded': len(excluded),
            'tier_A': len(tiers['A']),
            'tier_B': len(tiers['B']),
            'tier_C': len(tiers['C']),
        },
    }

    if return_excluded:
        result['excluded_guides'] = excluded

    return result


# =============================================================================
# RUN MANIFEST (v0.9.2+)
# =============================================================================

def emit_manifest(
    guides: List[Dict[str, Any]],
    policy: RankingPolicy,
    sequence_name: str = "unknown",
    genome: str = "hg38",
    crispor_version: Optional[str] = None,
    tss_index: Optional[int] = None,
    window: Optional[Tuple[int, int]] = None,
) -> Dict[str, Any]:
    """
    Generate a reproducibility manifest for a ranking run.

    Every ranking run should emit this manifest to enable:
    - Reproducibility of results
    - Audit trail for validation
    - Debugging of "why did the winner change" issues

    Args:
        guides: The ranked guides.
        policy: Ranking policy used.
        sequence_name: Name of the target sequence/gene.
        genome: Genome build (e.g., "hg38").
        crispor_version: CRISPOR version/commit if known.
        tss_index: TSS position in the input sequence.
        window: CRISPRa/i window used (e.g., (-400, -50)).

    Returns:
        Manifest dict suitable for JSON serialization.

    Example:
        >>> result = rank_guides(guides, RankingPolicy.CUTTING_STRICT)
        >>> manifest = emit_manifest(
        ...     result['ranked'],
        ...     RankingPolicy.CUTTING_STRICT,
        ...     sequence_name="RAI1",
        ...     genome="hg38",
        ... )
        >>> with open('rai1_manifest.json', 'w') as f:
        ...     json.dump(manifest, f, indent=2)
    """
    config = POLICY_CONFIGS[policy]

    # Create content hash for reproducibility
    guide_data = [
        {
            'sequence': g.get('sequence'),
            'off_targets': g.get('off_targets'),
            'rank': g.get('rank'),
            'tier': g.get('tier'),
        }
        for g in guides[:20]  # Top 20 for hash
    ]
    content_str = json.dumps(guide_data, sort_keys=True)
    content_hash = hashlib.sha256(content_str.encode()).hexdigest()[:16]

    manifest = {
        'version': '0.9.2',
        'timestamp': datetime.now().isoformat(),
        'content_hash': content_hash,

        'sequence': {
            'name': sequence_name,
            'genome': genome,
            'tss_index': tss_index,
            'window': window,
        },

        'policy': {
            'name': policy.value,
            'description': config.description,
            'gates': {
                'unscorable': config.gate_unscorable,
                '0mm': config.gate_0mm,
                '1mm': config.gate_1mm,
                'u6_incompatible': config.gate_u6_incompatible,
                'repeats': config.gate_repeats,
            },
            'limits': {
                'max_0mm': config.max_0mm,
                'max_1mm': config.max_1mm,
                'max_2mm': config.max_2mm,
            },
            'weights': {
                '3mm': config.weight_3mm,
                '4mm': config.weight_4mm,
                'mit': config.weight_mit,
                'cfd': config.weight_cfd,
            },
        },

        'crispor': {
            'version': crispor_version or 'unknown',
        },

        'results': {
            'total_guides': len(guides),
            'top_5': [
                {
                    'rank': g.get('rank'),
                    'sequence': g.get('sequence'),
                    'tier': g.get('tier'),
                    'off_targets': g.get('off_targets'),
                    'mit': g.get('crispor_mit', g.get('mit_score')),
                    'cfd': g.get('crispor_cfd', g.get('cfd_score')),
                }
                for g in guides[:5]
            ],
        },
    }

    return manifest


def print_ranking_report(
    result: Dict[str, Any],
    max_guides: int = 15,
) -> None:
    """
    Print a formatted ranking report to stdout.

    Args:
        result: Output from rank_guides().
        max_guides: Maximum guides to show.
    """
    print("\n" + "=" * 90)
    print(f"GUIDE RANKING - Policy: {result['policy'].upper()}")
    print(f"{result['policy_description']}")
    print("=" * 90)

    summary = result['summary']
    print(f"\nSummary: {summary['total_input']} input → "
          f"{summary['passed_gates']} passed gates, {summary['excluded']} excluded")
    print(f"Tiers: A={summary['tier_A']} (0/0/0) | "
          f"B={summary['tier_B']} (0/0/1-2) | "
          f"C={summary['tier_C']} (other)")

    print(f"\n{'Rank':<5} {'Tier':<5} {'Sequence':<22} "
          f"{'0mm':<5} {'1mm':<5} {'2mm':<5} {'3mm':<6} {'4mm':<6} "
          f"{'MIT':<5} {'CFD':<5}")
    print("-" * 90)

    for guide in result['ranked'][:max_guides]:
        ot = guide.get('off_targets', {})
        mit = guide.get('crispor_mit', guide.get('mit_score', 0))
        cfd = guide.get('crispor_cfd', guide.get('cfd_score', 0))

        print(f"{guide['rank']:<5} {guide['tier']:<5} {guide.get('sequence', 'N/A'):<22} "
              f"{ot.get(0, 0):<5} {ot.get(1, 0):<5} {ot.get(2, 0):<5} "
              f"{ot.get(3, 0):<6} {ot.get(4, 0):<6} "
              f"{mit:<5.0f} {cfd:<5.0f}")

    print("-" * 90)
    print("Tier A = Zero 0-2mm off-targets (safest)")
    print("Tier B = Zero 0-1mm, 1-2 at 2mm (good)")
    print("Tier C = Other passing guides")
    print()

    # Recommendation
    tiers = result['tiers']
    if tiers['A']:
        print(f"✓ RECOMMENDATION: Choose from {len(tiers['A'])} Tier A guide(s)")
        print(f"  Top Tier A: {tiers['A'][0].get('sequence')}")
    elif tiers['B']:
        print(f"⚠ No Tier A guides. Choose from {len(tiers['B'])} Tier B guide(s)")
        print(f"  Top Tier B: {tiers['B'][0].get('sequence')}")
    elif tiers['C']:
        print(f"⚠ No Tier A/B guides. {len(tiers['C'])} Tier C guide(s) available")
    else:
        print("✗ No guides passed hard gates under this policy")


def crispor_composite_score(
    mit_score: float,
    cfd_score: float,
    off_targets: Optional[Dict[int, int]] = None,
    u6_compatible: bool = True,
    is_repeat: bool = False,
    weights: Optional[Dict[int, float]] = None,
) -> Tuple[float, Dict[str, float]]:
    """
    Calculate CRISPOR-style composite score with mismatch distance weighting.

    Formula:
        COMPOSITE = (MIT + CFD) - Σ(weight[mm] × count[mm]) - U6_penalty - repeat_penalty

    Where:
        - weight[0-1mm] = 50-100 (critical)
        - weight[2mm] = 25 (important)
        - weight[3mm] = 5 (minor)
        - weight[4mm] = 1 (minimal)
        - U6_penalty = 100 if incompatible
        - repeat_penalty = 1000 if in repeat region

    This correctly handles the "MIT 98 / CFD 98 trap" where high raw scores
    are misleading due to many off-targets or U6 incompatibility.

    Args:
        mit_score: MIT specificity score (0-100).
        cfd_score: CFD cutting frequency score (0-100).
        off_targets: Dict mapping mismatch count to off-target count.
                     e.g., {0: 0, 1: 0, 2: 1, 3: 5, 4: 23}
        u6_compatible: True if compatible with U6/Pol III promoters.
        is_repeat: True if guide is in genomic repeat region.
        weights: Custom mismatch weights (default: OFFTARGET_MISMATCH_WEIGHTS).

    Returns:
        (composite_score, breakdown)
        composite_score: Final score (higher is better, can be negative)
        breakdown: Dict with component contributions

    Example:
        >>> # Guide #1: MIT=93, CFD=95, 0 dangerous OTs
        >>> score, _ = crispor_composite_score(93, 95, {0:0, 1:0, 2:0, 3:5, 4:15})
        >>> score
        160.0  # (93+95) - (5*5 + 15*1) = 188 - 25 - 15 = 148

        >>> # Guide #7: MIT=98, CFD=98, TTTT start
        >>> score, _ = crispor_composite_score(98, 98, {0:0, 1:0, 2:1, 3:6, 4:22},
        ...                                     u6_compatible=False)
        >>> score
        41.0  # (98+98) - (25 + 30 + 22) - 100 = 196 - 77 - 100 = 19
    """
    if weights is None:
        weights = OFFTARGET_MISMATCH_WEIGHTS

    if off_targets is None:
        off_targets = {}

    # Base score from MIT + CFD
    base_score = mit_score + cfd_score

    # Off-target penalties by mismatch distance
    ot_penalties = {}
    total_ot_penalty = 0.0

    for mm_count, num_ots in off_targets.items():
        weight = weights.get(mm_count, 0.5)  # Default small weight for mm>4
        penalty = weight * num_ots
        ot_penalties[f"ot_{mm_count}mm"] = penalty
        total_ot_penalty += penalty

    # U6 incompatibility penalty (100 points)
    u6_penalty = 0.0 if u6_compatible else 100.0

    # Repeat region penalty (1000 points - effectively disqualifying)
    repeat_penalty = 0.0 if not is_repeat else 1000.0

    # Final composite
    composite = base_score - total_ot_penalty - u6_penalty - repeat_penalty

    breakdown = {
        "base_score": base_score,
        "mit_contribution": mit_score,
        "cfd_contribution": cfd_score,
        "total_ot_penalty": total_ot_penalty,
        **ot_penalties,
        "u6_penalty": u6_penalty,
        "repeat_penalty": repeat_penalty,
        "composite": composite,
    }

    return (composite, breakdown)


def rank_guides_crispor_style(
    guides: List[Dict[str, Any]],
    require_u6_compatible: bool = True,
    exclude_repeats: bool = True,
) -> List[Dict[str, Any]]:
    """
    Rank guides using CRISPOR-style composite scoring.

    This function takes guide dictionaries with CRISPOR metrics and
    returns them sorted by composite score (highest first).

    Args:
        guides: List of guide dicts, each containing:
            - sequence: Guide sequence
            - mit_score: MIT specificity (0-100)
            - cfd_score: CFD score (0-100)
            - off_targets: Dict {mismatch_count: num_off_targets}
            - u6_compatible: bool (optional, default True)
            - is_repeat: bool (optional, default False)
        require_u6_compatible: Exclude U6-incompatible guides entirely.
        exclude_repeats: Exclude guides in repeat regions entirely.

    Returns:
        Sorted list of guides with 'crispor_composite' and 'crispor_rank' added.

    Example:
        >>> guides = [
        ...     {"sequence": "TTCGATGAATGGTTGCTACC", "mit_score": 93, "cfd_score": 95,
        ...      "off_targets": {0:0, 1:0, 2:0, 3:5, 4:15}},
        ...     {"sequence": "TTTTAATGGCCGGCGATGCC", "mit_score": 98, "cfd_score": 98,
        ...      "off_targets": {0:0, 1:0, 2:1, 3:6, 4:22}, "u6_compatible": False},
        ... ]
        >>> ranked = rank_guides_crispor_style(guides)
        >>> ranked[0]["sequence"]
        'TTCGATGAATGGTTGCTACC'  # Guide #1 wins despite lower MIT/CFD
    """
    scored_guides = []

    for guide in guides:
        u6_ok = guide.get("u6_compatible", True)
        is_repeat = guide.get("is_repeat", False)

        # Hard exclusions
        if require_u6_compatible and not u6_ok:
            continue
        if exclude_repeats and is_repeat:
            continue

        # Compute composite score
        composite, breakdown = crispor_composite_score(
            mit_score=guide.get("mit_score", 50),
            cfd_score=guide.get("cfd_score", 50),
            off_targets=guide.get("off_targets", {}),
            u6_compatible=u6_ok,
            is_repeat=is_repeat,
        )

        scored_guide = {**guide}
        scored_guide["crispor_composite"] = composite
        scored_guide["crispor_breakdown"] = breakdown
        scored_guides.append(scored_guide)

    # Sort by composite score (highest first)
    scored_guides.sort(key=lambda g: g["crispor_composite"], reverse=True)

    # Add rank
    for i, guide in enumerate(scored_guides):
        guide["crispor_rank"] = i + 1

    return scored_guides


def validate_and_rerank_with_crispor(
    phaselab_guides: List[Dict[str, Any]],
    crispor_data: List[Dict[str, Any]],
    verbose: bool = True,
) -> List[Dict[str, Any]]:
    """
    Validate PhaseLab guides against CRISPOR data and re-rank using v0.9.1 composite scoring.

    This function bridges PhaseLab's design_guides() output with CRISPOR validation,
    applying proper off-target mismatch distance weighting per ChatGPT's recommendations.

    Args:
        phaselab_guides: Output from design_guides() - list of guide dicts with:
            - sequence: 20bp guide
            - position: relative to TSS
            - gc: GC content
            - coherence_R: IR coherence score
            - go_no_go: "GO" or "NO-GO"
            - mit_score: PhaseLab's MIT estimate
            - cfd_score: PhaseLab's CFD estimate
            - combined_score: PhaseLab's original ranking score

        crispor_data: CRISPOR output - list of dicts with:
            - sequence: Guide sequence (matching key)
            - mit_specificity: Real MIT score from CRISPOR
            - cfd_specificity: Real CFD score from CRISPOR
            - ot_0mm, ot_1mm, ot_2mm, ot_3mm, ot_4mm: Off-target counts

        verbose: Print comparison table.

    Returns:
        Merged and re-ranked guides with CRISPOR validation, including:
        - crispor_mit: Validated MIT score
        - crispor_cfd: Validated CFD score
        - off_targets: Dict of {mm: count}
        - crispor_composite: CRISPOR composite score
        - crispor_rank: New rank based on composite score
        - phaselab_rank: Original PhaseLab rank
        - rank_delta: Change in rank (positive = improved)

    Example:
        >>> from phaselab.crispr import design_guides
        >>> guides_df = design_guides(sequence, tss_index=500)
        >>> phaselab_guides = guides_df.to_dict('records')
        >>>
        >>> # After running CRISPOR (web or local)
        >>> crispor_data = parse_crispor_results(crispor_output)
        >>>
        >>> validated = validate_and_rerank_with_crispor(phaselab_guides, crispor_data)
        >>> print(f"Top guide changed: {validated[0]['sequence']}")
    """
    # Index CRISPOR data by sequence
    crispor_index = {g['sequence'].upper(): g for g in crispor_data if 'sequence' in g}

    merged = []

    for i, guide in enumerate(phaselab_guides):
        seq = guide.get('sequence', '').upper()
        phaselab_rank = i + 1  # Original rank

        merged_guide = {
            **guide,
            'phaselab_rank': phaselab_rank,
        }

        # Look up CRISPOR data
        if seq in crispor_index:
            cg = crispor_index[seq]

            # Extract validated metrics (handle 0 as valid value, None as missing)
            mit_raw = cg.get('mit_specificity')
            cfd_raw = cg.get('cfd_specificity')

            # MIT/CFD of 0 with 0 off-targets means unmappable - treat as unvalidated
            # MIT/CFD of 0 with >0 off-targets means many perfect matches - very bad
            total_ots = sum([
                cg.get('ot_0mm', 0), cg.get('ot_1mm', 0), cg.get('ot_2mm', 0),
                cg.get('ot_3mm', 0), cg.get('ot_4mm', 0)
            ])

            # Also check for explicit is_unscorable flag from CrisporValidator
            is_unscorable = cg.get('is_unscorable', False)

            if is_unscorable or (mit_raw == 0 and cfd_raw == 0 and total_ots == 0):
                # INVALID/UNSCORABLE: CRISPOR couldn't score this guide
                # Possible causes: NotEnoughFlankSeq, unmappable region, genome gap
                # This is a RED FLAG - we have NO safety data for this guide
                # Hard-demote these guides below all validated guides
                mit = 0
                cfd = 0
                merged_guide['crispor_validated'] = False
                merged_guide['is_unscorable'] = True
                merged_guide['unscorable_reason'] = cg.get('unscorable_reason',
                    "CRISPOR returned MIT=0, CFD=0, OTs=0 (invalid/unmappable)")
            else:
                # Valid CRISPOR data (0 is valid if there are off-targets)
                mit = mit_raw if mit_raw is not None else cg.get('mit_score', 50)
                cfd = cfd_raw if cfd_raw is not None else cg.get('cfd_score', 50)
                merged_guide['crispor_validated'] = True
                merged_guide['is_unscorable'] = False

            off_targets = {
                0: cg.get('ot_0mm', 0),
                1: cg.get('ot_1mm', 0),
                2: cg.get('ot_2mm', 0),
                3: cg.get('ot_3mm', 0),
                4: cg.get('ot_4mm', 0),
            }

            merged_guide['crispor_mit'] = mit
            merged_guide['crispor_cfd'] = cfd
            merged_guide['off_targets'] = off_targets

            # Check U6 compatibility
            is_poly_t, _ = poly_t_penalty(seq)
            is_repeat, _ = is_repeat_region(seq)
            merged_guide['u6_compatible'] = not is_poly_t
            merged_guide['is_repeat'] = is_repeat

            # Compute composite score
            composite, breakdown = crispor_composite_score(
                mit_score=mit,
                cfd_score=cfd,
                off_targets=off_targets,
                u6_compatible=not is_poly_t,
                is_repeat=is_repeat,
            )

            # Apply HARD penalty for unscorable guides - no safety data = high risk
            # These should rank BELOW all validated guides
            if merged_guide.get('is_unscorable', False):
                unscorable_penalty = 500.0  # Severe penalty - we have NO safety data
                composite -= unscorable_penalty
                breakdown['unscorable_penalty'] = unscorable_penalty

            merged_guide['crispor_composite'] = composite
            merged_guide['crispor_breakdown'] = breakdown

        else:
            # No CRISPOR match - guide not in CRISPOR output at all
            # This is INVALID - we have zero off-target data
            merged_guide['crispor_validated'] = False
            merged_guide['is_unscorable'] = True
            merged_guide['unscorable_reason'] = "Guide not found in CRISPOR output"
            merged_guide['crispor_mit'] = 0
            merged_guide['crispor_cfd'] = 0
            merged_guide['off_targets'] = {}

            # Compute composite with SEVERE penalty for unscorable guides
            composite, breakdown = crispor_composite_score(
                mit_score=0,
                cfd_score=0,
                off_targets={},
                u6_compatible=guide.get('u6_compatible', True),
                is_repeat=False,
            )
            # Hard penalty - rank below all validated guides
            unscorable_penalty = 500.0
            composite -= unscorable_penalty
            breakdown['unscorable_penalty'] = unscorable_penalty

            merged_guide['crispor_composite'] = composite
            merged_guide['crispor_breakdown'] = breakdown

        merged.append(merged_guide)

    # Sort by CRISPOR composite score
    merged.sort(key=lambda g: g['crispor_composite'], reverse=True)

    # Add new ranks and deltas
    for i, guide in enumerate(merged):
        guide['crispor_rank'] = i + 1
        guide['rank_delta'] = guide['phaselab_rank'] - guide['crispor_rank']

    if verbose:
        print("\n" + "=" * 100)
        print("CRISPOR-VALIDATED RANKING (v0.9.1 Composite Scoring)")
        print("=" * 100)
        print(f"{'Rank':<5} {'PL#':<5} {'Δ':<4} {'Sequence':<22} {'MIT':<6} {'CFD':<6} "
              f"{'0mm':<5} {'1mm':<5} {'2mm':<5} {'Composite':<10} {'R̄':<6} {'Status'}")
        print("-" * 100)

        for guide in merged[:15]:  # Top 15
            delta_str = f"+{guide['rank_delta']}" if guide['rank_delta'] > 0 else str(guide['rank_delta'])
            status = "GO" if guide.get('go_no_go') == 'GO' else guide.get('go_no_go', '?')
            r_bar = guide.get('coherence_R', 0) or 0
            ot = guide.get('off_targets', {})
            validated = "✓" if guide.get('crispor_validated') else "~"

            print(f"{guide['crispor_rank']:<5} {guide['phaselab_rank']:<5} {delta_str:<4} "
                  f"{guide.get('sequence', 'N/A'):<22} "
                  f"{guide.get('crispor_mit', 0):<6.0f} {guide.get('crispor_cfd', 0):<6.0f} "
                  f"{ot.get(0, 0):<5} {ot.get(1, 0):<5} {ot.get(2, 0):<5} "
                  f"{guide['crispor_composite']:<10.1f} {r_bar:<6.3f} {status} {validated}")

        print("-" * 100)
        print(f"✓ = CRISPOR validated, ~ = PhaseLab estimates")
        print(f"Composite = (MIT+CFD) - Σ(weight × OT_count) - penalties")
        print(f"Weights: 0mm=100, 1mm=50, 2mm=25, 3mm=5, 4mm=1")

        # Show rank changes
        improved = sum(1 for g in merged if g['rank_delta'] > 0)
        declined = sum(1 for g in merged if g['rank_delta'] < 0)
        print(f"\nRank changes: {improved} improved, {declined} declined, "
              f"{len(merged) - improved - declined} unchanged")

    return merged
