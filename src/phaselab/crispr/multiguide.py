"""
PhaseLab Multi-Guide Synergy: Model interactions between multiple CRISPR guides.

Implements:
- Synergy prediction for CRISPRa guide combinations
- Steric clash detection between nearby guides
- Optimal spacing determination
- Combinatorial guide set design
- IR coherence for multi-guide reliability

References:
- CRISPR-M: Multi-view deep learning for predictions
- Combinatorial CRISPR screens (AI-driven prioritization)
- Multi-guide CRISPRa for complex promoters
- Experimental data on guide synergy (Konermann et al.)

Version: 0.5.0
"""

import numpy as np
from typing import Dict, List, Optional, Tuple, Union, Any
from dataclasses import dataclass, field
from itertools import combinations
import logging

from ..core.coherence import coherence_score, go_no_go
from ..core.constants import E_MINUS_2

logger = logging.getLogger(__name__)


# Empirical spacing parameters for guide interactions
# Based on published CRISPRa synergy studies
GUIDE_SPACING_PARAMS = {
    # Minimum spacing before steric clash (bp)
    'min_spacing': 30,

    # Optimal spacing range for synergy (bp)
    'optimal_spacing_min': 50,
    'optimal_spacing_max': 200,

    # Maximum spacing for significant interaction (bp)
    'max_interaction_distance': 500,

    # Steric clash radius (bp) - guides closer than this interfere
    'clash_radius': 25,
}

# Synergy coefficients by relative position
# Based on experimental CRISPRa data
POSITION_SYNERGY = {
    'upstream_upstream': 0.9,    # Both upstream of TSS
    'upstream_downstream': 0.7,  # One each side
    'downstream_downstream': 0.5, # Both downstream (less effective)
    'same_region': 0.8,          # Same promoter region
}


@dataclass
class GuideCandidate:
    """A guide candidate for multi-guide analysis."""
    sequence: str
    position: int              # Position relative to TSS
    strand: str                # '+' or '-'
    individual_score: float    # Score without considering other guides
    coherence_R: float = 0.0   # IR coherence
    name: str = ""

    # Computed properties
    efficiency: float = 0.0
    specificity: float = 0.0


@dataclass
class GuidePair:
    """Analysis of a pair of guides."""
    guide1: GuideCandidate
    guide2: GuideCandidate

    # Spacing analysis
    spacing: int
    has_steric_clash: bool
    in_optimal_range: bool

    # Synergy prediction
    synergy_score: float       # 0-2 scale (>1 = synergistic, <1 = antagonistic)
    combined_efficiency: float

    # Interaction type
    interaction_type: str      # 'synergistic', 'additive', 'antagonistic', 'clash'


@dataclass
class MultiGuideSet:
    """A set of guides designed to work together."""
    guides: List[GuideCandidate]
    pairs: List[GuidePair]

    # Overall scores
    combined_score: float
    synergy_factor: float      # Average synergy across pairs
    ensemble_coherence: float  # Combined IR coherence

    # Quality metrics
    has_clashes: bool
    min_spacing: int
    max_spacing: int

    # Recommendation
    recommendation: str
    go_no_go: str


@dataclass
class MultiGuideConfig:
    """Configuration for multi-guide analysis."""

    # Spacing parameters
    min_spacing: int = 30
    optimal_spacing_min: int = 50
    optimal_spacing_max: int = 200
    max_interaction_distance: int = 500
    clash_radius: int = 25

    # Synergy modeling
    enable_synergy: bool = True
    synergy_model: str = "empirical"  # 'empirical', 'learned', 'none'

    # Combination parameters
    max_guides: int = 4        # Maximum guides in a set
    min_guides: int = 2        # Minimum guides in a set

    # Scoring weights
    weight_individual: float = 0.4
    weight_synergy: float = 0.3
    weight_coherence: float = 0.3

    # IR coherence
    compute_coherence: bool = True


def calculate_spacing(guide1: GuideCandidate, guide2: GuideCandidate) -> int:
    """
    Calculate spacing between two guides.

    Args:
        guide1: First guide
        guide2: Second guide

    Returns:
        Absolute distance in bp
    """
    return abs(guide1.position - guide2.position)


def check_steric_clash(
    guide1: GuideCandidate,
    guide2: GuideCandidate,
    config: Optional[MultiGuideConfig] = None
) -> bool:
    """
    Check if two guides would have steric clash.

    Steric clash occurs when:
    - Guides are too close (<30bp typically)
    - dCas9 proteins would physically interfere
    - Same strand and overlapping

    Args:
        guide1: First guide
        guide2: Second guide
        config: MultiGuideConfig

    Returns:
        True if steric clash expected
    """
    if config is None:
        config = MultiGuideConfig()

    spacing = calculate_spacing(guide1, guide2)

    # Basic distance check
    if spacing < config.clash_radius:
        return True

    # Same strand amplifies clash
    if guide1.strand == guide2.strand and spacing < config.min_spacing:
        return True

    # Opposite strands can be slightly closer
    if guide1.strand != guide2.strand and spacing < config.clash_radius * 0.8:
        return True

    return False


def predict_pairwise_synergy(
    guide1: GuideCandidate,
    guide2: GuideCandidate,
    tss_position: int = 0,
    config: Optional[MultiGuideConfig] = None,
) -> float:
    """
    Predict synergy between two guides.

    Synergy > 1.0 means combined effect is greater than sum of parts.
    Synergy < 1.0 means interference or redundancy.

    Args:
        guide1: First guide
        guide2: Second guide
        tss_position: TSS position for relative calculations
        config: MultiGuideConfig

    Returns:
        Synergy coefficient (0-2, with 1.0 = additive)
    """
    if config is None:
        config = MultiGuideConfig()

    spacing = calculate_spacing(guide1, guide2)

    # Check for clash first
    if check_steric_clash(guide1, guide2, config):
        return 0.0  # Complete interference

    # Base synergy from spacing
    if config.optimal_spacing_min <= spacing <= config.optimal_spacing_max:
        # Optimal range - maximum synergy potential
        spacing_factor = 1.2
    elif spacing < config.optimal_spacing_min:
        # Too close - some interference
        spacing_factor = 0.7 + 0.3 * (spacing / config.optimal_spacing_min)
    elif spacing <= config.max_interaction_distance:
        # Beyond optimal but still interacting
        decay = (spacing - config.optimal_spacing_max) / (
            config.max_interaction_distance - config.optimal_spacing_max
        )
        spacing_factor = 1.2 - 0.5 * decay
    else:
        # Too far apart - independent action
        return 1.0  # Additive only

    # Position-based synergy (relative to TSS)
    pos1_rel = guide1.position - tss_position
    pos2_rel = guide2.position - tss_position

    # Both upstream of TSS
    if pos1_rel < 0 and pos2_rel < 0:
        position_factor = POSITION_SYNERGY['upstream_upstream']
    # One upstream, one downstream
    elif (pos1_rel < 0) != (pos2_rel < 0):
        position_factor = POSITION_SYNERGY['upstream_downstream']
    # Both downstream
    else:
        position_factor = POSITION_SYNERGY['downstream_downstream']

    # Strand consideration
    # Same strand can have cooperative loading
    if guide1.strand == guide2.strand:
        strand_factor = 1.1
    else:
        strand_factor = 1.0

    # Coherence contribution
    # High coherence guides synergize better
    if guide1.coherence_R > 0 and guide2.coherence_R > 0:
        coherence_factor = 0.5 * (guide1.coherence_R + guide2.coherence_R) + 0.5
    else:
        coherence_factor = 1.0

    # Combined synergy
    synergy = spacing_factor * position_factor * strand_factor * coherence_factor

    # Clamp to reasonable range
    return float(np.clip(synergy, 0, 2.0))


def analyze_guide_pair(
    guide1: GuideCandidate,
    guide2: GuideCandidate,
    tss_position: int = 0,
    config: Optional[MultiGuideConfig] = None,
) -> GuidePair:
    """
    Analyze interaction between two guides.

    Args:
        guide1: First guide
        guide2: Second guide
        tss_position: TSS position
        config: MultiGuideConfig

    Returns:
        GuidePair analysis object
    """
    if config is None:
        config = MultiGuideConfig()

    spacing = calculate_spacing(guide1, guide2)
    has_clash = check_steric_clash(guide1, guide2, config)

    in_optimal = (
        config.optimal_spacing_min <= spacing <= config.optimal_spacing_max
    )

    synergy = predict_pairwise_synergy(guide1, guide2, tss_position, config)

    # Determine interaction type
    if has_clash:
        interaction_type = "clash"
    elif synergy > 1.1:
        interaction_type = "synergistic"
    elif synergy >= 0.9:
        interaction_type = "additive"
    else:
        interaction_type = "antagonistic"

    # Combined efficiency
    # For synergistic: E_combined > E1 + E2
    # For additive: E_combined ≈ E1 + E2 - E1*E2 (diminishing returns)
    e1 = guide1.efficiency if guide1.efficiency > 0 else guide1.individual_score
    e2 = guide2.efficiency if guide2.efficiency > 0 else guide2.individual_score

    if has_clash:
        combined = max(e1, e2) * 0.5  # Severe reduction
    else:
        # Multiplicative model with synergy
        combined = 1.0 - (1.0 - e1) * (1.0 - e2) * (2.0 - synergy)
        combined = min(1.0, combined)

    return GuidePair(
        guide1=guide1,
        guide2=guide2,
        spacing=spacing,
        has_steric_clash=has_clash,
        in_optimal_range=in_optimal,
        synergy_score=synergy,
        combined_efficiency=combined,
        interaction_type=interaction_type,
    )


def design_multiguide_set(
    candidates: List[GuideCandidate],
    n_guides: int = 2,
    tss_position: int = 0,
    config: Optional[MultiGuideConfig] = None,
) -> List[MultiGuideSet]:
    """
    Design optimal multi-guide sets from candidates.

    Evaluates all combinations and ranks by combined effectiveness.

    Args:
        candidates: List of guide candidates
        n_guides: Number of guides in each set
        tss_position: TSS position
        config: MultiGuideConfig

    Returns:
        List of MultiGuideSet objects, sorted by score

    Example:
        >>> candidates = [GuideCandidate(...) for guide in top_guides]
        >>> sets = design_multiguide_set(candidates, n_guides=2)
        >>> print(f"Best set synergy: {sets[0].synergy_factor:.2f}")
    """
    if config is None:
        config = MultiGuideConfig()

    if len(candidates) < n_guides:
        raise ValueError(
            f"Need at least {n_guides} candidates, got {len(candidates)}"
        )

    sets = []

    # Evaluate all combinations
    for combo in combinations(range(len(candidates)), n_guides):
        guides = [candidates[i] for i in combo]

        # Analyze all pairs
        pairs = []
        for i, j in combinations(range(n_guides), 2):
            pair = analyze_guide_pair(guides[i], guides[j], tss_position, config)
            pairs.append(pair)

        # Check for any clashes
        has_clashes = any(p.has_steric_clash for p in pairs)

        # Calculate spacing stats
        spacings = [p.spacing for p in pairs]
        min_spacing = min(spacings) if spacings else 0
        max_spacing = max(spacings) if spacings else 0

        # Average synergy
        synergy_factor = np.mean([p.synergy_score for p in pairs]) if pairs else 1.0

        # Combined score
        individual_avg = np.mean([g.individual_score for g in guides])

        # Ensemble coherence (geometric mean)
        coherences = [g.coherence_R for g in guides if g.coherence_R > 0]
        if coherences:
            ensemble_coherence = float(np.exp(np.mean(np.log(coherences))))
        else:
            ensemble_coherence = 0.5

        # Combined efficiency (accounts for synergy and diminishing returns)
        combined_score = _compute_set_score(
            guides, pairs, synergy_factor, ensemble_coherence, config
        )

        # GO/NO-GO
        if has_clashes:
            go_status = "NO-GO"
            recommendation = "CLASH: Guides too close, will interfere"
        elif ensemble_coherence < E_MINUS_2:
            go_status = "NO-GO"
            recommendation = "LOW COHERENCE: Unreliable prediction"
        elif synergy_factor < 0.8:
            go_status = "CAUTION"
            recommendation = "ANTAGONISTIC: Guides may interfere"
        elif synergy_factor > 1.1:
            go_status = "GO"
            recommendation = f"SYNERGISTIC: {synergy_factor:.1f}x synergy expected"
        else:
            go_status = "GO"
            recommendation = "ADDITIVE: Normal combined effect expected"

        guide_set = MultiGuideSet(
            guides=guides,
            pairs=pairs,
            combined_score=combined_score,
            synergy_factor=float(synergy_factor),
            ensemble_coherence=ensemble_coherence,
            has_clashes=has_clashes,
            min_spacing=min_spacing,
            max_spacing=max_spacing,
            recommendation=recommendation,
            go_no_go=go_status,
        )
        sets.append(guide_set)

    # Sort by combined score (descending)
    sets.sort(key=lambda s: s.combined_score, reverse=True)

    return sets


def _compute_set_score(
    guides: List[GuideCandidate],
    pairs: List[GuidePair],
    synergy_factor: float,
    ensemble_coherence: float,
    config: MultiGuideConfig,
) -> float:
    """Compute combined score for a guide set."""
    # Individual contribution
    individual_avg = np.mean([g.individual_score for g in guides])

    # Synergy contribution (normalized)
    synergy_norm = synergy_factor / 1.2  # 1.2 is approximate max

    # Combined
    score = (
        config.weight_individual * individual_avg +
        config.weight_synergy * synergy_norm +
        config.weight_coherence * ensemble_coherence
    )

    # Penalty for clashes
    clash_count = sum(1 for p in pairs if p.has_steric_clash)
    if clash_count > 0:
        score *= 0.3 ** clash_count

    return float(score)


def optimize_guide_spacing(
    anchor_guide: GuideCandidate,
    candidates: List[GuideCandidate],
    tss_position: int = 0,
    config: Optional[MultiGuideConfig] = None,
) -> List[Tuple[GuideCandidate, GuidePair]]:
    """
    Find optimal second guide given an anchor guide.

    Useful when one guide is already selected and you need
    to find the best partner.

    Args:
        anchor_guide: The fixed guide
        candidates: Potential partner guides
        tss_position: TSS position
        config: MultiGuideConfig

    Returns:
        List of (guide, pair_analysis) tuples, sorted by synergy
    """
    if config is None:
        config = MultiGuideConfig()

    results = []

    for candidate in candidates:
        # Skip self-comparison
        if candidate.position == anchor_guide.position:
            continue

        pair = analyze_guide_pair(anchor_guide, candidate, tss_position, config)

        if not pair.has_steric_clash:
            results.append((candidate, pair))

    # Sort by synergy score
    results.sort(key=lambda x: x[1].synergy_score, reverse=True)

    return results


def predict_combinatorial_effect(
    guides: List[GuideCandidate],
    individual_effects: List[float],
    tss_position: int = 0,
    config: Optional[MultiGuideConfig] = None,
) -> Dict[str, Any]:
    """
    Predict combined effect of multiple guides.

    Uses multiplicative model with synergy corrections.

    Args:
        guides: List of guides
        individual_effects: Individual fold-change or efficiency for each
        tss_position: TSS position
        config: MultiGuideConfig

    Returns:
        Prediction dictionary with combined effect

    Example:
        >>> result = predict_combinatorial_effect(
        ...     guides=[g1, g2, g3],
        ...     individual_effects=[2.0, 1.8, 1.5]  # Fold changes
        ... )
        >>> print(f"Combined fold change: {result['combined_effect']:.1f}x")
    """
    if config is None:
        config = MultiGuideConfig()

    if len(guides) != len(individual_effects):
        raise ValueError("Must have same number of guides and effects")

    n = len(guides)

    # Analyze all pairs
    pair_synergies = []
    for i in range(n):
        for j in range(i + 1, n):
            synergy = predict_pairwise_synergy(
                guides[i], guides[j], tss_position, config
            )
            pair_synergies.append(synergy)

    avg_synergy = np.mean(pair_synergies) if pair_synergies else 1.0

    # Combined effect model
    # For fold changes: multiplicative with diminishing returns
    if all(e >= 1.0 for e in individual_effects):
        # Fold changes (multiplicative)
        base_product = np.prod(individual_effects)

        # Apply diminishing returns (log-sum model)
        log_effects = np.log(individual_effects)
        log_combined = np.sum(log_effects) * (0.6 + 0.4 * avg_synergy)
        combined = np.exp(log_combined)

        # Cap at reasonable maximum
        combined = min(combined, base_product * avg_synergy)
    else:
        # Efficiency values (0-1)
        # P(all active) model with synergy
        combined = 1.0
        for effect in individual_effects:
            combined *= (1.0 - (1.0 - effect) * (2.0 - avg_synergy) / 2.0)

    return {
        'combined_effect': float(combined),
        'individual_effects': individual_effects,
        'average_synergy': float(avg_synergy),
        'n_guides': n,
        'synergy_type': 'synergistic' if avg_synergy > 1.1 else (
            'additive' if avg_synergy >= 0.9 else 'antagonistic'
        ),
    }


def validate_guide_set(
    guides: List[GuideCandidate],
    config: Optional[MultiGuideConfig] = None,
) -> Dict[str, Any]:
    """
    Validate a proposed guide set for compatibility.

    Checks:
    - Steric clashes
    - Spacing appropriateness
    - Coherence quality
    - Overall synergy potential

    Args:
        guides: List of guides to validate
        config: MultiGuideConfig

    Returns:
        Validation results dictionary
    """
    if config is None:
        config = MultiGuideConfig()

    n = len(guides)
    issues = []
    warnings = []

    # Check pairs
    clashes = 0
    min_spacing = float('inf')
    max_spacing = 0
    synergies = []

    for i in range(n):
        for j in range(i + 1, n):
            spacing = calculate_spacing(guides[i], guides[j])
            min_spacing = min(min_spacing, spacing)
            max_spacing = max(max_spacing, spacing)

            if check_steric_clash(guides[i], guides[j], config):
                clashes += 1
                issues.append(
                    f"Steric clash between guide {i+1} and {j+1} "
                    f"(spacing: {spacing}bp)"
                )

            synergy = predict_pairwise_synergy(guides[i], guides[j], 0, config)
            synergies.append(synergy)
            if synergy < 0.8:
                warnings.append(
                    f"Low synergy ({synergy:.2f}) between guide {i+1} and {j+1}"
                )

    # Coherence check
    low_coherence = [
        i for i, g in enumerate(guides)
        if g.coherence_R > 0 and g.coherence_R < E_MINUS_2
    ]
    if low_coherence:
        issues.append(
            f"Low coherence guides (NO-GO): {[i+1 for i in low_coherence]}"
        )

    # Overall assessment
    valid = len(issues) == 0
    avg_synergy = np.mean(synergies) if synergies else 1.0

    return {
        'valid': valid,
        'issues': issues,
        'warnings': warnings,
        'n_clashes': clashes,
        'min_spacing': min_spacing if min_spacing != float('inf') else 0,
        'max_spacing': max_spacing,
        'average_synergy': float(avg_synergy),
        'recommendation': (
            "VALID" if valid else "INVALID: " + "; ".join(issues)
        ),
    }


# Experimental synergy data for validation
KNOWN_SYNERGY_EXAMPLES = {
    # From Konermann et al. and other CRISPRa studies
    'high_synergy': {
        'spacing_range': (80, 150),
        'expected_synergy': 1.3,
        'description': 'Optimal spacing shows ~30% synergy',
    },
    'close_spacing': {
        'spacing_range': (20, 40),
        'expected_synergy': 0.6,
        'description': 'Very close guides interfere',
    },
    'far_spacing': {
        'spacing_range': (400, 600),
        'expected_synergy': 1.0,
        'description': 'Distant guides act independently',
    },
}
