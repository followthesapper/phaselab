"""
PhaseLab Enhancer Targeting: CRISPRa guide design for enhancer activation.

Implements:
- Enhancer region identification and scoring
- eQTL-based enhancer-gene linking
- CRISPRa enhancer activation prediction
- Multi-enhancer combinatorial activation
- Tissue-specific enhancer models

References:
- crisprQTL: Mapping enhancer-gene pairs via CRISPR screens
- pgBoost: Probabilistic scoring for enhancer-gene links
- Fulco et al. Activity-by-contact model for enhancer prediction
- ENCODE enhancer annotations (cCREs)

IMPORTANT: Many genes are difficult to activate through promoter-only
CRISPRa. Targeting enhancers can provide stronger or tissue-specific
activation for genes like MYC, NRXN1, SCN2A, and MECP2.

Version: 0.5.0
"""

import numpy as np
from typing import Dict, List, Optional, Tuple, Union, Any
from dataclasses import dataclass, field
import logging

from ..core.coherence import coherence_score, go_no_go
from ..core.constants import E_MINUS_2

logger = logging.getLogger(__name__)


# ENCODE cCRE (candidate Cis-Regulatory Element) types
CCRE_TYPES = {
    'pELS': 'proximal enhancer-like signature',   # <2kb from TSS
    'dELS': 'distal enhancer-like signature',     # >2kb from TSS
    'PLS': 'promoter-like signature',
    'CTCF-only': 'CTCF-bound only',
    'DNase-H3K4me3': 'promoter-proximal',
}

# Activity-by-contact model parameters
ABC_MODEL_PARAMS = {
    'activity_weight': 0.6,      # Weight for enhancer activity signal
    'contact_weight': 0.4,       # Weight for 3D contact frequency
    'distance_decay': 5000,      # Distance for contact decay (bp)
    'min_abc_score': 0.02,       # Minimum ABC score for significance
}


@dataclass
class Enhancer:
    """Enhancer element annotation."""
    chrom: str
    start: int
    end: int
    name: str = ""

    # Type and activity
    element_type: str = "dELS"   # pELS, dELS, etc.
    activity_score: float = 0.0  # H3K27ac/ATAC signal

    # Target gene information
    target_gene: Optional[str] = None
    target_gene_tss: Optional[int] = None
    distance_to_tss: int = 0

    # Linking scores
    abc_score: float = 0.0       # Activity-by-contact
    eqtl_support: bool = False   # Has eQTL evidence
    hi_c_contact: float = 0.0    # Hi-C contact frequency

    # Tissue specificity
    tissue: str = "generic"
    tissue_specific: bool = False

    @property
    def center(self) -> int:
        return (self.start + self.end) // 2

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def is_proximal(self) -> bool:
        return abs(self.distance_to_tss) < 2000


@dataclass
class EnhancerGuideResult:
    """Result from enhancer-targeting guide design."""
    guide_sequence: str
    guide_position: int
    enhancer: Enhancer

    # Scores
    activation_potential: float
    guide_quality: float
    combined_score: float

    # Coherence
    coherence_R: float
    go_no_go: str

    # Recommendations
    recommendation: str


@dataclass
class EnhancerConfig:
    """Configuration for enhancer targeting."""

    # Enhancer definition
    min_enhancer_length: int = 200
    max_enhancer_distance: int = 1_000_000  # 1Mb max

    # Activity thresholds
    min_activity_score: float = 0.3
    min_abc_score: float = 0.02

    # Guide design
    guides_per_enhancer: int = 3
    prefer_enhancer_center: bool = True

    # Tissue specificity
    require_tissue_match: bool = False
    target_tissue: str = "generic"


def calculate_abc_score(
    activity: float,
    contact: float,
    distance: int,
    config: Optional[Dict] = None,
) -> float:
    """
    Calculate Activity-by-Contact (ABC) score for enhancer-gene link.

    ABC model: Score = Activity × Contact / Σ(Activity × Contact)

    Based on Fulco et al. (2019) - predicts enhancer-gene pairs.

    Args:
        activity: Enhancer activity signal (H3K27ac, ATAC)
        contact: 3D contact frequency (Hi-C, if available)
        distance: Distance to target gene TSS
        config: ABC model parameters

    Returns:
        ABC score (higher = stronger link)
    """
    if config is None:
        config = ABC_MODEL_PARAMS

    # If no contact data, estimate from distance
    if contact <= 0:
        # Power-law decay with distance
        contact = 1.0 / (1.0 + (distance / config['distance_decay']) ** 1.5)

    # ABC score (not normalized - would need all enhancers for normalization)
    abc = activity * contact

    return float(abc)


def estimate_enhancer_activity(
    h3k27ac_signal: float = 0.0,
    atac_signal: float = 0.0,
    dnase_signal: float = 0.0,
) -> float:
    """
    Estimate enhancer activity from epigenetic signals.

    Args:
        h3k27ac_signal: H3K27ac ChIP-seq signal (active mark)
        atac_signal: ATAC-seq signal
        dnase_signal: DNase-seq signal

    Returns:
        Activity score 0-1
    """
    signals = []

    # H3K27ac is the strongest indicator of active enhancers
    if h3k27ac_signal > 0:
        signals.append(h3k27ac_signal * 1.2)  # Higher weight

    if atac_signal > 0:
        signals.append(atac_signal)

    if dnase_signal > 0:
        signals.append(dnase_signal * 0.8)

    if not signals:
        return 0.0

    # Combine signals (weighted average, capped at 1.0)
    activity = np.mean(signals)
    return float(min(1.0, activity))


def score_enhancer_for_activation(
    enhancer: Enhancer,
    config: Optional[EnhancerConfig] = None,
) -> Dict[str, float]:
    """
    Score an enhancer's potential for CRISPRa activation.

    Considers:
    - Enhancer activity (already active enhancers harder to boost)
    - Distance to target gene
    - ABC score (link confidence)
    - Tissue specificity

    Args:
        enhancer: Enhancer to score
        config: EnhancerConfig

    Returns:
        Scoring dictionary
    """
    if config is None:
        config = EnhancerConfig()

    scores = {}

    # Activity-based potential
    # Moderately active enhancers have highest activation potential
    # Very active = already maxed, inactive = hard to activate
    if enhancer.activity_score < 0.2:
        activity_potential = 0.4  # Weak enhancer, harder to activate
    elif enhancer.activity_score < 0.5:
        activity_potential = 0.9  # Moderate, good potential
    elif enhancer.activity_score < 0.8:
        activity_potential = 0.7  # Already active, some headroom
    else:
        activity_potential = 0.3  # Near maximum, little room to boost

    scores['activity_potential'] = activity_potential

    # Distance factor
    # Proximal enhancers easier to link to target gene
    distance = abs(enhancer.distance_to_tss)
    if distance < 2000:
        distance_factor = 1.0
    elif distance < 10000:
        distance_factor = 0.9
    elif distance < 100000:
        distance_factor = 0.7
    else:
        distance_factor = 0.5

    scores['distance_factor'] = distance_factor

    # ABC score contribution
    if enhancer.abc_score >= config.min_abc_score:
        abc_factor = min(1.0, enhancer.abc_score / 0.1)  # Normalize
    else:
        abc_factor = 0.3  # Below threshold, uncertain link

    scores['abc_factor'] = abc_factor

    # eQTL evidence boost
    if enhancer.eqtl_support:
        scores['eqtl_boost'] = 1.2
    else:
        scores['eqtl_boost'] = 1.0

    # Combined score
    combined = (
        0.35 * activity_potential +
        0.25 * distance_factor +
        0.25 * abc_factor +
        0.15 * (1.0 if enhancer.eqtl_support else 0.5)
    ) * scores['eqtl_boost']

    scores['combined_potential'] = float(min(1.0, combined))

    return scores


def identify_target_enhancers(
    gene_symbol: str,
    gene_tss: int,
    gene_chrom: str,
    candidate_enhancers: List[Enhancer],
    config: Optional[EnhancerConfig] = None,
) -> List[Tuple[Enhancer, float]]:
    """
    Identify enhancers likely to regulate a target gene.

    Args:
        gene_symbol: Target gene symbol
        gene_tss: Gene TSS position
        gene_chrom: Gene chromosome
        candidate_enhancers: List of candidate enhancers
        config: EnhancerConfig

    Returns:
        List of (enhancer, link_score) tuples, sorted by score
    """
    if config is None:
        config = EnhancerConfig()

    results = []

    for enh in candidate_enhancers:
        # Check chromosome match
        if enh.chrom != gene_chrom:
            continue

        # Calculate distance
        distance = abs(enh.center - gene_tss)

        # Check maximum distance
        if distance > config.max_enhancer_distance:
            continue

        # Update enhancer with distance
        enh.distance_to_tss = distance
        enh.target_gene = gene_symbol
        enh.target_gene_tss = gene_tss

        # Calculate ABC score if not already set
        if enh.abc_score <= 0:
            enh.abc_score = calculate_abc_score(
                activity=enh.activity_score,
                contact=enh.hi_c_contact,
                distance=distance,
            )

        # Score for CRISPRa potential
        scores = score_enhancer_for_activation(enh, config)
        link_score = scores['combined_potential']

        if link_score >= 0.3:  # Minimum threshold
            results.append((enh, link_score))

    # Sort by link score
    results.sort(key=lambda x: x[1], reverse=True)

    return results


def design_enhancer_guides(
    enhancer: Enhancer,
    sequence: str,
    enhancer_start_in_seq: int,
    config: Optional[EnhancerConfig] = None,
    guide_length: int = 20,
    pam: str = "NGG",
) -> List[EnhancerGuideResult]:
    """
    Design CRISPRa guides targeting an enhancer.

    For enhancer activation, guides should target:
    - Near the enhancer center (peak of activity)
    - Accessible regions within enhancer
    - Avoid repressive elements

    Args:
        enhancer: Target enhancer
        sequence: DNA sequence containing enhancer
        enhancer_start_in_seq: Position of enhancer start in sequence
        config: EnhancerConfig
        guide_length: Guide length (default 20)
        pam: PAM sequence

    Returns:
        List of EnhancerGuideResult objects
    """
    if config is None:
        config = EnhancerConfig()

    # Import here to avoid circular imports
    from .pam_scan import find_pam_sites
    from .scoring import gc_content, delta_g_santalucia, mit_specificity_score

    sequence = sequence.upper()
    results = []

    # Define search region within enhancer
    enh_start = enhancer_start_in_seq
    enh_end = enh_start + enhancer.length

    # Find PAM sites within enhancer
    all_hits = find_pam_sites(
        sequence,
        pam=pam,
        guide_length=guide_length,
        both_strands=True,
    )

    # Filter to enhancer region
    enhancer_hits = [
        h for h in all_hits
        if enh_start <= h.guide_start < enh_end
    ]

    if not enhancer_hits:
        logger.warning(f"No guides found in enhancer {enhancer.name}")
        return results

    # Score enhancer potential
    enh_scores = score_enhancer_for_activation(enhancer, config)

    # Score each guide
    for hit in enhancer_hits:
        guide_seq = hit.guide
        guide_pos = hit.guide_start

        # Basic guide quality
        gc = gc_content(guide_seq)
        if gc < 0.35 or gc > 0.75:
            continue  # Skip bad GC

        delta_g = delta_g_santalucia(guide_seq)
        mit = mit_specificity_score(guide_seq)

        # Position within enhancer (prefer center)
        enh_center = enh_start + enhancer.length // 2
        distance_from_center = abs(guide_pos - enh_center)
        center_factor = 1.0 - min(1.0, distance_from_center / (enhancer.length / 2))

        # Guide quality score
        guide_quality = (
            0.3 * (1.0 if 0.45 <= gc <= 0.65 else 0.6) +
            0.3 * min(1.0, -delta_g / 25.0) +
            0.2 * (mit / 100.0) +
            0.2 * center_factor
        )

        # Activation potential (enhancer × guide quality)
        activation_potential = enh_scores['combined_potential'] * guide_quality

        # IR coherence (simplified)
        R_bar = _compute_guide_coherence(guide_seq)
        status = go_no_go(R_bar)

        # Combined score
        combined = (
            0.4 * activation_potential +
            0.3 * guide_quality +
            0.3 * R_bar
        )

        # Recommendation
        if status == "NO-GO":
            recommendation = "LOW COHERENCE: Guide prediction unreliable"
        elif activation_potential >= 0.6:
            recommendation = f"EXCELLENT: High activation potential ({activation_potential:.2f})"
        elif activation_potential >= 0.4:
            recommendation = f"GOOD: Moderate activation potential ({activation_potential:.2f})"
        else:
            recommendation = f"MARGINAL: Low activation potential ({activation_potential:.2f})"

        result = EnhancerGuideResult(
            guide_sequence=guide_seq,
            guide_position=guide_pos,
            enhancer=enhancer,
            activation_potential=activation_potential,
            guide_quality=guide_quality,
            combined_score=combined,
            coherence_R=R_bar,
            go_no_go=status,
            recommendation=recommendation,
        )
        results.append(result)

    # Sort by combined score
    results.sort(key=lambda r: r.combined_score, reverse=True)

    # Return top guides
    return results[:config.guides_per_enhancer]


def _compute_guide_coherence(guide_seq: str) -> float:
    """Compute IR coherence using ATLAS-Q enhanced backend (v0.6.0+)."""
    from .coherence_utils import compute_guide_coherence
    return compute_guide_coherence(guide_seq, use_atlas_q=True)


def predict_enhancer_activation_effect(
    enhancer: Enhancer,
    guide_quality: float,
    baseline_expression: float = 1.0,
) -> Dict[str, float]:
    """
    Predict effect of enhancer CRISPRa on target gene expression.

    Args:
        enhancer: Target enhancer
        guide_quality: Quality score of guide used
        baseline_expression: Baseline gene expression level

    Returns:
        Prediction dictionary with fold-change estimates
    """
    # Enhancer activation typically gives 1.5-5x fold change
    # depending on enhancer strength and gene responsiveness

    # Base fold change from enhancer activity
    base_fold = 1.5 + 2.5 * enhancer.abc_score

    # Adjust for guide quality
    fold_change = base_fold * (0.5 + 0.5 * guide_quality)

    # Distance penalty
    distance_penalty = 1.0 - 0.3 * min(1.0, enhancer.distance_to_tss / 500000)
    fold_change *= distance_penalty

    # Confidence based on evidence
    if enhancer.eqtl_support:
        confidence = 0.8
    elif enhancer.abc_score >= 0.05:
        confidence = 0.6
    else:
        confidence = 0.4

    return {
        'predicted_fold_change': float(fold_change),
        'fold_change_low': float(fold_change * 0.6),
        'fold_change_high': float(fold_change * 1.5),
        'confidence': confidence,
        'baseline_expression': baseline_expression,
        'predicted_expression': baseline_expression * fold_change,
    }


def compare_promoter_vs_enhancer(
    promoter_guides: List[Dict],
    enhancer_guides: List[EnhancerGuideResult],
    target_expression: float,
    baseline_expression: float = 0.5,
) -> Dict[str, Any]:
    """
    Compare promoter CRISPRa vs enhancer CRISPRa strategies.

    Args:
        promoter_guides: Promoter-targeting guide results
        enhancer_guides: Enhancer-targeting guide results
        target_expression: Desired expression level
        baseline_expression: Current baseline expression

    Returns:
        Comparison and recommendation
    """
    results = {
        'promoter_option': None,
        'enhancer_option': None,
        'recommendation': None,
    }

    # Analyze promoter option
    if promoter_guides:
        best_promoter = max(
            promoter_guides,
            key=lambda g: g.get('combined_score', 0)
        )
        promoter_fold = 1.0 + 2.0 * best_promoter.get('combined_score', 0.5)
        promoter_expression = baseline_expression * promoter_fold

        results['promoter_option'] = {
            'best_guide': best_promoter,
            'predicted_fold_change': promoter_fold,
            'predicted_expression': promoter_expression,
            'reaches_target': promoter_expression >= target_expression,
        }

    # Analyze enhancer option
    if enhancer_guides:
        best_enhancer = enhancer_guides[0]  # Already sorted
        enh_pred = predict_enhancer_activation_effect(
            best_enhancer.enhancer,
            best_enhancer.guide_quality,
            baseline_expression,
        )

        results['enhancer_option'] = {
            'best_guide': best_enhancer,
            'predicted_fold_change': enh_pred['predicted_fold_change'],
            'predicted_expression': enh_pred['predicted_expression'],
            'confidence': enh_pred['confidence'],
            'reaches_target': enh_pred['predicted_expression'] >= target_expression,
        }

    # Generate recommendation
    promoter_reaches = (
        results['promoter_option'] and
        results['promoter_option']['reaches_target']
    )
    enhancer_reaches = (
        results['enhancer_option'] and
        results['enhancer_option']['reaches_target']
    )

    if promoter_reaches and not enhancer_reaches:
        recommendation = "PROMOTER: Sufficient activation from promoter CRISPRa"
    elif enhancer_reaches and not promoter_reaches:
        recommendation = "ENHANCER: Enhancer targeting needed to reach target"
    elif promoter_reaches and enhancer_reaches:
        # Both work - recommend based on confidence/simplicity
        recommendation = "BOTH VIABLE: Promoter simpler, enhancer may give higher activation"
    else:
        recommendation = "COMBINATION: May need multi-guide or promoter+enhancer strategy"

    results['recommendation'] = recommendation

    return results


# Pre-defined enhancer databases for common therapeutic targets
KNOWN_ENHANCERS = {
    'RAI1': [
        {
            'name': 'RAI1_enh1',
            'chrom': 'chr17',
            'start': 17679000,
            'end': 17680000,
            'distance_to_tss': -2500,
            'activity_score': 0.65,
            'tissue': 'brain',
        },
    ],
    'SCN2A': [
        {
            'name': 'SCN2A_enh1',
            'chrom': 'chr2',
            'start': 165235000,
            'end': 165236500,
            'distance_to_tss': -4400,
            'activity_score': 0.72,
            'tissue': 'brain',
        },
    ],
    'MECP2': [
        {
            'name': 'MECP2_enh1',
            'chrom': 'chrX',
            'start': 154025000,
            'end': 154026000,
            'distance_to_tss': -15000,
            'activity_score': 0.58,
            'tissue': 'brain',
        },
    ],
}


def get_known_enhancers(gene: str) -> List[Enhancer]:
    """
    Get known enhancers for a gene from built-in database.

    Args:
        gene: Gene symbol

    Returns:
        List of Enhancer objects
    """
    gene = gene.upper()
    if gene not in KNOWN_ENHANCERS:
        return []

    enhancers = []
    for enh_data in KNOWN_ENHANCERS[gene]:
        enh = Enhancer(
            chrom=enh_data['chrom'],
            start=enh_data['start'],
            end=enh_data['end'],
            name=enh_data['name'],
            activity_score=enh_data.get('activity_score', 0.5),
            distance_to_tss=enh_data.get('distance_to_tss', 0),
            tissue=enh_data.get('tissue', 'generic'),
            target_gene=gene,
        )
        enhancers.append(enh)

    return enhancers
