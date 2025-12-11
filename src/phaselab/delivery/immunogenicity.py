"""
PhaseLab Immunogenicity: Predict immune responses to CRISPR components.

Implements:
- Cas9 T-cell epitope prediction
- Guide RNA innate immune stimulation
- TLR9 motif scoring
- Pre-existing immunity assessment
- Immunosilencing recommendations

References:
- Charlesworth et al. (2019): Pre-existing Cas9 immunity in humans
- Ferdosi et al. (2019): Engineered immunosilenced Cas9
- Kim et al. (2018): gRNA immunostimulation via TLRs
- Wagner et al. (2021): Cas9 epitope masking

IMPORTANT: Immune responses can reduce efficacy and cause safety issues.
This module helps predict and mitigate immunogenicity risks.

Version: 0.5.0
"""

import numpy as np
from typing import Dict, List, Optional, Tuple, Union, Any
from dataclasses import dataclass, field
import re
import logging

logger = logging.getLogger(__name__)


# TLR9 recognition motifs (CpG motifs that stimulate innate immunity)
# These can trigger immune responses when present in DNA
TLR9_STIMULATORY_MOTIFS = {
    'GACGTT': 1.0,   # Strong stimulatory (human-optimal CpG)
    'AACGTT': 0.9,   # Strong stimulatory
    'GACGTC': 0.8,   # Moderate stimulatory
    'GTCGTT': 0.85,  # Strong in mouse
    'AACGCC': 0.6,   # Weaker
    'GACGAT': 0.5,   # Weak
}

# Neutralizing/inhibitory motifs
TLR9_INHIBITORY_MOTIFS = {
    'GCGCGC': 0.8,   # Strong inhibitory
    'CCTGGA': 0.6,   # Moderate inhibitory
    'CCTGCA': 0.5,   # Weaker inhibitory
}

# Known Cas9 immunogenic epitopes (from published studies)
CAS9_IMMUNOGENIC_REGIONS = {
    'SpCas9': {
        # T-cell epitopes identified by Ferdosi et al.
        'epitope_1': {'position': (368, 376), 'sequence': 'NYTQKFSNM', 'HLA': 'A*02:01'},
        'epitope_2': {'position': (978, 986), 'sequence': 'GLFNEPVTR', 'HLA': 'A*02:01'},
        # B-cell epitopes
        'b_epitope_1': {'position': (1100, 1120), 'linear': True},
    },
    'SaCas9': {
        'epitope_1': {'position': (232, 240), 'HLA': 'A*02:01'},
        'epitope_2': {'position': (567, 575), 'HLA': 'A*02:01'},
    },
}

# Pre-existing immunity prevalence
PRE_EXISTING_IMMUNITY = {
    'SpCas9': {
        'antibody_prevalence': 0.58,  # ~58% have anti-SpCas9 antibodies
        'tcell_prevalence': 0.46,     # ~46% have T-cell responses
        'source': 'S. pyogenes commensal exposure',
    },
    'SaCas9': {
        'antibody_prevalence': 0.78,  # Higher due to S. aureus exposure
        'tcell_prevalence': 0.52,
        'source': 'S. aureus commensal exposure',
    },
    'AsCas12a': {
        'antibody_prevalence': 0.10,  # Lower for Acidaminococcus
        'tcell_prevalence': 0.08,
        'source': 'Rare environmental exposure',
    },
}


@dataclass
class ImmunogenicityScore:
    """Immunogenicity assessment result."""

    # Overall risk
    overall_risk: str  # LOW, MODERATE, HIGH
    risk_score: float  # 0-1

    # Component risks
    pre_existing_risk: float
    t_cell_risk: float
    b_cell_risk: float
    innate_risk: float

    # Specific findings
    tlr9_motif_count: int
    predicted_epitopes: int

    # Recommendations
    recommendations: List[str]

    # Cas variant used
    cas_variant: str = "SpCas9"


@dataclass
class ImmunogenicityConfig:
    """Configuration for immunogenicity analysis."""

    # Cas variant
    cas_variant: str = "SpCas9"

    # Risk weights
    weight_preexisting: float = 0.30
    weight_tcell: float = 0.25
    weight_bcell: float = 0.20
    weight_innate: float = 0.25

    # Thresholds
    high_risk_threshold: float = 0.6
    moderate_risk_threshold: float = 0.3

    # Analysis options
    check_tlr9: bool = True
    check_epitopes: bool = True


def count_cpg_motifs(sequence: str) -> int:
    """Count CpG dinucleotides in sequence."""
    sequence = sequence.upper()
    return len(re.findall(r'CG', sequence))


def score_tlr9_motifs(
    sequence: str,
    return_details: bool = False
) -> Union[float, Dict[str, Any]]:
    """
    Score sequence for TLR9 stimulatory potential.

    TLR9 recognizes unmethylated CpG motifs in DNA, triggering
    innate immune responses.

    Args:
        sequence: DNA sequence to analyze
        return_details: Return detailed motif analysis

    Returns:
        TLR9 stimulation score (0-1) or detailed dictionary

    Example:
        >>> score = score_tlr9_motifs("GACGTTCCGACGTT")
        >>> print(f"TLR9 score: {score:.2f}")
    """
    sequence = sequence.upper()

    stimulatory_score = 0.0
    inhibitory_score = 0.0
    found_motifs = []

    # Check stimulatory motifs
    for motif, weight in TLR9_STIMULATORY_MOTIFS.items():
        count = sequence.count(motif)
        if count > 0:
            stimulatory_score += count * weight
            found_motifs.append({'motif': motif, 'count': count, 'type': 'stimulatory'})

    # Check inhibitory motifs
    for motif, weight in TLR9_INHIBITORY_MOTIFS.items():
        count = sequence.count(motif)
        if count > 0:
            inhibitory_score += count * weight
            found_motifs.append({'motif': motif, 'count': count, 'type': 'inhibitory'})

    # Base CpG count
    cpg_count = count_cpg_motifs(sequence)
    base_cpg_score = min(1.0, cpg_count / 10)  # Normalize

    # Combined score
    net_score = stimulatory_score - 0.5 * inhibitory_score + 0.3 * base_cpg_score
    normalized_score = float(np.clip(net_score / 5.0, 0, 1))  # Scale to 0-1

    if return_details:
        return {
            'tlr9_score': normalized_score,
            'stimulatory_score': stimulatory_score,
            'inhibitory_score': inhibitory_score,
            'cpg_count': cpg_count,
            'found_motifs': found_motifs,
            'risk_level': 'HIGH' if normalized_score > 0.5 else (
                'MODERATE' if normalized_score > 0.2 else 'LOW'
            ),
        }

    return normalized_score


def predict_guide_immunogenicity(
    guide_sequence: str,
    synthesis_method: str = "chemical",
) -> Dict[str, Any]:
    """
    Predict immunogenicity of a guide RNA sequence.

    In vitro transcribed (IVT) gRNAs with 5'-triphosphate can
    trigger innate immune responses. Chemically synthesized
    gRNAs with 5'-OH are generally less immunogenic.

    Args:
        guide_sequence: Guide RNA sequence
        synthesis_method: "chemical" (5'-OH) or "ivt" (5'-triphosphate)

    Returns:
        Immunogenicity assessment dictionary
    """
    guide_upper = guide_sequence.upper().replace('U', 'T')

    result = {
        'sequence': guide_sequence,
        'synthesis_method': synthesis_method,
        'length': len(guide_sequence),
    }

    # TLR9 scoring (for DNA template)
    tlr9_result = score_tlr9_motifs(guide_upper, return_details=True)
    result['tlr9_analysis'] = tlr9_result

    # IVT produces 5'-triphosphate which triggers RIG-I
    if synthesis_method == "ivt":
        result['five_prime_risk'] = 'HIGH'
        result['five_prime_note'] = (
            "IVT gRNAs have 5'-triphosphate that activates RIG-I pathway. "
            "Consider phosphatase treatment or chemical synthesis."
        )
        base_innate_risk = 0.6
    else:
        result['five_prime_risk'] = 'LOW'
        result['five_prime_note'] = "Chemical synthesis produces 5'-OH, lower RIG-I activation."
        base_innate_risk = 0.2

    # GU-rich regions (TLR7/8 stimulation for RNA)
    gu_count = guide_upper.count('G') + guide_upper.count('U')
    gu_fraction = gu_count / len(guide_sequence) if guide_sequence else 0

    if gu_fraction > 0.5:
        result['gu_rich_risk'] = 'MODERATE'
        innate_adjustment = 0.1
    else:
        result['gu_rich_risk'] = 'LOW'
        innate_adjustment = 0.0

    # Combined innate risk
    innate_risk = base_innate_risk + 0.3 * tlr9_result['tlr9_score'] + innate_adjustment
    result['innate_immune_risk'] = min(1.0, innate_risk)

    # Overall assessment
    if result['innate_immune_risk'] > 0.5:
        result['overall_risk'] = 'HIGH'
        result['recommendations'] = [
            "Use chemically synthesized gRNA with 5'-OH",
            "Consider chemical modifications (2'-O-methyl, phosphorothioate)",
            "Avoid GU-rich sequences if possible",
        ]
    elif result['innate_immune_risk'] > 0.2:
        result['overall_risk'] = 'MODERATE'
        result['recommendations'] = [
            "Chemical synthesis recommended",
            "Monitor for inflammatory responses",
        ]
    else:
        result['overall_risk'] = 'LOW'
        result['recommendations'] = ["Standard precautions sufficient"]

    return result


def predict_cas9_immunogenicity(
    cas_variant: str = "SpCas9",
    delivery_method: str = "AAV",
    config: Optional[ImmunogenicityConfig] = None,
) -> ImmunogenicityScore:
    """
    Predict immunogenicity of Cas9 protein.

    Args:
        cas_variant: Cas9 variant (SpCas9, SaCas9, AsCas12a)
        delivery_method: Delivery method (AAV, RNP, mRNA, plasmid)
        config: ImmunogenicityConfig

    Returns:
        ImmunogenicityScore with assessment
    """
    if config is None:
        config = ImmunogenicityConfig()

    # Get pre-existing immunity data
    immunity_data = PRE_EXISTING_IMMUNITY.get(cas_variant, {
        'antibody_prevalence': 0.3,
        'tcell_prevalence': 0.2,
    })

    # Pre-existing risk
    pre_existing_risk = (
        immunity_data['antibody_prevalence'] * 0.6 +
        immunity_data['tcell_prevalence'] * 0.4
    )

    # T-cell epitope risk
    epitope_data = CAS9_IMMUNOGENIC_REGIONS.get(cas_variant, {})
    n_epitopes = len(epitope_data)
    t_cell_risk = min(1.0, n_epitopes * 0.3)  # Each epitope adds risk

    # B-cell risk (correlates with pre-existing antibodies)
    b_cell_risk = immunity_data['antibody_prevalence']

    # Delivery method affects risk
    delivery_factors = {
        'AAV': 1.0,      # Continuous expression
        'RNP': 0.5,      # Transient, lower risk
        'mRNA': 0.6,     # Transient
        'plasmid': 0.8,  # Some persistence
    }
    delivery_factor = delivery_factors.get(delivery_method, 1.0)

    # Innate immune risk (from delivery vector)
    if delivery_method == 'AAV':
        innate_risk = 0.4  # AAV capsid can trigger innate immunity
    elif delivery_method == 'RNP':
        innate_risk = 0.2  # Lower innate stimulation
    else:
        innate_risk = 0.3

    # Apply delivery factor
    t_cell_risk *= delivery_factor
    b_cell_risk *= delivery_factor

    # Combined risk score
    risk_score = (
        config.weight_preexisting * pre_existing_risk +
        config.weight_tcell * t_cell_risk +
        config.weight_bcell * b_cell_risk +
        config.weight_innate * innate_risk
    )

    # Classify risk
    if risk_score >= config.high_risk_threshold:
        overall_risk = "HIGH"
    elif risk_score >= config.moderate_risk_threshold:
        overall_risk = "MODERATE"
    else:
        overall_risk = "LOW"

    # Generate recommendations
    recommendations = []

    if pre_existing_risk > 0.5:
        recommendations.append(
            f"High pre-existing immunity to {cas_variant}. "
            "Consider: (1) Pre-screening patients, (2) Alternative Cas variant, "
            "(3) Immune modulation protocol."
        )

    if cas_variant == "SpCas9" and t_cell_risk > 0.3:
        recommendations.append(
            "Consider immunosilenced SpCas9 variants (Ferdosi et al.) "
            "with epitope modifications."
        )

    if delivery_method == "AAV":
        recommendations.append(
            "AAV delivery causes continuous Cas9 expression. "
            "Consider self-inactivating systems or RNP delivery."
        )

    if b_cell_risk > 0.5:
        recommendations.append(
            "High B-cell immunity risk. Pre-existing antibodies may "
            "neutralize therapeutic. Consider antibody screening."
        )

    if not recommendations:
        recommendations.append("Standard monitoring recommended.")

    return ImmunogenicityScore(
        overall_risk=overall_risk,
        risk_score=risk_score,
        pre_existing_risk=pre_existing_risk,
        t_cell_risk=t_cell_risk,
        b_cell_risk=b_cell_risk,
        innate_risk=innate_risk,
        tlr9_motif_count=0,  # Not applicable for protein
        predicted_epitopes=n_epitopes,
        recommendations=recommendations,
        cas_variant=cas_variant,
    )


def assess_immunogenic_risk(
    guide_sequence: str,
    cas_variant: str = "SpCas9",
    delivery_method: str = "AAV",
    synthesis_method: str = "chemical",
    config: Optional[ImmunogenicityConfig] = None,
) -> Dict[str, Any]:
    """
    Comprehensive immunogenicity assessment for CRISPR therapy.

    Combines:
    - Cas protein immunogenicity
    - Guide RNA immunogenicity
    - Delivery vector considerations
    - Pre-existing immunity

    Args:
        guide_sequence: Guide RNA sequence
        cas_variant: Cas9 variant
        delivery_method: Delivery method
        synthesis_method: gRNA synthesis method
        config: ImmunogenicityConfig

    Returns:
        Comprehensive immunogenicity assessment

    Example:
        >>> result = assess_immunogenic_risk(
        ...     guide_sequence="GCGACTGCTACATAGCCAGG",
        ...     cas_variant="SpCas9",
        ...     delivery_method="AAV"
        ... )
        >>> print(f"Overall risk: {result['overall_risk']}")
    """
    if config is None:
        config = ImmunogenicityConfig()

    result = {
        'cas_variant': cas_variant,
        'delivery_method': delivery_method,
        'synthesis_method': synthesis_method,
    }

    # Cas protein immunogenicity
    cas_immuno = predict_cas9_immunogenicity(cas_variant, delivery_method, config)
    result['cas_immunogenicity'] = {
        'risk': cas_immuno.overall_risk,
        'score': cas_immuno.risk_score,
        'pre_existing': cas_immuno.pre_existing_risk,
        'epitopes': cas_immuno.predicted_epitopes,
    }

    # Guide RNA immunogenicity
    guide_immuno = predict_guide_immunogenicity(guide_sequence, synthesis_method)
    result['guide_immunogenicity'] = {
        'risk': guide_immuno['overall_risk'],
        'score': guide_immuno['innate_immune_risk'],
        'tlr9_score': guide_immuno['tlr9_analysis']['tlr9_score'],
    }

    # Combined risk
    combined_score = (
        0.6 * cas_immuno.risk_score +
        0.4 * guide_immuno['innate_immune_risk']
    )

    if combined_score >= 0.6:
        overall_risk = "HIGH"
    elif combined_score >= 0.3:
        overall_risk = "MODERATE"
    else:
        overall_risk = "LOW"

    result['combined_risk'] = {
        'overall_risk': overall_risk,
        'score': combined_score,
    }

    # Aggregate recommendations
    all_recs = cas_immuno.recommendations + guide_immuno['recommendations']
    result['recommendations'] = list(set(all_recs))  # Deduplicate

    # GO/NO-GO
    if overall_risk == "HIGH":
        result['status'] = "CAUTION"
        result['status_note'] = (
            "High immunogenicity risk. Review recommendations and consider "
            "alternative approaches or immune monitoring protocols."
        )
    else:
        result['status'] = "PROCEED"
        result['status_note'] = "Acceptable immunogenicity risk with standard precautions."

    return result


def suggest_immunosilencing_strategy(
    cas_variant: str = "SpCas9",
    risk_assessment: Optional[ImmunogenicityScore] = None,
) -> List[Dict[str, str]]:
    """
    Suggest immunosilencing strategies for reducing immunogenicity.

    Args:
        cas_variant: Cas9 variant
        risk_assessment: Optional prior risk assessment

    Returns:
        List of strategy recommendations
    """
    strategies = []

    # Universal strategies
    strategies.append({
        'strategy': 'RNP delivery',
        'description': 'Use ribonucleoprotein delivery for transient expression',
        'effectiveness': 'HIGH',
        'tradeoff': 'Requires electroporation or lipofection, limited to ex vivo',
    })

    strategies.append({
        'strategy': 'Chemical gRNA modifications',
        'description': "Add 2'-O-methyl and phosphorothioate modifications to gRNA",
        'effectiveness': 'MODERATE',
        'tradeoff': 'Increased cost, may slightly reduce activity',
    })

    if cas_variant == "SpCas9":
        strategies.append({
            'strategy': 'Immunosilenced SpCas9 variants',
            'description': 'Use epitope-masked variants (Ferdosi et al., Wagner et al.)',
            'effectiveness': 'HIGH',
            'tradeoff': 'May have slightly reduced activity',
        })

        strategies.append({
            'strategy': 'Switch to AsCas12a',
            'description': 'Use AsCas12a which has lower pre-existing immunity',
            'effectiveness': 'HIGH',
            'tradeoff': 'Different PAM requirement (TTTV), different guide design',
        })

    strategies.append({
        'strategy': 'Immune modulation',
        'description': 'Co-administer tolerogenic agents or use Treg induction',
        'effectiveness': 'MODERATE',
        'tradeoff': 'Complex protocol, potential side effects',
    })

    strategies.append({
        'strategy': 'Patient screening',
        'description': 'Pre-screen for anti-Cas9 antibodies and T-cells',
        'effectiveness': 'MODERATE',
        'tradeoff': 'May exclude significant patient population',
    })

    return strategies


# Known immunogenic peptides for validation
KNOWN_IMMUNOGENIC_PEPTIDES = {
    # SpCas9 epitopes from literature
    'SpCas9_epitope_1': 'NYTQKFSNM',
    'SpCas9_epitope_2': 'GLFNEPVTR',
    # These should trigger T-cell responses in HLA-A*02:01 individuals
}
