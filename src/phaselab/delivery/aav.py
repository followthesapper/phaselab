"""
PhaseLab AAV Delivery: AAV serotype selection and tropism modeling.

Implements:
- AAV serotype profiles and tissue tropism
- Packaging size constraints
- Blood-brain barrier penetration scoring
- Liver de-targeting strategies
- Brain-specific serotype selection

References:
- AAV-PHP.eB: Enhanced CNS tropism (Chan et al.)
- AAV9: Broad tropism, crosses BBB
- Capsid engineering: BRC06 and other novel variants
- MDPI review: Natural and engineered AAV serotypes (2024)

IMPORTANT: Serotype selection is critical for gene therapy success.
Wrong serotype = poor tissue targeting, immunogenicity, or inefficacy.

Version: 0.5.0
"""

import numpy as np
from typing import Dict, List, Optional, Tuple, Union, Any
from dataclasses import dataclass, field
import logging

logger = logging.getLogger(__name__)


# AAV packaging limits
AAV_PACKAGING_LIMITS = {
    'standard': 4700,        # Standard AAV packaging limit (bp)
    'oversized': 5000,       # Functional but reduced efficiency
    'maximum': 5200,         # Absolute maximum, poor efficiency
    'dual_vector': 9400,     # Using split/dual vector approach
}

# Typical construct sizes
CONSTRUCT_SIZES = {
    'spcas9': 4200,          # SpCas9 (too large for single AAV)
    'sacas9': 3200,          # SaCas9 (fits in AAV)
    'dcas9_krab': 4500,      # dCas9-KRAB (borderline)
    'dcas9_vpr': 5100,       # dCas9-VPR (oversized)
    'guide_only': 500,       # gRNA expression cassette
    'crispr_mini': 3800,     # Compact Cas variants
}


@dataclass
class AAVSerotype:
    """AAV serotype profile with tropism data."""

    name: str
    full_name: str

    # Tissue tropism scores (0-1)
    tropism: Dict[str, float]

    # BBB penetration
    crosses_bbb: bool
    bbb_efficiency: float  # 0-1

    # Liver tropism (often want to minimize)
    liver_tropism: float

    # Immunogenicity
    immunogenicity: str  # 'low', 'moderate', 'high'
    pre_existing_immunity: float  # Fraction of population with antibodies

    # Species considerations
    human_validated: bool
    mouse_specific: bool = False  # Some like PHP.eB require specific receptors

    # Packaging
    packaging_capacity: int = 4700

    # Clinical status
    clinical_stage: str = "preclinical"  # preclinical, clinical, approved

    # Notes
    notes: str = ""


# Comprehensive serotype profiles based on literature
SEROTYPE_PROFILES = {
    'AAV1': AAVSerotype(
        name='AAV1',
        full_name='Adeno-associated virus serotype 1',
        tropism={
            'muscle': 0.85,
            'heart': 0.80,
            'brain': 0.40,
            'liver': 0.50,
            'lung': 0.60,
            'retina': 0.30,
        },
        crosses_bbb=False,
        bbb_efficiency=0.1,
        liver_tropism=0.50,
        immunogenicity='moderate',
        pre_existing_immunity=0.35,
        human_validated=True,
        clinical_stage='approved',  # Luxturna uses AAV2, Zolgensma uses AAV9
        notes='Good for muscle, limited CNS penetration',
    ),
    'AAV2': AAVSerotype(
        name='AAV2',
        full_name='Adeno-associated virus serotype 2',
        tropism={
            'muscle': 0.50,
            'liver': 0.60,
            'brain': 0.30,
            'retina': 0.90,
            'lung': 0.40,
            'heart': 0.40,
        },
        crosses_bbb=False,
        bbb_efficiency=0.05,
        liver_tropism=0.60,
        immunogenicity='moderate',
        pre_existing_immunity=0.60,  # Most common pre-existing immunity
        human_validated=True,
        clinical_stage='approved',
        notes='Classic serotype, good for retina, high pre-existing immunity',
    ),
    'AAV5': AAVSerotype(
        name='AAV5',
        full_name='Adeno-associated virus serotype 5',
        tropism={
            'lung': 0.85,
            'brain': 0.45,
            'liver': 0.40,
            'retina': 0.70,
            'muscle': 0.30,
        },
        crosses_bbb=False,
        bbb_efficiency=0.15,
        liver_tropism=0.40,
        immunogenicity='low',
        pre_existing_immunity=0.25,
        human_validated=True,
        clinical_stage='clinical',
        notes='Lower immunogenicity, good for lung',
    ),
    'AAV8': AAVSerotype(
        name='AAV8',
        full_name='Adeno-associated virus serotype 8',
        tropism={
            'liver': 0.95,
            'muscle': 0.70,
            'heart': 0.75,
            'brain': 0.30,
            'pancreas': 0.60,
        },
        crosses_bbb=False,
        bbb_efficiency=0.1,
        liver_tropism=0.95,
        immunogenicity='low',
        pre_existing_immunity=0.20,
        human_validated=True,
        clinical_stage='clinical',
        notes='Excellent liver transduction, low immunogenicity',
    ),
    'AAV9': AAVSerotype(
        name='AAV9',
        full_name='Adeno-associated virus serotype 9',
        tropism={
            'brain': 0.70,
            'spinal_cord': 0.80,
            'heart': 0.85,
            'muscle': 0.80,
            'liver': 0.75,
            'lung': 0.60,
        },
        crosses_bbb=True,
        bbb_efficiency=0.65,
        liver_tropism=0.75,
        immunogenicity='moderate',
        pre_existing_immunity=0.45,
        human_validated=True,
        clinical_stage='approved',  # Zolgensma
        notes='Crosses BBB, broad tropism, FDA approved for SMA',
    ),
    'AAV-PHP.eB': AAVSerotype(
        name='AAV-PHP.eB',
        full_name='PHP.eB engineered AAV variant',
        tropism={
            'brain': 0.95,
            'spinal_cord': 0.90,
            'liver': 0.30,
            'muscle': 0.40,
            'heart': 0.50,
        },
        crosses_bbb=True,
        bbb_efficiency=0.90,
        liver_tropism=0.30,
        immunogenicity='moderate',
        pre_existing_immunity=0.15,  # Novel variant, less immunity
        human_validated=False,  # Mouse-specific receptor (LY6A)
        mouse_specific=True,
        clinical_stage='preclinical',
        notes='Excellent CNS in C57BL/6 mice, requires LY6A receptor (not in humans)',
    ),
    'AAV-DJ': AAVSerotype(
        name='AAV-DJ',
        full_name='AAV-DJ shuffled variant',
        tropism={
            'liver': 0.85,
            'muscle': 0.65,
            'brain': 0.35,
            'lung': 0.50,
        },
        crosses_bbb=False,
        bbb_efficiency=0.1,
        liver_tropism=0.85,
        immunogenicity='low',
        pre_existing_immunity=0.20,
        human_validated=True,
        clinical_stage='clinical',
        notes='High efficiency in vitro, good for liver',
    ),
    'AAVrh10': AAVSerotype(
        name='AAVrh10',
        full_name='Rhesus macaque AAV serotype 10',
        tropism={
            'brain': 0.75,
            'spinal_cord': 0.80,
            'liver': 0.55,
            'muscle': 0.70,
            'lung': 0.60,
        },
        crosses_bbb=True,
        bbb_efficiency=0.55,
        liver_tropism=0.55,
        immunogenicity='low',
        pre_existing_immunity=0.10,  # Low in humans
        human_validated=True,
        clinical_stage='clinical',
        notes='Good CNS penetration, lower pre-existing immunity than AAV9',
    ),
}


@dataclass
class AAVConfig:
    """Configuration for AAV selection."""

    # Target tissue
    target_tissue: str = "brain"

    # Delivery route
    delivery_route: str = "IV"  # IV, ICV, intrathecal, intramuscular, etc.

    # Payload size
    payload_size: int = 4000  # Total construct size in bp

    # Preferences
    minimize_liver: bool = True
    require_bbb_crossing: bool = False
    avoid_high_immunity: bool = True
    require_human_validated: bool = True

    # Thresholds
    min_tropism_score: float = 0.5
    max_pre_existing_immunity: float = 0.5


def score_serotype_tropism(
    serotype: AAVSerotype,
    target_tissue: str,
    config: Optional[AAVConfig] = None,
) -> Dict[str, float]:
    """
    Score a serotype for targeting a specific tissue.

    Args:
        serotype: AAVSerotype to evaluate
        target_tissue: Target tissue name
        config: AAVConfig

    Returns:
        Scoring dictionary
    """
    if config is None:
        config = AAVConfig()

    scores = {}

    # Primary tropism
    target_tissue_lower = target_tissue.lower()
    tropism = serotype.tropism.get(target_tissue_lower, 0.0)
    scores['target_tropism'] = tropism

    # BBB consideration for CNS targets
    if target_tissue_lower in ['brain', 'spinal_cord', 'cns']:
        if config.require_bbb_crossing and not serotype.crosses_bbb:
            scores['bbb_penalty'] = 0.0
        else:
            scores['bbb_score'] = serotype.bbb_efficiency
    else:
        scores['bbb_score'] = 1.0  # Not relevant

    # Liver de-targeting
    if config.minimize_liver:
        liver_penalty = 1.0 - 0.5 * serotype.liver_tropism
        scores['liver_score'] = liver_penalty
    else:
        scores['liver_score'] = 1.0

    # Immunogenicity
    immuno_scores = {'low': 1.0, 'moderate': 0.7, 'high': 0.4}
    scores['immunogenicity_score'] = immuno_scores.get(
        serotype.immunogenicity, 0.5
    )

    # Pre-existing immunity
    if config.avoid_high_immunity:
        immunity_score = 1.0 - serotype.pre_existing_immunity
        scores['pre_existing_immunity_score'] = immunity_score
    else:
        scores['pre_existing_immunity_score'] = 1.0

    # Human validation requirement
    if config.require_human_validated and not serotype.human_validated:
        scores['human_validation'] = 0.0
    else:
        scores['human_validation'] = 1.0

    # Mouse-specific warning
    if serotype.mouse_specific:
        scores['mouse_specific_warning'] = True
        scores['human_validation'] *= 0.5

    # Combined score
    combined = (
        0.35 * scores['target_tropism'] +
        0.20 * scores.get('bbb_score', 1.0) +
        0.15 * scores['liver_score'] +
        0.15 * scores['immunogenicity_score'] +
        0.10 * scores['pre_existing_immunity_score'] +
        0.05 * scores['human_validation']
    )

    scores['combined_score'] = float(combined)

    return scores


def check_packaging_constraints(
    payload_size: int,
    serotype: Optional[AAVSerotype] = None,
) -> Dict[str, Any]:
    """
    Check if payload fits in AAV packaging constraints.

    Args:
        payload_size: Total construct size in bp
        serotype: Optional specific serotype (affects capacity)

    Returns:
        Packaging assessment dictionary
    """
    capacity = serotype.packaging_capacity if serotype else 4700

    result = {
        'payload_size': payload_size,
        'capacity': capacity,
        'fits': payload_size <= capacity,
        'utilization': payload_size / capacity,
    }

    if payload_size <= capacity:
        result['status'] = 'OK'
        result['recommendation'] = 'Payload fits in single AAV vector'
    elif payload_size <= AAV_PACKAGING_LIMITS['oversized']:
        result['status'] = 'WARNING'
        result['recommendation'] = (
            'Oversized payload. Reduced packaging efficiency expected. '
            'Consider compact Cas variants.'
        )
    elif payload_size <= AAV_PACKAGING_LIMITS['maximum']:
        result['status'] = 'CRITICAL'
        result['recommendation'] = (
            'At maximum capacity. Poor efficiency. '
            'Strongly recommend smaller construct.'
        )
    else:
        result['status'] = 'IMPOSSIBLE'
        result['recommendation'] = (
            'Payload too large for single AAV. '
            'Options: (1) Dual-vector approach, (2) Smaller Cas (SaCas9), '
            '(3) Mini-promoter, (4) Split-intein system'
        )

    return result


def select_optimal_serotype(
    target_tissue: str,
    payload_size: int,
    config: Optional[AAVConfig] = None,
    available_serotypes: Optional[List[str]] = None,
) -> List[Tuple[AAVSerotype, Dict[str, float]]]:
    """
    Select optimal AAV serotype for a therapeutic application.

    Args:
        target_tissue: Target tissue
        payload_size: Construct size
        config: AAVConfig
        available_serotypes: Optional list of serotypes to consider

    Returns:
        List of (serotype, scores) tuples, sorted by combined score

    Example:
        >>> results = select_optimal_serotype(
        ...     target_tissue="brain",
        ...     payload_size=3500,
        ...     config=AAVConfig(
        ...         minimize_liver=True,
        ...         require_bbb_crossing=True
        ...     )
        ... )
        >>> best = results[0]
        >>> print(f"Best serotype: {best[0].name} (score: {best[1]['combined_score']:.2f})")
    """
    if config is None:
        config = AAVConfig()

    if available_serotypes is None:
        serotypes = list(SEROTYPE_PROFILES.values())
    else:
        serotypes = [
            SEROTYPE_PROFILES[s] for s in available_serotypes
            if s in SEROTYPE_PROFILES
        ]

    results = []

    for serotype in serotypes:
        # Check packaging
        packaging = check_packaging_constraints(payload_size, serotype)
        if packaging['status'] == 'IMPOSSIBLE':
            continue  # Skip if can't package

        # Score tropism
        scores = score_serotype_tropism(serotype, target_tissue, config)

        # Apply packaging penalty
        if packaging['status'] == 'WARNING':
            scores['combined_score'] *= 0.85
        elif packaging['status'] == 'CRITICAL':
            scores['combined_score'] *= 0.6

        scores['packaging_status'] = packaging['status']

        # Filter by minimum tropism
        if scores['target_tropism'] < config.min_tropism_score:
            continue

        results.append((serotype, scores))

    # Sort by combined score
    results.sort(key=lambda x: x[1]['combined_score'], reverse=True)

    return results


def get_serotype_recommendations(
    target_tissue: str,
    application: str = "crispra",
) -> Dict[str, Any]:
    """
    Get serotype recommendations for common applications.

    Args:
        target_tissue: Target tissue
        application: 'crispra', 'knockout', 'base_editing', etc.

    Returns:
        Recommendation dictionary
    """
    tissue_lower = target_tissue.lower()

    # Application-specific size estimates
    app_sizes = {
        'crispra': 4500,  # dCas9-VPR
        'crispri': 4200,  # dCas9-KRAB
        'knockout': 4200,  # SpCas9
        'base_editing': 5200,  # Larger editors
        'prime_editing': 6000,  # Very large
        'sacas9': 3200,  # Smaller Cas
    }

    payload = app_sizes.get(application, 4000)

    # Tissue-specific recommendations
    recommendations = {
        'brain': {
            'primary': 'AAV9',
            'alternatives': ['AAVrh10', 'AAV5'],
            'notes': 'AAV9 is FDA-approved and crosses BBB. Consider AAVrh10 for lower immunity.',
            'avoid': 'AAV-PHP.eB (mouse-specific)',
        },
        'liver': {
            'primary': 'AAV8',
            'alternatives': ['AAV-DJ', 'AAV2'],
            'notes': 'AAV8 has excellent liver tropism with low immunogenicity.',
        },
        'muscle': {
            'primary': 'AAV9',
            'alternatives': ['AAV1', 'AAV8'],
            'notes': 'AAV9 or AAV1 for systemic muscle transduction.',
        },
        'heart': {
            'primary': 'AAV9',
            'alternatives': ['AAV1', 'AAV8'],
            'notes': 'AAV9 has strong cardiac tropism.',
        },
        'retina': {
            'primary': 'AAV2',
            'alternatives': ['AAV5', 'AAV8'],
            'notes': 'AAV2 is proven for retinal gene therapy (Luxturna).',
        },
    }

    rec = recommendations.get(tissue_lower, {
        'primary': 'AAV9',
        'alternatives': ['AAV8', 'AAV5'],
        'notes': 'AAV9 has broad tropism for unspecified tissues.',
    })

    rec['application'] = application
    rec['estimated_payload_size'] = payload

    # Check if payload fits
    if payload > 4700:
        rec['packaging_warning'] = (
            f'{application} construct (~{payload}bp) exceeds single AAV capacity. '
            'Consider SaCas9 variant or dual-vector approach.'
        )

    return rec


# Delivery route compatibility
DELIVERY_ROUTES = {
    'IV': {
        'description': 'Intravenous (systemic)',
        'best_for': ['liver', 'muscle', 'heart'],
        'cns_penetration': 'variable',
        'dose_range': '1e13 - 1e14 vg/kg',
    },
    'ICV': {
        'description': 'Intracerebroventricular',
        'best_for': ['brain', 'spinal_cord'],
        'cns_penetration': 'excellent',
        'dose_range': '1e10 - 1e12 vg total',
    },
    'IT': {
        'description': 'Intrathecal',
        'best_for': ['spinal_cord', 'brain'],
        'cns_penetration': 'good',
        'dose_range': '1e11 - 1e13 vg total',
    },
    'IM': {
        'description': 'Intramuscular',
        'best_for': ['muscle'],
        'cns_penetration': 'none',
        'dose_range': '1e12 - 1e13 vg/injection',
    },
    'subretinal': {
        'description': 'Subretinal injection',
        'best_for': ['retina'],
        'cns_penetration': 'none',
        'dose_range': '1e10 - 1e11 vg/eye',
    },
}


def recommend_delivery_route(
    target_tissue: str,
    serotype: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Recommend delivery route for target tissue.

    Args:
        target_tissue: Target tissue
        serotype: Optional serotype being used

    Returns:
        Delivery route recommendation
    """
    tissue_lower = target_tissue.lower()

    route_priority = {
        'brain': ['ICV', 'IT', 'IV'],
        'spinal_cord': ['IT', 'ICV', 'IV'],
        'liver': ['IV'],
        'muscle': ['IM', 'IV'],
        'heart': ['IV'],
        'retina': ['subretinal'],
        'lung': ['IV', 'intratracheal'],
    }

    recommended = route_priority.get(tissue_lower, ['IV'])
    primary_route = recommended[0]

    result = {
        'primary_route': primary_route,
        'alternative_routes': recommended[1:] if len(recommended) > 1 else [],
        'route_details': DELIVERY_ROUTES.get(primary_route, {}),
    }

    # Serotype-specific considerations
    if serotype and serotype in SEROTYPE_PROFILES:
        stype = SEROTYPE_PROFILES[serotype]
        if tissue_lower in ['brain', 'spinal_cord'] and not stype.crosses_bbb:
            result['warning'] = (
                f'{serotype} has limited BBB penetration. '
                f'Consider direct CNS delivery (ICV/IT) or BBB-crossing serotype.'
            )

    return result
