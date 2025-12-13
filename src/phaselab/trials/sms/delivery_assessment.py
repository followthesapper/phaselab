"""
Delivery Feasibility Assessment for Smith-Magenis Syndrome.

This trial evaluates whether top CRISPRa strategies can be delivered
to CNS tissues using AAV vectors.

Key considerations:
1. Payload size (CRISPRa requires dCas9-VPR/VP64, ~4.5-5.2kb)
2. Serotype selection for CNS tropism
3. Blood-brain barrier penetration
4. Immunogenicity and pre-existing immunity
5. Delivery route (IV vs ICV/intrathecal)

The assessment provides a GO/NO-GO for delivery feasibility.
"""

import numpy as np
from typing import Dict, Any, List, Optional
import logging

from .core import (
    SMSTrialResult,
    SMSTrialConfig,
    TrialType,
    TrialStatus,
)

logger = logging.getLogger(__name__)


# Payload sizes for different CRISPR modalities
PAYLOAD_SIZES = {
    'CRISPRa_VPR': 5100,       # dCas9-VPR - oversized
    'CRISPRa_VP64': 4500,      # dCas9-VP64 - borderline
    'CRISPRa_mini': 3800,      # miniCas + activator
    'CRISPRi_KRAB': 4200,      # dCas9-KRAB
    'knockout_SpCas9': 4200,   # SpCas9
    'knockout_SaCas9': 3200,   # SaCas9 (fits easily)
    'base_editing': 5200,      # ABE/CBE - large
    'prime_editing': 6000,     # PE2 - requires split
    'guide_only': 500,         # Guide expression cassette
}


def run_delivery_assessment(
    modality: str = "CRISPRa_VP64",
    target_tissue: str = "brain",
    config: Optional[SMSTrialConfig] = None,
) -> SMSTrialResult:
    """
    Assess delivery feasibility for SMS gene therapy.

    Evaluates whether the therapeutic payload can be delivered
    to CNS tissues using current AAV technology.

    Args:
        modality: CRISPR modality being used.
        target_tissue: Target tissue for delivery.
        config: Trial configuration.

    Returns:
        SMSTrialResult with delivery feasibility assessment.

    Example:
        >>> from phaselab.trials.sms import run_delivery_assessment
        >>> result = run_delivery_assessment(modality="CRISPRa_VP64")
        >>> print(f"Feasibility: {result.metrics['delivery_feasibility']}")
        >>> print(f"Recommended serotype: {result.best_candidate['recommended_serotype']}")
    """
    if config is None:
        config = SMSTrialConfig()

    payload_size = PAYLOAD_SIZES.get(modality, 4500)

    if config.verbose:
        print("=" * 60)
        print("SMS Delivery Feasibility Assessment")
        print("=" * 60)
        print(f"Modality: {modality}")
        print(f"Payload size: {payload_size}bp")
        print(f"Target tissue: {target_tissue}")
        print(f"Delivery route: {config.delivery_route}")
        print()

    warnings = []
    errors = []

    try:
        # Import delivery functions
        from phaselab.delivery.aav import (
            select_optimal_serotype,
            check_packaging_constraints,
            recommend_delivery_route,
            get_serotype_recommendations,
            AAVConfig,
            SEROTYPE_PROFILES,
        )
        from phaselab.fusion import ClaimLevel

        # Check packaging constraints
        packaging = check_packaging_constraints(payload_size)

        if config.verbose:
            print(f"Packaging status: {packaging['status']}")
            print(f"  Utilization: {packaging['utilization']:.0%}")
            print(f"  {packaging['recommendation']}")

        # Configure AAV selection
        aav_config = AAVConfig(
            target_tissue=target_tissue,
            payload_size=payload_size,
            minimize_liver=True,
            require_bbb_crossing=(target_tissue.lower() in ['brain', 'cns', 'hypothalamus']),
            require_human_validated=True,
            avoid_high_immunity=True,
        )

        # Select optimal serotypes
        serotype_results = select_optimal_serotype(
            target_tissue=target_tissue,
            payload_size=payload_size,
            config=aav_config,
        )

        if config.verbose:
            print(f"\nSerotype ranking ({len(serotype_results)} options):")
            for i, (serotype, scores) in enumerate(serotype_results[:5]):
                print(f"  {i+1}. {serotype.name}: {scores['combined_score']:.2f}")

        # Get delivery route recommendation
        route_rec = recommend_delivery_route(
            target_tissue=target_tissue,
            serotype=serotype_results[0][0].name if serotype_results else None,
        )

        # Get general recommendations
        general_rec = get_serotype_recommendations(
            target_tissue=target_tissue,
            application=modality.lower().split('_')[0],
        )

        # Build candidates (serotype options)
        candidates = []

        for serotype, scores in serotype_results[:5]:
            candidate = {
                'serotype_name': serotype.name,
                'full_name': serotype.full_name,
                'combined_score': scores['combined_score'],
                'target_tropism': scores['target_tropism'],
                'bbb_efficiency': serotype.bbb_efficiency,
                'crosses_bbb': serotype.crosses_bbb,
                'liver_tropism': serotype.liver_tropism,
                'immunogenicity': serotype.immunogenicity,
                'pre_existing_immunity': serotype.pre_existing_immunity,
                'clinical_stage': serotype.clinical_stage,
                'human_validated': serotype.human_validated,
                'packaging_status': packaging['status'],
                'notes': serotype.notes,
            }

            # Add warnings for specific issues
            candidate_warnings = []
            if serotype.mouse_specific:
                candidate_warnings.append("Mouse-specific receptor (not translatable)")
            if serotype.pre_existing_immunity > 0.4:
                candidate_warnings.append(f"High pre-existing immunity ({serotype.pre_existing_immunity:.0%})")
            if not serotype.crosses_bbb and target_tissue.lower() in ['brain', 'cns']:
                candidate_warnings.append("Does not efficiently cross BBB")

            candidate['candidate_warnings'] = candidate_warnings
            candidates.append(candidate)

        # Determine best serotype
        best_serotype = candidates[0] if candidates else None

        # Determine overall feasibility
        if packaging['status'] == 'IMPOSSIBLE':
            delivery_feasibility = 'NOT_FEASIBLE'
            feasibility_reason = (
                f"Payload ({payload_size}bp) too large for single AAV. "
                f"Requires dual-vector approach or smaller cargo."
            )
        elif packaging['status'] == 'CRITICAL':
            delivery_feasibility = 'CHALLENGING'
            feasibility_reason = (
                f"Payload at maximum AAV capacity. Packaging efficiency will be poor. "
                f"Consider compact Cas variants."
            )
        elif not serotype_results:
            delivery_feasibility = 'NOT_FEASIBLE'
            feasibility_reason = "No suitable serotypes found for target tissue."
        elif best_serotype and best_serotype['combined_score'] > 0.6:
            delivery_feasibility = 'FEASIBLE'
            feasibility_reason = (
                f"Delivery appears feasible with {best_serotype['serotype_name']}. "
                f"Good CNS tropism and acceptable immunogenicity profile."
            )
        elif best_serotype and best_serotype['combined_score'] > 0.4:
            delivery_feasibility = 'CONDITIONALLY_FEASIBLE'
            feasibility_reason = (
                f"Delivery possible but with caveats. "
                f"May need to optimize route or consider pre-screening for immunity."
            )
        else:
            delivery_feasibility = 'CHALLENGING'
            feasibility_reason = "Delivery feasible but suboptimal. Consider alternative approaches."

        if config.verbose:
            print(f"\nFeasibility: {delivery_feasibility}")
            print(f"  {feasibility_reason}")

        # Build the delivery recommendation
        if best_serotype:
            delivery_recommendation = {
                'recommended_serotype': best_serotype['serotype_name'],
                'recommended_route': route_rec['primary_route'],
                'alternative_routes': route_rec.get('alternative_routes', []),
                'feasibility': delivery_feasibility,
                'payload_size': payload_size,
                'packaging_status': packaging['status'],
                'dose_range': route_rec.get('route_details', {}).get('dose_range', 'Unknown'),
            }
        else:
            delivery_recommendation = {
                'recommended_serotype': None,
                'feasibility': delivery_feasibility,
                'payload_size': payload_size,
                'packaging_status': packaging['status'],
            }

        # Determine overall claim level
        if delivery_feasibility == 'FEASIBLE':
            overall_claim = ClaimLevel.CONTEXT_DEPENDENT
            claim_desc = (
                f"Context-dependent evidence for delivery feasibility. "
                f"{best_serotype['serotype_name']} with {route_rec['primary_route']} route "
                f"is a reasonable approach for SMS CNS delivery."
            )
        elif delivery_feasibility == 'CONDITIONALLY_FEASIBLE':
            overall_claim = ClaimLevel.EXPLORATORY
            claim_desc = (
                f"Exploratory evidence for delivery. Feasible but with caveats. "
                f"Pre-clinical validation strongly recommended."
            )
        elif delivery_feasibility == 'CHALLENGING':
            overall_claim = ClaimLevel.EXPLORATORY
            claim_desc = (
                f"Exploratory evidence only. Delivery faces significant challenges. "
                f"Alternative strategies may be needed."
            )
        else:
            overall_claim = ClaimLevel.UNKNOWN
            claim_desc = (
                f"Delivery not feasible with current approach. "
                f"Recommend payload size reduction or dual-vector system."
            )

        # Add standard warnings
        if best_serotype and best_serotype['pre_existing_immunity'] > 0.3:
            warnings.append(
                f"{best_serotype['serotype_name']} has {best_serotype['pre_existing_immunity']:.0%} "
                f"pre-existing immunity in population. Patient screening may be needed."
            )

        if packaging['status'] in ['WARNING', 'CRITICAL']:
            warnings.append(packaging['recommendation'])

        if target_tissue.lower() in ['brain', 'cns', 'hypothalamus']:
            warnings.append(
                "CNS delivery requires careful dose optimization to avoid toxicity. "
                "Consider intrathecal or ICV routes for more controlled exposure."
            )

        # Build metrics
        metrics = {
            'modality': modality,
            'payload_size': payload_size,
            'target_tissue': target_tissue,
            'delivery_feasibility': delivery_feasibility,
            'packaging_status': packaging['status'],
            'packaging_utilization': packaging['utilization'],
            'n_viable_serotypes': len(candidates),
            'recommended_serotype': best_serotype['serotype_name'] if best_serotype else None,
            'recommended_route': route_rec['primary_route'],
        }

        if best_serotype:
            metrics['serotype_score'] = best_serotype['combined_score']
            metrics['target_tropism'] = best_serotype['target_tropism']
            metrics['bbb_crossing'] = best_serotype['crosses_bbb']

        # Build summary
        if delivery_feasibility == 'FEASIBLE':
            summary = (
                f"Delivery Assessment: FEASIBLE. "
                f"Recommended: {best_serotype['serotype_name']} via {route_rec['primary_route']} route. "
                f"Payload ({payload_size}bp) fits single AAV ({packaging['utilization']:.0%} capacity). "
                f"Claim level: {overall_claim.value}."
            )
        elif delivery_feasibility == 'NOT_FEASIBLE':
            summary = (
                f"Delivery Assessment: NOT FEASIBLE. "
                f"Payload ({payload_size}bp) exceeds single AAV capacity. "
                f"Recommend: dual-vector system, SaCas9, or mini-activators."
            )
        else:
            summary = (
                f"Delivery Assessment: {delivery_feasibility}. "
                f"{feasibility_reason}"
            )

        return SMSTrialResult(
            trial_type=TrialType.DELIVERY_ASSESSMENT,
            status=TrialStatus.COMPLETED,
            summary=summary,
            candidates=candidates,
            best_candidate=delivery_recommendation,
            metrics=metrics,
            claim_level=overall_claim.value,
            claim_description=claim_desc,
            warnings=warnings,
            errors=errors,
            metadata={
                'modality': modality,
                'payload_sizes': PAYLOAD_SIZES,
                'packaging_details': packaging,
                'route_recommendation': route_rec,
                'general_recommendations': general_rec,
            },
        )

    except Exception as e:
        logger.error(f"Delivery assessment failed: {e}")
        return SMSTrialResult(
            trial_type=TrialType.DELIVERY_ASSESSMENT,
            status=TrialStatus.FAILED,
            summary=f"Delivery Assessment failed: {str(e)}",
            errors=[str(e)],
            claim_level="unknown",
            claim_description="Assessment execution failed.",
        )
