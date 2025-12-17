"""
CRISPRi Modifier Suppression Trial for Smith-Magenis Syndrome.

This trial identifies CRISPRi guides to suppress circadian-disrupting
modifier genes, potentially improving SMS phenotype alongside RAI1 activation.

Target genes include PER1/2, CRY1, and CLOCK family members.

CAUTION: These are essential circadian clock genes. Suppression must be
carefully calibrated to avoid complete loss of circadian rhythmicity.
"""

import numpy as np
from typing import Dict, Any, List, Optional
import logging

from .core import (
    SMSTrialResult,
    SMSTrialConfig,
    TrialType,
    TrialStatus,
    SMS_MODIFIER_GENES,
)

logger = logging.getLogger(__name__)


# Example modifier gene sequences (truncated for demonstration)
MODIFIER_SEQUENCES = {
    'PER1': (
        "CGCGATCGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGATCGATCG"
        "CGATCGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGATCGATCGCG"
        "ATCGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGATCGATCGCGAT"
        "CGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGATCGATCGCGATCG"
    ),
    'PER2': (
        "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
        "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
        "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
        "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
    ),
    'CRY1': (
        "TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA"
        "CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG"
        "TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA"
        "CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG"
    ),
}


def run_sms_crispri_trial(
    target_gene: str = "PER1",
    promoter_sequence: Optional[str] = None,
    config: Optional[SMSTrialConfig] = None,
) -> SMSTrialResult:
    """
    Run CRISPRi trial for modifier gene suppression in SMS.

    Identifies guides to suppress circadian modifier genes that may
    contribute to SMS phenotype when RAI1 is deficient.

    Args:
        target_gene: Modifier gene to target (PER1, PER2, CRY1, CLOCK).
        promoter_sequence: Gene promoter sequence. Uses example if None.
        config: Trial configuration.

    Returns:
        SMSTrialResult with ranked guide candidates.

    Example:
        >>> from phaselab.trials.sms import run_sms_crispri_trial
        >>> result = run_sms_crispri_trial(target_gene="PER1")
        >>> print(f"Found {result.n_candidates} CRISPRi guides for PER1")
    """
    if config is None:
        config = SMSTrialConfig()

    if promoter_sequence is None:
        promoter_sequence = MODIFIER_SEQUENCES.get(target_gene, MODIFIER_SEQUENCES['PER1'])

    # Get gene info
    gene_info = SMS_MODIFIER_GENES.get(target_gene, {
        'role': 'Unknown modifier gene',
        'rationale': 'Potential circadian modifier',
        'caution': 'Use with care',
    })

    if config.verbose:
        print("=" * 60)
        print(f"SMS CRISPRi Trial: {target_gene} Suppression")
        print("=" * 60)
        print(f"Target: {target_gene} ({gene_info['role']})")
        print(f"Rationale: {gene_info['rationale']}")
        print(f"CAUTION: {gene_info['caution']}")
        print()

    warnings = []
    errors = []

    # Safety check for BMAL1
    if target_gene == 'BMAL1':
        warnings.append(
            "BMAL1 suppression is NOT recommended for SMS. "
            "BMAL1 is already dysregulated and further suppression may worsen phenotype."
        )

    try:
        # Import design functions
        from phaselab.crispr import design_crispri_guides, CRISPRiConfig
        from phaselab.crispr.enhanced_pipeline import (
            design_enhanced_guides,
            EnhancedGuideConfig,
            Modality,
        )
        from phaselab.fusion import ClaimLevel

        # Design CRISPRi guides
        if config.use_virtual_assay:
            guide_config = EnhancedGuideConfig(
                modality=Modality.CRISPRI,
                top_n=config.top_n_guides * 2,
            )

            if config.verbose:
                print("Using Virtual Assay Stack (v0.7.0+)...")

            design_result = design_enhanced_guides(
                sequence=promoter_sequence,
                tss_index=len(promoter_sequence) // 2,
                config=guide_config,
            )

            raw_guides = [g.__dict__ for g in design_result.guides]

        else:
            crispri_config = CRISPRiConfig(
                top_n=config.top_n_guides * 2,
                repressor="KRAB",
            )

            raw_df = design_crispri_guides(
                sequence=promoter_sequence,
                tss_index=len(promoter_sequence) // 2,
                config=crispri_config,
            )

            raw_guides = raw_df.to_dict('records') if len(raw_df) > 0 else []

        if config.verbose:
            print(f"Initial candidates: {len(raw_guides)}")

        # Evaluate guides
        candidates = []

        # For modifier suppression, we want MODERATE repression (30-60%)
        # Complete knockout would be too severe
        target_repression_min = 0.30
        target_repression_max = 0.60

        for guide in raw_guides:
            coherence = guide.get('coherence_R') or guide.get('coherence', 0.5)
            go_status = guide.get('go_no_go', 'UNKNOWN')

            # Filter by GO status if required
            if config.require_go_status and go_status != 'GO':
                continue

            # Get repression efficiency
            repression_eff = guide.get('repression_efficiency', 0.5)
            steric = guide.get('steric_hindrance', 0.5)

            # Calculate expected suppression level
            # Higher repression efficiency = more suppression
            expected_suppression = repression_eff * 0.7  # Scale factor for CRISPRi reality

            # Penalize if suppression would be too strong
            if expected_suppression > 0.70:
                # Add warning but don't exclude
                guide_warning = f"Strong suppression ({expected_suppression:.0%}) - use lower dose"
            else:
                guide_warning = None

            candidate = {
                'sequence': guide.get('sequence', ''),
                'pam': guide.get('pam', 'NGG'),
                'position': guide.get('position', 0),
                'strand': guide.get('strand', '+'),

                # Scores
                'coherence_R': coherence,
                'go_no_go': go_status,
                'repression_efficiency': repression_eff,
                'steric_hindrance': steric,
                'fused_score': guide.get('fused_score', guide.get('combined_score', 0)),
                'delta_g': guide.get('delta_g', -15.0),

                # Suppression assessment
                'expected_suppression': expected_suppression,
                'in_target_range': target_repression_min <= expected_suppression <= target_repression_max,

                # Claim level
                'claim_level': guide.get('claim_level', 'exploratory'),
                'claim_description': guide.get('claim_description', ''),

                # Safety
                'guide_warning': guide_warning,
            }

            candidates.append(candidate)

        # Sort by suitability for moderate suppression
        candidates.sort(
            key=lambda x: (
                x['in_target_range'],  # In target range first
                -abs(x['expected_suppression'] - 0.45),  # Closest to 45% suppression
                x['coherence_R'],  # Higher coherence
            ),
            reverse=True,
        )

        # Trim to top_n
        candidates = candidates[:config.top_n_guides]

        if config.verbose:
            print(f"Viable candidates: {len(candidates)}")
            if candidates:
                best = candidates[0]
                print(f"\nTop candidate:")
                print(f"  Sequence: {best['sequence']}")
                print(f"  Expected suppression: {best['expected_suppression']:.0%}")
                print(f"  In target range: {best['in_target_range']}")

        # Determine overall claim level
        if len(candidates) == 0:
            overall_claim = ClaimLevel.UNKNOWN
            claim_desc = f"No viable CRISPRi guides found for {target_gene} suppression."
        elif candidates[0]['in_target_range']:
            if candidates[0]['coherence_R'] > 0.5:
                overall_claim = ClaimLevel.CONTEXT_DEPENDENT
                claim_desc = (
                    f"Context-dependent evidence for {target_gene} modulation. "
                    f"Guides predicted to achieve moderate suppression suitable for SMS."
                )
            else:
                overall_claim = ClaimLevel.EXPLORATORY
                claim_desc = (
                    f"Exploratory evidence for {target_gene} suppression. "
                    f"Guides found but coherence is moderate."
                )
        else:
            overall_claim = ClaimLevel.EXPLORATORY
            claim_desc = (
                f"Exploratory evidence only. Guides available but suppression "
                f"levels may not be optimal for SMS modifier strategy."
            )

        # Add standard warnings for clock gene manipulation
        warnings.append(
            f"CRISPRi targeting {target_gene} will affect circadian rhythms. "
            f"Monitor sleep-wake cycles and melatonin profiles after treatment."
        )

        # Build metrics
        metrics = {
            'n_initial_guides': len(raw_guides),
            'n_viable_candidates': len(candidates),
            'n_in_target_range': sum(1 for c in candidates if c['in_target_range']),
            'target_gene': target_gene,
            'target_suppression_range': (target_repression_min, target_repression_max),
            'gene_role': gene_info['role'],
        }

        if candidates:
            metrics['best_expected_suppression'] = candidates[0]['expected_suppression']

        # Build summary
        if len(candidates) > 0:
            best = candidates[0]
            summary = (
                f"CRISPRi {target_gene} Trial: Found {len(candidates)} guides. "
                f"Top candidate expected to suppress {target_gene} by {best['expected_suppression']:.0%}. "
                f"Role: {gene_info['role']}. Claim level: {overall_claim.value}."
            )
            status = TrialStatus.COMPLETED
        else:
            summary = (
                f"CRISPRi {target_gene} Trial: No viable guides found. "
                f"Consider alternative modifier genes or targeting approaches."
            )
            status = TrialStatus.COMPLETED

        return SMSTrialResult(
            trial_type=TrialType.CRISPRI_MODIFIER,
            status=status,
            summary=summary,
            candidates=candidates,
            best_candidate=candidates[0] if candidates else None,
            metrics=metrics,
            claim_level=overall_claim.value,
            claim_description=claim_desc,
            warnings=warnings,
            errors=errors,
            metadata={
                'target_gene': target_gene,
                'gene_info': gene_info,
                'modality': 'CRISPRi',
                'promoter_length': len(promoter_sequence),
            },
        )

    except Exception as e:
        logger.error(f"CRISPRi trial failed: {e}")
        return SMSTrialResult(
            trial_type=TrialType.CRISPRI_MODIFIER,
            status=TrialStatus.FAILED,
            summary=f"CRISPRi {target_gene} Trial failed: {str(e)}",
            errors=[str(e)],
            claim_level="unknown",
            claim_description="Trial execution failed.",
        )
