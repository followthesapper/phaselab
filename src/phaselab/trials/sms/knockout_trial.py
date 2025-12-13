"""
Knockout Validation Trial for Smith-Magenis Syndrome Models.

This trial generates clean knockout guides for RAI1 modifier genes
to validate disease models and study gene function - NOT for therapy.

Purpose: Create knockout cell lines or animal models to:
1. Validate modifier gene contributions to SMS phenotype
2. Study circadian pathway interactions
3. Develop preclinical models for drug screening
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


# Example gene CDS sequences for knockout targeting
KNOCKOUT_SEQUENCES = {
    'RAI1': (
        "ATGCGATCGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGATCGATCG"
        "CGATCGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGATCGATCGCG"
        "ATCGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGATCGATCGCGAT"
        "CGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGATCGATCGCGATCG"
    ),
    'PER1': (
        "ATGGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA"
        "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG"
        "CATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA"
        "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG"
    ),
    'PER2': (
        "ATGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"
        "GTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"
        "GTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    ),
}


def run_sms_knockout_trial(
    target_gene: str = "PER1",
    cds_sequence: Optional[str] = None,
    config: Optional[SMSTrialConfig] = None,
) -> SMSTrialResult:
    """
    Run knockout trial for SMS model validation.

    Generates knockout guides for creating disease models or
    validating modifier gene function - NOT for therapeutic use.

    Args:
        target_gene: Gene to knock out (for model creation).
        cds_sequence: Coding sequence. Uses example if None.
        config: Trial configuration.

    Returns:
        SMSTrialResult with knockout guide candidates.

    Example:
        >>> from phaselab.trials.sms import run_sms_knockout_trial
        >>> result = run_sms_knockout_trial(target_gene="PER1")
        >>> print(f"Found {result.n_candidates} knockout guides")
        >>> if result.best_candidate:
        ...     print(f"Frameshift probability: {result.best_candidate['frameshift_prob']:.0%}")
    """
    if config is None:
        config = SMSTrialConfig()

    if cds_sequence is None:
        cds_sequence = KNOCKOUT_SEQUENCES.get(target_gene, KNOCKOUT_SEQUENCES['PER1'])

    if config.verbose:
        print("=" * 60)
        print(f"SMS Knockout Trial: {target_gene} Model Validation")
        print("=" * 60)
        print(f"Target: {target_gene} coding sequence ({len(cds_sequence)} bp)")
        print("PURPOSE: Model creation/validation, NOT therapy")
        print()

    warnings = []
    errors = []

    # Add clear warning about research use
    warnings.append(
        f"RESEARCH USE ONLY: {target_gene} knockout guides are for model "
        f"validation, not therapeutic application."
    )

    try:
        # Import design functions
        from phaselab.crispr import design_knockout_guides, KnockoutConfig
        from phaselab.fusion import ClaimLevel

        # Design knockout guides
        ko_config = KnockoutConfig(
            top_n=config.top_n_guides * 2,
            min_cut_efficiency=0.4,  # Higher threshold for knockout
        )

        raw_df = design_knockout_guides(
            sequence=cds_sequence,
            cds_start=0,  # Sequence starts at CDS
            config=ko_config,
        )

        raw_guides = raw_df.to_dict('records') if len(raw_df) > 0 else []

        if config.verbose:
            print(f"Initial candidates: {len(raw_guides)}")

        # Evaluate guides
        candidates = []

        for guide in raw_guides:
            coherence = guide.get('coherence_R') or guide.get('coherence', 0.5)
            go_status = guide.get('go_no_go', 'UNKNOWN')

            # Filter by GO status if required
            if config.require_go_status and go_status != 'GO':
                continue

            # Get knockout-specific scores
            cut_efficiency = guide.get('cut_efficiency', 0.5)
            frameshift_prob = guide.get('frameshift_prob', 0.5)
            nhej_prob = guide.get('nhej_prob', 0.8)

            # For knockout, we want high cut efficiency and frameshift
            ko_quality = cut_efficiency * frameshift_prob * nhej_prob

            candidate = {
                'sequence': guide.get('sequence', ''),
                'pam': guide.get('pam', 'NGG'),
                'cds_position': guide.get('cds_position', 0),
                'strand': guide.get('strand', '+'),

                # Scores
                'coherence_R': coherence,
                'go_no_go': go_status,
                'cut_efficiency': cut_efficiency,
                'frameshift_prob': frameshift_prob,
                'nhej_prob': nhej_prob,
                'repair_pathway': guide.get('repair_pathway', 'NHEJ'),

                # Quality assessment
                'ko_quality_score': ko_quality,
                'is_high_quality': ko_quality > 0.4,

                # Binding
                'delta_g': guide.get('delta_g', -15.0),
                'mit_score': guide.get('mit_score', 50),
                'cfd_score': guide.get('cfd_score', 50),

                # Claim level
                'claim_level': 'context_dependent',  # Knockout is well-established
            }

            candidates.append(candidate)

        # Sort by knockout quality
        candidates.sort(
            key=lambda x: (
                x['is_high_quality'],
                x['ko_quality_score'],
                x['coherence_R'],
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
                print(f"  Cut efficiency: {best['cut_efficiency']:.0%}")
                print(f"  Frameshift probability: {best['frameshift_prob']:.0%}")
                print(f"  CDS position: {best['cds_position']}")

        # Determine overall claim level
        if len(candidates) == 0:
            overall_claim = ClaimLevel.UNKNOWN
            claim_desc = f"No viable knockout guides found for {target_gene}."
        elif candidates[0]['is_high_quality']:
            overall_claim = ClaimLevel.STRONG_COMPUTATIONAL
            claim_desc = (
                f"Strong computational evidence for {target_gene} knockout. "
                f"High-quality guides with excellent cut efficiency and frameshift probability."
            )
        else:
            overall_claim = ClaimLevel.CONTEXT_DEPENDENT
            claim_desc = (
                f"Context-dependent evidence for {target_gene} knockout. "
                f"Guides available but efficiency may be moderate."
            )

        # Build metrics
        metrics = {
            'n_initial_guides': len(raw_guides),
            'n_viable_candidates': len(candidates),
            'n_high_quality': sum(1 for c in candidates if c['is_high_quality']),
            'target_gene': target_gene,
            'purpose': 'model_validation',
        }

        if candidates:
            metrics['best_cut_efficiency'] = candidates[0]['cut_efficiency']
            metrics['best_frameshift_prob'] = candidates[0]['frameshift_prob']

        # Build summary
        if len(candidates) > 0:
            best = candidates[0]
            summary = (
                f"Knockout {target_gene} Trial: Found {len(candidates)} guides. "
                f"Top guide has {best['cut_efficiency']:.0%} cut efficiency and "
                f"{best['frameshift_prob']:.0%} frameshift probability. "
                f"FOR MODEL VALIDATION ONLY."
            )
            status = TrialStatus.COMPLETED
        else:
            summary = (
                f"Knockout {target_gene} Trial: No viable guides found. "
                f"Consider alternative targeting regions or exons."
            )
            status = TrialStatus.COMPLETED

        return SMSTrialResult(
            trial_type=TrialType.KNOCKOUT_VALIDATION,
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
                'modality': 'Knockout',
                'purpose': 'model_validation',
                'cds_length': len(cds_sequence),
            },
        )

    except Exception as e:
        logger.error(f"Knockout trial failed: {e}")
        return SMSTrialResult(
            trial_type=TrialType.KNOCKOUT_VALIDATION,
            status=TrialStatus.FAILED,
            summary=f"Knockout {target_gene} Trial failed: {str(e)}",
            errors=[str(e)],
            claim_level="unknown",
            claim_description="Trial execution failed.",
        )
