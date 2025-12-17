"""
CRISPRa RAI1 Activation Trial for Smith-Magenis Syndrome.

This is the core therapeutic strategy for SMS: identify CRISPRa guides
that safely upregulate RAI1 expression into the therapeutic window (70-110%).

The trial:
1. Designs CRISPRa guides targeting the RAI1 promoter
2. Uses the Virtual Assay Stack for multi-layer scoring
3. Validates guides against therapeutic window constraints
4. Simulates expected expression boost
5. Returns ranked candidates with claim levels
"""

import numpy as np
from typing import Dict, Any, List, Optional
import logging

from .core import (
    SMSTrialResult,
    SMSTrialConfig,
    TrialType,
    TrialStatus,
    RAI1_INFO,
    EXAMPLE_RAI1_PROMOTER,
)

logger = logging.getLogger(__name__)


def run_sms_crispra_trial(
    promoter_sequence: Optional[str] = None,
    config: Optional[SMSTrialConfig] = None,
) -> SMSTrialResult:
    """
    Run CRISPRa trial for RAI1 activation in SMS.

    This trial designs and evaluates CRISPRa guides to upregulate RAI1
    from the haploinsufficient baseline (~50%) to the therapeutic window
    (70-110% of normal expression).

    Args:
        promoter_sequence: RAI1 promoter sequence. Uses example if None.
        config: Trial configuration.

    Returns:
        SMSTrialResult with ranked guide candidates.

    Example:
        >>> from phaselab.trials.sms import run_sms_crispra_trial
        >>> result = run_sms_crispra_trial()
        >>> print(f"Found {result.n_candidates} candidates")
        >>> if result.best_candidate:
        ...     print(f"Best guide: {result.best_candidate['sequence']}")
        ...     print(f"Expected boost: {result.best_candidate['expected_fold_change']:.2f}x")
    """
    if config is None:
        config = SMSTrialConfig()

    if promoter_sequence is None:
        promoter_sequence = EXAMPLE_RAI1_PROMOTER

    if config.verbose:
        print("=" * 60)
        print("SMS CRISPRa Trial: RAI1 Activation")
        print("=" * 60)
        print(f"Target: RAI1 promoter ({len(promoter_sequence)} bp)")
        print(f"Therapeutic window: {config.therapeutic_window[0]:.0%} - {config.therapeutic_window[1]:.0%}")
        print(f"Baseline expression: {config.baseline_expression:.0%}")
        print()

    warnings = []
    errors = []

    try:
        # Import design functions
        from phaselab.crispr import design_guides, GuideDesignConfig
        from phaselab.crispr.enhanced_pipeline import (
            design_enhanced_guides,
            EnhancedGuideConfig,
            Modality,
        )
        from phaselab.therapy.dosage import (
            TherapeuticWindow,
            estimate_expression_change,
            validate_therapeutic_level,
        )
        from phaselab.fusion import ClaimLevel

        # Set up therapeutic window
        therapeutic_window = TherapeuticWindow(
            baseline_expression=config.baseline_expression,
            therapeutic_min=config.therapeutic_window[0],
            therapeutic_max=config.therapeutic_window[1],
            optimal_expression=config.optimal_expression,
            disease_name="Smith-Magenis Syndrome",
            gene_symbol="RAI1",
        )

        # Design guides using enhanced pipeline if available
        if config.use_virtual_assay:
            guide_config = EnhancedGuideConfig(
                modality=Modality.CRISPRA,
                top_n=config.top_n_guides * 2,  # Get extra for filtering
            )

            if config.verbose:
                print("Using Virtual Assay Stack (v0.7.0+)...")

            design_result = design_enhanced_guides(
                sequence=promoter_sequence,
                tss_index=len(promoter_sequence) // 2,  # Assume TSS in middle
                config=guide_config,
            )

            raw_guides = [g.__dict__ for g in design_result.guides]

        else:
            # Fallback to basic design
            guide_config = GuideDesignConfig(
                top_n=config.top_n_guides * 2,
            )

            raw_df = design_guides(
                sequence=promoter_sequence,
                tss_index=len(promoter_sequence) // 2,
                config=guide_config,
            )

            raw_guides = raw_df.to_dict('records') if len(raw_df) > 0 else []

        if config.verbose:
            print(f"Initial candidates: {len(raw_guides)}")

        # Evaluate each guide for therapeutic potential
        candidates = []

        for guide in raw_guides:
            # Get coherence and binding energy
            coherence = guide.get('coherence_R') or guide.get('coherence', 0.5)
            delta_g = guide.get('delta_g', -15.0)

            # Estimate expression change
            expression_est = estimate_expression_change(
                guide_coherence=coherence,
                binding_energy=delta_g,
            )

            expected_fold = expression_est['estimated_fold_change']
            achieved_expression = config.baseline_expression * expected_fold

            # Validate against therapeutic window
            validation = validate_therapeutic_level(achieved_expression, therapeutic_window)

            # Check GO/NO-GO status
            go_status = guide.get('go_no_go', 'UNKNOWN')
            if config.require_go_status and go_status != 'GO':
                continue

            # Check claim level
            claim_level = guide.get('claim_level', 'unknown')
            if _claim_level_rank(claim_level) < _claim_level_rank(config.min_claim_level):
                continue

            # Build candidate record
            candidate = {
                'sequence': guide.get('sequence', ''),
                'pam': guide.get('pam', 'NGG'),
                'position': guide.get('position', 0),
                'strand': guide.get('strand', '+'),

                # Scores
                'coherence_R': coherence,
                'go_no_go': go_status,
                'fused_score': guide.get('fused_score', guide.get('combined_score', 0)),
                'delta_g': delta_g,

                # Therapeutic assessment
                'expected_fold_change': expected_fold,
                'expected_expression': achieved_expression,
                'in_therapeutic_window': validation['in_therapeutic_window'],
                'therapeutic_score': validation['therapeutic_score'],
                'therapeutic_classification': validation['classification'],

                # Reliability
                'fold_change_low': expression_est['fold_change_low'],
                'fold_change_high': expression_est['fold_change_high'],
                'reliability': expression_est['reliability'],

                # Claim level
                'claim_level': claim_level,
                'claim_description': guide.get('claim_description', ''),

                # Safety
                'safety_margin': validation['safety_margin'],
                'distance_to_optimal': validation['distance_to_optimal'],
            }

            # Add to candidates if in or near therapeutic window
            if validation['in_therapeutic_window'] or achieved_expression < config.therapeutic_window[1]:
                candidates.append(candidate)

        # Sort by therapeutic suitability
        candidates.sort(
            key=lambda x: (
                x['in_therapeutic_window'],  # In window first
                x['therapeutic_score'],       # Higher score
                -x['distance_to_optimal'],    # Closer to optimal
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
                print(f"  Expected expression: {best['expected_expression']:.0%}")
                print(f"  In therapeutic window: {best['in_therapeutic_window']}")
                print(f"  Coherence R̄: {best['coherence_R']:.4f}")

        # Determine overall claim level
        if len(candidates) == 0:
            overall_claim = ClaimLevel.UNKNOWN
            claim_desc = "No viable CRISPRa guides found for RAI1 activation."
        elif all(c['in_therapeutic_window'] for c in candidates[:3]):
            # Multiple candidates in window
            avg_reliability = np.mean([c['reliability'] for c in candidates[:3]])
            if avg_reliability > 0.7:
                overall_claim = ClaimLevel.STRONG_COMPUTATIONAL
                claim_desc = "Strong computational evidence for therapeutic RAI1 activation. Multiple high-reliability guides predicted to achieve therapeutic expression."
            else:
                overall_claim = ClaimLevel.CONTEXT_DEPENDENT
                claim_desc = "Context-dependent evidence for therapeutic RAI1 activation. Guides predicted to achieve therapeutic expression, but reliability varies."
        elif candidates[0]['in_therapeutic_window']:
            overall_claim = ClaimLevel.CONTEXT_DEPENDENT
            claim_desc = "Context-dependent evidence. Top guide(s) predicted in therapeutic window, but limited alternatives."
        else:
            overall_claim = ClaimLevel.EXPLORATORY
            claim_desc = "Exploratory evidence only. No guides confidently predicted to achieve therapeutic window."

        # Check for warnings
        if len(candidates) > 0:
            best = candidates[0]
            if best['expected_expression'] > 1.0:
                warnings.append(
                    f"Top candidate may cause overexpression ({best['expected_expression']:.0%}). "
                    "Consider lower-efficiency guides or reduced dosing."
                )
            if best['reliability'] < 0.5:
                warnings.append(
                    f"Top candidate has low prediction reliability ({best['reliability']:.0%}). "
                    "Experimental validation strongly recommended."
                )

        # Build metrics
        metrics = {
            'n_initial_guides': len(raw_guides),
            'n_viable_candidates': len(candidates),
            'n_in_therapeutic_window': sum(1 for c in candidates if c['in_therapeutic_window']),
            'therapeutic_window': config.therapeutic_window,
            'baseline_expression': config.baseline_expression,
            'target_gene': 'RAI1',
        }

        if candidates:
            metrics['best_expected_expression'] = candidates[0]['expected_expression']
            metrics['best_fold_change'] = candidates[0]['expected_fold_change']
            metrics['avg_reliability'] = np.mean([c['reliability'] for c in candidates])

        # Build summary
        if len(candidates) > 0:
            best = candidates[0]
            summary = (
                f"CRISPRa RAI1 Trial: Found {len(candidates)} viable guides. "
                f"Top candidate expected to boost expression to {best['expected_expression']:.0%} "
                f"(from {config.baseline_expression:.0%} baseline). "
                f"Claim level: {overall_claim.value}."
            )
            status = TrialStatus.COMPLETED
        else:
            summary = (
                f"CRISPRa RAI1 Trial: No viable guides found that meet criteria. "
                f"Consider relaxing constraints or trying alternative targeting regions."
            )
            status = TrialStatus.COMPLETED

        return SMSTrialResult(
            trial_type=TrialType.CRISPRA_RAI1,
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
                'gene': 'RAI1',
                'modality': 'CRISPRa',
                'promoter_length': len(promoter_sequence),
                'config': {
                    'therapeutic_window': config.therapeutic_window,
                    'coherence_mode': config.coherence_mode,
                    'use_virtual_assay': config.use_virtual_assay,
                },
            },
        )

    except Exception as e:
        logger.error(f"CRISPRa trial failed: {e}")
        return SMSTrialResult(
            trial_type=TrialType.CRISPRA_RAI1,
            status=TrialStatus.FAILED,
            summary=f"CRISPRa RAI1 Trial failed: {str(e)}",
            errors=[str(e)],
            claim_level="unknown",
            claim_description="Trial execution failed.",
        )


def _claim_level_rank(level: str) -> int:
    """Convert claim level to numeric rank for comparison."""
    ranks = {
        'unknown': 0,
        'exploratory': 1,
        'context_dependent': 2,
        'strong_computational': 3,
    }
    return ranks.get(level.lower(), 0)
