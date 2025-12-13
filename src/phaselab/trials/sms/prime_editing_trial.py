"""
Prime Editing Trial for Smith-Magenis Syndrome.

This trial identifies pegRNA designs for precise correction of
RAI1 regulatory motifs or small insertions/deletions.

Use case: Patients with small indels or regulatory mutations that
cannot be corrected by single-base editors.

Note: Prime editing has larger payload requirements and may need
dual-AAV delivery, which affects CNS feasibility.
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
)

logger = logging.getLogger(__name__)


# Example edits that could be addressed by prime editing
RAI1_PRIME_EDIT_TARGETS = {
    'regulatory_insertion': {
        'description': 'Insert enhancer binding motif',
        'edit_position': 100,
        'edit_from': '',
        'edit_to': 'CACGTG',  # E-box motif
        'edit_type': 'insertion',
        'rationale': 'Enhance transcription factor binding',
    },
    'codon_correction': {
        'description': 'Correct premature stop codon',
        'edit_position': 150,
        'edit_from': 'TAG',
        'edit_to': 'CAG',  # Stop → Gln
        'edit_type': 'substitution',
        'rationale': 'Restore full-length protein',
    },
    'small_deletion_repair': {
        'description': 'Repair 3bp deletion causing frameshift',
        'edit_position': 200,
        'edit_from': '',
        'edit_to': 'GCT',
        'edit_type': 'insertion',
        'rationale': 'Restore reading frame',
    },
}

# Example sequence for prime editing targeting
EXAMPLE_PE_SEQUENCE = (
    "GCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
    "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
)


def run_sms_prime_editing_trial(
    edit_target: str = "codon_correction",
    target_sequence: Optional[str] = None,
    config: Optional[SMSTrialConfig] = None,
) -> SMSTrialResult:
    """
    Run prime editing trial for RAI1 precise correction.

    Identifies pegRNA designs for precise edits that cannot be
    accomplished with base editors alone.

    Args:
        edit_target: Target edit identifier.
        target_sequence: Sequence containing edit site. Uses example if None.
        config: Trial configuration.

    Returns:
        SMSTrialResult with pegRNA candidates.

    Example:
        >>> from phaselab.trials.sms import run_sms_prime_editing_trial
        >>> result = run_sms_prime_editing_trial(edit_target="codon_correction")
        >>> if result.best_candidate:
        ...     print(f"PBS length: {result.best_candidate['pbs_length']}")
        ...     print(f"RT length: {result.best_candidate['rt_length']}")
    """
    if config is None:
        config = SMSTrialConfig()

    # Get edit target info
    target_info = RAI1_PRIME_EDIT_TARGETS.get(edit_target, {
        'description': 'Custom edit',
        'edit_position': 100,
        'edit_from': 'A',
        'edit_to': 'G',
        'edit_type': 'substitution',
        'rationale': 'Unknown',
    })

    if target_sequence is None:
        target_sequence = EXAMPLE_PE_SEQUENCE

    edit_position = target_info['edit_position']
    edit_from = target_info['edit_from']
    edit_to = target_info['edit_to']
    edit_type = target_info['edit_type']

    if config.verbose:
        print("=" * 60)
        print(f"SMS Prime Editing Trial: {edit_target}")
        print("=" * 60)
        print(f"Description: {target_info['description']}")
        print(f"Edit type: {edit_type}")
        print(f"Edit: {edit_from or '[none]'} → {edit_to}")
        print(f"Rationale: {target_info['rationale']}")
        print()

    warnings = []
    errors = []

    # Delivery warning for prime editing
    warnings.append(
        "Prime editing requires larger payload (~6kb for PE2). "
        "May require dual-AAV delivery for CNS applications, "
        "which reduces efficiency."
    )

    try:
        # Import design functions
        from phaselab.crispr import design_prime_edit, PrimeEditConfig
        from phaselab.fusion import ClaimLevel

        # Design prime editing guides
        pe_config = PrimeEditConfig(
            top_n=config.top_n_guides * 2,
            edit_type=edit_type,
            check_secondary_structure=True,
        )

        raw_df = design_prime_edit(
            sequence=target_sequence,
            edit_position=edit_position,
            edit_from=edit_from,
            edit_to=edit_to,
            config=pe_config,
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

            # Get prime editing-specific scores
            pbs_length = guide.get('pbs_length', 13)
            rt_length = guide.get('rt_length', 15)
            pbs_gc = guide.get('pbs_gc', 0.5)
            rt_gc = guide.get('rt_gc', 0.5)
            nick_distance = guide.get('nick_to_edit_distance', 20)
            ss_dg = guide.get('extension_ss_dg', 0)

            # Quality assessment
            # Optimal: PBS 13-15bp, RT 10-16bp, nick 10-50bp from edit
            pbs_optimal = 0.8 if 13 <= pbs_length <= 15 else 0.5
            rt_optimal = 0.8 if 10 <= rt_length <= 16 else 0.5
            nick_optimal = 1.0 if 10 <= nick_distance <= 50 else 0.5
            ss_optimal = 1.0 if ss_dg > -5 else (0.6 if ss_dg > -10 else 0.3)

            pe_quality = pbs_optimal * rt_optimal * nick_optimal * ss_optimal

            # Classify quality
            quality_class = 'HIGH' if pe_quality > 0.5 else ('MODERATE' if pe_quality > 0.25 else 'LOW')

            candidate = {
                'spacer': guide.get('spacer', ''),
                'pam': guide.get('pam', 'NGG'),
                'strand': guide.get('strand', '+'),
                'nick_position': guide.get('nick_position', 0),
                'nick_to_edit_distance': nick_distance,

                # PBS info
                'pbs_sequence': guide.get('pbs_sequence', ''),
                'pbs_length': pbs_length,
                'pbs_gc': pbs_gc,

                # RT template info
                'rt_sequence': guide.get('rt_sequence', ''),
                'rt_length': rt_length,
                'rt_gc': rt_gc,

                # Secondary structure
                'extension_ss_dg': ss_dg,

                # Scores
                'coherence_R': coherence,
                'go_no_go': go_status,
                'combined_score': guide.get('combined_score', 0),

                # Quality assessment
                'pe_quality_score': pe_quality,
                'quality_class': quality_class,

                # Binding
                'delta_g': guide.get('delta_g', -15.0),
                'mit_score': guide.get('mit_score', 50),

                # Claim level
                'claim_level': 'exploratory',  # PE is newer technology
            }

            candidates.append(candidate)

        # Sort by PE quality
        candidates.sort(
            key=lambda x: (
                x['quality_class'] == 'HIGH',
                x['pe_quality_score'],
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
                print(f"  Spacer: {best['spacer']}")
                print(f"  PBS length: {best['pbs_length']}bp (optimal: 13-15)")
                print(f"  RT length: {best['rt_length']}bp (optimal: 10-16)")
                print(f"  Nick-to-edit: {best['nick_to_edit_distance']}bp")
                print(f"  Quality: {best['quality_class']}")

        # Determine overall claim level
        # Prime editing is newer, so we're more conservative
        if len(candidates) == 0:
            overall_claim = ClaimLevel.UNKNOWN
            claim_desc = f"No viable pegRNA designs found for {edit_target}."
        elif candidates[0]['quality_class'] == 'HIGH':
            overall_claim = ClaimLevel.EXPLORATORY
            claim_desc = (
                f"Exploratory evidence for {edit_target} prime edit correction. "
                f"High-quality pegRNA design identified, but prime editing efficiency "
                f"in vivo is variable and delivery-dependent."
            )
        else:
            overall_claim = ClaimLevel.EXPLORATORY
            claim_desc = (
                f"Exploratory evidence only. pegRNA designs available but "
                f"quality is moderate. Prime editing efficiency may be limited."
            )

        # Add technology maturity warning
        warnings.append(
            "Prime editing is less mature than Cas9 knockout or base editing. "
            "In vivo efficiency data is limited, especially for CNS targets."
        )

        # Build metrics
        metrics = {
            'n_initial_guides': len(raw_guides),
            'n_viable_candidates': len(candidates),
            'n_high_quality': sum(1 for c in candidates if c['quality_class'] == 'HIGH'),
            'edit_target': edit_target,
            'edit_type': edit_type,
            'edit_size': abs(len(edit_to) - len(edit_from)),
        }

        if candidates:
            metrics['best_pbs_length'] = candidates[0]['pbs_length']
            metrics['best_rt_length'] = candidates[0]['rt_length']

        # Build summary
        if len(candidates) > 0:
            best = candidates[0]
            summary = (
                f"Prime Editing {edit_target} Trial: Found {len(candidates)} pegRNA designs. "
                f"Best has PBS={best['pbs_length']}bp, RT={best['rt_length']}bp, "
                f"quality: {best['quality_class']}. "
                f"Note: Delivery may require dual-AAV approach."
            )
            status = TrialStatus.COMPLETED
        else:
            summary = (
                f"Prime Editing {edit_target} Trial: No viable pegRNA designs found. "
                f"Edit may be too far from available PAM sites or have structural issues."
            )
            status = TrialStatus.COMPLETED

        return SMSTrialResult(
            trial_type=TrialType.PRIME_EDITING,
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
                'edit_target': edit_target,
                'target_info': target_info,
                'modality': 'PrimeEditing',
                'technology_maturity': 'emerging',
                'payload_size': '~6kb (PE2)',
            },
        )

    except Exception as e:
        logger.error(f"Prime editing trial failed: {e}")
        return SMSTrialResult(
            trial_type=TrialType.PRIME_EDITING,
            status=TrialStatus.FAILED,
            summary=f"Prime Editing {edit_target} Trial failed: {str(e)}",
            errors=[str(e)],
            claim_level="unknown",
            claim_description="Trial execution failed.",
        )
