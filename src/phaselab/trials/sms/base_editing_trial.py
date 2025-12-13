"""
Base Editing Trial for Smith-Magenis Syndrome Variant Correction.

This trial identifies ABE or CBE guides to correct hypomorphic
RAI1 variants that cause partial loss of function.

Use case: Patients with RAI1 point mutations (not deletions) may
benefit from targeted base correction to restore RAI1 function.

Note: This is only applicable to specific SMS patients with
correctable point mutations, not the 17p11.2 deletion majority.
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


# Known RAI1 pathogenic variants that might be correctable with base editing
# In practice, these would come from ClinVar/HGMD databases
RAI1_CORRECTABLE_VARIANTS = {
    'p.R1217Q': {
        'cdna': 'c.3650G>A',
        'consequence': 'missense',
        'position_in_cds': 3650,
        'ref_base': 'G',
        'alt_base': 'A',
        'correction_editor': 'ABE',  # A→G to restore
        'clinical_significance': 'pathogenic',
        'frequency': 'rare',
    },
    'p.Q1562R': {
        'cdna': 'c.4685A>G',
        'consequence': 'missense',
        'position_in_cds': 4685,
        'ref_base': 'A',
        'alt_base': 'G',
        'correction_editor': 'CBE',  # Can't directly fix; example only
        'clinical_significance': 'pathogenic',
        'frequency': 'rare',
    },
    'p.S1808N': {
        'cdna': 'c.5423G>A',
        'consequence': 'missense',
        'position_in_cds': 5423,
        'ref_base': 'G',
        'alt_base': 'A',
        'correction_editor': 'ABE',  # A→G to restore
        'clinical_significance': 'pathogenic',
        'frequency': 'rare',
    },
}

# Example sequence around a variant site
EXAMPLE_VARIANT_CONTEXT = (
    "GCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
)


def run_sms_base_editing_trial(
    variant_id: str = "p.R1217Q",
    variant_context_sequence: Optional[str] = None,
    variant_position_in_seq: Optional[int] = None,
    config: Optional[SMSTrialConfig] = None,
) -> SMSTrialResult:
    """
    Run base editing trial for RAI1 variant correction.

    Identifies ABE or CBE guides to correct specific pathogenic
    RAI1 variants in patients with point mutations.

    Args:
        variant_id: Variant identifier (e.g., "p.R1217Q").
        variant_context_sequence: Sequence around variant. Uses example if None.
        variant_position_in_seq: Position of variant base in sequence.
        config: Trial configuration.

    Returns:
        SMSTrialResult with base editing guide candidates.

    Example:
        >>> from phaselab.trials.sms import run_sms_base_editing_trial
        >>> result = run_sms_base_editing_trial(variant_id="p.R1217Q")
        >>> if result.best_candidate:
        ...     print(f"Editor: {result.metadata['editor']}")
        ...     print(f"Position efficiency: {result.best_candidate['position_efficiency']:.0%}")
    """
    if config is None:
        config = SMSTrialConfig()

    # Get variant info
    variant_info = RAI1_CORRECTABLE_VARIANTS.get(variant_id, {
        'cdna': 'unknown',
        'consequence': 'unknown',
        'position_in_cds': 0,
        'ref_base': 'N',
        'alt_base': 'N',
        'correction_editor': 'ABE',
        'clinical_significance': 'unknown',
    })

    editor_type = variant_info.get('correction_editor', 'ABE')
    target_base = 'A' if editor_type == 'ABE' else 'C'

    if variant_context_sequence is None:
        variant_context_sequence = EXAMPLE_VARIANT_CONTEXT

    if variant_position_in_seq is None:
        variant_position_in_seq = len(variant_context_sequence) // 2

    if config.verbose:
        print("=" * 60)
        print(f"SMS Base Editing Trial: {variant_id}")
        print("=" * 60)
        print(f"Variant: {variant_info['cdna']} ({variant_info['consequence']})")
        print(f"Editor: {editor_type} ({target_base}→{'G' if target_base == 'A' else 'T'})")
        print(f"Clinical significance: {variant_info['clinical_significance']}")
        print()

    warnings = []
    errors = []

    # Applicability warning
    warnings.append(
        f"Base editing is only applicable to SMS patients with {variant_id} mutation. "
        f"Not suitable for 17p11.2 deletion cases (~90% of SMS patients)."
    )

    try:
        # Import design functions
        from phaselab.crispr import (
            design_base_edit_guides,
            BaseEditConfig,
            design_abe_guides,
            design_cbe_guides,
        )
        from phaselab.fusion import ClaimLevel

        # Design base editing guides
        be_config = BaseEditConfig(
            editor="ABE8e" if editor_type == "ABE" else "BE4",
            target_base=target_base,
            top_n=config.top_n_guides * 2,
            check_bystanders=True,
            max_bystanders_in_window=2,
        )

        raw_df = design_base_edit_guides(
            sequence=variant_context_sequence,
            target_position=variant_position_in_seq,
            target_base=target_base,
            config=be_config,
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

            # Get base editing-specific scores
            position_eff = guide.get('position_efficiency', 0.5)
            context_score = guide.get('context_score', 1.0)
            combined_eff = guide.get('combined_efficiency', position_eff * context_score)
            n_bystanders = guide.get('n_bystanders', 0)

            # Quality assessment
            be_quality = combined_eff * (1 - 0.1 * n_bystanders)  # Penalize bystanders

            # Correction potential assessment
            # High efficiency in activity window = likely to correct variant
            correction_potential = 'HIGH' if combined_eff > 0.7 else ('MODERATE' if combined_eff > 0.4 else 'LOW')

            candidate = {
                'sequence': guide.get('sequence', ''),
                'pam': guide.get('pam', 'NGG'),
                'strand': guide.get('strand', '+'),
                'target_in_window_pos': guide.get('target_in_window_pos', 5),

                # Scores
                'coherence_R': coherence,
                'go_no_go': go_status,
                'position_efficiency': position_eff,
                'context_score': context_score,
                'combined_efficiency': combined_eff,

                # Bystander info
                'n_bystanders': n_bystanders,
                'bystander_positions': guide.get('bystander_positions', []),

                # Quality assessment
                'be_quality_score': be_quality,
                'correction_potential': correction_potential,

                # Binding
                'delta_g': guide.get('delta_g', -15.0),
                'mit_score': guide.get('mit_score', 50),

                # Claim level
                'claim_level': 'context_dependent',
            }

            candidates.append(candidate)

        # Sort by base editing quality
        candidates.sort(
            key=lambda x: (
                x['correction_potential'] == 'HIGH',
                x['be_quality_score'],
                -x['n_bystanders'],  # Fewer bystanders preferred
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
                print(f"  Target position in window: {best['target_in_window_pos']}")
                print(f"  Combined efficiency: {best['combined_efficiency']:.0%}")
                print(f"  Bystanders: {best['n_bystanders']}")
                print(f"  Correction potential: {best['correction_potential']}")

        # Determine overall claim level
        if len(candidates) == 0:
            overall_claim = ClaimLevel.UNKNOWN
            claim_desc = f"No viable {editor_type} guides found for {variant_id} correction."
        elif candidates[0]['correction_potential'] == 'HIGH' and candidates[0]['n_bystanders'] <= 1:
            overall_claim = ClaimLevel.CONTEXT_DEPENDENT
            claim_desc = (
                f"Context-dependent evidence for {variant_id} correction with {editor_type}. "
                f"High-efficiency guides with minimal bystander editing."
            )
        elif candidates[0]['correction_potential'] in ['HIGH', 'MODERATE']:
            overall_claim = ClaimLevel.EXPLORATORY
            claim_desc = (
                f"Exploratory evidence for {variant_id} correction. "
                f"Guides available but bystander editing or efficiency concerns."
            )
        else:
            overall_claim = ClaimLevel.EXPLORATORY
            claim_desc = (
                f"Exploratory evidence only. {editor_type} correction possible but "
                f"efficiency may be limited for {variant_id}."
            )

        # Bystander warning
        if candidates and candidates[0]['n_bystanders'] > 0:
            bystander_warning = (
                f"Top guide has {candidates[0]['n_bystanders']} potential bystander edit(s) "
                f"at positions {candidates[0]['bystander_positions']}. "
                f"Verify these don't affect RAI1 function."
            )
            warnings.append(bystander_warning)

        # Build metrics
        metrics = {
            'n_initial_guides': len(raw_guides),
            'n_viable_candidates': len(candidates),
            'n_high_potential': sum(1 for c in candidates if c['correction_potential'] == 'HIGH'),
            'variant_id': variant_id,
            'editor_type': editor_type,
            'target_base': target_base,
        }

        if candidates:
            metrics['best_combined_efficiency'] = candidates[0]['combined_efficiency']
            metrics['best_correction_potential'] = candidates[0]['correction_potential']

        # Build summary
        if len(candidates) > 0:
            best = candidates[0]
            summary = (
                f"Base Editing {variant_id} Trial: Found {len(candidates)} {editor_type} guides. "
                f"Top guide has {best['combined_efficiency']:.0%} efficiency with "
                f"{best['n_bystanders']} bystander(s). "
                f"Correction potential: {best['correction_potential']}."
            )
            status = TrialStatus.COMPLETED
        else:
            summary = (
                f"Base Editing {variant_id} Trial: No viable {editor_type} guides found. "
                f"Variant may not be in accessible activity window or PAM context."
            )
            status = TrialStatus.COMPLETED

        return SMSTrialResult(
            trial_type=TrialType.BASE_EDITING,
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
                'variant_id': variant_id,
                'variant_info': variant_info,
                'editor': editor_type,
                'modality': 'BaseEditing',
                'applicability': 'point_mutations_only',
            },
        )

    except Exception as e:
        logger.error(f"Base editing trial failed: {e}")
        return SMSTrialResult(
            trial_type=TrialType.BASE_EDITING,
            status=TrialStatus.FAILED,
            summary=f"Base Editing {variant_id} Trial failed: {str(e)}",
            errors=[str(e)],
            claim_level="unknown",
            claim_description="Trial execution failed.",
        )
