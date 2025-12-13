"""
SMS Pipeline Orchestrator: Complete preclinical decision engine.

This orchestrator runs all SMS trials in sequence, integrating results
across modalities to produce a comprehensive therapeutic assessment.

Pipeline stages:
1. CRISPRa RAI1 activation trial (core therapy)
2. Circadian rescue simulation (phenotype prediction)
3. Delivery feasibility assessment
4. (Optional) CRISPRi modifier trials
5. (Optional) Knockout validation trials
6. (Optional) Base/prime editing for specific variants
7. Validation study design recommendations

The pipeline produces:
- Ranked therapeutic strategies
- GO/NO-GO decision with claim levels
- Falsification test recommendations
- Next steps for wet lab validation
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional
from datetime import datetime
import logging

from .core import (
    SMSTrial,
    SMSTrialResult,
    SMSTrialConfig,
    TrialType,
    TrialStatus,
    RAI1_INFO,
    EXAMPLE_RAI1_PROMOTER,
)
from .crispra_trial import run_sms_crispra_trial
from .crispri_trial import run_sms_crispri_trial
from .knockout_trial import run_sms_knockout_trial
from .base_editing_trial import run_sms_base_editing_trial
from .prime_editing_trial import run_sms_prime_editing_trial
from .circadian_simulation import run_circadian_rescue_simulation
from .delivery_assessment import run_delivery_assessment

logger = logging.getLogger(__name__)


@dataclass
class SMSPipelineResult:
    """Complete SMS pipeline result."""

    # Status
    status: str  # "complete", "partial", "failed"
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())

    # Trial results
    crispra_result: Optional[SMSTrialResult] = None
    circadian_result: Optional[SMSTrialResult] = None
    delivery_result: Optional[SMSTrialResult] = None
    crispri_results: List[SMSTrialResult] = field(default_factory=list)
    knockout_results: List[SMSTrialResult] = field(default_factory=list)
    base_editing_results: List[SMSTrialResult] = field(default_factory=list)
    prime_editing_results: List[SMSTrialResult] = field(default_factory=list)

    # Integrated assessment
    overall_go_nogo: str = "UNKNOWN"
    overall_claim_level: str = "unknown"
    therapeutic_recommendation: str = ""

    # Falsification tests
    falsification_tests: List[Dict[str, Any]] = field(default_factory=list)

    # Next steps
    validation_priorities: List[str] = field(default_factory=list)
    wet_lab_recommendations: List[str] = field(default_factory=list)

    # Warnings and errors
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'status': self.status,
            'timestamp': self.timestamp,
            'overall_go_nogo': self.overall_go_nogo,
            'overall_claim_level': self.overall_claim_level,
            'therapeutic_recommendation': self.therapeutic_recommendation,
            'crispra_summary': self.crispra_result.summary if self.crispra_result else None,
            'circadian_summary': self.circadian_result.summary if self.circadian_result else None,
            'delivery_summary': self.delivery_result.summary if self.delivery_result else None,
            'falsification_tests': self.falsification_tests,
            'validation_priorities': self.validation_priorities,
            'wet_lab_recommendations': self.wet_lab_recommendations,
            'warnings': self.warnings,
            'errors': self.errors,
        }


class SMSPipeline:
    """
    Complete SMS preclinical decision engine.

    Orchestrates all trials and produces integrated therapeutic assessment.

    Example:
        >>> from phaselab.trials.sms import SMSPipeline
        >>> pipeline = SMSPipeline()
        >>> result = pipeline.run_full_pipeline()
        >>> print(f"GO/NO-GO: {result.overall_go_nogo}")
        >>> print(f"Claim level: {result.overall_claim_level}")
        >>> for rec in result.wet_lab_recommendations:
        ...     print(f"  - {rec}")
    """

    def __init__(self, config: Optional[SMSTrialConfig] = None):
        """Initialize pipeline with configuration."""
        self.config = config or SMSTrialConfig()
        self.results: Dict[str, SMSTrialResult] = {}

    def run_full_pipeline(
        self,
        promoter_sequence: Optional[str] = None,
        include_modifiers: bool = False,
        include_editing: bool = False,
    ) -> SMSPipelineResult:
        """
        Run complete SMS pipeline.

        Args:
            promoter_sequence: RAI1 promoter sequence.
            include_modifiers: Run CRISPRi modifier trials.
            include_editing: Run base/prime editing trials.

        Returns:
            Complete pipeline result.
        """
        if self.config.verbose:
            print("=" * 70)
            print("  SMS PRECLINICAL DECISION ENGINE")
            print("  PhaseLab v0.9.0")
            print("=" * 70)
            print()

        pipeline_result = SMSPipelineResult(status="running")

        # Stage 1: CRISPRa RAI1 activation
        if self.config.verbose:
            print("[1/7] Running CRISPRa RAI1 activation trial...")

        crispra_result = run_sms_crispra_trial(
            promoter_sequence=promoter_sequence,
            config=self.config,
        )
        pipeline_result.crispra_result = crispra_result

        # Stage 2: Circadian rescue simulation
        if crispra_result.is_successful and crispra_result.best_candidate:
            predicted_expression = crispra_result.best_candidate.get(
                'expected_expression', 0.80
            )
        else:
            predicted_expression = 0.80  # Default

        if self.config.verbose:
            print(f"[2/7] Running circadian rescue simulation (RAI1={predicted_expression:.0%})...")

        circadian_result = run_circadian_rescue_simulation(
            predicted_rai1_expression=predicted_expression,
            config=self.config,
        )
        pipeline_result.circadian_result = circadian_result

        # Stage 3: Delivery assessment
        if self.config.verbose:
            print("[3/7] Running delivery feasibility assessment...")

        delivery_result = run_delivery_assessment(
            modality="CRISPRa_VP64",
            target_tissue="brain",
            config=self.config,
        )
        pipeline_result.delivery_result = delivery_result

        # Stage 4: Optional CRISPRi modifier trials
        if include_modifiers:
            if self.config.verbose:
                print("[4/7] Running CRISPRi modifier trials...")

            for gene in ['PER1', 'CRY1']:
                result = run_sms_crispri_trial(
                    target_gene=gene,
                    config=self.config,
                )
                pipeline_result.crispri_results.append(result)
        else:
            if self.config.verbose:
                print("[4/7] Skipping CRISPRi modifier trials")

        # Stage 5: Optional knockout validation
        if include_modifiers:
            if self.config.verbose:
                print("[5/7] Running knockout validation trials...")

            ko_result = run_sms_knockout_trial(
                target_gene="PER1",
                config=self.config,
            )
            pipeline_result.knockout_results.append(ko_result)
        else:
            if self.config.verbose:
                print("[5/7] Skipping knockout validation trials")

        # Stage 6: Optional editing trials
        if include_editing:
            if self.config.verbose:
                print("[6/7] Running base/prime editing trials...")

            be_result = run_sms_base_editing_trial(
                variant_id="p.R1217Q",
                config=self.config,
            )
            pipeline_result.base_editing_results.append(be_result)

            pe_result = run_sms_prime_editing_trial(
                edit_target="codon_correction",
                config=self.config,
            )
            pipeline_result.prime_editing_results.append(pe_result)
        else:
            if self.config.verbose:
                print("[6/7] Skipping base/prime editing trials")

        # Stage 7: Integration and assessment
        if self.config.verbose:
            print("[7/7] Integrating results and generating recommendations...")

        self._integrate_results(pipeline_result)

        pipeline_result.status = "complete"

        if self.config.verbose:
            print()
            print("=" * 70)
            print("  PIPELINE COMPLETE")
            print("=" * 70)
            print(f"  Overall GO/NO-GO: {pipeline_result.overall_go_nogo}")
            print(f"  Claim Level: {pipeline_result.overall_claim_level}")
            print()
            print("  Recommendation:")
            print(f"    {pipeline_result.therapeutic_recommendation}")
            print()

        return pipeline_result

    def _integrate_results(self, result: SMSPipelineResult):
        """Integrate trial results into overall assessment."""

        # Determine overall GO/NO-GO
        go_criteria = []

        # CRISPRa success
        if result.crispra_result and result.crispra_result.has_viable_candidates:
            best = result.crispra_result.best_candidate
            if best.get('in_therapeutic_window', False):
                go_criteria.append(('crispra', 'GO', 'Guides in therapeutic window'))
            else:
                go_criteria.append(('crispra', 'CAUTION', 'Guides may not achieve therapeutic window'))
        else:
            go_criteria.append(('crispra', 'NO-GO', 'No viable CRISPRa guides'))

        # Circadian rescue
        if result.circadian_result:
            rescue_status = result.circadian_result.metrics.get('rescue_status', 'UNKNOWN')
            if rescue_status in ['FULL_RESCUE', 'PARTIAL_RESCUE']:
                go_criteria.append(('circadian', 'GO', f'{rescue_status} predicted'))
            else:
                go_criteria.append(('circadian', 'CAUTION', f'{rescue_status}'))

        # Delivery feasibility
        if result.delivery_result:
            feasibility = result.delivery_result.metrics.get('delivery_feasibility', 'UNKNOWN')
            if feasibility == 'FEASIBLE':
                go_criteria.append(('delivery', 'GO', 'Delivery feasible'))
            elif feasibility == 'CONDITIONALLY_FEASIBLE':
                go_criteria.append(('delivery', 'CAUTION', 'Delivery conditionally feasible'))
            else:
                go_criteria.append(('delivery', 'NO-GO', f'Delivery {feasibility}'))

        # Determine overall status
        no_gos = [c for c in go_criteria if c[1] == 'NO-GO']
        cautions = [c for c in go_criteria if c[1] == 'CAUTION']

        if no_gos:
            result.overall_go_nogo = 'NO-GO'
            result.warnings.extend([f"{c[0]}: {c[2]}" for c in no_gos])
        elif len(cautions) > 1:
            result.overall_go_nogo = 'CONDITIONAL'
            result.warnings.extend([f"{c[0]}: {c[2]}" for c in cautions])
        else:
            result.overall_go_nogo = 'GO'

        # Determine overall claim level
        claim_levels = []
        if result.crispra_result:
            claim_levels.append(result.crispra_result.claim_level)
        if result.circadian_result:
            claim_levels.append(result.circadian_result.claim_level)
        if result.delivery_result:
            claim_levels.append(result.delivery_result.claim_level)

        # Most conservative claim level
        level_order = ['unknown', 'exploratory', 'context_dependent', 'strong_computational']
        if claim_levels:
            result.overall_claim_level = min(claim_levels, key=lambda x: level_order.index(x))

        # Generate therapeutic recommendation
        if result.overall_go_nogo == 'GO':
            result.therapeutic_recommendation = self._generate_go_recommendation(result)
        elif result.overall_go_nogo == 'CONDITIONAL':
            result.therapeutic_recommendation = self._generate_conditional_recommendation(result)
        else:
            result.therapeutic_recommendation = self._generate_nogo_recommendation(result)

        # Generate falsification tests
        result.falsification_tests = self._generate_falsification_tests(result)

        # Generate validation priorities
        result.validation_priorities = self._generate_validation_priorities(result)

        # Generate wet lab recommendations
        result.wet_lab_recommendations = self._generate_wet_lab_recommendations(result)

    def _generate_go_recommendation(self, result: SMSPipelineResult) -> str:
        """Generate GO recommendation."""
        best_guide = result.crispra_result.best_candidate if result.crispra_result else None
        serotype = result.delivery_result.best_candidate.get('recommended_serotype') if result.delivery_result else 'AAV9'

        parts = [
            f"PROCEED with CRISPRa-based RAI1 activation therapy.",
        ]

        if best_guide:
            parts.append(
                f"Top guide ({best_guide['sequence'][:10]}...) predicted to boost RAI1 to "
                f"{best_guide.get('expected_expression', 0.8):.0%}, within therapeutic window."
            )

        parts.append(f"Recommended delivery: {serotype} via ICV/intrathecal route.")
        parts.append(f"Claim level: {result.overall_claim_level} - proceed with validation.")

        return " ".join(parts)

    def _generate_conditional_recommendation(self, result: SMSPipelineResult) -> str:
        """Generate conditional recommendation."""
        return (
            f"CONDITIONAL GO for SMS CRISPRa therapy. "
            f"Multiple cautions identified: {', '.join(result.warnings[:3])}. "
            f"Address these concerns before advancing. "
            f"Claim level: {result.overall_claim_level}."
        )

    def _generate_nogo_recommendation(self, result: SMSPipelineResult) -> str:
        """Generate NO-GO recommendation."""
        blockers = [w for w in result.warnings if 'NO-GO' in w or 'not' in w.lower()]
        return (
            f"NO-GO for current approach. Critical issues: {', '.join(blockers[:2])}. "
            f"Consider: (1) alternative targeting strategies, (2) dual-AAV delivery, "
            f"(3) smaller CRISPR variants (SaCas9), or (4) non-CRISPRa approaches."
        )

    def _generate_falsification_tests(self, result: SMSPipelineResult) -> List[Dict[str, Any]]:
        """Generate falsification tests per ChatGPT analysis."""
        tests = []

        # Test A: Ranking validity
        tests.append({
            'id': 'A',
            'name': 'Ranking Validity',
            'description': (
                'Compare top 10 PhaseLab-ranked guides against 10 randomly selected '
                'CRISPOR-acceptable guides in the same cell model.'
            ),
            'failure_condition': (
                'PhaseLab guides do NOT outperform random controls in '
                'activation efficiency, safety, or reproducibility.'
            ),
            'required_for': 'any_claim',
        })

        # Test B: Risk prediction validity
        if result.crispra_result and result.crispra_result.candidates:
            tests.append({
                'id': 'B',
                'name': 'Risk Prediction Validity',
                'description': (
                    'Test guides PhaseLab flags as high-risk (high risk mass, exonic tail risk) '
                    'in wet lab assays.'
                ),
                'failure_condition': (
                    'Flagged guides show no elevated off-target activity or toxicity.'
                ),
                'required_for': 'strong_computational_claim',
            })

        # Test C: Dosage prediction validity
        if result.circadian_result:
            predicted_window = result.circadian_result.metrics.get('therapeutic_window')
            tests.append({
                'id': 'C',
                'name': 'Dosage Prediction Validity',
                'description': (
                    f'Verify that RAI1 activation within predicted range '
                    f'{predicted_window if predicted_window else "(70-110%)"} '
                    f'rescues circadian phenotype without toxicity.'
                ),
                'failure_condition': (
                    'Either no phenotype rescue in predicted range, '
                    'or toxicity within predicted safe window.'
                ),
                'required_for': 'therapeutic_claim',
            })

        # Test D: UNKNOWN bucket sanity
        tests.append({
            'id': 'D',
            'name': 'UNKNOWN Bucket Sanity',
            'description': (
                'Test guides marked as UNKNOWN claim level anyway. '
                'Expect high variance and inconsistent behavior.'
            ),
            'failure_condition': (
                'UNKNOWN guides perform consistently well, indicating '
                'uncertainty logic is flawed.'
            ),
            'required_for': 'methodology_validation',
        })

        return tests

    def _generate_validation_priorities(self, result: SMSPipelineResult) -> List[str]:
        """Generate validation priorities."""
        priorities = []

        if result.crispra_result and result.crispra_result.best_candidate:
            priorities.append(
                f"1. CRITICAL: Validate top CRISPRa guide "
                f"({result.crispra_result.best_candidate['sequence'][:10]}...) "
                f"in iPSC-derived neurons"
            )

        priorities.append(
            "2. HIGH: Measure actual RAI1 upregulation vs predicted "
            f"({result.crispra_result.best_candidate.get('expected_fold_change', 1.7):.1f}x)"
            if result.crispra_result and result.crispra_result.best_candidate else
            "2. HIGH: Measure RAI1 upregulation in cell model"
        )

        priorities.append(
            "3. HIGH: Run off-target analysis (GUIDE-seq or similar)"
        )

        if result.circadian_result:
            priorities.append(
                "4. MEDIUM: Validate circadian phenotype rescue in cell model"
            )

        priorities.append(
            "5. MEDIUM: Run falsification test A (ranking validity)"
        )

        return priorities

    def _generate_wet_lab_recommendations(self, result: SMSPipelineResult) -> List[str]:
        """Generate wet lab recommendations."""
        recs = []

        # Cell model
        recs.append(
            "Use iPSC-derived neurons from SMS patients (preferable) "
            "or HEK293T with RAI1 reporter as screening model"
        )

        # Guide testing
        if result.crispra_result:
            n_guides = min(10, result.crispra_result.n_candidates)
            recs.append(
                f"Test top {n_guides} PhaseLab guides + {n_guides} random controls "
                f"(total: {n_guides * 2} conditions)"
            )

        # Readouts
        recs.append(
            "Primary readout: RAI1 mRNA (qPCR) and protein (Western) at 48-72h"
        )

        recs.append(
            "Secondary readout: Circadian reporters (PER2::Luc) if available"
        )

        # Safety
        recs.append(
            "Safety assessment: Off-target editing (GUIDE-seq), cell viability, "
            "transcriptome profiling for unexpected changes"
        )

        # Controls
        recs.append(
            "Include: non-targeting control, dCas9-only control, "
            "positive control (known active guide if available)"
        )

        return recs

    def run_crispra_only(
        self,
        promoter_sequence: Optional[str] = None,
    ) -> SMSTrialResult:
        """Run only CRISPRa trial (quick assessment)."""
        return run_sms_crispra_trial(
            promoter_sequence=promoter_sequence,
            config=self.config,
        )

    def run_circadian_only(
        self,
        predicted_expression: float = 0.80,
    ) -> SMSTrialResult:
        """Run only circadian simulation."""
        return run_circadian_rescue_simulation(
            predicted_rai1_expression=predicted_expression,
            config=self.config,
        )

    def run_delivery_only(
        self,
        modality: str = "CRISPRa_VP64",
    ) -> SMSTrialResult:
        """Run only delivery assessment."""
        return run_delivery_assessment(
            modality=modality,
            target_tissue="brain",
            config=self.config,
        )
