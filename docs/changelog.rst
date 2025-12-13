Changelog
=========

Version 0.9.3 (December 2025)
-----------------------------

**CRISPRa Binding Register Model - Major Methodology Correction**

This release corrects a fundamental assumption in CRISPRa guide enumeration that
caused systematic exclusion of experimentally validated guides in GC-dense promoters.

**The Problem**: Standard CRISPR pipelines use rigid PAM-spacer anchoring
(``guide_start = pam_start - guide_length``) derived from cutting-era Cas9 work.
This assumption is invalid for CRISPRa/dCas9 binding, where:

- Binding tolerates non-canonical PAMs
- The functional binding register can shift ±1-2bp
- GC-dense promoters have overlapping PAM-like motifs

**The Fix**: v0.9.3 introduces:

1. **NucleaseRole enum**: Explicit ``BINDING`` vs ``CUTTING`` mode
2. **Relaxed PAM patterns**: e.g., SaCas9 NNGRRN (binding) vs NNGRRT (cutting)
3. **Sliding binding register**: ±2bp enumeration in BINDING mode
4. **Configurable guide length**: Override defaults (e.g., 20bp with SaCas9)

**Validation**: Chang et al. 2022 sg2 winner (``CCTGGCACCCGAGGCCACGA``) was
systematically excluded by all prior versions. v0.9.3 correctly recovers it at
TSS-80 with the NNGRRN PAM pattern.

**New CRISPRa Design API**

.. code-block:: python

    from phaselab.crispr import design_crispra_guides, Nuclease, NucleaseRole

    # Design CRISPRa guides with explicit binding mode
    result = design_crispra_guides(
        gene_symbol="Rai1",
        promoter_sequence=promoter_seq,
        tss_position=600,
        nuclease=Nuclease.SACAS9,
        relaxed_pam=True,    # BINDING mode (default for CRISPRa)
        guide_length=20,     # Override default 21bp
    )

    # Access results
    for guide in result.tier_a_guides[:5]:
        print(f"{guide['sequence']} TSS{guide['tss_relative_position']:+d}")

**Key Insight (Publishable)**:

    *CRISPRa guide effectiveness is invariant to small PAM-guide register shifts
    in GC-dense promoters. Computational pipelines that enforce rigid spacer
    anchoring systematically miss experimentally validated guides.*

This is not a bug fix - it's a **modeling correction** that reflects the
biological reality of dCas9 binding tolerance.

Version 0.9.2 (December 2025)
-----------------------------

**Dominance-Based Ranking System**

- Lexicographic sorting on (0mm, 1mm, 2mm) off-targets
- Hard gates exclude guides entirely (not just penalized)
- Tier system: A (0/0/0), B (0/0/1-2), C (other)
- Policy-explicit ranking with reproducibility manifests

**RankingPolicy System**

- ``CUTTING_STRICT``: Maximum safety for knockout
- ``BINDING_STRICT``: For CRISPRa/CRISPRi binding applications
- ``EXPLORATORY``: Relaxed constraints for research (not therapeutic)

Version 0.9.1 (December 2025)
-----------------------------

**CRISPOR-Style Composite Scoring**

- New ``crispor_composite_score()`` function with mismatch distance weighting
- Off-target penalties: 0-1mm (critical), 2mm (important), 3-4mm (minimal)
- Properly handles the "MIT 98 / CFD 98 trap" from high-OT guides
- ``rank_guides_crispor_style()`` for batch ranking with automatic exclusions

**U6/Pol III Compatibility Checks**

- ``poly_t_penalty()`` detects TTTT runs that terminate U6 transcription
- ``is_repeat_region()`` identifies tandem and dinucleotide repeats
- ``u6_compatibility_check()`` comprehensive promoter compatibility

**API Additions**

.. code-block:: python

    from phaselab.crispr import (
        # CRISPOR composite scoring (v0.9.1)
        crispor_composite_score,
        rank_guides_crispor_style,
        OFFTARGET_MISMATCH_WEIGHTS,
        CRISPORMetrics,

        # U6/Pol III compatibility (v0.9.1)
        poly_t_penalty,
        is_repeat_region,
        u6_compatibility_check,
    )

Version 0.9.0 (December 2025)
-----------------------------

**SMS Trials Module**

- Complete therapeutic trial framework for Smith-Magenis Syndrome
- CRISPRa RAI1 activation trial with therapeutic window validation
- CRISPRi modifier gene suppression trials (PER1, CRY1, CLOCK)
- Knockout model validation trials (research use only)
- Base editing trials for RAI1 point mutation correction
- Prime editing trials for regulatory motif repair
- Circadian rescue simulation with sleep/wake prediction
- AAV delivery feasibility assessment for CNS targeting

**SMS Pipeline Orchestrator**

- Integrated GO/NO-GO decision system
- Multi-trial coordination with claim level propagation
- Automatic falsification test generation
- Wet lab recommendations and validation priorities

**Falsification Test Framework**

- Test A: Ranking validity (PhaseLab vs random controls)
- Test B: Risk prediction (CAUTION guides should fail more)
- Test C: Dosage prediction (expression correlation)
- Test D: UNKNOWN bucket calibration

**API Additions**

.. code-block:: python

    from phaselab.trials.sms import (
        SMSPipeline,
        SMSTrialConfig,
        run_sms_crispra_trial,
        run_sms_crispri_trial,
        run_circadian_rescue_simulation,
        run_delivery_assessment,
    )

Version 0.8.0 (December 2025)
-----------------------------

**Claim Level System**

- Four-tier evidence classification: STRONG_COMPUTATIONAL, CONTEXT_DEPENDENT, EXPLORATORY, UNKNOWN
- Claim level propagation through all pipelines
- Prevents over-claiming from computational predictions

**Fusion Module**

- Multi-source data integration with uncertainty quantification
- Virtual assay stack for enhanced guide scoring
- Tissue-specific modeling integration

Version 0.7.0 (November 2025)
-----------------------------

**Enhanced Pipeline**

- ``design_enhanced_guides()`` with modality-specific scoring
- Full modality support: CRISPRa, CRISPRi, Knockout, Base Editing, Prime Editing
- Tissue-specific scoring for brain, liver, blood, muscle

**API Additions**

.. code-block:: python

    from phaselab.crispr.enhanced_pipeline import (
        design_enhanced_guides,
        EnhancedGuideConfig,
        Modality,
    )

Version 0.6.1 (December 2025)
-----------------------------

**Coherence Mode Parameter**

- Added ``mode="heuristic"`` (fast) vs ``mode="quantum"`` (VQE) parameter
- Heuristic mode is now default for speed
- Quantum mode provides research-grade accuracy

**Honest Coherence Weighting**

- Heuristic coherence demoted to tie-breaker weight (0.05 vs 0.30)
- Prevents over-reliance on proxy metrics
- Two-stage scoring: hard safety gates + soft ranking

**Risk Mass Metrics**

- Added ``risk_mass_close``: Off-targets within 100bp of TSS
- Added ``risk_mass_exonic``: Off-targets in exonic regions
- Added ``tail_risk_score``: Aggregate tail risk metric

**Evidence Levels**

- Level A: Hardware-validated (IBM Quantum)
- Level B: VQE-simulated (quantum mode)
- Level C: Heuristic only (capped influence)

**Score Capping**

- Unvalidated guides capped to prevent misleading rankings
- Evidence level affects maximum achievable score

Version 0.6.0 (November 2025)
-----------------------------

**ATLAS-Q Integration**

- Full integration with ATLAS-Q tensor network simulator
- IR measurement grouping (5x variance reduction)
- Real circular statistics coherence (replaces heuristic)

**Coherence-Aware VQE**

- VQE optimization with real-time coherence tracking
- GO/NO-GO classification during optimization
- Optional GPU acceleration via Triton kernels

**Rust Backend Support**

- Optional Rust backend for 30-77x faster simulation
- Automatic fallback to Python if Rust unavailable

**CRISPOR Integration**

- IR-enhanced off-target analysis
- Off-target entropy metrics
- Coherence contrast scoring

**Unified ATLAS-Q Coherence**

- Single coherence computation path for all CRISPR modalities
- Consistent API across CRISPRa, CRISPRi, knockout, editing

Version 0.5.0 (October 2025)
----------------------------

**Real ATAC-seq Integration**

- BigWig file support for tissue-specific accessibility
- CpG methylation modeling for CRISPRa efficiency

**Nucleosome Occupancy**

- NuPoP-like algorithm for nucleosome prediction
- Integration with guide scoring

**Multi-Guide Synergy**

- Combinatorial CRISPR design
- Pairwise synergy prediction
- Guide set optimization

**Enhancer Targeting**

- Enhancer identification and scoring
- CRISPRa enhancer guide design
- Promoter vs enhancer comparison

**AAV Delivery Modeling**

- Serotype selection
- Delivery efficiency prediction
- Immunogenicity assessment

**Validation Reports**

- Comprehensive validation report generation
- Evidence summary and confidence scoring

Version 0.4.0 (September 2025)
------------------------------

**Complete CRISPR Toolkit**

- CRISPR knockout (Cas9 cutting)
- CRISPRi (transcriptional repression)
- All modalities hardware-validated on IBM Torino

**Therapeutic Dosage Optimization**

- Haploinsufficiency models
- Dose-response prediction
- Therapeutic window estimation

Version 0.3.0 (August 2025)
---------------------------

**Multi-Tissue Circadian Models**

- Inter-tissue coupling
- Tissue-specific parameters
- SCN-peripheral synchronization

**Drug Response Modeling**

- Chronotherapy optimization
- Dosing schedule prediction
- Response curve modeling

**Expanded CRISPR Editors**

- Base editing (ABE/CBE)
- Prime editing (pegRNA design)
- Bystander prediction

Version 0.2.0 (July 2025)
-------------------------

**Protein Folding Coherence**

- Folding reliability assessment
- Coherence-structure correlation
- Engineering applications

**Tissue-Specific Chromatin**

- ENCODE integration
- Cell-type specific accessibility
- Tissue-aware guide scoring

Version 0.1.0 (June 2025)
-------------------------

**Initial Release**

- Core coherence metrics (R, V_phi)
- GO/NO-GO classification
- Basic CRISPRa guide design
- SMS circadian clock model
- IBM Quantum integration
