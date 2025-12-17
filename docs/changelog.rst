Changelog
=========

Version 1.0.0 (December 2025)
-----------------------------

**Spatial Coherence Paradigm - Major Methodology Shift**

This release represents a fundamental paradigm shift based on experimental validation:

**E200-E211 showed guide-sequence coherence DOES NOT WORK (r ≈ 0)**

Computing coherence from guide sequences (GC content, thermodynamic properties,
structural features) has no predictive value for experimental outcomes. This
approach has been **deprecated** in v1.0.0.

**E213-E216 validated SPATIAL coherence DOES WORK**

Measuring spatial coherence of response landscapes predicts perturbation reliability:

- Correlation: r = -0.24 to -0.50 with outcome variance
- Variance reduction: 32-49% when selecting from stable regions
- Validation scale: 115,251 sgRNAs across 6 genes

**The key insight**: "The guide is the probe, not the structure." IR coherence
measures the SYSTEM'S response consistency, not properties of the perturbation itself.

**New Modules**

1. **phaselab.landscapes** - Core perturbation-response data structures

   - ``ResponseLandscape`` - Generic position → response mapping
   - ``CoherenceProfile`` - Per-position coherence values with validation
   - ``StabilityClass`` - STABLE, MIXED, AMPLIFYING, IRRELEVANT
   - ``classify_regions()`` - Region classification algorithm

2. **phaselab.spatial** - E213-validated tiling screen analysis

   - ``analyze_tiling_coherence()`` - Full coherence analysis pipeline
   - ``load_tiling_screen()`` - Data loading utilities
   - ``TilingResult`` - Structured result with stable/amplifying regions

3. **phaselab.surf** - CRISPR-SURF integration

   - ``parse_surf_output()`` - Parse SURF deconvolution output
   - ``compute_surf_coherence()`` - Coherence on deconvolved data
   - ``SURFPipeline`` - End-to-end SURF + coherence pipeline
   - ``compare_raw_vs_surf()`` - Raw vs deconvolved comparison

4. **phaselab.omics** - Genomics assay coherence

   - ``analyze_atac_coherence()`` - ATAC-seq stable accessibility
   - ``analyze_chip_coherence()`` - ChIP-seq stable binding
   - ``analyze_expression_coherence()`` - RNA-seq reliable changes

5. **phaselab.microbio** - Microbial screen analysis

   - ``analyze_tnseq_coherence()`` - TnSeq essential domains
   - ``analyze_crispri_coherence()`` - Bacterial CRISPRi screens
   - ``analyze_drug_coherence()`` - Drug dose-response stability

6. **phaselab.chem** - Chemical/biochemical systems

   - ``analyze_binding_coherence()`` - Stable binding hot spots
   - ``analyze_reaction_coherence()`` - Stable reaction conditions
   - ``analyze_screening_coherence()`` - HTS reliable hits

7. **phaselab.protein.mutscan** - Mutational scanning analysis

   - ``analyze_mutscan_coherence()`` - Functional domain identification
   - ``local_coherence_profile()`` - Per-residue coherence
   - ``map_coherence_to_structure()`` - PDB B-factor mapping

**Quantum Mode Configuration**

New quantum execution modes:

.. code-block:: python

    from phaselab.quantum import QuantumMode, set_quantum_mode

    set_quantum_mode(QuantumMode.OFF)       # Classical only (fastest)
    set_quantum_mode(QuantumMode.AUDIT)     # Classical + quantum validation
    set_quantum_mode(QuantumMode.REQUIRED)  # Quantum mandatory (slowest)

**Breaking Changes**

- ``compute_coherence=True`` in CRISPR pipeline is **deprecated** (no effect)
- ``weight_coherence`` defaults to 0.0 (was 1.0)
- Guide-sequence coherence functions emit deprecation warnings
- SMS trials config now uses spatial coherence by default

**Upgrade Guide**

Replace guide-sequence coherence:

.. code-block:: python

    # OLD (deprecated):
    from phaselab.crispr import design_guides
    guides = design_guides(seq, tss, compute_coherence=True)

    # NEW (v1.0.0):
    from phaselab.spatial import analyze_tiling_coherence
    result = analyze_tiling_coherence(tiling_landscape)
    stable_positions = [r['start'] for r in result.stable_regions]


Version 0.9.5 (December 2025)
-----------------------------

**Quantum Discriminator for Late-Stage Guide Selection**

This release adds a quantum chemistry module for discriminating between elite
CRISPRa guides that are classically indistinguishable. Uses IBM Quantum hardware
or simulation to resolve binding energy differences.

**Core Claim**:

    *IR-enhanced quantum VQE on current IBM hardware can resolve binding energy
    differences between CRISPRa guides that are indistinguishable under classical
    scoring, providing a physically grounded late-stage discriminator for
    therapeutic guide selection.*

**New Components**:

1. **Effective Binding Hamiltonian**: H = H_HB + H_stack + H_charge + H_constraint

   - Watson-Crick hydrogen bonding (G-C: -0.18 eV, A-T: -0.12 eV)
   - π-π stacking stabilization
   - Backbone electrostatics with screening
   - 12-qubit seed region encoding

2. **Quantum VQE Execution**:

   - EfficientSU2 ansatz (2 reps, linear entanglement)
   - COBYLA optimizer with 30 max iterations
   - 1000 shots per measurement
   - Hardware support: ibm_torino, ibm_brisbane, etc.

3. **GO/NO-GO Threshold**: R̄ > e⁻² ≈ 0.135 for execution quality

**New API**:

.. code-block:: python

    from phaselab.crispr import (
        run_quantum_discriminator,
        design_guides_with_quantum_discriminator,
        DiscriminatorStatus,
        DISCRIMINATOR_GATES,
    )

    # Run discriminator on degenerate guides
    result = run_quantum_discriminator(
        guides=elite_guides,
        dna_context=promoter_sequence,
        use_hardware=False,  # True for IBM Quantum
    )

    print(result.summary())
    # Quantum-resolved ranking with energy separations

**Pre-Quantum Gates**:

- min_mit_score: 50
- max_exonic_ot: 0
- min_delta_r: 0.30
- min_phase_coherence: 0.90

**Status Codes**: QUANTUM_SUCCESS, NO_DEGENERACY, INSUFFICIENT_GUIDES, QUANTUM_FAILED


Version 0.9.4 (December 2025)
-----------------------------

**Three Breakthrough Paths for CRISPRa Guide Ranking**

This release adds three complementary scoring paths for CRISPRa guide selection:

1. **Path A: Binding Energy Landscape** - Quantum chemistry for relative binding
   energetics using effective Hamiltonians with Watson-Crick base pairing.

2. **Path B: Transcriptional Phase Alignment** - IR dynamics for phase perturbation
   modeling using vectorized Kuramoto oscillator simulation.

3. **Path C: Off-Target Landscape Geometry** - Coherence contrast between on-target
   and off-target binding for specificity scoring.

**Multi-Evidence Fusion**: Combines all three paths using weighted geometric mean
for unified guide ranking.

.. code-block:: python

    from phaselab.crispr import (
        compute_binding_energy,
        compute_phase_alignment,
        compute_offtarget_geometry,
        compute_multi_evidence_score,
    )

    # Combined scoring
    result = compute_multi_evidence_score(
        guide_sequence=guide,
        promoter_sequence=promoter,
        tss_position=tss,
        guide_position=guide_pos,
    )
    print(f"Combined score: {result.combined_score:.3f}")


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
