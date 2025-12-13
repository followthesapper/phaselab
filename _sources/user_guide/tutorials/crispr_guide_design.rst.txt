CRISPR Guide Design
===================

This tutorial covers designing guide RNAs with phase-coherence validation.

Overview
--------

PhaseLab provides a complete CRISPR toolkit:

- **CRISPRa**: Transcriptional activation
- **CRISPRi**: Transcriptional interference/repression
- **Knockout**: Gene disruption via DSB
- **Prime editing**: Precise edits via pegRNA
- **Base editing**: Single nucleotide changes (ABE/CBE)

CRISPRa Guide Design
--------------------

Basic Usage
^^^^^^^^^^^

.. code-block:: python

   from phaselab.crispr import design_guides

   # Target sequence
   sequence = """
   ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
   AGGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
   """

   tss_index = 50  # Transcription start site

   # Design guides
   guides = design_guides(sequence, tss_index)

   # View results
   print(guides[['sequence', 'position', 'gc_content', 'combined_score', 'go_no_go']])

Custom Configuration
^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from phaselab.crispr import design_guides, GuideDesignConfig

   config = GuideDesignConfig(
       pam="NGG",                    # SpCas9 PAM
       guide_length=20,              # 20bp guide
       window_upstream=500,          # Search 500bp upstream of TSS
       window_downstream=100,        # Search 100bp downstream
       min_gc=0.40,                  # Minimum GC content
       max_gc=0.70,                  # Maximum GC content
       coherence_mode="heuristic",   # Fast screening
   )

   guides = design_guides(sequence, tss_index, config=config)

CRISPRi Guide Design
--------------------

For gene repression:

.. code-block:: python

   from phaselab.crispr import design_crispri_guides, CRISPRiConfig

   config = CRISPRiConfig(
       pam="NGG",
       guide_length=20,
       window_tss=-50,     # Optimal: 50bp downstream of TSS
       window_end=300,     # To +300bp
       min_gc=0.40,
       max_gc=0.70,
   )

   crispri_guides = design_crispri_guides(sequence, tss_index, config=config)

   # Includes repression efficiency and steric hindrance scores
   print(crispri_guides[['sequence', 'position', 'repression_score', 'steric_score']])

Knockout Guide Design
---------------------

For gene disruption:

.. code-block:: python

   from phaselab.crispr import design_knockout_guides, KnockoutConfig

   config = KnockoutConfig(
       pam="NGG",
       guide_length=20,
       target_exon=True,           # Target exonic regions
       avoid_first_exon=True,      # Skip first exon (may have alternative start)
       frameshift_required=True,   # Prioritize frameshift-inducing sites
   )

   ko_guides = design_knockout_guides(
       sequence,
       exon_start=100,
       exon_end=250,
       config=config
   )

   # Includes cut efficiency and frameshift probability
   print(ko_guides[['sequence', 'cut_efficiency', 'frameshift_prob', 'repair_pathway']])

Prime Editing
-------------

For precise edits:

.. code-block:: python

   from phaselab.crispr import design_prime_edit, PrimeEditConfig

   config = PrimeEditConfig(
       pbs_length_range=(10, 15),      # Primer binding site length
       rt_template_range=(10, 20),     # RT template length
       nick_offset_range=(-3, 3),      # Nick site offset
   )

   # Design pegRNA for a substitution
   pegRNA = design_prime_edit(
       sequence,
       edit_position=150,
       edit_type="substitution",
       new_base="G",
       config=config
   )

   print(f"Spacer: {pegRNA['spacer']}")
   print(f"PBS: {pegRNA['pbs']}")
   print(f"RT template: {pegRNA['rt_template']}")
   print(f"Full pegRNA: {pegRNA['full_sequence']}")

Base Editing
------------

For single nucleotide changes:

.. code-block:: python

   from phaselab.crispr import design_base_edit_guides, BaseEditConfig

   # ABE (A-to-G editing)
   abe_config = BaseEditConfig(
       editor_type="ABE",              # Adenine base editor
       activity_window=(4, 8),         # Positions 4-8 from PAM-distal end
       bystander_tolerance=0,          # No bystanders allowed
   )

   abe_guides = design_base_edit_guides(
       sequence,
       target_position=150,
       target_base="A",
       config=abe_config
   )

   # CBE (C-to-T editing)
   cbe_config = BaseEditConfig(
       editor_type="CBE",              # Cytosine base editor
       activity_window=(4, 8),
   )

   cbe_guides = design_base_edit_guides(
       sequence,
       target_position=150,
       target_base="C",
       config=cbe_config
   )

   # Check bystanders
   from phaselab.crispr import find_bystanders
   bystanders = find_bystanders(abe_guides.iloc[0], sequence)
   print(f"Bystander edits: {bystanders}")

Multi-Guide Sets
----------------

For combinatorial CRISPR:

.. code-block:: python

   from phaselab.crispr import design_multiguide_set, MultiGuideConfig

   config = MultiGuideConfig(
       n_guides=3,                     # Number of guides
       min_spacing=50,                 # Minimum bp between guides
       max_spacing=500,                # Maximum bp between guides
       require_synergy=True,           # Prioritize synergistic pairs
   )

   multiguide = design_multiguide_set(sequence, tss_index, config=config)

   print(f"Guides: {[g['sequence'] for g in multiguide['guides']]}")
   print(f"Synergy score: {multiguide['synergy_score']:.3f}")
   print(f"Combined efficacy: {multiguide['combined_efficacy']:.3f}")

Enhancer Targeting
------------------

For enhancer-based activation:

.. code-block:: python

   from phaselab.crispr import design_enhancer_guides, EnhancerConfig

   config = EnhancerConfig(
       enhancer_regions=[(1000, 1500), (2000, 2500)],  # Known enhancer coordinates
       min_h3k27ac=0.5,               # Minimum H3K27ac signal
       prefer_super_enhancers=True,   # Prioritize super-enhancers
   )

   enhancer_guides = design_enhancer_guides(sequence, config=config)

   # Compare promoter vs enhancer targeting
   from phaselab.crispr import compare_promoter_vs_enhancer
   comparison = compare_promoter_vs_enhancer(
       promoter_guides=guides,
       enhancer_guides=enhancer_guides
   )
   print(f"Promoter efficacy: {comparison['promoter_efficacy']:.3f}")
   print(f"Enhancer efficacy: {comparison['enhancer_efficacy']:.3f}")

Scoring Functions
-----------------

Individual scoring components:

.. code-block:: python

   from phaselab.crispr import (
       gc_content,
       delta_g_santalucia,
       mit_specificity_score,
       cfd_score,
       max_homopolymer_run,
       chromatin_accessibility_score,
   )

   guide = "GAAGGAGAGCAAGAGCGCGA"

   # GC content (optimal: 40-70%)
   gc = gc_content(guide)
   print(f"GC: {gc:.2%}")

   # Thermodynamic stability (SantaLucia)
   dg = delta_g_santalucia(guide)
   print(f"deltaG: {dg:.1f} kcal/mol")

   # MIT specificity score
   mit = mit_specificity_score(guide, offtargets=[...])
   print(f"MIT: {mit:.3f}")

   # CFD score
   cfd = cfd_score(guide, offtarget="GAAGGAGAGCAAGAGCGCGG")
   print(f"CFD: {cfd:.3f}")

   # Homopolymer runs (avoid >4)
   homo = max_homopolymer_run(guide)
   print(f"Max homopolymer: {homo}")

   # Chromatin accessibility
   access = chromatin_accessibility_score(position=150, dnase_peaks=[(100, 200)])
   print(f"Accessibility: {access:.3f}")

Complete Pipeline Example
-------------------------

.. code-block:: python

   from phaselab.crispr import (
       design_guides,
       GuideDesignConfig,
       compute_coherence_batch,
   )
   from phaselab.integrations.crispor import analyze_offtarget_landscape

   # RAI1 promoter sequence
   rai1_promoter = """GCGCGCTCGCGCGCTCGCGCGAAGGAGAGCAAGAGCGCGACGGCTA..."""

   # Step 1: Design with heuristic coherence (fast)
   config = GuideDesignConfig(
       pam="NGG",
       window_upstream=400,
       window_downstream=50,
       coherence_mode="heuristic",
   )

   guides = design_guides(rai1_promoter, tss_index=len(rai1_promoter)//2, config=config)
   print(f"Total candidates: {len(guides)}")

   # Step 2: Filter by basic criteria
   filtered = guides[
       (guides['go_no_go'] == 'GO') &
       (guides['gc_content'].between(0.4, 0.7)) &
       (guides['homopolymer'] <= 4)
   ]
   print(f"After filtering: {len(filtered)}")

   # Step 3: Quantum validation of top 10
   top_10 = filtered.sort_values('combined_score', ascending=False).head(10)

   quantum_r = compute_coherence_batch(
       top_10['sequence'].tolist(),
       mode="quantum"
   )
   top_10['quantum_coherence'] = quantum_r

   # Step 4: Off-target analysis
   for idx, row in top_10.iterrows():
       result = analyze_offtarget_landscape(
           guide_seq=row['sequence'],
           coherence_mode="quantum"
       )
       top_10.loc[idx, 'ir_score'] = result.ir_enhanced_score
       top_10.loc[idx, 'evidence_level'] = result.evidence_level

   # Step 5: Final ranking
   final = top_10.sort_values('ir_score', ascending=False)
   print(final[['sequence', 'quantum_coherence', 'ir_score', 'evidence_level']])

Next Steps
----------

- :doc:`coherence_modes` - Understand heuristic vs quantum
- :doc:`quantum_validation` - Hardware validation
