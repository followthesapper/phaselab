Examples
========

Working code examples for common use cases.

.. toctree::
   :maxdepth: 1

   sms_gene_therapy
   cancer_crispra
   chronotherapy

SMS Gene Therapy Pipeline
-------------------------

Complete pipeline for Smith-Magenis Syndrome RAI1 gene therapy design.

.. code-block:: python

   from phaselab.crispr import design_guides, compute_coherence_batch
   from phaselab.circadian import simulate_sms_clock

   # Design RAI1-targeting guides
   rai1_promoter = """GCGCGCTCGCGCGCTCGCGCGAAGGAGAGCAAGAGCGCGACGGCTA..."""

   guides = design_guides(rai1_promoter, tss_index=200)
   go_guides = guides[guides['go_no_go'] == 'GO']

   # Quantum validation
   quantum_r = compute_coherence_batch(
       go_guides.head(10)['sequence'].tolist(),
       mode="quantum"
   )

   # Predict therapeutic outcomes
   for guide, r in zip(go_guides.head(10)['sequence'], quantum_r):
       predicted_rai1 = 0.5 + 0.35 * r  # Efficacy model
       sim = simulate_sms_clock(rai1_level=predicted_rai1)
       print(f"{guide[:10]}... R={r:.3f} -> Period={sim['period']:.1f}h")

Cancer CRISPRa Workflow
-----------------------

Tumor suppressor reactivation via CRISPRa.

.. code-block:: python

   from phaselab.crispr import design_guides, design_enhancer_guides

   # TP53 reactivation
   tp53_promoter = """..."""

   # Promoter-targeting guides
   promoter_guides = design_guides(tp53_promoter, tss_index=500)

   # Enhancer-targeting guides
   enhancer_guides = design_enhancer_guides(
       tp53_promoter,
       enhancer_regions=[(1000, 1500)],
   )

   # Compare approaches
   print(f"Promoter guides: {len(promoter_guides)}")
   print(f"Enhancer guides: {len(enhancer_guides)}")

Chronotherapy Optimization
--------------------------

Optimize drug dosing schedules using circadian modeling.

.. code-block:: python

   from phaselab.circadian import optimize_dosing_schedule

   schedule = optimize_dosing_schedule(
       target_gene='RAI1',
       delivery_method='AAV',
       tissue_targets=['SCN', 'liver'],
   )

   print(f"Optimal times: {schedule['dose_times']}")
   print(f"Predicted coherence: {schedule['predicted_coherence']:.4f}")
