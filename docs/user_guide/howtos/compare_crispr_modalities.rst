Compare CRISPR Modalities
=========================

This guide helps you choose between CRISPR modalities for your application.

Modality Overview
-----------------

============  ============================  ================  ==============
Modality      Use Case                      Reversible?       Complexity
============  ============================  ================  ==============
CRISPRa       Gene activation               Yes               Low
CRISPRi       Gene repression               Yes               Low
Knockout      Gene disruption               No                Low
Prime edit    Precise edits                 No                High
Base edit     Single nucleotide changes     No                Medium
============  ============================  ================  ==============

Quick Decision Guide
--------------------

**Want to increase gene expression?**
  Use CRISPRa

**Want to decrease gene expression?**
  Use CRISPRi (reversible) or Knockout (permanent)

**Need a specific mutation?**
  - Single base change: Base editing
  - Larger change: Prime editing

**Gene therapy for haploinsufficiency?**
  CRISPRa (restore expression from remaining allele)

Comparing CRISPRa vs CRISPRi
----------------------------

.. code-block:: python

   from phaselab.crispr import design_guides, design_crispri_guides

   # Same target, different approaches
   sequence = "..."  # Target gene promoter
   tss_index = 200

   # CRISPRa: Activate gene
   crispra = design_guides(sequence, tss_index)
   print(f"CRISPRa candidates: {len(crispra)}")
   print(f"Best CRISPRa: {crispra.iloc[0]['combined_score']:.3f}")

   # CRISPRi: Repress gene
   crispri = design_crispri_guides(sequence, tss_index)
   print(f"CRISPRi candidates: {len(crispri)}")
   print(f"Best CRISPRi: {crispri.iloc[0]['combined_score']:.3f}")

Knockout vs CRISPRi
-------------------

.. code-block:: python

   from phaselab.crispr import design_crispri_guides, design_knockout_guides

   # CRISPRi: Reversible repression
   crispri = design_crispri_guides(sequence, tss_index)

   # Knockout: Permanent disruption
   knockout = design_knockout_guides(
       sequence,
       exon_start=500,
       exon_end=800
   )

   print("CRISPRi advantages:")
   print("- Reversible (remove dCas9 to restore)")
   print("- No DNA damage")
   print("- Tunable (VP64/KRAB fusion strength)")

   print("Knockout advantages:")
   print("- Complete loss-of-function")
   print("- Permanent effect")
   print("- Well-established protocols")

Base Editing vs Prime Editing
-----------------------------

.. code-block:: python

   from phaselab.crispr import design_base_edit_guides, design_prime_edit

   target_position = 150
   sequence = "..."

   # Base editing: Single base change (A->G or C->T)
   # Check if target base is editable
   target_base = sequence[target_position]

   if target_base == 'A':
       # ABE can convert A->G
       abe_guides = design_base_edit_guides(
           sequence,
           target_position=target_position,
           target_base='A',
           editor_type='ABE'
       )
       print(f"ABE candidates: {len(abe_guides)}")

   elif target_base == 'C':
       # CBE can convert C->T
       cbe_guides = design_base_edit_guides(
           sequence,
           target_position=target_position,
           target_base='C',
           editor_type='CBE'
       )
       print(f"CBE candidates: {len(cbe_guides)}")

   else:
       print("Base editing not applicable for this base")

   # Prime editing: Any edit type
   prime = design_prime_edit(
       sequence,
       edit_position=target_position,
       edit_type='substitution',
       new_base='G'
   )
   print(f"Prime editing pegRNA designed")

Decision Matrix
---------------

.. code-block:: python

   def recommend_modality(
       goal: str,
       reversible: bool = False,
       specific_mutation: str = None
   ) -> str:
       """Recommend CRISPR modality based on requirements."""

       if goal == "increase_expression":
           return "CRISPRa"

       elif goal == "decrease_expression":
           return "CRISPRi" if reversible else "Knockout"

       elif goal == "specific_mutation":
           if specific_mutation in ['A>G', 'T>C']:
               return "ABE (base editing)"
           elif specific_mutation in ['C>T', 'G>A']:
               return "CBE (base editing)"
           else:
               return "Prime editing"

       elif goal == "gene_disruption":
           return "Knockout"

       else:
           return "Consult literature for your specific case"

   # Examples
   print(recommend_modality("increase_expression"))  # CRISPRa
   print(recommend_modality("decrease_expression", reversible=True))  # CRISPRi
   print(recommend_modality("specific_mutation", specific_mutation="A>G"))  # ABE

Full Comparison Pipeline
------------------------

.. code-block:: python

   from phaselab.crispr import (
       design_guides,
       design_crispri_guides,
       design_knockout_guides,
       compute_coherence_batch,
   )
   import pandas as pd

   # Design for all modalities
   crispra = design_guides(sequence, tss_index).head(5)
   crispra['modality'] = 'CRISPRa'

   crispri = design_crispri_guides(sequence, tss_index).head(5)
   crispri['modality'] = 'CRISPRi'

   knockout = design_knockout_guides(sequence, 500, 800).head(5)
   knockout['modality'] = 'Knockout'

   # Combine and compare
   all_guides = pd.concat([crispra, crispri, knockout], ignore_index=True)

   # Add quantum coherence
   all_guides['quantum_coherence'] = compute_coherence_batch(
       all_guides['sequence'].tolist(),
       mode="quantum"
   )

   # Compare by modality
   comparison = all_guides.groupby('modality').agg({
       'combined_score': 'mean',
       'quantum_coherence': 'mean',
       'gc_content': 'mean'
   })
   print(comparison)

See Also
--------

- :doc:`filter_guides_by_coherence` - Filter by coherence
- :doc:`/user_guide/tutorials/crispr_guide_design` - Full design tutorial
