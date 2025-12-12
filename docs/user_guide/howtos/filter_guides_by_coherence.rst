Filter Guides by Coherence
==========================

This guide shows how to filter CRISPR guides using coherence metrics and evidence levels.

Basic Coherence Filtering
-------------------------

Filter guides that pass the GO/NO-GO threshold:

.. code-block:: python

   from phaselab.crispr import design_guides

   # Design guides
   guides = design_guides(sequence, tss_index)

   # Filter by GO status
   go_guides = guides[guides['go_no_go'] == 'GO']

   print(f"Total: {len(guides)}, GO: {len(go_guides)}")

Filtering with Evidence Levels
------------------------------

v0.6.1 assigns evidence levels based on validation:

.. code-block:: python

   from phaselab.crispr import design_guides, compute_coherence_batch

   # Step 1: Design with heuristic (Level C)
   guides = design_guides(sequence, tss_index)
   go_guides = guides[guides['go_no_go'] == 'GO']

   # Step 2: Upgrade top candidates to Level B (quantum simulation)
   top_candidates = go_guides.head(20)
   quantum_r = compute_coherence_batch(
       top_candidates['sequence'].tolist(),
       mode="quantum"
   )
   top_candidates['quantum_coherence'] = quantum_r
   top_candidates['evidence_level'] = 'B'

   # Step 3: Filter by evidence level
   level_b_guides = top_candidates[top_candidates['evidence_level'] == 'B']

Multi-Criteria Filtering
------------------------

Combine coherence with other criteria:

.. code-block:: python

   from phaselab.crispr import design_guides
   from phaselab.core.constants import E_MINUS_2

   guides = design_guides(sequence, tss_index)

   # Multi-criteria filter
   filtered = guides[
       # Coherence criteria
       (guides['go_no_go'] == 'GO') &
       (guides['coherence'] > 0.5) &  # Above moderate threshold

       # Specificity criteria
       (guides['offtarget_score'] > 0.7) &

       # Sequence criteria
       (guides['gc_content'].between(0.40, 0.65)) &
       (guides['homopolymer'] <= 4)
   ]

   # Sort by coherence
   filtered = filtered.sort_values('coherence', ascending=False)

Coherence Z-Score Filtering
---------------------------

Use z-scores for relative ranking within a locus:

.. code-block:: python

   from phaselab.crispr import compute_coherence_with_zscore

   # Get coherence with z-scores
   results = compute_coherence_with_zscore(
       guides['sequence'].tolist(),
       mode="heuristic"
   )

   # Add to dataframe
   guides['coherence'], guides['z_score'] = zip(*results)

   # Filter by z-score (above average coherence)
   above_avg = guides[guides['z_score'] > 0]

   # Top performers (>1 std above mean)
   top_performers = guides[guides['z_score'] > 1.0]

Tiered Filtering Pipeline
-------------------------

Use a tiered approach for large candidate sets:

.. code-block:: python

   from phaselab.crispr import (
       design_guides,
       compute_coherence_batch,
       is_guide_coherent,
   )

   # Tier 1: Initial design (thousands of candidates)
   all_guides = design_guides(sequence, tss_index)
   print(f"Tier 1: {len(all_guides)} candidates")

   # Tier 2: Basic GO filter
   tier2 = all_guides[all_guides['go_no_go'] == 'GO']
   print(f"Tier 2 (GO): {len(tier2)} candidates")

   # Tier 3: Top 100 by combined score
   tier3 = tier2.nlargest(100, 'combined_score')
   print(f"Tier 3 (top 100): {len(tier3)} candidates")

   # Tier 4: Quantum validation of top 20
   tier4 = tier3.head(20).copy()
   tier4['quantum_coherence'] = compute_coherence_batch(
       tier4['sequence'].tolist(),
       mode="quantum"
   )
   print(f"Tier 4 (quantum validated): {len(tier4)} candidates")

   # Final: Best 5 by quantum coherence
   final = tier4.nlargest(5, 'quantum_coherence')
   print(f"Final candidates: {len(final)}")

Custom Threshold Filtering
--------------------------

Use custom coherence thresholds:

.. code-block:: python

   from phaselab.crispr import design_guides, is_guide_coherent

   guides = design_guides(sequence, tss_index)

   # Custom high threshold
   HIGH_THRESHOLD = 0.8

   # Check each guide
   guides['high_coherence'] = [
       is_guide_coherent(seq, threshold=HIGH_THRESHOLD)
       for seq in guides['sequence']
   ]

   # Filter
   high_coherence_guides = guides[guides['high_coherence']]

See Also
--------

- :doc:`compare_crispr_modalities` - Compare modality approaches
- :doc:`/user_guide/tutorials/coherence_modes` - Coherence mode details
