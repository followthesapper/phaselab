Coherence Modes (v0.6.1)
========================

PhaseLab v0.6.1 introduces honest coherence computation with two distinct modes. Understanding when to use each mode is critical for reliable results.

Overview
--------

============  =============  ===========  ================  ===================
Mode          Speed          R Range      ATLAS-Q Benefit   Use Case
============  =============  ===========  ================  ===================
heuristic     ~0.1ms/guide   0.68-0.69    None              Screening, tie-breaking
quantum       ~100-500ms     0.84-0.97    Yes               Research, validation
============  =============  ===========  ================  ===================

Heuristic Mode (Default)
------------------------

The heuristic mode uses Hamiltonian coefficient variance as a proxy for coherence. It is **fast but not true quantum coherence**.

.. code-block:: python

   from phaselab.crispr import compute_guide_coherence

   guide = "ATCGATCGATCGATCGATCG"

   # Heuristic mode (default)
   r_bar = compute_guide_coherence(guide, mode="heuristic")
   print(f"Heuristic R: {r_bar:.4f}")  # ~0.68-0.69

**When to use heuristic mode:**

- Initial screening of thousands of guides
- Tie-breaking between guides with similar specificity
- When speed is critical
- Development and testing

**Important limitations:**

- Does NOT benefit from ATLAS-Q acceleration
- R values are compressed (~0.68-0.69 range)
- Cannot distinguish truly excellent guides from merely good ones
- Should NOT be used as primary ranking signal

Quantum Mode
------------

The quantum mode runs actual VQE simulation on the gRNA Hamiltonian. Results match IBM hardware validation.

.. code-block:: python

   from phaselab.crispr import compute_guide_coherence

   guide = "ATCGATCGATCGATCGATCG"

   # Quantum mode
   r_bar = compute_guide_coherence(guide, mode="quantum")
   print(f"Quantum R: {r_bar:.4f}")  # ~0.84-0.97

**When to use quantum mode:**

- Final validation of top candidates
- Research-grade analysis
- When accuracy matters more than speed
- Comparing with hardware validation results

**Benefits:**

- ATLAS-Q acceleration (if available)
- Full R range with good discrimination
- Results match IBM Quantum hardware
- Suitable for publication-quality results

Checking Method Details
-----------------------

Use ``get_coherence_eligibility_info()`` to understand what will be computed:

.. code-block:: python

   from phaselab.crispr import get_coherence_eligibility_info

   # Check heuristic mode
   info_h = get_coherence_eligibility_info(mode="heuristic")
   print(f"Mode: {info_h['mode']}")
   print(f"Method: {info_h['method']}")
   print(f"ATLAS-Q active: {info_h['acceleration_active']}")
   print(f"Expected R range: {info_h['expected_r_bar_range']}")
   print(f"Time per guide: {info_h['expected_time_per_guide']}")

   # Check quantum mode
   info_q = get_coherence_eligibility_info(mode="quantum")
   print(f"Mode: {info_q['mode']}")
   print(f"ATLAS-Q active: {info_q['acceleration_active']}")

Full Results with Details
-------------------------

Get complete coherence results including the method used:

.. code-block:: python

   from phaselab.crispr import compute_guide_coherence_with_details

   guide = "ATCGATCGATCGATCGATCG"

   # Get full details
   r_bar, v_phi, is_go, method = compute_guide_coherence_with_details(
       guide, mode="quantum"
   )

   print(f"R: {r_bar:.4f}")
   print(f"V_phi: {v_phi:.4f}")
   print(f"GO status: {is_go}")
   print(f"Method: {method}")  # "quantum_atlas_q" or "quantum_native"

Batch Processing
----------------

Process multiple guides efficiently:

.. code-block:: python

   from phaselab.crispr import compute_coherence_batch

   guides = [
       "ATCGATCGATCGATCGATCG",
       "GCTAGCTAGCTAGCTAGCTA",
       "TACGATCGATCGATCGATCG",
   ]

   # Heuristic batch (fast)
   r_bars_h = compute_coherence_batch(guides, mode="heuristic")

   # Quantum batch (slow, use for final candidates only)
   r_bars_q = compute_coherence_batch(guides, mode="quantum")

   for guide, r_h, r_q in zip(guides, r_bars_h, r_bars_q):
       print(f"{guide[:10]}... H:{r_h:.3f} Q:{r_q:.3f}")

Z-Score for Relative Ranking
----------------------------

When guides have similar R values, use z-scores for ranking:

.. code-block:: python

   from phaselab.crispr import compute_coherence_with_zscore

   guides = [
       "ATCGATCGATCGATCGATCG",
       "GCTAGCTAGCTAGCTAGCTA",
       "TACGATCGATCGATCGATCG",
       "AGTCAGTCAGTCAGTCAGTC",
   ]

   # Get R and z-score
   results = compute_coherence_with_zscore(guides, mode="heuristic")

   for guide, (r_bar, zscore) in zip(guides, results):
       print(f"{guide[:10]}... R:{r_bar:.4f} z:{zscore:+.2f}")

Z-scores are useful when:

- All guides have similar absolute R values
- You need relative ranking within a locus
- GC-rich regions where R values cluster

Evidence Levels
---------------

PhaseLab assigns evidence levels based on validation:

=========  =============================  ================
Level      Description                    Score Weight
=========  =============================  ================
A          IBM hardware validated         Full weight
B          VQE quantum simulation         Full weight
C          Heuristic only                 Capped (0.05)
=========  =============================  ================

.. code-block:: python

   # Level C guides have capped influence
   # This prevents heuristic-only guides from ranking
   # above properly validated candidates

Best Practices
--------------

**Screening Pipeline:**

.. code-block:: python

   from phaselab.crispr import (
       design_guides,
       compute_coherence_batch,
   )

   # Step 1: Initial design with heuristic mode (fast)
   guides = design_guides(sequence, tss_index)

   # Step 2: Filter by safety and specificity
   safe_guides = guides[
       (guides['offtarget_score'] > 0.5) &
       (guides['gc_content'].between(0.4, 0.7))
   ]

   # Step 3: Quantum validation of top candidates only
   top_10 = safe_guides.head(10)
   quantum_r = compute_coherence_batch(
       top_10['sequence'].tolist(),
       mode="quantum"
   )

   # Step 4: Final ranking with quantum coherence
   top_10['quantum_coherence'] = quantum_r
   final = top_10.sort_values('quantum_coherence', ascending=False)

**Publication-Quality Analysis:**

.. code-block:: python

   from phaselab.crispr import compute_guide_coherence

   # Always use quantum mode for publication
   result = compute_guide_coherence(
       guide,
       mode="quantum",
       return_full_result=True
   )

   print(f"Method: {result.method}")  # Verify method for reproducibility

Next Steps
----------

- :doc:`quantum_validation` - Hardware validation on IBM Quantum
- :doc:`crispr_guide_design` - Complete guide design workflow
