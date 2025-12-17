Select Guides from Stable Regions
==================================

This guide shows how to use spatial coherence to select CRISPR guides
from reliable targeting regions.

Problem
-------

You have tiling screen data and want to design guides that will give
reproducible results.

Solution
--------

1. Analyze spatial coherence to identify stable regions
2. Design guides restricted to those regions
3. Annotate guides with local coherence

Step 1: Analyze Your Tiling Data
--------------------------------

.. code-block:: python

   from phaselab.spatial import load_tiling_screen, analyze_tiling_coherence

   # Load your tiling screen
   landscape = load_tiling_screen(
       'tiling_screen.tsv',
       position_col='position',
       response_col='log2fc',
   )

   # Analyze coherence
   result = analyze_tiling_coherence(
       landscape,
       window=50,
       stable_threshold=0.7,
   )

   # Check if validated
   if not result.is_validated:
       print("WARNING: Coherence not validated (p > 0.05)")
       print("Results may not be reliable")

Step 2: Get Stable Region Boundaries
------------------------------------

.. code-block:: python

   # Extract stable regions
   stable_positions = set()
   for region in result.stable_regions:
       for pos in range(region['start'], region['end'] + 1):
           stable_positions.add(pos)

   print(f"Found {len(result.stable_regions)} stable regions")
   print(f"Total stable positions: {len(stable_positions)}")

Step 3: Design Guides in Stable Regions
---------------------------------------

.. code-block:: python

   from phaselab.crispr import design_crispra_guides

   # Design guides restricted to stable regions
   guides = design_crispra_guides(
       promoter_sequence=promoter_seq,
       tss_position=tss,
   )

   # Filter to stable regions
   stable_guides = []
   for guide in guides.candidates:
       pos = guide['tss_relative_position']
       # Check if guide position is in a stable region
       if any(r['start'] <= pos <= r['end'] for r in result.stable_regions):
           guide['region_status'] = 'STABLE'
           guide['local_coherence'] = result.get_coherence_at(pos)
           stable_guides.append(guide)

   print(f"Guides in stable regions: {len(stable_guides)}")

Step 4: Rank by Coherence-Weighted Score
----------------------------------------

.. code-block:: python

   # Sort by local coherence (secondary to standard scoring)
   stable_guides.sort(
       key=lambda g: (g['combined_score'], g['local_coherence']),
       reverse=True
   )

   # View top candidates
   print("\nTop guides from stable regions:")
   for i, guide in enumerate(stable_guides[:5], 1):
       print(f"{i}. {guide['sequence']}")
       print(f"   Position: TSS{guide['tss_relative_position']:+d}")
       print(f"   Score: {guide['combined_score']:.3f}")
       print(f"   Local coherence: {guide['local_coherence']:.3f}")

Alternative: Direct Position Filtering
--------------------------------------

.. code-block:: python

   from phaselab.crispr import design_crispra_guides

   # Design with position filter
   guides = design_crispra_guides(
       promoter_sequence=promoter_seq,
       tss_position=tss,
       position_filter=lambda pos: any(
           r['start'] <= pos <= r['end']
           for r in result.stable_regions
       ),
   )

Complete Example
----------------

.. code-block:: python

   """Complete workflow: tiling data → stable region guides"""

   from phaselab.spatial import load_tiling_screen, analyze_tiling_coherence
   from phaselab.crispr import design_crispra_guides

   # 1. Load and analyze tiling data
   landscape = load_tiling_screen('my_screen.tsv', 'position', 'log2fc')
   coherence = analyze_tiling_coherence(landscape)

   print(f"Stable regions: {len(coherence.stable_regions)}")
   print(f"Amplifying regions: {len(coherence.amplifying_regions)}")

   # 2. Design guides
   guides = design_crispra_guides(promoter_seq, tss=500)

   # 3. Filter and annotate
   final_guides = []
   for guide in guides.candidates:
       pos = guide['tss_relative_position']
       coh = coherence.get_coherence_at(pos)

       if coh >= 0.7:
           guide['coherence'] = coh
           guide['stability'] = 'STABLE'
           final_guides.append(guide)
       elif coh >= 0.4:
           guide['coherence'] = coh
           guide['stability'] = 'MIXED'
           # Include but flag as caution
           final_guides.append(guide)

   # 4. Output
   print(f"\nGuides from stable/mixed regions: {len(final_guides)}")
   for g in final_guides[:3]:
       print(f"  {g['sequence']} [{g['stability']}] coh={g['coherence']:.2f}")

See Also
--------

- :doc:`/user_guide/tutorials/spatial_coherence` - Understanding spatial coherence
- :doc:`filter_guides_by_coherence` - Additional filtering options

