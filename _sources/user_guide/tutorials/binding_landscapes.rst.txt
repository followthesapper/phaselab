Binding Landscape Analysis
==========================

Analyze protein-ligand and protein-protein binding data to identify
stable binding determinants and hot spots.

Overview
--------

Binding landscapes map mutations to affinity changes:

- **Position**: Residue or modification site
- **Response**: ΔΔG, IC50 change, or binding score

Spatial coherence identifies regions where binding effects are consistent,
revealing true binding hot spots vs. noise.

Basic Example
-------------

.. code-block:: python

   from phaselab.chem import (
       BindingLandscape,
       analyze_binding_coherence,
   )

   # Your binding data (e.g., alanine scan)
   landscape = BindingLandscape(
       positions=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
       affinities=[0.1, 0.2, 2.5, 2.8, 3.0, 2.7, 0.3, 0.1, -0.1, 0.0],
       target="ABL1",
       ligand="Imatinib",
       affinity_type="ddG",  # kcal/mol
   )

   # Analyze coherence
   result = analyze_binding_coherence(landscape)

   # View results
   print(f"Hot spots identified: {len(result.hot_spots)}")
   print(f"Validated: {result.is_validated}")

Finding Binding Hot Spots
-------------------------

.. code-block:: python

   # Hot spots: large ΔΔG + high coherence
   for hotspot in result.hot_spots:
       print(f"Hot spot at position {hotspot['position']}")
       print(f"  ΔΔG: {hotspot['effect']:.2f} kcal/mol")
       print(f"  Coherence: {hotspot['coherence']:.3f}")
       print(f"  Reliable: {'YES' if hotspot['coherence'] > 0.7 else 'NO'}")

   # Stable binding regions
   for region in result.stable_regions:
       print(f"Stable binding region: {region['start']}-{region['end']}")

IC50 and Other Metrics
----------------------

.. code-block:: python

   # Works with any affinity metric
   landscape = BindingLandscape(
       positions=mutation_positions,
       affinities=fold_change_ic50,  # Log2 fold-change
       target="EGFR",
       ligand="Compound_X",
       affinity_type="log2_ic50_fc",
   )

   result = analyze_binding_coherence(
       landscape,
       effect_threshold=1.0,  # >2-fold change
   )

Drug Discovery Application
--------------------------

.. code-block:: python

   from phaselab.chem import (
       load_binding_data,
       identify_druggable_sites,
   )

   # Load binding data
   landscape = load_binding_data(
       'binding_scan.csv',
       position_col='residue',
       binding_col='delta_delta_g',
   )

   # Analyze
   result = analyze_binding_coherence(landscape)

   # Find druggable regions
   druggable = identify_druggable_sites(
       result,
       min_coherence=0.7,      # High confidence
       min_effect=2.0,          # Significant binding contribution
   )

   for site in druggable:
       print(f"Druggable site: {site['region']}")
       print(f"  Binding contribution: {site['total_ddg']:.1f} kcal/mol")
       print(f"  Confidence: {site['confidence']}")

Reaction Optimization
---------------------

The chemistry module also supports reaction condition optimization:

.. code-block:: python

   from phaselab.chem import (
       ReactionLandscape,
       analyze_reaction_coherence,
   )

   # Reaction screen data
   landscape = ReactionLandscape(
       conditions=temperatures,  # Or pH, concentration, etc.
       yields=yield_values,
       reaction="Suzuki_coupling",
   )

   result = analyze_reaction_coherence(landscape)

   # Find stable operating window
   for window in result.stable_windows:
       print(f"Stable range: {window['min']}-{window['max']}")
       print(f"  Expected yield: {window['mean_yield']:.0f}%")
       print(f"  Variance: ±{window['std']:.0f}%")

HTS Screening Analysis
----------------------

Identify reliable hits from high-throughput screens:

.. code-block:: python

   from phaselab.chem import (
       ScreeningLandscape,
       analyze_screening_coherence,
   )

   # HTS data
   landscape = ScreeningLandscape(
       compounds=compound_ids,
       activities=inhibition_values,
       assay="kinase_inhibition",
   )

   result = analyze_screening_coherence(landscape)

   # Get reliable hits (active + coherent)
   reliable_hits = result.reliable_hits
   print(f"Reliable hits: {len(reliable_hits)} / {len(result.all_hits)}")

   for hit in reliable_hits[:10]:
       print(f"{hit['compound']}: {hit['activity']:.1f}% (coh={hit['coherence']:.2f})")

Interpretation
--------------

.. list-table:: Binding Coherence Interpretation
   :header-rows: 1

   * - Coherence
     - Effect Size
     - Interpretation
   * - High (>0.7)
     - Large ΔΔG
     - Reliable binding hot spot
   * - High (>0.7)
     - Small ΔΔG
     - Stable but minor contribution
   * - Low (<0.4)
     - Large ΔΔG
     - Artifact or context-dependent
   * - Low (<0.4)
     - Small ΔΔG
     - Not a binding determinant

Next Steps
----------

- :doc:`/user_guide/howtos/integrate_with_rosetta` - Combine with structure-based scoring
- :doc:`/user_guide/explanations/claim_levels` - Understanding confidence levels

