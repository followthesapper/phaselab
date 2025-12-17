Spatial Coherence Analysis (v1.0.0)
===================================

This tutorial covers the core v1.0.0 methodology: **spatial coherence of response landscapes**.

.. note::

   This is the validated methodology. Guide-sequence coherence (pre-v1.0.0) has been
   **deprecated** based on experimental validation showing r ≈ 0 correlation with outcomes.

The Core Insight
----------------

When you perturb a biological system at different positions, you get a **response landscape** -
a mapping from perturbation position to measured effect.

**Spatial coherence** measures how smoothly effects change across positions:

- **High coherence**: Nearby positions give similar effects → reliable, reproducible
- **Low coherence**: Effects jump unpredictably → unreliable, high variance

This was validated across 115,251 sgRNAs (6 genes):

- Correlation with outcome variance: r = -0.24 to -0.50
- Variance reduction in stable regions: 32-49%

Basic Example
-------------

.. code-block:: python

   from phaselab.landscapes import ResponseLandscape
   from phaselab.spatial import analyze_tiling_coherence

   # Your tiling screen data: positions and measured effects
   positions = [100, 110, 120, 130, 140, 150, 160, 170, 180, 190]
   effects = [0.8, 0.75, 0.82, 0.78, 0.1, 0.95, 0.12, 0.88, 0.85, 0.9]

   # Create a response landscape
   landscape = ResponseLandscape(
       coords=positions,
       responses=effects,
       metadata={'gene': 'MYC', 'assay': 'CRISPRa'}
   )

   # Analyze spatial coherence
   result = analyze_tiling_coherence(landscape)

   # Print summary
   print(result.summary())

Understanding the Output
------------------------

The analysis returns a ``TilingResult`` with:

.. code-block:: python

   # Validation metrics
   print(f"Correlation: {result.profile.correlation:.3f}")
   print(f"P-value: {result.profile.p_value:.4f}")
   print(f"Validated: {result.is_validated}")

   # Stable regions (good for targeting)
   for region in result.stable_regions:
       print(f"Stable: {region['start']}-{region['end']}")
       print(f"  Coherence: {region['coherence']:.3f}")
       print(f"  Mean effect: {region['mean_response']:.3f}")

   # Amplifying regions (avoid these)
   for region in result.amplifying_regions:
       print(f"Amplifying: {region['start']}-{region['end']}")
       print(f"  WARNING: High variance zone")

Stability Classes
-----------------

PhaseLab classifies regions into four categories:

.. list-table:: Stability Classification
   :header-rows: 1

   * - Class
     - Coherence
     - Meaning
     - Recommendation
   * - **STABLE**
     - > 0.7
     - Low variance, reproducible
     - Safe to target
   * - **MIXED**
     - 0.4 - 0.7
     - Moderate variance
     - Use with caution
   * - **AMPLIFYING**
     - < 0.4
     - High variance, unpredictable
     - Avoid
   * - **IRRELEVANT**
     - N/A
     - Below response threshold
     - Skip (no measurable effect)

Loading Real Data
-----------------

.. code-block:: python

   from phaselab.spatial import load_tiling_screen

   # From TSV file
   landscape = load_tiling_screen(
       'my_tiling_screen.tsv',
       position_col='tss_distance',      # Column with positions
       response_col='log2fc',            # Column with responses
       gene_symbol='RAI1',               # Optional metadata
   )

   # From pandas DataFrame
   import pandas as pd
   df = pd.read_csv('my_data.csv')

   landscape = ResponseLandscape(
       coords=df['position'].values,
       responses=df['effect'].values,
   )

Configuring Analysis
--------------------

.. code-block:: python

   result = analyze_tiling_coherence(
       landscape,
       window=50,                  # Positions for local coherence (default: 10)
       stable_threshold=0.7,      # Coherence threshold for STABLE (default: 0.7)
       amplifying_threshold=0.4,  # Coherence threshold for AMPLIFYING (default: 0.4)
       min_region_size=5,         # Minimum positions per region (default: 3)
       response_threshold=0.1,    # Minimum effect to consider (default: 0.0)
   )

The Key Insight
---------------

.. epigraph::

   "The guide is the probe, not the structure."

   -- PhaseLab v1.0.0 Paradigm Shift

Coherence doesn't measure properties of your perturbation (guide sequence, thermodynamics).
It measures the **system's response consistency**.

Two identical guides at different positions can have completely different outcomes because
the underlying regulatory landscape is different - and spatial coherence captures this.

Next Steps
----------

- :doc:`protein_mutscan` - Apply spatial coherence to protein mutational scanning
- :doc:`binding_landscapes` - Analyze protein-ligand binding data
- :doc:`/user_guide/howtos/select_guides_from_stable_regions` - Practical guide selection

