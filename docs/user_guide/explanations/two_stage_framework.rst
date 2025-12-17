Two-Stage Framework
===================

PhaseLab uses a two-stage approach to reliability assessment:
feasibility inference (pre-tiling) and landscape resolution (with tiling).

Overview
--------

.. list-table::
   :header-rows: 1

   * - Stage
     - Data Required
     - Output
     - Claim Level
   * - Stage I
     - Structural priors only
     - GO/NO-GO feasibility
     - EXPLORATORY
   * - Stage II
     - 16-20 perturbations
     - Resolved landscape
     - CONTEXT_DEPENDENT

Stage I: Feasibility (Pre-Tiling)
---------------------------------

**Purpose**: Should we invest in a tiling screen?

**Inputs**:

- Gene identity and promoter region
- ENCODE/Roadmap chromatin accessibility
- Known structural features

**Method**:

.. code-block:: python

   from phaselab.spatial import estimate_feasibility

   feasibility = estimate_feasibility(
       gene_symbol='RAI1',
       promoter_coords=('chr17', 17555000, 17558000),
       cell_type='K562',
   )

   print(f"Predicted stable fraction: {feasibility.stable_fraction:.0%}")
   print(f"Recommendation: {feasibility.recommendation}")

**Outputs**:

1. **Predicted stable regions** (candidate positions)
2. **Expected signal-to-noise** (will tiling be informative?)
3. **GO/NO-GO** for tiling investment

**What Stage I CANNOT Do**:

- Definitively identify stable regions
- Guarantee guide selection success
- Replace empirical validation

Stage II: Resolution (With Tiling)
----------------------------------

**Purpose**: Resolve the actual response landscape

**Inputs**:

- Tiling screen data (16-20+ perturbations minimum)
- Position and response values

**Method**:

.. code-block:: python

   from phaselab.spatial import load_tiling_screen, analyze_tiling_coherence

   # Load tiling data
   landscape = load_tiling_screen(
       'tiling_results.tsv',
       position_col='position',
       response_col='log2fc',
   )

   # Analyze
   result = analyze_tiling_coherence(landscape)

   # Check validation
   if result.is_validated:
       print(f"Validated correlation: r={result.profile.correlation:.3f}")
       print(f"Stable regions: {len(result.stable_regions)}")
   else:
       print("Insufficient data for validation")

**Outputs**:

1. **Validated coherence profile** (empirically tested)
2. **Definitive stable/amplifying regions**
3. **Quantified variance reduction**

Why Two Stages?
---------------

**Cost-benefit optimization**:

.. list-table::
   :header-rows: 1

   * - Cost Type
     - Stage I
     - Stage II
   * - Computational
     - Minutes
     - Minutes
   * - Experimental
     - $0
     - ~$2,000 (tiling)
   * - Time
     - Hours
     - Weeks

**Risk reduction**:

Stage I filters out cases where tiling is unlikely to help:

- Uniformly accessible promoters (no landscape structure)
- Uniformly repressed regions (no signal)
- Insufficient prior knowledge

Minimum Viable Tiling
---------------------

Stage II requires a minimum dataset:

.. list-table::
   :header-rows: 1

   * - Metric
     - Minimum
     - Recommended
   * - Perturbations
     - 16
     - 20-30
   * - Spacing
     - Even across region
     - 25-50bp intervals
   * - Replicates
     - 2
     - 3+
   * - Coverage
     - Full window
     - ±500bp of target

Why 16-20? This provides sufficient power to:

1. Detect coherence-outcome correlation (r > 0.2)
2. Identify at least 2 stable regions
3. Estimate variance reduction

Transition Criteria
-------------------

**Stage I → Stage II** when:

.. code-block:: python

   # Automatic recommendation
   if feasibility.recommendation == "PROCEED_TO_TILING":
       # Run tiling screen
       pass

   # Manual override (with justification)
   if feasibility.stable_fraction < 0.3 and must_proceed:
       print("Warning: Low expected yield from tiling")

**Stage II validation** requires:

1. Correlation p-value < 0.05
2. At least one stable region identified
3. Variance reduction > 20%

Practical Workflow
------------------

.. code-block:: python

   """Complete two-stage workflow"""

   from phaselab.spatial import (
       estimate_feasibility,
       load_tiling_screen,
       analyze_tiling_coherence,
   )

   # === STAGE I ===
   print("Stage I: Feasibility Assessment")
   print("=" * 40)

   feasibility = estimate_feasibility(
       gene_symbol='MYC',
       promoter_coords=('chr8', 128747680, 128750680),
       cell_type='K562',
   )

   print(f"Predicted stable fraction: {feasibility.stable_fraction:.0%}")
   print(f"Recommendation: {feasibility.recommendation}")

   if feasibility.recommendation != "PROCEED_TO_TILING":
       print("Stage I indicates tiling unlikely to help")
       print("Consider alternative approaches")
       exit()

   # === RUN TILING SCREEN (wet lab) ===
   # ... weeks pass ...

   # === STAGE II ===
   print("\nStage II: Landscape Resolution")
   print("=" * 40)

   landscape = load_tiling_screen('myc_tiling.tsv', 'position', 'log2fc')
   result = analyze_tiling_coherence(landscape)

   print(f"Validated: {result.is_validated}")
   print(f"Correlation: r={result.profile.correlation:.3f}")
   print(f"Stable regions: {len(result.stable_regions)}")
   print(f"Claim level: {result.claim_level}")

Comparison Table
----------------

.. list-table::
   :header-rows: 1

   * - Aspect
     - Stage I
     - Stage II
   * - Data source
     - ENCODE, structure
     - Your tiling data
   * - Confidence
     - EXPLORATORY
     - CONTEXT_DEPENDENT
   * - Stable regions
     - Predicted
     - Validated
   * - Can select guides?
     - Not reliably
     - Yes
   * - When to use
     - Before experiment
     - After tiling

Summary
-------

1. **Stage I** tells you IF tiling will help (cheap, fast)
2. **Stage II** tells you WHERE to target (requires data)
3. **Never skip Stage I** - it saves wasted experiments
4. **Never skip Stage II** - Stage I predictions need validation

See Also
--------

- :doc:`/user_guide/tutorials/spatial_coherence` - Core methodology
- :doc:`claim_levels` - Understanding confidence levels

