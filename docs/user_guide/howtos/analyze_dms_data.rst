Analyze Deep Mutational Scanning Data
=====================================

This guide shows how to use PhaseLab to analyze deep mutational scanning (DMS)
data from MaveDB or custom experiments.

Problem
-------

You have DMS fitness scores for a protein and want to:

- Identify functional domains
- Find essential residues
- Distinguish real effects from noise

Solution
--------

Use spatial coherence analysis to identify regions where mutational effects
are consistent and reproducible.

Step 1: Load Your Data
----------------------

From MaveDB format:

.. code-block:: python

   from phaselab.protein.mutscan import load_mavedb_scores

   # Load MaveDB score set
   landscape = load_mavedb_scores(
       'urn:mavedb:00000001-a-1',  # MaveDB accession
       cache_dir='~/.phaselab/cache',
   )

From CSV/TSV:

.. code-block:: python

   from phaselab.protein.mutscan import MutScanLandscape
   import pandas as pd

   # Load your data
   df = pd.read_csv('my_dms.csv')

   # Create landscape
   landscape = MutScanLandscape(
       positions=df['position'].values,
       effects=df['fitness_score'].values,
       protein_id='MY_PROTEIN',
       protein_name='My Target Protein',
       effect_type='fitness',  # or 'binding', 'stability'
   )

Step 2: Analyze Coherence
-------------------------

.. code-block:: python

   from phaselab.protein.mutscan import analyze_mutscan_coherence

   result = analyze_mutscan_coherence(
       landscape,
       window=5,                  # Residues for local coherence
       essential_threshold=-1.0,  # Score threshold for essential
       tolerant_threshold=-0.3,   # Score threshold for tolerant
   )

   # Summary
   print(f"Total residues: {len(landscape)}")
   print(f"Essential domains: {len(result.essential_regions)}")
   print(f"Variable regions: {len(result.variable_regions)}")
   print(f"Mean coherence: {result.mean_coherence:.3f}")

Step 3: Extract Functional Regions
----------------------------------

.. code-block:: python

   # Essential regions (low fitness + high coherence)
   print("\nEssential Domains:")
   for domain in result.essential_regions:
       print(f"  Residues {domain.start}-{domain.end}")
       print(f"    Mean fitness: {domain.mean_effect:.2f}")
       print(f"    Coherence: {domain.coherence:.3f}")

   # Variable regions (can be mutated)
   print("\nTolerant Regions (engineering targets):")
   for region in result.variable_regions:
       print(f"  Residues {region.start}-{region.end}")
       print(f"    Mean fitness: {region.mean_effect:.2f}")

Step 4: Per-Residue Analysis
----------------------------

.. code-block:: python

   # Get coherence profile
   profile = result.coherence_profile

   # Identify specific residues
   for i, (pos, coh, eff) in enumerate(zip(
       profile.positions,
       profile.coherence,
       profile.effects
   )):
       if coh > 0.7 and eff < -1.5:
           print(f"Critical residue {pos}: fitness={eff:.2f}, coh={coh:.3f}")

Step 5: Export Results
----------------------

.. code-block:: python

   # To DataFrame
   df_results = result.to_dataframe()
   df_results.to_csv('dms_coherence_analysis.csv')

   # To PDB B-factors (for visualization)
   result.to_pdb_bfactors(
       input_pdb='structure.pdb',
       output_pdb='structure_coherence.pdb',
       chain='A',
       value='coherence',  # or 'effect'
   )

Complete Example
----------------

.. code-block:: python

   """Analyze TEM-1 β-lactamase DMS data"""

   from phaselab.protein.mutscan import (
       load_mavedb_scores,
       analyze_mutscan_coherence,
   )

   # 1. Load data
   landscape = load_mavedb_scores('urn:mavedb:00000097-a-1')
   print(f"Loaded {len(landscape)} positions")

   # 2. Analyze
   result = analyze_mutscan_coherence(
       landscape,
       window=5,
       essential_threshold=-1.0,
   )

   # 3. Report
   print(f"\n{'='*50}")
   print("TEM-1 β-LACTAMASE COHERENCE ANALYSIS")
   print(f"{'='*50}")
   print(f"Essential domains: {len(result.essential_regions)}")
   print(f"Variable regions: {len(result.variable_regions)}")
   print(f"Overall coherence: {result.mean_coherence:.3f}")

   # 4. Key findings
   print("\nKey Essential Regions:")
   for d in sorted(result.essential_regions, key=lambda x: x.mean_effect)[:3]:
       print(f"  {d.start}-{d.end}: fitness={d.mean_effect:.2f}")

   print("\nGood for Engineering:")
   for r in result.variable_regions[:3]:
       print(f"  {r.start}-{r.end}: fitness≈{r.mean_effect:.2f}")

Tips
----

1. **Window size**: Use 3-5 for small proteins, 5-10 for large proteins
2. **Thresholds**: Adjust based on your assay's dynamic range
3. **Missing data**: Coherence handles gaps gracefully
4. **Multiple conditions**: Analyze each condition separately

See Also
--------

- :doc:`/user_guide/tutorials/protein_mutscan` - Full tutorial
- :doc:`integrate_with_alphafold` - Combine with structure prediction

