Protein Mutational Scanning
===========================

Apply spatial coherence to deep mutational scanning (DMS) data to identify
functional domains and essential residues.

Overview
--------

Deep mutational scanning generates fitness landscapes across proteins:

- **Position**: Residue number
- **Response**: Fitness effect (growth rate, binding, activity)

Spatial coherence identifies regions where mutations have consistent effects,
distinguishing functional cores from tolerant regions.

Basic Example
-------------

.. code-block:: python

   from phaselab.protein.mutscan import (
       MutScanLandscape,
       analyze_mutscan_coherence,
   )

   # Your DMS data
   landscape = MutScanLandscape(
       positions=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
       effects=[-0.1, -0.2, -2.5, -2.8, -3.0, -2.7, -0.3, -0.1, 0.0, 0.1],
       protein_id="TEM1",
       protein_name="TEM-1 β-Lactamase",
       effect_type="fitness",
   )

   # Analyze coherence
   result = analyze_mutscan_coherence(
       landscape,
       window=5,                    # Residues for local coherence
       essential_threshold=-1.0,   # Effect threshold for essential
   )

   # View results
   print(f"Functional domains: {len(result.functional_domains)}")
   print(f"Essential regions: {len(result.essential_regions)}")
   print(f"Variable regions: {len(result.variable_regions)}")

Finding Essential Domains
-------------------------

.. code-block:: python

   # Essential domains: low fitness + high coherence
   for domain in result.essential_regions:
       print(f"Essential: residues {domain.start}-{domain.end}")
       print(f"  Mean fitness effect: {domain.mean_effect:.2f}")
       print(f"  Coherence: {domain.coherence:.3f}")
       print(f"  Likely function: {domain.annotation}")

   # Variable regions: mutations tolerated
   for region in result.variable_regions:
       print(f"Tolerant: residues {region.start}-{region.end}")
       print(f"  Good for engineering: mutations here may be neutral")

Per-Residue Coherence Profile
-----------------------------

.. code-block:: python

   from phaselab.protein.mutscan import local_coherence_profile

   # Get coherence at each position
   profile = local_coherence_profile(landscape, window=5)

   for pos, coh in zip(profile.positions, profile.coherence):
       status = "essential" if coh > 0.7 else "variable"
       print(f"Residue {pos}: coherence={coh:.3f} ({status})")

Loading Real DMS Data
---------------------

.. code-block:: python

   from phaselab.protein.mutscan import load_mutscan_data

   # From MaveDB-style file
   landscape = load_mutscan_data(
       'dms_scores.csv',
       position_col='position',
       effect_col='score',
       protein_id='GFP',
       protein_name='Green Fluorescent Protein',
   )

   # From pandas DataFrame
   import pandas as pd
   df = pd.read_csv('my_dms.csv')

   landscape = MutScanLandscape(
       positions=df['residue'].values,
       effects=df['fitness'].values,
       protein_id='MYC',
   )

Mapping to Structure
--------------------

Map coherence values to a PDB file for visualization:

.. code-block:: python

   from phaselab.protein.mutscan import map_coherence_to_structure

   # Map to B-factor column
   map_coherence_to_structure(
       result,
       pdb_file='protein.pdb',
       output_file='protein_coherence.pdb',
       chain='A',
   )

   # Visualize in PyMOL:
   # spectrum b, blue_white_red
   # Blue = low coherence (variable)
   # Red = high coherence (essential)

Interpretation Guide
--------------------

.. list-table:: DMS Coherence Interpretation
   :header-rows: 1

   * - Coherence
     - Fitness Effect
     - Interpretation
   * - High (>0.7)
     - Strong negative
     - Essential functional domain
   * - High (>0.7)
     - Near zero
     - Structural (tolerates mutation but constrained)
   * - Low (<0.4)
     - Variable
     - Flexible region, good for engineering
   * - Low (<0.4)
     - Strong negative
     - Context-dependent (mutation effect varies)

Use Cases
---------

1. **Protein Engineering**: Target low-coherence regions for mutations
2. **Drug Discovery**: High-coherence essential domains are drug targets
3. **Variant Interpretation**: High-coherence positions are likely pathogenic
4. **Structure Validation**: Compare coherence to structural features

Next Steps
----------

- :doc:`binding_landscapes` - Protein-ligand binding analysis
- :doc:`/user_guide/howtos/integrate_with_alphafold` - Combine with structure prediction

