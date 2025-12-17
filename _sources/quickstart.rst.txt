Quick Start
===========

This guide covers the essential PhaseLab workflows in 5 minutes.

Core Concept: Spatial Coherence
-------------------------------

PhaseLab's foundation is **spatial coherence of response landscapes**:

.. code-block:: python

   from phaselab.landscapes import ResponseLandscape
   from phaselab.spatial import analyze_tiling_coherence

   # Your perturbation-response data
   landscape = ResponseLandscape(
       coords=[100, 150, 200, 250, 300],  # Perturbation positions
       responses=[0.8, 0.75, 0.1, 0.82, 0.78],  # Response values
   )

   # Analyze spatial coherence
   result = analyze_tiling_coherence(landscape)

   # Check if model is validated
   print(f"Correlation: {result.profile.correlation:.3f}")
   print(f"Validated: {'YES' if result.is_validated else 'NO'}")

   # Find stable regions
   for region in result.stable_regions:
       print(f"Stable: {region['start']}-{region['end']}")
       print(f"  Coherence: {region['coherence']:.3f}")

Stability Classification
------------------------

PhaseLab classifies perturbation regions into four categories:

- **STABLE** (coherence > 0.7): Low variance, reliable predictions
- **MIXED** (0.4-0.7): Moderate variance, context-dependent
- **AMPLIFYING** (< 0.4): High variance, unreliable
- **IRRELEVANT**: Below response threshold

CRISPR Tiling Analysis
----------------------

Analyze CRISPR tiling screen data:

.. code-block:: python

   from phaselab.spatial import load_tiling_screen, analyze_tiling_coherence

   # Load tiling screen data
   landscape = load_tiling_screen(
       "tiling_screen.tsv",
       position_col="tss_distance",
       response_col="log2fc",
   )

   # Analyze coherence
   result = analyze_tiling_coherence(
       landscape,
       window=5,
       stable_threshold=0.7,
   )

   # Select guides from stable regions
   for region in result.stable_regions:
       print(f"Select guides from {region['start']} to {region['end']}")

Protein Mutational Scanning
---------------------------

Analyze deep mutational scanning (DMS) data:

.. code-block:: python

   from phaselab.protein.mutscan import MutScanLandscape, analyze_mutscan_coherence

   # Create landscape from DMS data
   landscape = MutScanLandscape(
       positions=residue_positions,
       effects=fitness_effects,
       protein_id="TEM1",
       protein_name="TEM-1 β-Lactamase",
       effect_type="fitness",
   )

   # Analyze coherence
   result = analyze_mutscan_coherence(
       landscape,
       window=10,
       essential_threshold=-1.0,
   )

   # Find essential domains
   print(f"Essential regions: {len(result.essential_regions)}")
   for domain in result.essential_regions[:3]:
       print(f"  Residues {domain.start}-{domain.end}: effect={domain.mean_effect:.2f}")

Binding Landscape Analysis
--------------------------

Analyze protein-ligand binding data:

.. code-block:: python

   from phaselab.chem import BindingLandscape, analyze_binding_coherence

   # Create binding landscape
   landscape = BindingLandscape(
       positions=mutation_positions,
       affinities=delta_delta_g_values,
       target="ABL1",
       ligand="Imatinib",
       affinity_type="ddG",
   )

   # Analyze coherence
   result = analyze_binding_coherence(landscape)

   # Find reliable hot spots
   for hotspot in result.hot_spots[:5]:
       print(f"Position {hotspot['position']}: ΔΔG={hotspot['effect']:.2f}")

Claim Levels
------------

PhaseLab reports honest uncertainty:

.. code-block:: python

   from phaselab.fusion import ClaimLevel

   # Claim levels propagate through analyses
   print("UNKNOWN - Cannot assess reliability")
   print("EXPLORATORY - Structural priors only")
   print("CONTEXT_DEPENDENT - Single tiling dataset")
   print("STRONG_COMPUTATIONAL - Cross-validated")

Quantum Mode
------------

For most use cases, quantum mode should be OFF:

.. code-block:: python

   from phaselab.quantum import set_quantum_mode, QuantumMode

   # Default: Classical only (fastest)
   set_quantum_mode(QuantumMode.OFF)

   # Validation: Classical + quantum subset
   set_quantum_mode(QuantumMode.AUDIT)

   # Research: Quantum mandatory
   set_quantum_mode(QuantumMode.REQUIRED)

**Rule:** If a classical experiment can falsify a claim, quantum is optional.

SMS Therapeutic Pipeline
------------------------

Run the complete SMS gene therapy assessment:

.. code-block:: python

   from phaselab.trials.sms import SMSPipeline, SMSTrialConfig

   # Configure pipeline
   config = SMSTrialConfig(
       therapeutic_window=(0.70, 1.10),
       verbose=True,
   )

   # Run assessment
   pipeline = SMSPipeline(config=config)
   result = pipeline.run_full_pipeline()

   # Get GO/NO-GO decision
   print(f"Overall: {result.overall_go_nogo}")
   print(f"Claim level: {result.overall_claim_level}")

   # Get falsification tests
   for test in result.falsification_tests:
       print(f"Test {test['id']}: {test['name']}")

Next Steps
----------

- :doc:`user_guide/tutorials/index` - In-depth tutorials
- :doc:`reference/index` - Complete API reference
- :doc:`examples/index` - Working code examples
- :doc:`research/index` - Research applications
