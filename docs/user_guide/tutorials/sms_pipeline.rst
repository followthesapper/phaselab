SMS Therapeutic Pipeline
========================

This tutorial covers the complete Smith-Magenis Syndrome (SMS) gene therapy
assessment pipeline in PhaseLab.

Overview
--------

Smith-Magenis Syndrome is caused by RAI1 haploinsufficiency (~50% expression).
The SMS trials module provides:

1. **CRISPRa RAI1 activation** - Boost remaining allele
2. **CRISPRi modifier suppression** - Target PER1, CRY1, CLOCK
3. **Circadian rescue simulation** - Predict sleep/wake improvement
4. **Delivery assessment** - AAV feasibility for CNS

Quick Start
-----------

.. code-block:: python

   from phaselab.trials.sms import SMSPipeline, SMSTrialConfig

   # Configure pipeline
   config = SMSTrialConfig(
       therapeutic_window=(0.70, 1.10),  # 70-110% of normal RAI1
       optimal_expression=0.80,           # Target 80%
       verbose=True,
   )

   # Run full assessment
   pipeline = SMSPipeline(config=config)
   result = pipeline.run_full_pipeline()

   # Get GO/NO-GO decision
   print(f"Overall: {result.overall_go_nogo}")
   print(f"Claim level: {result.overall_claim_level}")

Understanding Results
---------------------

.. code-block:: python

   # Individual trial results
   print(f"CRISPRa candidates: {result.crispra_result.n_candidates}")
   print(f"Circadian rescue: {result.circadian_result.metrics['rescue_status']}")
   print(f"Delivery feasible: {result.delivery_result.metrics['delivery_feasibility']}")

   # Get falsification tests for wet lab
   for test in result.falsification_tests:
       print(f"Test {test['id']}: {test['name']}")
       print(f"  Failure condition: {test['failure_condition']}")

Individual Trials
-----------------

CRISPRa RAI1 Activation
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from phaselab.trials.sms import run_sms_crispra_trial

   result = run_sms_crispra_trial(
       promoter_sequence=rai1_promoter,  # Optional: uses default if not provided
       config=config,
   )

   print(f"Status: {result.status}")
   print(f"Candidates: {result.n_candidates}")

   if result.best_candidate:
       print(f"Best guide: {result.best_candidate['sequence']}")
       print(f"Position: TSS{result.best_candidate['position']:+d}")

CRISPRi Modifier Suppression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from phaselab.trials.sms import run_sms_crispri_trial

   # Target circadian modifier genes
   for target in ["PER1", "CRY1", "CLOCK"]:
       result = run_sms_crispri_trial(
           target_gene=target,
           config=config,
       )

       print(f"{target}: {result.n_candidates} candidates")
       if result.best_candidate:
           supp = result.best_candidate.get('expected_suppression', 'N/A')
           print(f"  Expected suppression: {supp}")

Circadian Rescue Simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from phaselab.trials.sms import run_circadian_rescue_simulation

   result = run_circadian_rescue_simulation(
       predicted_rai1_expression=0.80,  # After CRISPRa
       config=config,
   )

   print(f"Rescue status: {result.metrics['rescue_status']}")
   print(f"Final R̄: {result.metrics['final_R_bar']:.3f}")
   print(f"Sleep quality: {result.metrics['sleep_quality_prediction']}")

AAV Delivery Assessment
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from phaselab.trials.sms import run_delivery_assessment

   result = run_delivery_assessment(
       modality="CRISPRa_VP64",
       target_tissue="brain",
       config=config,
   )

   print(f"Feasibility: {result.metrics['delivery_feasibility']}")
   print(f"Payload size: {result.metrics['payload_size']}bp")
   print(f"Serotype: {result.best_candidate['recommended_serotype']}")

Falsification Tests
-------------------

PhaseLab generates specific tests that would prove predictions wrong:

.. list-table:: Falsification Tests
   :header-rows: 1

   * - Test
     - Name
     - Failure Condition
   * - A
     - Ranking Validity
     - PhaseLab guides don't outperform random CRISPOR controls
   * - B
     - Risk Prediction
     - CAUTION guides don't fail more than SAFE guides
   * - C
     - Dosage Prediction
     - Expression correlation < 0.6
   * - D
     - UNKNOWN Calibration
     - UNKNOWN guides don't fail at random rate

These tests are designed to be run in the wet lab to validate or falsify
PhaseLab's predictions.

Claim Levels
------------

The pipeline reports uncertainty honestly:

- **UNKNOWN**: Cannot assess reliability
- **EXPLORATORY**: Structural priors only (no tiling data)
- **CONTEXT_DEPENDENT**: Single dataset validation
- **STRONG_COMPUTATIONAL**: Cross-validated predictions

.. important::

   All SMS therapeutic claims require wet lab validation before clinical use.
   PhaseLab provides computational guidance, not clinical recommendations.

Configuration Options
---------------------

.. code-block:: python

   config = SMSTrialConfig(
       # Therapeutic parameters
       therapeutic_window=(0.70, 1.10),
       optimal_expression=0.80,
       baseline_expression=0.50,  # SMS patients

       # Pipeline settings
       use_virtual_assay=True,     # Enhanced scoring
       coherence_mode="heuristic", # Fast mode (or "quantum")

       # Simulation parameters
       simulation_hours=120.0,     # 5 days
       n_circadian_trials=5,       # Replicates

       # Output
       verbose=True,
   )

Full Example
------------

.. code-block:: python

   from phaselab.trials.sms import (
       SMSPipeline,
       SMSTrialConfig,
   )

   # Configure
   config = SMSTrialConfig(
       therapeutic_window=(0.70, 1.10),
       verbose=True,
   )

   # Run pipeline
   pipeline = SMSPipeline(config=config)
   result = pipeline.run_full_pipeline(
       include_modifiers=True,   # CRISPRi trials
       include_editing=True,     # Base/prime editing
   )

   # Summary
   print("=" * 60)
   print("SMS GENE THERAPY ASSESSMENT")
   print("=" * 60)
   print(f"Overall GO/NO-GO: {result.overall_go_nogo}")
   print(f"Claim level: {result.overall_claim_level}")
   print()

   # CRISPRa results
   print(f"CRISPRa candidates: {result.crispra_result.n_candidates}")
   if result.crispra_result.best_candidate:
       guide = result.crispra_result.best_candidate
       print(f"  Lead guide: {guide['sequence'][:20]}...")
       print(f"  Position: TSS{guide['position']:+d}")

   # Circadian prediction
   print(f"Circadian rescue: {result.circadian_result.metrics['rescue_status']}")

   # Wet lab recommendations
   print("\nWet Lab Recommendations:")
   for rec in result.wet_lab_recommendations[:5]:
       print(f"  - {rec}")

   # Falsification tests
   print("\nFalsification Tests:")
   for test in result.falsification_tests:
       print(f"  {test['id']}: {test['name']}")

Next Steps
----------

- :doc:`/reference/trials` - Full API reference
- :doc:`/user_guide/explanations/claim_levels` - Understanding claim levels
- :doc:`/user_guide/howtos/export_for_wetlab` - Export guides for ordering

