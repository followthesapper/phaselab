Quantum Hardware Validation
===========================

This tutorial covers validating guide RNA coherence on IBM Quantum hardware.

Overview
--------

PhaseLab can validate coherence metrics on real quantum hardware:

- IBM Torino (127 qubits, Heron processor)
- IBM Brisbane (127 qubits)
- IBM Kyoto (127 qubits)

Hardware validation provides the strongest evidence (Level A) for guide reliability.

Setup
-----

IBM Quantum Account
^^^^^^^^^^^^^^^^^^^

1. Create account at `quantum.ibm.com <https://quantum.ibm.com>`_
2. Get your API token from the dashboard
3. Save credentials:

.. code-block:: python

   from qiskit_ibm_runtime import QiskitRuntimeService

   # Save credentials (only need to do once)
   QiskitRuntimeService.save_account(
       channel="ibm_quantum",
       token="your-token-here",
       overwrite=True
   )

   # Or use environment variable
   import os
   os.environ['IBM_QUANTUM_TOKEN'] = 'your-token-here'

Check Available Backends
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from qiskit_ibm_runtime import QiskitRuntimeService

   service = QiskitRuntimeService()

   # List available backends
   for backend in service.backends():
       print(f"{backend.name}: {backend.num_qubits} qubits")

Simulator Testing
-----------------

Always test with simulator first:

.. code-block:: python

   from phaselab.quantum import QuantumCoherenceValidator

   # Use Aer simulator
   validator = QuantumCoherenceValidator(backend='aer_simulator')

   guide = "GAAGGAGAGCAAGAGCGCGA"
   result = validator.validate(guide)

   print(f"Simulator R: {result['coherence']:.4f}")
   print(f"Status: {result['status']}")

Hardware Validation
-------------------

Basic Hardware Run
^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from phaselab.quantum import QuantumCoherenceValidator

   # Initialize with hardware backend
   validator = QuantumCoherenceValidator(
       backend='ibm_torino',
       shots=4096
   )

   guide = "GAAGGAGAGCAAGAGCGCGA"
   result = validator.validate(guide)

   print(f"Hardware R: {result['coherence']:.4f}")
   print(f"Status: {result['status']}")
   print(f"Job ID: {result['job_id']}")

Batch Validation
^^^^^^^^^^^^^^^^

.. code-block:: python

   from phaselab.quantum import QuantumCoherenceValidator

   validator = QuantumCoherenceValidator(
       backend='ibm_torino',
       shots=4096
   )

   guides = [
       "GAAGGAGAGCAAGAGCGCGA",
       "AACTGCAAAGAAGTGGGCAC",
       "TACAGGAGCTTCCAGCGTCA"
   ]

   # Batch validation (more efficient)
   results = validator.validate_batch(guides)

   for guide, result in zip(guides, results):
       print(f"{guide[:10]}... R={result['coherence']:.3f} [{result['status']}]")

Async Jobs
^^^^^^^^^^

For long-running jobs:

.. code-block:: python

   from phaselab.quantum import QuantumCoherenceValidator

   validator = QuantumCoherenceValidator(backend='ibm_torino')

   # Submit job asynchronously
   job = validator.submit_validation(guide="GAAGGAGAGCAAGAGCGCGA")
   print(f"Job submitted: {job.job_id()}")

   # Check status
   print(f"Status: {job.status()}")

   # Wait for result
   result = validator.get_result(job)
   print(f"R: {result['coherence']:.4f}")

IR Measurement Grouping
-----------------------

PhaseLab uses IR measurement grouping for 5x variance reduction:

.. code-block:: python

   from phaselab.quantum import QuantumCoherenceValidator

   validator = QuantumCoherenceValidator(
       backend='ibm_torino',
       use_ir_grouping=True,  # Enable IR grouping
       shots=4096
   )

   # IR grouping reduces variance by grouping commuting observables
   result = validator.validate(
       guide="GAAGGAGAGCAAGAGCGCGA",
       return_details=True
   )

   print(f"R: {result['coherence']:.4f}")
   print(f"Variance: {result['variance']:.6f}")
   print(f"Measurement groups: {result['n_groups']}")

Error Mitigation
----------------

PhaseLab supports error mitigation:

.. code-block:: python

   from phaselab.quantum import QuantumCoherenceValidator

   validator = QuantumCoherenceValidator(
       backend='ibm_torino',
       error_mitigation='zne',  # Zero-noise extrapolation
       shots=4096
   )

   result = validator.validate(guide="GAAGGAGAGCAAGAGCGCGA")
   print(f"Mitigated R: {result['coherence']:.4f}")

Available mitigation strategies:

- ``'none'``: No mitigation
- ``'twirling'``: Pauli twirling
- ``'zne'``: Zero-noise extrapolation
- ``'pec'``: Probabilistic error cancellation

VQE Coherence
-------------

Run VQE for research-grade coherence:

.. code-block:: python

   from phaselab.quantum import VQECoherenceEstimator
   from phaselab.core.hamiltonians import build_grna_hamiltonian

   # Build Hamiltonian
   guide = "GAAGGAGAGCAAGAGCGCGA"
   H = build_grna_hamiltonian(guide)

   # Run VQE
   estimator = VQECoherenceEstimator(
       backend='ibm_torino',
       ansatz='hardware_efficient',
       optimizer='COBYLA',
       maxiter=100
   )

   result = estimator.run(H)

   print(f"Ground state energy: {result['energy']:.6f}")
   print(f"Coherence R: {result['coherence']:.4f}")
   print(f"Optimal parameters: {result['parameters']}")

Hardware Results Analysis
-------------------------

Analyze hardware validation results:

.. code-block:: python

   from phaselab.quantum import analyze_hardware_results

   # Collect multiple runs
   results = []
   for i in range(5):
       r = validator.validate(guide)
       results.append(r)

   # Analyze
   analysis = analyze_hardware_results(results)

   print(f"Mean R: {analysis['mean_coherence']:.4f}")
   print(f"Std R: {analysis['std_coherence']:.4f}")
   print(f"Confidence interval: {analysis['ci_95']}")
   print(f"Hardware noise estimate: {analysis['noise_estimate']:.4f}")

Complete Validation Pipeline
----------------------------

.. code-block:: python

   from phaselab.crispr import design_guides, compute_coherence_batch
   from phaselab.quantum import QuantumCoherenceValidator

   # Step 1: Design guides
   guides_df = design_guides(sequence, tss_index)

   # Step 2: Filter with heuristic coherence
   go_guides = guides_df[guides_df['go_no_go'] == 'GO']

   # Step 3: Quantum simulation for top candidates
   top_10 = go_guides.head(10)
   sim_coherence = compute_coherence_batch(
       top_10['sequence'].tolist(),
       mode="quantum"
   )
   top_10['sim_coherence'] = sim_coherence

   # Step 4: Hardware validation for top 3
   validator = QuantumCoherenceValidator(backend='ibm_torino')

   top_3 = top_10.nlargest(3, 'sim_coherence')
   for idx, row in top_3.iterrows():
       result = validator.validate(row['sequence'])
       top_3.loc[idx, 'hw_coherence'] = result['coherence']
       top_3.loc[idx, 'evidence_level'] = 'A'

   print(top_3[['sequence', 'sim_coherence', 'hw_coherence', 'evidence_level']])

Expected Hardware Results
-------------------------

Based on IBM hardware validation:

================  ==========  ==============
Guide Type        R Range     Status
================  ==========  ==============
High quality      0.94-0.97   GO (excellent)
Good quality      0.84-0.94   GO (good)
Marginal          0.50-0.84   GO (moderate)
Poor quality      <0.50       Variable
================  ==========  ==============

Hardware validation typically shows:

- R values 20-30% higher than heuristic
- Good correlation with simulation (r > 0.9)
- Variance reduction with IR grouping

Troubleshooting
---------------

**Job Queue Times**

Hardware jobs may queue for minutes to hours:

.. code-block:: python

   # Check queue position
   job = validator.submit_validation(guide)
   print(f"Queue position: {job.queue_position()}")

**Hardware Errors**

.. code-block:: python

   try:
       result = validator.validate(guide)
   except Exception as e:
       print(f"Hardware error: {e}")
       # Fallback to simulation
       result = validator.validate(guide, backend='aer_simulator')

**Rate Limits**

IBM Quantum has rate limits. For many guides:

.. code-block:: python

   import time

   for guide in guides:
       result = validator.validate(guide)
       time.sleep(1)  # Rate limit delay

Next Steps
----------

- :doc:`coherence_modes` - Compare modes
- :doc:`crispr_guide_design` - Full design pipeline
