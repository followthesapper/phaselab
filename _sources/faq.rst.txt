Frequently Asked Questions
==========================

General Questions
-----------------

What is PhaseLab?
^^^^^^^^^^^^^^^^^

PhaseLab is a phase-coherence analysis framework for quantum, biological, and dynamical systems. It provides tools for CRISPR guide RNA design, circadian clock modeling, and quantum coherence validation.

What is the GO/NO-GO threshold?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The GO/NO-GO threshold is e^-2 (approximately 0.135). Systems with coherence R > e^-2 are classified as "GO" (coherent, reliable), while those below are "NO-GO" (decoherent, unreliable).

This threshold emerges from information-theoretic considerations in the Informational Relativity framework.

What's the difference between heuristic and quantum coherence modes?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Heuristic mode** (default):

- Fast (~0.1ms per guide)
- Uses Hamiltonian coefficient variance as proxy
- R values cluster around 0.68-0.69
- Good for screening, use as tie-breaker

**Quantum mode**:

- Slow (~100-500ms per guide)
- Runs actual VQE simulation
- R values range 0.84-0.97 (matches hardware)
- Use for research-grade analysis

CRISPR Questions
----------------

Which CRISPR modality should I use?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

============  =============================  ================
Modality      Use Case                       Typical R Range
============  =============================  ================
CRISPRa       Gene activation                0.84-0.97
CRISPRi       Gene repression                0.84-0.97
Knockout      Gene disruption                0.80-0.95
Prime edit    Precise edits                  0.75-0.92
Base edit     Single nucleotide changes      0.78-0.94
============  =============================  ================

How do I interpret the combined_score?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The combined score (0-1) integrates:

- GC content (40-70% optimal)
- Thermodynamic stability
- Position relative to TSS
- Chromatin accessibility
- Off-target specificity
- Coherence (heuristic or quantum)

Higher scores indicate better candidates.

What are evidence levels?
^^^^^^^^^^^^^^^^^^^^^^^^^

- **Level A**: Hardware-validated on IBM Quantum (strongest)
- **Level B**: VQE-simulated with quantum mode
- **Level C**: Heuristic only (weakest)

Prefer guides with Level A or B evidence for experimental validation.

Quantum Questions
-----------------

Do I need IBM Quantum access?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

No. PhaseLab works fully with:

- Heuristic mode (fast, no quantum needed)
- Quantum simulation mode (local VQE)
- ATLAS-Q backend (if installed)

IBM Quantum is only needed for hardware validation (Level A evidence).

What is ATLAS-Q?
^^^^^^^^^^^^^^^^

ATLAS-Q is a GPU-accelerated quantum tensor network simulator. When available, it provides:

- 5x variance reduction via IR measurement grouping
- Faster coherence computation
- Better accuracy

Install with: ``pip install atlas-quantum``

How do I set up IBM Quantum?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from qiskit_ibm_runtime import QiskitRuntimeService

   QiskitRuntimeService.save_account(
       channel="ibm_quantum",
       token="your-token-here"
   )

Get your token at `quantum.ibm.com <https://quantum.ibm.com>`_.

Troubleshooting
---------------

ImportError: No module named 'phaselab'
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Ensure you're in the correct Python environment:

.. code-block:: bash

   python -m pip install phaselab
   python -c "import phaselab; print(phaselab.__version__)"

ATLAS-Q not detected
^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from phaselab.quantum import is_atlas_q_available
   print(is_atlas_q_available())  # Should be True

If False, install ATLAS-Q:

.. code-block:: bash

   pip install atlas-quantum

Quantum mode is slow
^^^^^^^^^^^^^^^^^^^^

Quantum mode runs VQE simulation (~100-500ms per guide). For screening:

1. Use heuristic mode for initial filtering
2. Apply quantum mode only to top candidates

.. code-block:: python

   # Screen with heuristic
   guides = design_guides(seq, tss, config=Config(coherence_mode="heuristic"))

   # Validate top 10 with quantum
   top_10 = guides.head(10)
   quantum_r = compute_coherence_batch(top_10['sequence'].tolist(), mode="quantum")

Getting Help
------------

- GitHub Issues: `github.com/followthesapper/phaselab/issues <https://github.com/followthesapper/phaselab/issues>`_
- Documentation: `phaselab.readthedocs.io <https://phaselab.readthedocs.io>`_
- Email: Contact via GitHub
