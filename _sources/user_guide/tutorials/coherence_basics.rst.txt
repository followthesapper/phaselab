Coherence Basics
================

This tutorial covers the fundamentals of Informational Relativity (IR) coherence metrics in PhaseLab.

What is Coherence?
------------------

In the IR framework, coherence (R) measures phase synchronization in a system. It is defined as:

.. math::

   \bar{R} = e^{-V_\phi / 2}

where :math:`V_\phi` is the phase variance.

- R = 1: Perfect coherence (all phases aligned)
- R = 0: Complete decoherence (random phases)
- R = e^-2 ~ 0.135: Universal GO/NO-GO threshold

The GO/NO-GO Threshold
----------------------

The threshold e^-2 emerges from information-theoretic considerations:

.. code-block:: python

   from phaselab.core.constants import E_MINUS_2, FOUR_PI_SQUARED

   print(f"GO/NO-GO threshold: {E_MINUS_2:.6f}")
   print(f"4*pi^2 constant: {FOUR_PI_SQUARED:.6f}")

Systems with R > e^-2 are classified as **GO** (coherent, reliable). Systems below are **NO-GO** (decoherent, unreliable).

Computing Coherence
-------------------

From Phase Data
^^^^^^^^^^^^^^^

.. code-block:: python

   import phaselab as pl
   import numpy as np

   # Phase measurements (radians)
   phases = [0.1, 0.15, 0.12, 0.08, 0.11]

   # Compute coherence
   R_bar = pl.coherence_score(phases)
   print(f"R = {R_bar:.4f}")

   # Phase variance
   V_phi = pl.phase_variance(phases)
   print(f"V_phi = {V_phi:.6f}")

   # Verify relationship
   import math
   R_from_V = math.exp(-V_phi / 2)
   print(f"R from V_phi: {R_from_V:.4f}")

From Quantum Expectations
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from phaselab.quantum import compute_coherence_from_expectations
   import numpy as np

   # Pauli expectation values from quantum measurement
   expectations = np.array([0.85, 0.82, 0.88, 0.79, 0.91])

   result = compute_coherence_from_expectations(expectations)
   print(f"R = {result.R_bar:.4f}")
   print(f"V_phi = {result.V_phi:.4f}")
   print(f"GO status: {result.is_go}")

Classification
--------------

Binary Classification
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   import phaselab as pl

   # GO/NO-GO classification
   print(pl.go_no_go(0.85))   # "GO"
   print(pl.go_no_go(0.10))   # "NO-GO"

   # Custom threshold
   print(pl.go_no_go(0.45, threshold=0.5))  # "NO-GO"

Detailed Categories
^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from phaselab.core.coherence import classify_coherence

   # Detailed classification
   print(classify_coherence(0.95))  # "EXCELLENT"
   print(classify_coherence(0.65))  # "GOOD"
   print(classify_coherence(0.25))  # "MODERATE"
   print(classify_coherence(0.10))  # "SEVERE"
   print(classify_coherence(0.02))  # "CRITICAL"

Classification thresholds:

============  ==============
R Range       Category
============  ==============
> 0.8         EXCELLENT
0.5 - 0.8     GOOD
e^-2 - 0.5    MODERATE
0.05 - e^-2   SEVERE
< 0.05        CRITICAL
============  ==============

Circular Statistics
-------------------

PhaseLab uses circular statistics for proper phase handling:

.. code-block:: python

   import numpy as np
   from phaselab.core.coherence import circular_mean, circular_variance

   # Phases near 0/2*pi boundary
   phases = np.array([0.1, 0.05, 6.2, 6.25, 0.02])

   # Linear statistics would be wrong!
   linear_mean = np.mean(phases)
   print(f"Linear mean: {linear_mean:.4f}")  # ~2.5 (wrong!)

   # Circular statistics handle wrap-around correctly
   circ_mean = circular_mean(phases)
   print(f"Circular mean: {circ_mean:.4f}")  # ~0.02 (correct!)

   # Circular variance
   circ_var = circular_variance(phases)
   print(f"Circular variance: {circ_var:.6f}")

Physical Interpretation
-----------------------

The coherence metric has direct physical meaning:

1. **Quantum Systems**: R measures how well quantum phases are synchronized across measurements. High R indicates reliable quantum state preparation.

2. **Biological Clocks**: R measures circadian rhythm regularity. SMS patients have reduced R due to RAI1 haploinsufficiency.

3. **CRISPR Design**: R indicates guide reliability. Guides with low R have inconsistent targeting behavior.

Example: Quantum vs Classical Noise
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   import numpy as np
   import phaselab as pl

   # Simulate quantum-coherent system
   n_measurements = 100
   true_phase = 0.5
   quantum_noise = 0.05

   quantum_phases = true_phase + np.random.normal(0, quantum_noise, n_measurements)
   R_quantum = pl.coherence_score(quantum_phases)
   print(f"Quantum coherent R: {R_quantum:.4f}")  # ~0.99

   # Simulate decoherent system
   classical_noise = 1.0
   classical_phases = true_phase + np.random.normal(0, classical_noise, n_measurements)
   R_classical = pl.coherence_score(classical_phases)
   print(f"Decoherent R: {R_classical:.4f}")  # ~0.6

   # Simulate completely random phases
   random_phases = np.random.uniform(0, 2*np.pi, n_measurements)
   R_random = pl.coherence_score(random_phases)
   print(f"Random R: {R_random:.4f}")  # ~0.1

Next Steps
----------

- :doc:`crispr_guide_design` - Apply coherence to CRISPR design
- :doc:`coherence_modes` - Learn about heuristic vs quantum modes
- :doc:`quantum_validation` - Validate on IBM Quantum hardware
