Coherence Theory
================

This document explains the mathematical foundations of PhaseLab's coherence calculations.

Circular Statistics
-------------------

Phase data requires circular statistics because phases are periodic (0 = 2π).

**Circular Mean**

The circular mean handles wrap-around correctly:

.. math::

   \bar{\theta} = \text{atan2}\left(\frac{1}{N}\sum_{i=1}^N \sin(\theta_i), \frac{1}{N}\sum_{i=1}^N \cos(\theta_i)\right)

Standard (linear) mean fails for phases near 0/2π boundary.

**Circular Variance**

Circular variance measures phase spread:

.. math::

   V_\phi = 1 - R

where R is the mean resultant length:

.. math::

   R = \frac{1}{N}\sqrt{\left(\sum_i \cos(\theta_i)\right)^2 + \left(\sum_i \sin(\theta_i)\right)^2}

**Properties**

- V_φ = 0: All phases identical (perfect alignment)
- V_φ = 1: Phases uniformly distributed (complete disorder)

Phase Extraction
----------------

Different data sources require different phase extraction methods.

**From Quantum Expectations**

Pauli expectation values E ∈ [-1, 1] map to phases:

.. math::

   \theta = \arccos(E)

This maps:

- E = 1 → θ = 0
- E = 0 → θ = π/2
- E = -1 → θ = π

**From Hamiltonian Coefficients**

Hamiltonian terms have coefficients that encode structure:

.. math::

   \theta_i = \arctan2(\text{Im}(c_i), \text{Re}(c_i))

For real coefficients:

.. math::

   \theta_i = \pi \cdot \mathbb{1}[c_i < 0]

Coherence Computation
---------------------

**From Phase Variance**

The coherence R-bar is:

.. math::

   \bar{R} = e^{-V_\phi / 2}

**From Expectation Values (Direct)**

For N expectation values:

.. math::

   \bar{R} = \frac{1}{N}\sqrt{\left(\sum_i \cos(\arccos(E_i))\right)^2 + \left(\sum_i \sin(\arccos(E_i))\right)^2}

Simplifies to:

.. math::

   \bar{R} = \frac{1}{N}\sqrt{\left(\sum_i E_i\right)^2 + \left(\sum_i \sqrt{1-E_i^2}\right)^2}

**Heuristic Mode**

For fast screening, PhaseLab uses coefficient variance as a proxy:

.. math::

   V_{coeff} = \text{Var}(|c_1|, |c_2|, ..., |c_n|)

This correlates with but is not identical to true phase variance.

ATLAS-Q Acceleration
--------------------

ATLAS-Q provides accelerated coherence computation through:

1. **IR Measurement Grouping**: Groups commuting Pauli terms to reduce measurements by 5x

2. **Batch Processing**: Vectorized computation over multiple circuits

3. **GPU Kernels**: Custom Triton kernels for phase statistics

The acceleration applies to quantum mode, not heuristic mode.

Variance Reduction
------------------

IR measurement grouping reduces variance:

**Without Grouping**

Each Pauli term measured independently:

.. math::

   \sigma^2_{total} = \sum_i \sigma^2_i / N_{shots}

**With Grouping**

Commuting terms measured together:

.. math::

   \sigma^2_{grouped} = \sum_g \sigma^2_g / N_{shots}

Since |groups| < |terms|, variance decreases.

Typical improvement: 5x variance reduction for molecular Hamiltonians.

Statistical Considerations
--------------------------

**Sample Size**

More measurements improve coherence estimates:

.. math::

   \sigma_{\bar{R}} \propto \frac{1}{\sqrt{N}}

Minimum recommended: N > 100 measurements

**Confidence Intervals**

For large N, R-bar is approximately normal:

.. math::

   \bar{R} \pm 1.96 \cdot \frac{\sigma_R}{\sqrt{N}}

**Bias Correction**

Small samples overestimate R. Correction factor:

.. math::

   \bar{R}_{corrected} = \bar{R} - \frac{1-\bar{R}^2}{2N}

Implementation Details
----------------------

PhaseLab implements coherence in ``phaselab/quantum/coherence.py``:

.. code-block:: python

   def compute_coherence_from_expectations(expectations, use_atlas_q=True):
       # Convert to phases
       phases = np.arccos(np.clip(expectations, -1, 1))

       # Circular statistics
       cos_sum = np.sum(np.cos(phases))
       sin_sum = np.sum(np.sin(phases))

       # Mean resultant length
       R = np.sqrt(cos_sum**2 + sin_sum**2) / len(phases)

       # Phase variance
       V_phi = 1 - R

       # Coherence
       R_bar = np.exp(-V_phi / 2)

       return CoherenceResult(R_bar=R_bar, V_phi=V_phi, ...)

See Also
--------

- :doc:`informational_relativity` - Theoretical foundation
- :doc:`go_no_go_threshold` - Threshold derivation
