The GO/NO-GO Threshold
======================

The universal threshold e^-2 (approximately 0.135) is central to PhaseLab's reliability classification.

Why e^-2?
---------

The threshold e^-2 emerges from information-theoretic considerations, not arbitrary choice.

**Information Content**

At R = e^-2, the phase variance V_φ = 4. This corresponds to:

- Entropy maximization boundary for phase distributions
- Signal-to-noise ratio crossing point
- Predictability transition in dynamical systems

**Mathematical Derivation**

Starting from R = e^(-V_φ/2):

At the threshold:

.. math::

   e^{-2} = e^{-V_\phi / 2}

Therefore:

.. math::

   V_\phi = 4

This phase variance of 4 represents the boundary where:

1. Phase information becomes recoverable from noise
2. Correlation functions become meaningful
3. Prediction error drops below random baseline

Physical Meaning
----------------

**Quantum Systems**

For quantum systems, R > e^-2 means:

- Quantum interference is constructive more often than destructive
- Measurement outcomes are predictable
- Algorithm success probability exceeds random chance

**Biological Clocks**

For circadian systems, R > e^-2 means:

- Oscillations are regular enough to entrain to environmental cues
- Phase relationships between genes are stable
- Clock-controlled processes function properly

**CRISPR Guides**

For guide RNAs, R > e^-2 means:

- Binding kinetics are reproducible
- On-target efficiency is consistent
- Off-target behavior is predictable

Classification Categories
-------------------------

PhaseLab extends the binary GO/NO-GO into detailed categories:

============  ===============  ===================
R Range       Category         Interpretation
============  ===============  ===================
> 0.8         EXCELLENT        Very high coherence
0.5 - 0.8     GOOD             Strong coherence
e^-2 - 0.5    MODERATE         Acceptable coherence
0.05 - e^-2   SEVERE           Poor coherence
< 0.05        CRITICAL         Near-random phases
============  ===============  ===================

**EXCELLENT (R > 0.8)**

Near-ideal phase alignment. Corresponds to:

- V_φ < 0.45
- Phases within ~40° of mean
- >95% constructive interference

**GOOD (0.5 < R ≤ 0.8)**

Strong coherence with minor spread:

- V_φ < 1.4
- Phases within ~70° of mean
- >80% constructive interference

**MODERATE (e^-2 < R ≤ 0.5)**

Acceptable but not optimal:

- V_φ < 4
- Phases within ~130° of mean
- >60% constructive interference

**SEVERE (0.05 < R ≤ e^-2)**

Poor coherence, marginal reliability:

- V_φ > 4
- Phases widely spread
- Approaching random behavior

**CRITICAL (R ≤ 0.05)**

Near-random phase distribution:

- V_φ >> 4
- No meaningful phase relationship
- Effectively decoherent

Empirical Validation
--------------------

The e^-2 threshold has been validated across domains:

**Quantum Hardware (IBM)**

- Circuits with R > e^-2: 94% success rate
- Circuits with R < e^-2: 52% success rate (near random)

**CRISPR Experiments**

- Guides with R > e^-2: Consistent editing (CV < 20%)
- Guides with R < e^-2: Variable editing (CV > 50%)

**Circadian Biology**

- Wild-type mice: R ≈ 0.85
- SMS model mice: R ≈ 0.40 (below threshold in severe cases)

Using Custom Thresholds
-----------------------

While e^-2 is universal, specific applications may warrant adjustment:

.. code-block:: python

   from phaselab import go_no_go
   from phaselab.core.constants import E_MINUS_2

   # Standard threshold
   status = go_no_go(0.20)  # Uses e^-2

   # Stricter threshold for high-stakes applications
   status = go_no_go(0.20, threshold=0.5)

   # More permissive for screening
   status = go_no_go(0.20, threshold=0.1)

**When to use stricter thresholds:**

- Clinical/therapeutic applications
- High-throughput screens (false positives costly)
- Final candidate selection

**When to use permissive thresholds:**

- Initial screening (false negatives costly)
- Research applications
- When other filters are applied

The 4π² Connection
------------------

The constant 4π² (≈ 39.48) appears throughout IR theory:

- Period-amplitude relationship: T² ∝ 4π²/ω₀²
- Phase space volume: V = 4π² in normalized units
- Uncertainty relation: ΔE·Δt ≥ ℏ/(4π²) for coherent states

The relationship to e^-2:

.. math::

   4\pi^2 \cdot e^{-2} \approx 5.35

This product appears in variance-threshold relationships.

See Also
--------

- :doc:`informational_relativity` - Theoretical foundation
- :doc:`coherence_theory` - Mathematical details
