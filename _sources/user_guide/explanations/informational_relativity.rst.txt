Informational Relativity
========================

Informational Relativity (IR) is the theoretical framework underlying PhaseLab's coherence metrics.

Core Concept
------------

IR proposes that information processing systems exhibit phase coherence that can be quantified using circular statistics. The central insight is that phase synchronization is a universal indicator of system reliability.

The coherence metric R (R-bar) measures how well phases are aligned:

.. math::

   \bar{R} = e^{-V_\phi / 2}

where :math:`V_\phi` is the phase variance computed using circular statistics.

Physical Interpretation
-----------------------

**Quantum Systems**

In quantum mechanics, phases determine interference patterns. When phases are coherent (aligned), constructive interference produces reliable, repeatable outcomes. When phases are random, destructive interference leads to noise.

**Biological Clocks**

Circadian rhythms depend on synchronized oscillations across cells and tissues. High coherence indicates robust, predictable rhythms. Low coherence indicates disrupted timing, as seen in Smith-Magenis Syndrome.

**Guide RNA Binding**

CRISPR guide RNAs with coherent binding kinetics show consistent on-target activity. Guides with phase-disordered binding show variable efficiency and higher off-target rates.

The Universal Threshold
-----------------------

The GO/NO-GO threshold e^-2 (approximately 0.135) emerges from information-theoretic considerations:

- Systems with R > e^-2 have sufficient phase coherence for reliable operation
- Systems with R < e^-2 are dominated by phase noise and unreliable

This threshold is not arbitrary. It represents the boundary where:

1. Signal exceeds noise in phase measurements
2. Predictive power of the system becomes meaningful
3. Information transfer remains stable

Mathematical Foundation
-----------------------

**Circular Statistics**

Phase data lives on a circle (0 to 2π), requiring circular statistics:

.. math::

   \bar{\theta} = \text{atan2}\left(\sum_i \sin(\theta_i), \sum_i \cos(\theta_i)\right)

The circular variance is:

.. math::

   V = 1 - \frac{1}{N}\sqrt{\left(\sum_i \cos(\theta_i)\right)^2 + \left(\sum_i \sin(\theta_i)\right)^2}

**Phase Variance to Coherence**

The transformation from phase variance to coherence:

.. math::

   \bar{R} = e^{-V_\phi / 2}

This exponential relationship ensures:

- R = 1 when V_φ = 0 (perfect coherence)
- R decreases smoothly as V_φ increases
- R approaches 0 for large V_φ (complete decoherence)

Connection to Quantum Mechanics
-------------------------------

In quantum systems, the coherence metric relates to:

1. **Expectation Values**: Pauli measurements yield expectation values that encode phase information

2. **Density Matrix Off-Diagonals**: Coherence R relates to the magnitude of off-diagonal elements in the density matrix

3. **Decoherence Time**: Systems with higher R maintain coherence longer under environmental noise

Historical Context
------------------

IR builds on:

- **Random Matrix Theory**: Universal spectral statistics in complex systems
- **Circular Statistics**: Proper handling of angular/phase data
- **Quantum Information Theory**: Entropy and coherence measures
- **Synchronization Theory**: Kuramoto oscillators and coupled systems

The framework was developed to provide a unified language for coherence across quantum, biological, and engineered systems.

Practical Implications
----------------------

**For CRISPR Design**

- Guides with high R show consistent editing efficiency
- Low R guides have unpredictable performance
- The e^-2 threshold filters out unreliable candidates

**For Gene Therapy**

- Therapeutic interventions should restore R above threshold
- Dose-response curves should track R, not just expression levels
- Multi-tissue coherence indicates systemic restoration

**For Quantum Computing**

- Algorithm reliability correlates with R
- Circuit optimization should maximize R
- Hardware validation requires R measurements

See Also
--------

- :doc:`coherence_theory` - Mathematical details of coherence
- :doc:`go_no_go_threshold` - Why e^-2 is the threshold
