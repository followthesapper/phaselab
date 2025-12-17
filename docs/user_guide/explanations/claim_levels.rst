Claim Levels
============

PhaseLab uses a hierarchical system to report uncertainty honestly.
This document explains each claim level and when they apply.

Why Claim Levels?
-----------------

Most prediction tools report point estimates without uncertainty.
PhaseLab acknowledges that reliability depends on:

1. **Data availability**: Have we seen similar cases?
2. **Validation depth**: Single dataset or cross-validated?
3. **Experimental support**: Structural priors or tiling data?

Rather than hiding this, we report it explicitly.

The Four Levels
---------------

UNKNOWN
~~~~~~~

**Meaning**: Cannot assess reliability

**Requirements**: Insufficient data or context

**When it applies**:

- No tiling data for the region
- Novel target with no structural homologs
- Missing critical information

**Interpretation**: Don't make predictions. Gather more data.

.. code-block:: python

   result.claim_level == ClaimLevel.UNKNOWN
   # "We can't tell you anything reliable"

EXPLORATORY
~~~~~~~~~~~

**Meaning**: Preliminary assessment only

**Requirements**: Structural priors available

**When it applies**:

- Before tiling screen (Stage I)
- Using ENCODE/Roadmap chromatin data
- Structural predictions (AlphaFold, etc.)

**Interpretation**: Directional guidance only. Not for decisions.

.. code-block:: python

   result.claim_level == ClaimLevel.EXPLORATORY
   # "This is our best guess, but unvalidated"

CONTEXT_DEPENDENT
~~~~~~~~~~~~~~~~~

**Meaning**: Valid within a specific context

**Requirements**: Single tiling dataset validation

**When it applies**:

- One tiling screen analyzed
- Spatial coherence validated (p < 0.05)
- Same cell type/conditions as tiling

**Interpretation**: Reliable IF conditions match the validation context.

.. code-block:: python

   result.claim_level == ClaimLevel.CONTEXT_DEPENDENT
   # "Validated in K562, may differ in other cells"

STRONG_COMPUTATIONAL
~~~~~~~~~~~~~~~~~~~~

**Meaning**: High computational confidence

**Requirements**: Cross-validated across datasets

**When it applies**:

- Multiple independent tiling screens
- Validated across cell types
- Consistent coherence-outcome correlation

**Interpretation**: Strong computational evidence. Still needs wet lab.

.. code-block:: python

   result.claim_level == ClaimLevel.STRONG_COMPUTATIONAL
   # "Validated across 6 genes, 115k guides"

How Claims Propagate
--------------------

When combining analyses, the overall claim level is the **minimum**
of component levels:

.. code-block:: python

   # Example: Combining CRISPRa and circadian simulations
   crispra_claim = ClaimLevel.CONTEXT_DEPENDENT
   circadian_claim = ClaimLevel.EXPLORATORY

   overall_claim = min(crispra_claim, circadian_claim)
   # Result: EXPLORATORY (weakest link)

This ensures you never overstate confidence.

The Claim Level Hierarchy
-------------------------

.. code-block:: text

   UNKNOWN < EXPLORATORY < CONTEXT_DEPENDENT < STRONG_COMPUTATIONAL

   More data / validation required ←──────────→ More confidence

Upgrading Claims
----------------

Claims only upgrade with **additional evidence**:

.. list-table:: Claim Upgrades
   :header-rows: 1

   * - From
     - To
     - Required Evidence
   * - UNKNOWN
     - EXPLORATORY
     - Structural priors or homolog data
   * - EXPLORATORY
     - CONTEXT_DEPENDENT
     - Tiling screen validation (p < 0.05)
   * - CONTEXT_DEPENDENT
     - STRONG_COMPUTATIONAL
     - Cross-validation across datasets

.. important::

   Claims **never** upgrade automatically. PhaseLab requires you to
   provide additional evidence to increase confidence.

Using Claim Levels in Code
--------------------------

.. code-block:: python

   from phaselab.fusion import ClaimLevel

   # Check claim level
   if result.claim_level >= ClaimLevel.CONTEXT_DEPENDENT:
       # Safe to make recommendations
       print("Validated prediction")
   else:
       # Need more data
       print("Preliminary only - validation required")

   # Filter by claim level
   reliable_regions = [
       r for r in result.regions
       if r.claim_level >= ClaimLevel.CONTEXT_DEPENDENT
   ]

Decision Guide
--------------

.. list-table:: What to Do at Each Level
   :header-rows: 1

   * - Claim Level
     - Action
     - Use For
   * - UNKNOWN
     - Gather data first
     - Nothing
   * - EXPLORATORY
     - Run tiling screen
     - Experimental planning
   * - CONTEXT_DEPENDENT
     - Proceed with caution
     - Guide selection, pilot studies
   * - STRONG_COMPUTATIONAL
     - Proceed confidently
     - Large-scale experiments

Therapeutic Claims
------------------

For therapeutic applications (e.g., SMS pipeline), claim levels are
particularly important:

.. warning::

   Even STRONG_COMPUTATIONAL claims require **wet lab validation**
   before any clinical application. PhaseLab provides computational
   guidance, not clinical recommendations.

The SMS pipeline generates falsification tests specifically designed
to validate or invalidate computational predictions.

See Also
--------

- :doc:`/user_guide/tutorials/sms_pipeline` - SMS therapeutic assessment
- :doc:`evidence_levels` - A/B/C evidence classification

