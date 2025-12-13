CRISPR Scoring
==============

This document explains how PhaseLab scores CRISPR guide RNAs.

Scoring Components
------------------

Each guide receives scores from multiple components:

**Sequence Features**

1. **GC Content** (0-1)
   - Optimal range: 40-65%
   - Penalized outside this range
   - Affects binding stability

2. **Homopolymer Runs** (penalty)
   - Runs of 4+ identical bases penalized
   - Causes synthesis and binding issues

3. **Thermodynamic Stability** (ΔG)
   - SantaLucia nearest-neighbor model
   - Optimal: -30 to -45 kcal/mol
   - Too strong: slow turnover
   - Too weak: poor binding

**Position Features**

4. **Distance from TSS** (0-1)
   - CRISPRa: -400 to +50 optimal
   - CRISPRi: +50 to +300 optimal
   - Knockout: Exonic regions

5. **Chromatin Accessibility** (0-1)
   - DNase peaks boost score
   - ATAC-seq integration available
   - Closed chromatin penalized

**Specificity Features**

6. **Off-Target Score** (0-1)
   - MIT specificity algorithm
   - CFD score for mismatches
   - Lower = more off-targets

7. **Seed Region** (penalty)
   - Mismatches in seed (positions 1-12) heavily weighted
   - PAM-proximal region critical

**Coherence Features (v0.6.0+)**

8. **Coherence Score** (0-1)
   - IR framework metric
   - Two modes: heuristic vs quantum
   - GO/NO-GO classification

Score Combination
-----------------

**v0.6.0 Weighting**

.. code-block:: python

   combined_score = (
       0.20 * gc_score +
       0.15 * thermo_score +
       0.15 * position_score +
       0.15 * accessibility_score +
       0.35 * specificity_score
   )

   # Coherence as secondary filter
   if not is_go:
       combined_score *= 0.5  # Penalty for NO-GO

**v0.6.1 Two-Stage Scoring**

Stage 1: Hard safety gates (must pass all):

- Off-target count < threshold
- No exact matches in critical regions
- GC content in acceptable range

Stage 2: Soft ranking (weighted sum):

.. code-block:: python

   soft_score = (
       0.25 * gc_score +
       0.20 * thermo_score +
       0.15 * position_score +
       0.10 * accessibility_score +
       0.25 * specificity_score +
       weight * coherence_score  # 0.30 for quantum, 0.05 for heuristic
   )

Coherence Weighting (v0.6.1)
----------------------------

**The Problem**

v0.6.0 weighted heuristic coherence at 0.30, but:

- Heuristic R clusters around 0.68-0.69
- Poor discrimination between guides
- Overinflated influence on ranking

**The Solution**

v0.6.1 adjusts weights by mode:

============  ========  =============
Mode          Weight    Rationale
============  ========  =============
heuristic     0.05      Tie-breaker only
quantum       0.30      Research-grade
============  ========  =============

Evidence Levels
---------------

v0.6.1 assigns evidence levels affecting final scores:

**Level A: Hardware-Validated**

- Validated on IBM Quantum hardware
- Full score weight
- Strongest evidence

**Level B: VQE-Simulated**

- Quantum mode coherence
- Full score weight
- Good evidence

**Level C: Heuristic Only**

- Fast proxy metric
- Capped score influence
- Weaker evidence

.. code-block:: python

   if evidence_level == 'C':
       # Cap heuristic-only guides
       combined_score = min(combined_score, 0.85)

Risk Mass Metrics (v0.6.1)
--------------------------

New metrics for off-target risk:

**risk_mass_close**

Off-targets within 100bp of any TSS:

.. math::

   \text{risk\_mass\_close} = \sum_{ot \in \text{close}} \text{CFD}(ot)

**risk_mass_exonic**

Off-targets in exonic regions:

.. math::

   \text{risk\_mass\_exonic} = \sum_{ot \in \text{exonic}} \text{CFD}(ot)

**tail_risk_score**

Aggregate tail risk:

.. math::

   \text{tail\_risk} = \frac{\sum_{i > 90\%} \text{CFD}_i}{\sum_i \text{CFD}_i}

Modality-Specific Scoring
-------------------------

**CRISPRa**

Additional factors:

- VP64/VPR fusion compatibility
- Synergistic activation domain distance
- Enhancer proximity bonus

**CRISPRi**

Additional factors:

- KRAB domain compatibility
- Steric hindrance score
- Repression efficiency model

**Knockout**

Additional factors:

- Frameshift probability
- Repair pathway prediction (NHEJ vs HDR)
- Essential exon targeting

**Base Editing**

Additional factors:

- Activity window position (4-8)
- Bystander edit count
- C or A presence at target

**Prime Editing**

Additional factors:

- PBS binding strength
- RT template efficiency
- Hairpin formation risk

Score Interpretation
--------------------

============  =================
Score Range   Interpretation
============  =================
> 0.8         Excellent candidate
0.6 - 0.8     Good candidate
0.4 - 0.6     Acceptable
0.2 - 0.4     Marginal
< 0.2         Poor candidate
============  =================

Recommended workflow:

1. Filter by GO status
2. Filter by score > 0.6
3. Rank by combined_score
4. Validate top 3-5 with quantum coherence

See Also
--------

- :doc:`evidence_levels` - Evidence classification
- :doc:`/user_guide/tutorials/crispr_guide_design` - Design tutorial
