Evidence Levels
===============

PhaseLab v0.6.1 introduces evidence levels to classify guide validation strength.

The Three Levels
----------------

**Level A: Hardware-Validated**

Guides validated on real quantum hardware (IBM Quantum).

- Coherence measured on actual qubits
- Includes hardware noise effects
- Strongest evidence of reliability
- R values typically 0.84-0.97

**Level B: VQE-Simulated**

Guides validated using quantum simulation (VQE).

- Full quantum algorithm execution
- Simulated noise models
- Good evidence of reliability
- R values match hardware validation

**Level C: Heuristic Only**

Guides scored using fast heuristic proxy.

- Hamiltonian coefficient variance
- No quantum simulation
- Weak evidence of reliability
- R values cluster around 0.68-0.69

Why Evidence Levels Matter
--------------------------

**The Problem with Heuristic-Only**

Before v0.6.1, all coherence scores were weighted equally. This caused issues:

1. **False Confidence**: Heuristic R ≈ 0.68 looks good but doesn't reflect true coherence

2. **Poor Discrimination**: All guides score similarly in heuristic mode

3. **Misleading Rankings**: Heuristic coherence dominated ranking despite being a proxy

**The Solution**

Evidence levels ensure:

1. **Honest Reporting**: Users know what validation was performed

2. **Appropriate Weighting**: Heuristic scores have reduced influence

3. **Clear Recommendations**: Experimental validation prioritizes Level B+ guides

How Levels Affect Scoring
-------------------------

**Score Weight by Level**

============  ===============  ==================
Level         Coherence Weight Score Cap
============  ===============  ==================
A             0.30             None
B             0.30             None
C             0.05             0.85
============  ===============  ==================

**Example Impact**

Consider two guides with identical non-coherence scores of 0.70:

.. code-block:: python

   # Guide 1: Level B (quantum R = 0.92)
   score_1 = 0.70 * 0.70 + 0.30 * 0.92 = 0.49 + 0.28 = 0.77

   # Guide 2: Level C (heuristic R = 0.68)
   score_2 = 0.70 * 0.95 + 0.05 * 0.68 = 0.67 + 0.03 = 0.70
   # Capped at 0.85 max

Guide 1 (Level B) appropriately ranks higher.

Upgrading Evidence Levels
-------------------------

**From C to B**

Run quantum simulation:

.. code-block:: python

   from phaselab.crispr import compute_guide_coherence

   # Upgrade to Level B
   r_bar = compute_guide_coherence(guide_seq, mode="quantum")
   guide['evidence_level'] = 'B'
   guide['quantum_coherence'] = r_bar

**From B to A**

Validate on IBM Quantum hardware:

.. code-block:: python

   from phaselab.quantum import QuantumCoherenceValidator

   validator = QuantumCoherenceValidator(backend='ibm_torino')
   result = validator.validate(guide_seq)

   guide['evidence_level'] = 'A'
   guide['hardware_coherence'] = result['coherence']

Recommended Workflow
--------------------

1. **Initial Design** (Level C)

   - Fast screening of thousands of candidates
   - Use heuristic mode
   - Filter by GO status and basic criteria

2. **Candidate Refinement** (Level B)

   - Top 10-20 candidates
   - Run quantum simulation
   - Re-rank by quantum coherence

3. **Final Validation** (Level A)

   - Top 3-5 candidates
   - Hardware validation if available
   - Confirm experimental priorities

.. code-block:: python

   from phaselab.crispr import design_guides, compute_coherence_batch

   # Level C: Initial design
   guides = design_guides(sequence, tss_index)
   guides['evidence_level'] = 'C'

   # Filter
   candidates = guides[guides['go_no_go'] == 'GO'].head(20)

   # Level B: Quantum validation
   quantum_r = compute_coherence_batch(
       candidates['sequence'].tolist(),
       mode="quantum"
   )
   candidates['quantum_coherence'] = quantum_r
   candidates['evidence_level'] = 'B'

   # Level A: Hardware (if available)
   # ... IBM Quantum validation ...

Displaying Evidence Levels
--------------------------

Include evidence level in reports:

.. code-block:: python

   for idx, row in guides.iterrows():
       level = row.get('evidence_level', 'C')
       r_bar = row.get('quantum_coherence', row.get('coherence', 'N/A'))

       level_desc = {
           'A': 'Hardware-validated',
           'B': 'VQE-simulated',
           'C': 'Heuristic'
       }

       print(f"Guide: {row['sequence'][:15]}...")
       print(f"  Evidence: Level {level} ({level_desc[level]})")
       print(f"  Coherence: {r_bar}")
       print(f"  Score: {row['combined_score']:.3f}")

Best Practices
--------------

1. **Always report evidence level** in publications and reports

2. **Prioritize Level B+** for experimental validation

3. **Use Level C** only for initial screening

4. **Document validation method** for reproducibility

5. **Consider upgrading** critical candidates to Level A

See Also
--------

- :doc:`crispr_scoring` - Full scoring details
- :doc:`/user_guide/tutorials/coherence_modes` - Coherence mode guide
