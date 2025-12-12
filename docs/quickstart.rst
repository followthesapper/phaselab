Quick Start
===========

This guide covers the essential PhaseLab workflows in 5 minutes.

Core Coherence
--------------

PhaseLab's foundation is the Informational Relativity coherence metric:

.. code-block:: python

   import phaselab as pl

   # Compute coherence from phase data
   phases = [0.1, 0.15, 0.12, 0.08, 0.11]
   R_bar = pl.coherence_score(phases)
   print(f"Coherence R: {R_bar:.4f}")

   # GO/NO-GO classification (threshold: e^-2 ~ 0.135)
   status = pl.go_no_go(R_bar)
   print(f"Status: {status}")  # "GO" or "NO-GO"

   # Phase variance (V_phi)
   V_phi = pl.phase_variance(phases)
   print(f"Phase variance: {V_phi:.6f}")

The universal GO/NO-GO threshold is e^-2 (approximately 0.1353). Systems with R > e^-2 are considered coherent ("GO"), while those below are unreliable ("NO-GO").

CRISPR Guide Design
-------------------

Design guide RNAs with coherence validation:

.. code-block:: python

   from phaselab.crispr import design_guides, GuideDesignConfig

   # Target sequence with TSS
   sequence = """
   ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
   AGGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
   """
   tss_index = 50  # Transcription start site

   # Design guides
   guides = design_guides(sequence, tss_index)

   # View top candidates
   print(guides[['sequence', 'position', 'gc_content', 'combined_score', 'go_no_go']].head())

Coherence Modes (v0.6.1)
------------------------

PhaseLab v0.6.1 introduces honest coherence computation with two modes:

.. code-block:: python

   from phaselab.crispr import compute_guide_coherence, get_coherence_eligibility_info

   guide = "ATCGATCGATCGATCGATCG"

   # HEURISTIC mode (default, fast ~0.1ms)
   # Uses Hamiltonian coefficient variance as proxy
   # R ~ 0.68-0.69, use as tie-breaker only
   r_heuristic = compute_guide_coherence(guide, mode="heuristic")
   print(f"Heuristic R: {r_heuristic:.4f}")

   # QUANTUM mode (slow ~100-500ms, research-grade)
   # Runs actual VQE simulation
   # R ~ 0.84-0.97, matches hardware validation
   r_quantum = compute_guide_coherence(guide, mode="quantum")
   print(f"Quantum R: {r_quantum:.4f}")

   # Check what method will be used
   info = get_coherence_eligibility_info(mode="quantum")
   print(f"Method: {info['method']}")
   print(f"ATLAS-Q active: {info['acceleration_active']}")

Evidence Levels
---------------

v0.6.1 classifies guides by validation evidence:

- **Level A**: Hardware-validated on IBM Quantum (strongest evidence)
- **Level B**: VQE-simulated with quantum mode (good evidence)
- **Level C**: Heuristic only (use as tie-breaker)

.. code-block:: python

   from phaselab.integrations.crispor import analyze_offtarget_landscape

   # Get detailed analysis with evidence levels
   result = analyze_offtarget_landscape(
       guide_seq="ATCGATCGATCGATCGATCG",
       offtarget_sites=[...],
       coherence_mode="quantum"
   )

   print(f"Evidence level: {result.evidence_level}")
   print(f"Coherence: {result.coherence:.4f}")

Different CRISPR Modalities
---------------------------

PhaseLab supports the complete CRISPR toolkit:

.. code-block:: python

   from phaselab.crispr import (
       design_guides,           # CRISPRa (activation)
       design_crispri_guides,   # CRISPRi (interference)
       design_knockout_guides,  # Knockout
       design_prime_edit,       # Prime editing
       design_base_edit_guides, # Base editing (ABE/CBE)
   )

   # CRISPRa for gene activation
   crispra_guides = design_guides(sequence, tss_index)

   # CRISPRi for gene repression
   crispri_guides = design_crispri_guides(sequence, tss_index)

   # Knockout for gene disruption
   knockout_guides = design_knockout_guides(sequence, exon_start=100, exon_end=200)

   # Prime editing for precise edits
   pegRNA = design_prime_edit(
       sequence,
       edit_position=150,
       edit_type="substitution",
       new_base="G"
   )

   # Base editing for single nucleotide changes
   abe_guides = design_base_edit_guides(
       sequence,
       target_position=150,
       editor_type="ABE"  # or "CBE"
   )

Circadian Clock Modeling
------------------------

Model SMS circadian disruption for gene therapy research:

.. code-block:: python

   from phaselab.circadian import simulate_sms_clock

   # Simulate SMS condition (50% RAI1)
   baseline = simulate_sms_clock(rai1_level=0.5)
   print(f"SMS Period: {baseline['period']:.2f} hours")
   print(f"SMS Coherence: {baseline['coherence']:.4f}")

   # Simulate with CRISPRa restoration (85% RAI1)
   treated = simulate_sms_clock(rai1_level=0.85)
   print(f"Treated Period: {treated['period']:.2f} hours")
   print(f"Treated Coherence: {treated['coherence']:.4f}")

Quantum Coherence from Hamiltonians
-----------------------------------

For advanced quantum applications:

.. code-block:: python

   from phaselab.quantum import (
       compute_coherence_from_hamiltonian,
       compute_coherence_from_expectations,
       is_atlas_q_available,
   )
   import numpy as np

   # Check ATLAS-Q availability
   print(f"ATLAS-Q: {is_atlas_q_available()}")

   # Compute coherence from Hamiltonian coefficients
   coefficients = np.array([0.5, 0.3, 0.2, 0.1])
   result = compute_coherence_from_hamiltonian(
       coefficients,
       use_atlas_q=True
   )
   print(f"R: {result.R_bar:.4f}, V_phi: {result.V_phi:.4f}")
   print(f"Status: {'GO' if result.is_go else 'NO-GO'}")

Next Steps
----------

- :doc:`user_guide/tutorials/index` - In-depth tutorials
- :doc:`reference/index` - Complete API reference
- :doc:`examples/index` - Working code examples
- :doc:`research/index` - Research applications
