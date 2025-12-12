PhaseLab Documentation
======================

PhaseLab is a comprehensive phase-coherence analysis framework for quantum, biological, and dynamical systems. Built on the Informational Relativity (IR) framework, PhaseLab provides quantum coherence metrics, CRISPR guide RNA design pipelines, circadian clock modeling, and therapeutic dosage optimization.

Version 0.6.1 (December 2025)

**Phase-Coherence Meets Precision Medicine**

PhaseLab bridges theoretical physics and biotechnology by applying IR coherence metrics to guide RNA design and gene therapy optimization. Every guide RNA candidate is validated using the universal GO/NO-GO threshold (R > e^-2), ensuring only reliable candidates proceed to experimental validation.

Key Capabilities
----------------

**Quantum Coherence (v0.6.0+)**

- ATLAS-Q integration for advanced quantum simulation
- IR measurement grouping with 5x variance reduction
- Real circular statistics coherence (replaces heuristic)
- Coherence-aware VQE optimization
- Optional GPU acceleration via Triton kernels

**CRISPR Design (v0.6.1)**

- **Coherence Mode Parameter**: ``mode="heuristic"`` (fast) vs ``mode="quantum"`` (VQE)
- **Honest Coherence Weighting**: Heuristic demoted to tie-breaker (0.05 vs 0.30)
- **Two-Stage Scoring**: Hard safety gates + soft ranking
- **Risk Mass Metrics**: ``risk_mass_close``, ``risk_mass_exonic``, ``tail_risk_score``
- **Evidence Levels**: A/B/C classification for validation status
- **Score Capping**: Unvalidated guides capped to prevent misleading rankings

**Complete CRISPR Toolkit**

- CRISPRa guide design for transcriptional activation
- CRISPRi guide design for transcriptional interference
- CRISPR knockout guide design for gene disruption
- Prime editing pegRNA design for precise edits
- Base editing guide design (ABE/CBE) for single-nucleotide changes
- Multi-guide synergy modeling for combinatorial CRISPR
- Enhancer targeting for CRISPRa

**Biological Modeling**

- Circadian clock modeling (SMS gene therapy research)
- Multi-tissue circadian models with inter-tissue coupling
- Drug response modeling and chronotherapy optimization
- Protein folding coherence assessment
- Tissue-specific chromatin accessibility models

Quick Example
-------------

.. code-block:: python

   import phaselab as pl
   from phaselab.crispr import design_guides

   # Compute coherence from phase data
   phases = [0.1, 0.15, 0.12, 0.08, 0.11]
   R_bar = pl.coherence_score(phases)
   print(f"Coherence: {R_bar:.4f}")
   print(f"Status: {pl.go_no_go(R_bar)}")

   # Design CRISPR guides with coherence validation
   sequence = "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCG..."
   guides = design_guides(sequence, tss_index=100)

   # v0.6.1: Use coherence mode for honest validation
   from phaselab.crispr import compute_guide_coherence
   r_bar = compute_guide_coherence("ATCGATCGATCGATCGATCG", mode="quantum")

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   user_guide/tutorials/index
   user_guide/howtos/index
   user_guide/explanations/index

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   reference/index

.. toctree::
   :maxdepth: 1
   :caption: Research Applications

   research/index

.. toctree::
   :maxdepth: 1
   :caption: Developer Documentation

   developer/index

.. toctree::
   :maxdepth: 1
   :caption: Additional Resources

   examples/index
   faq
   changelog
   citing

Quick Links
-----------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
* `GitHub Repository <https://github.com/followthesapper/phaselab>`_
* `Issue Tracker <https://github.com/followthesapper/phaselab/issues>`_
* `PyPI Package <https://pypi.org/project/phaselab/>`_

Citation
--------

.. code-block:: bibtex

   @software{phaselab2025,
     title={PhaseLab: Phase-Coherence Analysis Framework},
     author={Vaca, Dylan},
     year={2025},
     url={https://github.com/followthesapper/phaselab},
     version={0.6.1}
   }

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
