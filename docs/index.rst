PhaseLab Documentation
======================

PhaseLab is a comprehensive phase-coherence analysis framework for quantum, biological, and dynamical systems. Built on the Informational Relativity (IR) framework, PhaseLab provides quantum coherence metrics, CRISPR guide RNA design pipelines, circadian clock modeling, and therapeutic dosage optimization.

Version 0.9.3 (December 2025)

**Phase-Coherence Meets Precision Medicine**

PhaseLab bridges theoretical physics and biotechnology by applying IR coherence metrics to guide RNA design and gene therapy optimization. Every guide RNA candidate is validated using the universal GO/NO-GO threshold (R > e^-2), ensuring only reliable candidates proceed to experimental validation.

Key Capabilities
----------------

**CRISPRa Binding Register Model (v0.9.3)**

- **NucleaseRole enum**: Explicit ``BINDING`` vs ``CUTTING`` mode for PAM stringency
- **Relaxed PAM patterns**: SaCas9 NNGRRN (binding) vs NNGRRT (cutting)
- **Sliding binding register**: ±2bp enumeration captures guides rigid anchoring misses
- **Configurable guide length**: Override defaults for literature reproduction
- **Validated**: Chang et al. 2022 sg2 winner correctly recovered at TSS-80

**Guide Enumeration & Policy System (v0.9.2)**

- Region declaration with multi-TSS support
- PAM scanning for SpCas9, SaCas9, Cas12a
- Dominance-based ranking: lexicographic sort on (0mm, 1mm, 2mm) off-targets
- Policy system: ``CUTTING_STRICT``, ``BINDING_STRICT``, ``EXPLORATORY``
- Tier system: A (0/0/0), B (0/0/1-2), C (other)
- Reproducibility manifests for audit trails

**SMS Trials Module (v0.9.0)**

- Complete therapeutic trial framework for Smith-Magenis Syndrome
- CRISPRa RAI1 activation with therapeutic window validation (70-110%)
- CRISPRi modifier gene suppression (PER1, CRY1, CLOCK)
- Base/prime editing trials for mutation correction
- GO/NO-GO decision system with claim level propagation

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

   from phaselab.crispr import design_crispra_guides, Nuclease

   # v0.9.3: Design CRISPRa guides with binding mode
   result = design_crispra_guides(
       gene_symbol="Rai1",
       promoter_sequence=promoter_seq,
       tss_position=600,
       nuclease=Nuclease.SACAS9,
       relaxed_pam=True,    # BINDING mode (default)
       guide_length=20,     # Override default 21bp
   )

   # Access ranked guides by tier
   for guide in result.tier_a_guides[:5]:
       print(f"{guide['sequence']} TSS{guide['tss_relative_position']:+d}")

   # Coherence validation
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
     version={0.9.3}
   }

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
