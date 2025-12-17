PhaseLab Documentation
======================

**A reliability layer for perturbation science.**

PhaseLab assesses whether experimental results from perturbation experiments (CRISPR, mutagenesis, chemical modification) will be reproducible—before you run the experiment.

Version 1.0.0 (December 2025)

The Core Insight
----------------

**Spatial coherence of response landscapes predicts perturbation reliability.**

Validated across 115,251 sgRNAs (6 genes):

- Correlation with outcome variance: r = -0.24 to -0.50
- Variance reduction in stable regions: 32-49%

Key Capabilities
----------------

**Spatial Coherence Paradigm (v1.0.0)**

This release represents a fundamental paradigm shift based on experimental validation:

- **E200-E211**: Guide-sequence coherence does NOT work (r ≈ 0)
- **E213-E216**: Spatial coherence of response landscapes DOES work

The key insight: *"The guide is the probe, not the structure."* Coherence measures the system's response consistency, not properties of the perturbation itself.

**New Domain Modules**

- ``phaselab.landscapes`` - Core perturbation-response data structures
- ``phaselab.spatial`` - E213-validated tiling screen analysis
- ``phaselab.protein.mutscan`` - Deep mutational scanning coherence
- ``phaselab.protein.folding`` - Structure prediction quality control
- ``phaselab.chem`` - Binding landscape and reaction optimization
- ``phaselab.omics`` - ATAC-seq, ChIP-seq, RNA-seq reliability
- ``phaselab.microbio`` - TnSeq and bacterial CRISPRi essentiality

**Two-Stage Framework**

- **Stage I**: Feasibility inference (pre-tiling) using structural priors
- **Stage II**: Landscape resolution with minimum viable tiling (16-20 perturbations)

**Claim Levels**

Honest uncertainty reporting:

- ``UNKNOWN`` - Cannot assess reliability
- ``EXPLORATORY`` - Structural priors only
- ``CONTEXT_DEPENDENT`` - Single tiling dataset
- ``STRONG_COMPUTATIONAL`` - Cross-validated

Quick Example
-------------

.. code-block:: python

   from phaselab.spatial import analyze_tiling_coherence
   from phaselab.landscapes import ResponseLandscape

   # Your tiling screen data
   landscape = ResponseLandscape(
       coords=[100, 150, 200, 250, 300],
       responses=[0.8, 0.75, 0.1, 0.82, 0.78],
   )

   # Analyze spatial coherence
   result = analyze_tiling_coherence(landscape)

   # Find stable regions
   for region in result.stable_regions:
       print(f"Stable: {region['start']}-{region['end']}")
       print(f"  Coherence: {region['coherence']:.3f}")

Quantum Mode
------------

For most use cases, quantum mode should be OFF (default):

.. code-block:: python

   from phaselab.quantum import set_quantum_mode, QuantumMode

   set_quantum_mode(QuantumMode.OFF)       # Classical only (fastest)
   set_quantum_mode(QuantumMode.AUDIT)     # Validation subset
   set_quantum_mode(QuantumMode.REQUIRED)  # Research only

**Rule:** If a classical experiment can falsify a claim, quantum is optional.

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
     title={PhaseLab: A Reliability Layer for Perturbation Science},
     author={Vaca, Dylan},
     year={2025},
     url={https://github.com/followthesapper/phaselab},
     version={1.0.0}
   }

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
