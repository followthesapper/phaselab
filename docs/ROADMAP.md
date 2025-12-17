# PhaseLab Development Roadmap

This document outlines the planned development trajectory for PhaseLab, from the current release through full clinical integration.

## Version History

| Version | Status | Release Date | Summary |
|---------|--------|--------------|---------|
| 0.1.0 | Released | Jun 2025 | Core coherence, CRISPR pipeline, SMS model |
| 0.2.0 | Released | Jul 2025 | Protein folding, tissue-specific chromatin |
| 0.3.0 | Released | Aug 2025 | Multi-tissue circadian, base/prime editing |
| 0.4.0 | Released | Sep 2025 | Complete CRISPR toolkit, therapy dosage |
| 0.5.0 | Released | Oct 2025 | ATAC-seq integration, enhancer targeting |
| 0.6.0 | Released | Nov 2025 | ATLAS-Q integration, CRISPOR integration |
| 0.7.0 | Released | Nov 2025 | Enhanced pipeline, modality support |
| 0.8.0 | Released | Dec 2025 | Claim levels, fusion module |
| 0.9.0 | Released | Dec 2025 | SMS Trials Module, falsification tests |
| **1.0.0** | **Current** | Dec 2025 | **Spatial coherence paradigm, multi-domain modules** |
| 1.1.0 | Planned | 2026 | Full wet lab integration, clinical validation |

---

## v1.0.0 - Spatial Coherence Paradigm (Current Release)

**Status**: Released

### Major Paradigm Shift

**Guide-sequence coherence has been deprecated.** Experiments E200-E211 showed that computing coherence from guide sequences (GC content, thermodynamic properties) has no predictive value for experimental outcomes (r ≈ 0).

**Spatial coherence is the validated approach.** Experiments E213-E216 validated that measuring spatial coherence of response landscapes predicts perturbation reliability:
- r = -0.24 to -0.50 correlation with outcome variance
- 32-49% variance reduction when selecting from stable regions
- 115,251 sgRNAs tested across 6 genes

### New Modules

- **phaselab.landscapes** - Core perturbation-response data structures
  - `ResponseLandscape` - Generic position → response mapping
  - `CoherenceProfile` - Per-position coherence values
  - `StabilityClass` - STABLE, MIXED, AMPLIFYING, IRRELEVANT

- **phaselab.spatial** - E213-validated tiling screen analysis
  - `analyze_tiling_coherence()` - Full coherence analysis
  - `load_tiling_screen()` - Data loading utilities
  - Region classification and stable region identification

- **phaselab.surf** - CRISPR-SURF integration
  - `parse_surf_output()` - Parse SURF deconvolution output
  - `compute_surf_coherence()` - Coherence on deconvolved data
  - `SURFPipeline` - End-to-end analysis pipeline

- **phaselab.omics** - Genomics assay coherence
  - ATAC-seq: `analyze_atac_coherence()` - Stable accessibility regions
  - ChIP-seq: `analyze_chip_coherence()` - Stable binding sites
  - RNA-seq: `analyze_expression_coherence()` - Reliable expression changes

- **phaselab.microbio** - Microbial screen analysis
  - TnSeq: `analyze_tnseq_coherence()` - Essential domain identification
  - CRISPRi: `analyze_crispri_coherence()` - Bacterial screen analysis
  - Drug: `analyze_drug_coherence()` - Dose-response stability

- **phaselab.chem** - Chemical systems
  - Binding: `analyze_binding_coherence()` - Stable binding hot spots
  - Reaction: `analyze_reaction_coherence()` - Stable operating windows
  - HTS: `analyze_screening_coherence()` - Reliable hit identification

- **phaselab.protein.mutscan** - Mutational scanning
  - `analyze_mutscan_coherence()` - Functional domain identification
  - `local_coherence_profile()` - Per-residue coherence
  - `map_coherence_to_structure()` - PDB B-factor mapping

### Quantum Mode Configuration

New quantum execution modes:
- `QuantumMode.OFF` - Classical only (default, fastest)
- `QuantumMode.AUDIT` - Classical + quantum validation
- `QuantumMode.REQUIRED` - Quantum mandatory

```python
from phaselab.quantum import set_quantum_mode, QuantumMode
set_quantum_mode(QuantumMode.AUDIT)
```

### Breaking Changes

- `compute_coherence=True` in CRISPR pipeline is deprecated (no effect)
- `weight_coherence` defaults to 0.0 (was 1.0)
- Guide-sequence coherence functions emit deprecation warnings
- SMS trials now use spatial coherence by default

### API Quick Start

```python
from phaselab.landscapes import ResponseLandscape
from phaselab.spatial import analyze_tiling_coherence

# Create response landscape
landscape = ResponseLandscape(
    coords=[100, 200, 300, 400, 500],
    responses=[0.8, 0.75, 0.1, 0.82, 0.78],
)

# Analyze coherence
result = analyze_tiling_coherence(landscape)

# Find stable regions
for region in result.stable_regions:
    print(f"Stable: {region['start']}-{region['end']}, coh={region['coherence']:.3f}")
```

---

## v0.9.0 - SMS Trials Module (Released)

**Status**: Released

### Features
- **Complete SMS Therapeutic Trial Framework**
  - CRISPRa RAI1 activation trial with therapeutic window validation
  - CRISPRi modifier gene suppression trials (PER1, CRY1, CLOCK)
  - Knockout model validation trials (research use)
  - Base editing trials for point mutation correction
  - Prime editing trials for regulatory motif repair
  - Circadian rescue simulation with sleep/wake prediction
  - AAV delivery feasibility assessment for CNS
- **SMS Pipeline Orchestrator** - Integrated GO/NO-GO decision system
- **Falsification Test Framework** - Automated wet-lab validation experiment generation
- **Claim Level System** - STRONG_COMPUTATIONAL, CONTEXT_DEPENDENT, EXPLORATORY, UNKNOWN

### API Additions

```python
from phaselab.trials.sms import (
    SMSPipeline,
    SMSTrialConfig,
    run_sms_crispra_trial,
    run_sms_crispri_trial,
    run_sms_knockout_trial,
    run_sms_base_editing_trial,
    run_sms_prime_editing_trial,
    run_circadian_rescue_simulation,
    run_delivery_assessment,
)

# Run full SMS therapeutic pipeline
pipeline = SMSPipeline(config=SMSTrialConfig())
result = pipeline.run_full_pipeline()
print(f"GO/NO-GO: {result.overall_go_nogo}")
print(f"Falsification tests: {len(result.falsification_tests)}")
```

---

## v0.8.0 - Claim Levels & Fusion (Released)

**Status**: Released

### Features
- **Claim Level System**: Four-tier evidence classification
- **Fusion Module**: Multi-source data integration with uncertainty quantification
- **Virtual Assay Stack**: Enhanced guide scoring with tissue-specific modeling

---

## v0.7.0 - Enhanced Pipeline (Released)

**Status**: Released

### Features
- **Enhanced Pipeline**: `design_enhanced_guides()` with modality-specific scoring
- **Modality Support**: CRISPRa, CRISPRi, Knockout, Base Editing, Prime Editing
- **Tissue-Specific Scoring**: Brain, liver, blood, and muscle accessibility models

---

## v0.6.x - ATLAS-Q Integration (Released)

**Status**: Released

### Features
- ATLAS-Q tensor network backend integration
- CRISPOR integration with IR-enhanced off-target analysis
- Coherence mode parameter (heuristic vs quantum)
- Evidence levels (A/B/C classification)

---

## v0.2.0-0.5.0 - Foundation Releases (Released)

**Status**: Released

### Features Summary
- **v0.2.0**: Protein folding coherence, tissue-specific chromatin models
- **v0.3.0**: Multi-tissue circadian models, base/prime editing support
- **v0.4.0**: Complete CRISPR toolkit (knockout, CRISPRi), therapy dosage optimization
- **v0.5.0**: ATAC-seq integration, nucleosome occupancy, enhancer targeting, AAV delivery

---

## v1.0.0 - Clinical Integration (Planned)

**Status**: Planned
**Target**: 2026

### Goals
- Full wet lab validation pipeline
- Clinical trial support features
- Regulatory documentation helpers
- Production-ready reliability

### Planned Features

#### Wet Lab Integration
- Automated qPCR analysis
- Western blot quantification
- Flow cytometry integration
- Microscopy image analysis

#### Clinical Validation
- Patient stratification tools
- Outcome prediction models
- Safety monitoring dashboards
- Efficacy tracking

#### Regulatory Support
- GLP-compliant logging
- Audit trail generation
- Data provenance tracking
- Report generation

### Requirements for 1.0
- [ ] Published peer-reviewed validation
- [ ] At least one wet lab confirmation
- [ ] IRB-approved study integration
- [ ] Partner institution collaboration

---

## Future Directions (Post-1.0)

### Potential v1.x Features
- Additional gene targets (SHANK3, CHD8, NRXN1)
- In vivo delivery optimization
- AAV serotype selection
- Blood-brain barrier modeling

### Research Collaborations
- Academic partnerships for validation
- Pharmaceutical industry pilots
- Patient foundation collaborations

### Platform Extensions
- Cloud-based analysis service
- GUI application
- Jupyter notebook integrations
- Clinical decision support

---

## Contributing to the Roadmap

We welcome community input on the development roadmap. To suggest features:

1. Open an issue on GitHub with `[ROADMAP]` prefix
2. Describe the use case and proposed solution
3. Include any relevant references or prior art

Priority is given to features that:
- Have clear therapeutic applications
- Can be validated experimentally
- Extend the IR framework coherently
- Benefit multiple user groups

---

## Version Numbering

PhaseLab follows semantic versioning:

- **MAJOR** (1.x.x): Breaking API changes, clinical-grade release
- **MINOR** (x.1.x): New features, backward compatible
- **PATCH** (x.x.1): Bug fixes, documentation updates

Development versions use `-dev` suffix (e.g., `0.2.0-dev`).

---

*Last updated: December 2025*
*Maintainer: Dylan Vaca*
