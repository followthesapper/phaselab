# PhaseLab Development Roadmap

This document outlines the planned development trajectory for PhaseLab, from the current release through full clinical integration.

## Version History

| Version | Status | Release Date | Summary |
|---------|--------|--------------|---------|
| 0.1.0 | Released | Dec 2025 | Core coherence, CRISPR pipeline, SMS model |
| 0.2.0 | In Development | Q1 2026 | Protein folding, expanded chromatin models |
| 0.3.0 | Planned | Q2 2026 | Multi-tissue circadian, drug response modeling |
| 1.0.0 | Planned | 2026 | Full wet lab integration, clinical validation |

---

## v0.1.0 - Foundation (Current Release)

**Status**: Released

### Features
- Core IR coherence metrics (R̄, V_φ)
- GO/NO-GO classification with e⁻² threshold
- SpCas9 gRNA design pipeline
- MIT and CFD scoring
- Basic chromatin accessibility model
- SMS circadian clock model (Kuramoto-based)
- Gene target library (RAI1, SCN2A)

### Validated Targets
| Target | Disease | Hardware | Status |
|--------|---------|----------|--------|
| RAI1 | Smith-Magenis Syndrome | IBM Torino | Validated |
| SCN2A | Autism-linked NDD | IBM Torino | Validated |

---

## v0.2.0 - Protein & Chromatin (In Development)

**Status**: In Development
**Branch**: `feature/v0.2.0-protein-chromatin`

### New Features

#### Protein Folding Coherence (`phaselab.protein`)
- Ramachandran angle coherence (φ, ψ as literal phases)
- Contact map consistency scoring
- Secondary structure coherence
- AlphaFold pLDDT integration
- GO/NO-GO for structure predictions

#### Molecular Dynamics Coherence
- RMSD trajectory stability
- Ensemble convergence assessment
- Block averaging for sampling quality
- RMSF-based flexibility coherence

#### Binding Prediction Coherence
- Docking pose clustering
- Affinity prediction confidence
- Interaction fingerprint consistency

#### Expanded Chromatin Models (`phaselab.chromatin`)
- Tissue-specific accessibility models:
  - Brain (cortex, hippocampus)
  - Liver (hepatocyte)
  - Blood (PBMC, T cells)
  - Muscle, heart, lung, kidney
- ENCODE data integration framework
- PsychENCODE brain-specific profiles
- Chromatin-weighted CRISPR scoring

### API Additions

```python
# Protein folding
from phaselab.protein import compute_structure_coherence, alphafold_coherence
result = alphafold_coherence("prediction.pdb")
print(f"Structure R̄: {result.R_bar}, Status: {result.go_no_go}")

# Chromatin accessibility
from phaselab.chromatin import brain_accessibility, combine_with_crispr
score = brain_accessibility("chr2", 165239414, region="cortex")
print(f"Accessibility: {score.accessibility}, State: {score.state}")
```

---

## v0.3.0 - Multi-Tissue & Drug Response (Planned)

**Status**: Planned
**Target**: Q2 2026

### Planned Features

#### Multi-Tissue Circadian Models
- Tissue-specific clock gene expression
- Inter-tissue phase coupling
- Chronotherapy optimization
- Shift work / jet lag modeling

#### Drug Response Modeling
- Pharmacokinetic phase dynamics
- Drug-clock interactions
- Dosing schedule optimization
- Individual variability modeling

#### Expanded CRISPR Systems
- SaCas9 (NNGRRT PAM)
- Cas12a (TTTV PAM)
- Base editors (ABE, CBE)
- Prime editing support

#### Enhanced Chromatin
- Single-cell ATAC-seq integration
- Cell-state specific models
- Dynamic chromatin (time-varying)

### Planned API

```python
# Multi-tissue circadian
from phaselab.circadian import MultiTissueClock
clock = MultiTissueClock(tissues=["liver", "brain", "muscle"])
result = clock.simulate(drug_timing="08:00")

# Drug response
from phaselab.drug import DoseOptimizer
optimizer = DoseOptimizer(target="RAI1", tissue="brain")
schedule = optimizer.optimize_timing()
```

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
