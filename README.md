# PhaseLab

**Phase-coherence analysis framework for quantum, biological, and dynamical systems.**

[![PyPI version](https://badge.fury.io/py/phaselab.svg)](https://badge.fury.io/py/phaselab)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

PhaseLab implements the **Informational Relativity (IR) framework** for assessing simulation reliability across domains. It provides:

- **Quantum coherence metrics** (R̄, V_φ) validated on IBM Quantum hardware
- **CRISPR/CRISPRa guide RNA design** with multi-layer validation
- **Circadian clock modeling** for gene therapy dosage optimization

## Quick Start

### Installation

```bash
pip install phaselab

# With quantum computing support
pip install phaselab[quantum]

# Full installation
pip install phaselab[all]
```

### Design CRISPR Guides in 5 Lines

```python
from phaselab.crispr import design_guides

# Your promoter sequence (1kb upstream of TSS)
promoter = "ATGC..."  # Full sequence here

# Design and rank guide RNAs
guides = design_guides(
    sequence=promoter,
    tss_index=500,  # TSS position in sequence
)

# View top candidates
print(guides[['sequence', 'position', 'mit_score', 'coherence_R', 'go_no_go']])
```

### Simulate Circadian Clock (SMS Model)

```python
from phaselab.circadian import simulate_sms_clock, therapeutic_scan

# Simulate SMS patient (50% RAI1)
result = simulate_sms_clock(rai1_level=0.5)
print(f"Synchronization: {result['final_R_bar']:.3f}")
print(f"Classification: {result['classification']}")

# Find therapeutic window
scan = therapeutic_scan()
print(f"Optimal RAI1: {scan['optimal_level']}")
print(f"Required boost: +{scan['required_boost']*100:.0f}%")
```

### Compute Coherence Metrics

```python
from phaselab.core import coherence_score, go_no_go
import numpy as np

# From phase data (Kuramoto order parameter)
phases = np.array([0.1, 0.15, 0.12, 0.18])
R_bar = coherence_score(phases, mode='phases')
print(f"R̄ = {R_bar:.4f}, Status: {go_no_go(R_bar)}")

# From variance (IR formula: R̄ = exp(-V_φ/2))
V_phi = 0.5
R_bar = coherence_score(V_phi, mode='variance')
print(f"R̄ = {R_bar:.4f}")
```

## The IR Framework

PhaseLab is built on **Informational Relativity**, a framework that provides:

### Core Equation
```
R̄ = exp(-V_φ/2)
```

Where:
- **R̄** (R-bar): Coherence/order parameter [0, 1]
- **V_φ** (V-phi): Phase variance

### GO/NO-GO Threshold
```
R̄ > e⁻² ≈ 0.135 → GO (reliable)
R̄ < e⁻² ≈ 0.135 → NO-GO (unreliable)
```

This threshold has been validated on:
- IBM Quantum hardware (H₂ VQE: R̄ = 0.891)
- gRNA binding simulations (R̄ = 0.84)
- Circadian oscillator models

## CRISPR Pipeline Features

The `design_guides()` function provides:

| Layer | Method | Purpose |
|-------|--------|---------|
| **PAM Scanning** | NGG, NNGRRT, TTTV | Find Cas binding sites |
| **GC Content** | 40-70% filter | Optimal binding |
| **Thermodynamics** | SantaLucia ΔG | Binding energy |
| **MIT Score** | Position-weighted | Off-target specificity |
| **CFD Score** | Mismatch penalty | Cutting frequency |
| **Chromatin** | DNase HS model | Accessibility |
| **IR Coherence** | R̄ metric | Simulation reliability |

## Circadian Model Features

The SMS model includes:

- **Kuramoto oscillator base** - Phase coupling dynamics
- **PER delayed feedback** - Realistic negative loop
- **REV-ERBα/RORα modulation** - BMAL1 regulation
- **RAI1 dosage effects** - SMS-specific coupling
- **Therapeutic window analysis** - Find optimal boost

## Example: Smith-Magenis Syndrome Gene Therapy

This framework was developed to design CRISPRa guides for RAI1 upregulation in SMS:

```python
from phaselab.crispr import design_guides
from phaselab.circadian import therapeutic_scan
from phaselab.io import export_crispor_batch

# 1. Design guides for RAI1 promoter
rai1_promoter = """TGTCTCTTCCCACCAGGATGCC..."""  # 1kb sequence
guides = design_guides(rai1_promoter, tss_index=500)

# 2. Export for CRISPOR validation
export_crispor_batch(rai1_promoter, "crispor_input.fa")

# 3. Predict therapeutic window
scan = therapeutic_scan()
print(f"Target RAI1 boost: +{scan['required_boost']*100:.0f}%")

# 4. Top candidates
for i, row in guides.head(3).iterrows():
    print(f"{row['sequence']} | pos {row['position']} | R̄={row['coherence_R']:.3f} | {row['go_no_go']}")
```

**Result**: Hardware-validated gRNA `TACAGGAGCTTCCAGCGTCA` with:
- MIT specificity: 83
- CFD score: 93
- Zero off-targets ≤2 mismatches
- IBM Torino coherence: R̄ = 0.839

## Documentation

- **[API Guide](docs/API_GUIDE.md)** - Complete API reference with detailed examples
- **[Examples](docs/EXAMPLES.md)** - Practical code examples for common use cases
- **[SMS Gene Therapy Research](docs/SMS_GENE_THERAPY_RESEARCH.md)** - Full research paper on IBM Quantum-validated CRISPRa design for Smith-Magenis Syndrome

## Research Papers

Three publishable papers establishing PhaseLab and its applications:

| Paper | Title | Target Journals |
|-------|-------|-----------------|
| **[Paper 1](docs/papers/PAPER_1_PHASELAB_FRAMEWORK.md)** | PhaseLab: A Generalized Coherence Framework for Quantum-Biological Simulation | *Nature Computational Science*, *NPJ Quantum Information* |
| **[Paper 2](docs/papers/PAPER_2_CRISPRA_GRNA_DESIGN.md)** | Quantum-Informed CRISPRa gRNA Design for RAI1 Activation in SMS | *Nature Biotechnology*, *Nucleic Acids Research* |
| **[Paper 3](docs/papers/PAPER_3_CIRCADIAN_MODELING.md)** | Phase-Based Modeling of Circadian Dysregulation in SMS | *Cell Systems*, *eLife*, *Journal of Biological Rhythms* |

## API Reference

### Core (`phaselab.core`)

```python
from phaselab.core import (
    coherence_score,      # Compute R̄ from various inputs
    go_no_go,             # GO/NO-GO classification
    phase_variance,       # Compute V_φ from phases
    E_MINUS_2,            # e⁻² threshold constant
    build_pauli_hamiltonian,  # Generic Hamiltonian builder
)
```

### CRISPR (`phaselab.crispr`)

```python
from phaselab.crispr import (
    design_guides,        # Main pipeline
    GuideDesignConfig,    # Configuration dataclass
    find_pam_sites,       # PAM scanner
    gc_content,           # GC calculation
    delta_g_santalucia,   # Thermodynamic ΔG
    mit_specificity_score,  # MIT off-target score
)
```

### Circadian (`phaselab.circadian`)

```python
from phaselab.circadian import (
    simulate_sms_clock,   # SMS circadian model
    SMSClockParams,       # Model parameters
    therapeutic_scan,     # RAI1 level sweep
    classify_synchronization,  # R̄ to class
    kuramoto_order_parameter,  # Base Kuramoto R̄
)
```

## Validation

PhaseLab metrics have been validated against:

| System | Method | R̄ | Hardware |
|--------|--------|-----|----------|
| H₂ molecule | VQE | 0.891 | IBM Brisbane |
| gRNA binding | Hamiltonian sim | 0.839-0.854 | IBM Torino |
| Circadian clock | Kuramoto ODE | 0.73-0.99 | Classical |

## Citation

If you use PhaseLab in research, please cite:

```bibtex
@software{phaselab2025,
  author = {Vaca, Dylan},
  title = {PhaseLab: Phase-coherence analysis for quantum and biological systems},
  year = {2025},
  url = {https://github.com/followthesapper/phaselab}
}
```

## Contributing

Contributions welcome! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

MIT License - see [LICENSE](LICENSE) for details.

---

*Developed as part of the Informational Relativity research program.*
*Hardware validation: IBM Torino, December 2025.*
