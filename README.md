# PhaseLab

**Phase-coherence analysis framework for quantum, biological, and dynamical systems.**

[![PyPI version](https://badge.fury.io/py/phaselab.svg)](https://badge.fury.io/py/phaselab)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

PhaseLab implements the **Informational Relativity (IR) framework** for assessing simulation reliability across domains. It provides:

- **Quantum coherence metrics** (R̄, V_φ) validated on IBM Quantum hardware
- **Comprehensive CRISPR toolkit**:
  - **CRISPRa** - Transcriptional activation guide design
  - **CRISPRi** - Transcriptional interference/repression
  - **Knockout** - Cas9 cutting for gene disruption
  - **Prime editing** - pegRNA design for precise edits
  - **Base editing** - ABE (A→G) and CBE (C→T) single-nucleotide changes
- **Therapeutic dosage optimization** for haploinsufficiency disorders
- **Circadian clock modeling** for gene therapy timing
- **Gene target library** for disorders (RAI1, SCN2A, SHANK3, CHD8)

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

### Design CRISPR Knockout Guides

```python
from phaselab.crispr import design_knockout_guides

# Design guides to disrupt a gene
guides = design_knockout_guides(
    sequence=gene_sequence,
    cds_start=200,  # Start of coding sequence
)

# View top candidates with cutting efficiency and frameshift probability
print(guides[['sequence', 'cut_efficiency', 'frameshift_prob', 'go_no_go']])
```

### Design CRISPRi Guides (Gene Repression)

```python
from phaselab.crispr import design_crispri_guides

# Design guides for transcriptional repression
guides = design_crispri_guides(
    sequence=promoter_seq,
    tss_index=500,  # Transcription start site
)

# View candidates with repression efficiency
print(guides[['sequence', 'position', 'repression_efficiency', 'steric_hindrance']])
```

### Design Prime Editing pegRNAs

```python
from phaselab.crispr import design_prime_edit

# Design pegRNAs for precise A→G edit
pegrnas = design_prime_edit(
    sequence=target_region,
    edit_position=150,
    edit_from="A",
    edit_to="G",
)

print(pegrnas[['spacer', 'pbs_length', 'rt_length', 'combined_score']])
```

### Design Base Editing Guides

```python
from phaselab.crispr import design_abe_guides, design_cbe_guides

# ABE: A→G editing
abe_guides = design_abe_guides(sequence, target_position=100)

# CBE: C→T editing
cbe_guides = design_cbe_guides(sequence, target_position=100)

print(abe_guides[['sequence', 'target_in_window_pos', 'combined_efficiency']])
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

## What PhaseLab Is and Is Not

### What PhaseLab IS:
- **A prioritization tool** that helps rank guide candidates for experimental validation
- **A structural analysis layer** that examines off-target risk topology, not just counts
- **A complement to CRISPOR** that adds IR coherence metrics to standard specificity scores
- **A reliability assessment framework** using phase coherence to identify trustworthy simulations

### What PhaseLab is NOT:
- **NOT a predictor** of fold-change, editing efficiency, or wet-lab outcomes
- **NOT a replacement for CRISPOR** - we use CRISPOR's MIT/CFD scores as ground truth for off-targets
- **NOT a guarantee** - all top candidates require experimental validation
- **NOT specialized to any gene** - SMS/RAI1 is a case study, the framework is general

### The Core Value Proposition:
PhaseLab reduces experimental entropy by providing **structure-aware prioritization**. Instead of testing 50 guides blindly, PhaseLab helps you identify which 5-10 are worth testing first based on:
1. IR coherence (simulation reliability)
2. Off-target risk distribution (entropy)
3. Coherence contrast (ΔR̄)
4. Exonic risk flagging

## What's New in v0.6.1

- **Coherence Mode Parameter**: Choose between `mode="heuristic"` (fast, default) and `mode="quantum"` (VQE simulation)
- **Honest Coherence Weighting**: Heuristic coherence demoted to tie-breaker (0.05 weight vs 0.30 for quantum)
- **Two-Stage Scoring**: Hard safety gates (Stage 1) + soft ranking (Stage 2)
- **Risk Mass Metrics**: `risk_mass_close`, `risk_mass_exonic`, `tail_risk_score`
- **Concentration Measures**: Gini coefficient and Herfindahl-Hirschman Index for off-target distribution
- **Evidence Levels**: A/B/C classification prevents over-trust in unvalidated guides
- **Score Capping**: Level C (unvalidated) guides capped at 0.20 to prevent misleading rankings

## What's New in v0.6.0

- **ATLAS-Q Integration**: All CRISPR modules now use unified coherence computation via ATLAS-Q backend
- **IR-Enhanced Off-Target Analysis**: CRISPOR integration with off-target entropy, coherence contrast (ΔR̄), and energy spectrum
- **Off-Target Clustering**: Detection of dangerous off-target "families" indicating systematic vulnerabilities
- **Region Difficulty Index**: "Soup index" for GC-rich/repetitive regions to normalize expectations
- **Real Circular Statistics**: Proper phase variance calculation replacing heuristic estimates
- **Exonic Risk Flagging**: Automatic detection of dangerous off-targets in coding regions
- **Score Adjustment System**: IR-based penalties/bonuses for off-target landscape quality

## CRISPR Toolkit (v0.4.0+)

PhaseLab provides a complete genome engineering toolkit:

| Module | Function | Use Case |
|--------|----------|----------|
| **CRISPRa** | `design_guides()` | Transcriptional activation |
| **CRISPRi** | `design_crispri_guides()` | Transcriptional repression |
| **Knockout** | `design_knockout_guides()` | Gene disruption via DSB |
| **Prime Editing** | `design_prime_edit()` | Precise insertions/deletions |
| **Base Editing** | `design_abe_guides()`, `design_cbe_guides()` | Single-nucleotide changes |

### Scoring Layers

All CRISPR modules include multi-layer scoring:

| Layer | Method | Purpose |
|-------|--------|---------|
| **PAM Scanning** | NGG, NNGRRT, TTTV | Find Cas binding sites |
| **GC Content** | 40-70% filter | Optimal binding |
| **Thermodynamics** | SantaLucia ΔG | Binding energy |
| **MIT Score** | Position-weighted | Off-target specificity |
| **CFD Score** | Mismatch penalty | Cutting frequency |
| **Chromatin** | DNase HS model | Accessibility |
| **IR Coherence** | R̄ metric | Simulation reliability |
| **Off-target Entropy** | Shannon entropy | Risk distribution (v0.6.0) |
| **Coherence Contrast** | ΔR̄ | Off-target competitiveness (v0.6.0) |

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

## Gene Targets

PhaseLab includes pre-configured targets for haploinsufficiency disorders:

| Target | Disease | Therapeutic Window | Status |
|--------|---------|-------------------|--------|
| **RAI1** | Smith-Magenis Syndrome | 70-110% | Hardware validated |
| **SCN2A** | Autism-linked NDD, epilepsy | 65-115% | Hardware validated |
| **SHANK3** | Phelan-McDermid Syndrome | 60-110% | Hardware validated |
| **CHD8** | CHD8-related ASD | 65-115% | Hardware validated |

```python
from phaselab.targets import load_target_config, list_available_targets

# List all targets
print(list_available_targets())  # ['RAI1', 'SCN2A']

# Load SCN2A configuration
scn2a = load_target_config("SCN2A")
print(f"Gene: {scn2a.gene_symbol}")
print(f"Disease: {scn2a.disease}")
print(f"TSS: chr{scn2a.chrom}:{scn2a.tss_genomic}")
```

See [Target Library Documentation](docs/TARGETS.md) for adding new targets.

## Documentation

- **[API Guide](docs/API_GUIDE.md)** - Complete API reference with detailed examples
- **[Examples](docs/EXAMPLES.md)** - Practical code examples for common use cases
- **[Target Library](docs/TARGETS.md)** - Gene target configurations for CRISPRa experiments
- **[SMS Gene Therapy Research](docs/SMS_GENE_THERAPY_RESEARCH.md)** - IBM Quantum-validated CRISPRa design for Smith-Magenis Syndrome (RAI1)
- **[SCN2A Gene Therapy Research](docs/SCN2A_GENE_THERAPY_RESEARCH.md)** - IBM Quantum-validated CRISPRa design for Autism-linked NDD (SCN2A)

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
    # CRISPRa (activation)
    design_guides,
    GuideDesignConfig,

    # CRISPRi (repression)
    design_crispri_guides,
    CRISPRiConfig,

    # Knockout
    design_knockout_guides,
    KnockoutConfig,
    cut_efficiency_score,
    frameshift_probability,

    # Prime editing
    design_prime_edit,
    PrimeEditConfig,
    design_pbs,
    design_rt_template,

    # Base editing
    design_abe_guides,
    design_cbe_guides,
    BaseEditConfig,
    editing_efficiency_at_position,

    # Scoring utilities
    find_pam_sites,
    gc_content,
    delta_g_santalucia,
    mit_specificity_score,

    # ATLAS-Q Coherence (v0.6.0)
    compute_guide_coherence,
    compute_guide_coherence_with_details,
    is_guide_coherent,
    get_coherence_method,
)
```

### CRISPOR Integration (`phaselab.integrations.crispor`)

```python
from phaselab.integrations.crispor import (
    # CRISPOR Client
    CrisporClient,
    CrisporConfig,

    # Pipeline
    design_guides_with_crispor,
    GuideCandidate,
    ScoringWeights,

    # IR-Enhanced Off-Target Analysis (v0.6.0)
    OffTargetSite,
    OffTargetIRAnalysis,
    analyze_offtarget_landscape,
    compute_offtarget_entropy,
    compute_coherence_contrast,
    compute_ir_enhanced_score,
)
```

### Therapy (`phaselab.therapy`)

```python
from phaselab.therapy import (
    TherapeuticWindow,
    optimize_dosage,
    validate_therapeutic_level,
    estimate_expression_change,
    predict_therapeutic_efficacy,
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

PhaseLab metrics have been validated against IBM Quantum hardware:

| Module | System | R̄ Range | Hardware | Status |
|--------|--------|---------|----------|--------|
| Core | H₂ molecule VQE | 0.891 | IBM Brisbane | ✓ GO |
| CRISPRa | RAI1 gRNA (SMS) | 0.839 | IBM Torino | ✓ GO |
| CRISPRa | SCN2A gRNA (Autism) | 0.970 | IBM Torino | ✓ GO |
| **Knockout** | Cut efficiency | 0.44-0.63 | IBM Torino | ✓ GO |
| **CRISPRi** | Repression scoring | 0.62-0.89 | IBM Torino | ✓ GO |
| **Prime Editing** | pegRNA design | 0.81-0.93 | IBM Torino | ✓ GO |
| **Base Editing** | ABE/CBE guides | 0.62-0.94 | IBM Torino | ✓ GO |
| **Therapy** | Dosage optimization | 0.65-0.98 | IBM Torino | ✓ GO |
| Circadian | Kuramoto ODE | 0.73-0.99 | Classical | ✓ GO |

## Known Limitations

PhaseLab has inherent limitations that users should understand:

### Biological Limitations
- **No efficacy prediction**: We cannot predict fold-change in gene expression or editing efficiency
- **Context-dependent effects**: Chromatin state, cell type, and delivery method all affect real-world outcomes
- **GC-rich regions remain hard**: High-GC promoters (like RAI1) have many off-targets regardless of tool used

### Technical Limitations
- **Heuristic coherence**: Guide coherence uses Hamiltonian coefficient variance as a proxy, not actual VQE measurement
- **Hardware R̄ differs from heuristic R̄**: IBM Quantum validation gave R̄ = 0.84-0.97, heuristic gives ~0.68-0.69
- **CRISPOR dependency**: Off-target counts and MIT/CFD scores come from CRISPOR's empirical models

### Methodological Limitations
- **Dual leaderboard required**: Unvalidated guides must not be compared directly to validated guides
- **Experimental validation required**: All top candidates must be tested in wet lab
- **Region-specific thresholds needed**: A "good" score in a GC-soup region differs from a balanced region

### What This Means in Practice
PhaseLab is best used as a **prioritization filter**, not an oracle:
1. Generate candidates with `design_guides()`
2. Validate top candidates on CRISPOR
3. Use IR-enhanced analysis to identify structural risks
4. Select 5-10 guides for wet-lab testing
5. Iterate based on experimental results

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
*Version 0.6.0: ATLAS-Q integration with unified coherence, IR-enhanced off-target analysis, and CRISPOR integration.*
