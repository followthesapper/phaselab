# PhaseLab API Guide

**Phase-coherence analysis framework for quantum, biological, and dynamical systems**

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Core Module](#core-module)
   - [Coherence Score](#coherence-score)
   - [Phase Variance](#phase-variance)
   - [GO/NO-GO Classification](#gono-go-classification)
4. [CRISPR Module](#crispr-module)
   - [PAM Scanning](#pam-scanning)
   - [Guide Scoring](#guide-scoring)
   - [Design Pipeline](#design-pipeline)
5. [Circadian Module](#circadian-module)
   - [SMS Clock Model](#sms-clock-model)
   - [Kuramoto Oscillators](#kuramoto-oscillators)
6. [Quantum Integration](#quantum-integration)
7. [Complete Examples](#complete-examples)
8. [API Reference](#api-reference)

---

## Installation

### Basic Installation

```bash
pip install phaselab
```

### With Quantum Support (IBM Quantum)

```bash
pip install phaselab[quantum]
```

### With Plotting Support

```bash
pip install phaselab[plotting]
```

### Full Installation

```bash
pip install phaselab[all]
```

### Development Installation

```bash
git clone https://github.com/dylanvaca/phaselab.git
cd phaselab
pip install -e ".[dev]"
```

---

## Quick Start

```python
import phaselab as pl

# Compute coherence from phase data
phases = [0.1, 0.15, 0.12, 0.08, 0.11]
R_bar = pl.coherence_score(phases)
print(f"Coherence: {R_bar:.4f}")

# GO/NO-GO classification
status = pl.go_no_go(R_bar)
print(f"Status: {status}")

# Design CRISPR guides
from phaselab.crispr import design_guides

sequence = "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCG..."
guides = design_guides(sequence, tss_index=100)
print(guides[['sequence', 'position', 'gc_content', 'combined_score']])
```

---

## Core Module

The core module implements the Informational Relativity coherence framework.

### Coherence Score

The coherence score R̄ quantifies phase synchronization:

$$\bar{R} = e^{-V_\phi / 2}$$

where $V_\phi$ is the phase variance.

```python
from phaselab import coherence_score

# From raw phase values (radians)
phases = [0.1, 0.15, 0.12, 0.08, 0.11]
R_bar = coherence_score(phases, mode='phases')
# Returns: 0.9995

# From quantum expectation values
expectations = [0.85, 0.82, 0.88, 0.79]
R_bar = coherence_score(expectations, mode='expectations')
# Returns: 0.9823

# From direct variance value
R_bar = coherence_score(0.05, mode='variance')
# Returns: 0.9753

# Auto-detect mode (default)
R_bar = coherence_score(phases)  # Automatically detects list of phases
```

**Parameters:**
- `data`: Input data (phases, expectations, or variance)
- `mode`: One of `'auto'`, `'phases'`, `'expectations'`, `'variance'`

**Returns:** Float between 0 and 1

### Phase Variance

Compute phase variance directly:

```python
from phaselab import phase_variance

phases = [0.1, 0.15, 0.12, 0.08, 0.11]
V_phi = phase_variance(phases)
# Returns: 0.00052
```

### GO/NO-GO Classification

Binary classification based on the e⁻² threshold (≈0.1353):

```python
from phaselab import go_no_go
from phaselab.core.constants import E_MINUS_2

# Basic classification
status = go_no_go(0.85)  # Returns: "GO"
status = go_no_go(0.10)  # Returns: "NO-GO"

# Custom threshold
status = go_no_go(0.45, threshold=0.5)  # Returns: "NO-GO"

# Get detailed classification
from phaselab.core.coherence import classify_coherence

category = classify_coherence(0.85)  # Returns: "EXCELLENT"
category = classify_coherence(0.45)  # Returns: "MODERATE"
category = classify_coherence(0.10)  # Returns: "CRITICAL"
```

**Classification Thresholds:**
| R̄ Range | Category |
|---------|----------|
| > 0.8 | EXCELLENT |
| 0.5 - 0.8 | GOOD |
| e⁻² - 0.5 | MODERATE |
| 0.05 - e⁻² | SEVERE |
| < 0.05 | CRITICAL |

---

## CRISPR Module

The CRISPR module provides tools for CRISPRa guide RNA design with phase-coherence validation.

### PAM Scanning

Find all PAM sites in a sequence:

```python
from phaselab.crispr import find_pam_sites

sequence = "ATGCGATCGATCGAGGCGATCGATCGATCGATCGATCGATCG"

# Find NGG PAM sites (default for SpCas9)
sites = find_pam_sites(sequence, pam="NGG")
# Returns: [(position, strand, guide_sequence), ...]

# Find multiple PAM types
sites_ngg = find_pam_sites(sequence, pam="NGG")   # SpCas9
sites_nag = find_pam_sites(sequence, pam="NAG")   # SpCas9 alternate
sites_tttn = find_pam_sites(sequence, pam="TTTN") # Cas12a
```

**Parameters:**
- `sequence`: DNA sequence string
- `pam`: PAM sequence pattern (default: `"NGG"`)
- `guide_length`: Length of guide sequence (default: `20`)

### Guide Scoring

Score individual guide sequences:

```python
from phaselab.crispr import score_guide, gc_content, santalucia_delta_g

guide = "GAAGGAGAGCAAGAGCGCGA"

# GC content (optimal: 40-70%)
gc = gc_content(guide)
# Returns: 0.60

# Thermodynamic stability (SantaLucia ΔG)
delta_g = santalucia_delta_g(guide)
# Returns: -45.2 (kcal/mol)

# Combined scoring
score = score_guide(
    guide,
    gc_weight=0.3,
    thermo_weight=0.3,
    position_weight=0.2,
    accessibility_weight=0.2
)
# Returns: 0.78
```

### Design Pipeline

The high-level pipeline for guide RNA design:

```python
from phaselab.crispr import design_guides

# Basic usage
sequence = """
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
AGGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
"""
tss_index = 50  # Transcription start site position

guides = design_guides(sequence, tss_index)
print(guides.head())
```

**Output DataFrame columns:**
- `sequence`: 20bp guide sequence
- `position`: Distance from TSS (negative = upstream)
- `strand`: `+` or `-`
- `gc_content`: GC fraction (0-1)
- `delta_g`: Thermodynamic stability
- `combined_score`: Weighted composite score
- `go_no_go`: Classification status

**Advanced usage with custom configuration:**

```python
from phaselab.crispr import design_guides, CRISPRConfig

# Custom configuration
config = CRISPRConfig(
    pam="NGG",
    guide_length=20,
    window_upstream=500,    # Search 500bp upstream of TSS
    window_downstream=100,  # Search 100bp downstream
    min_gc=0.40,
    max_gc=0.70,
    weights={
        'gc': 0.25,
        'thermo': 0.25,
        'position': 0.25,
        'accessibility': 0.25
    }
)

# With DNase accessibility peaks
dnase_peaks = [(45, 85), (120, 180)]  # Accessible regions

guides = design_guides(
    sequence,
    tss_index,
    config=config,
    dnase_peaks=dnase_peaks,
    verbose=True
)

# Filter top candidates
top_guides = guides[guides['go_no_go'] == 'GO'].head(10)
```

---

## Circadian Module

The circadian module models biological clock systems for Smith-Magenis Syndrome research.

### SMS Clock Model

Simulate the SMS circadian clock with RAI1 haploinsufficiency:

```python
from phaselab.circadian import simulate_sms_clock, SMSClockParams

# Basic simulation with 50% RAI1 (SMS condition)
results = simulate_sms_clock(rai1_level=0.5)

print(f"Period: {results['period']:.2f} hours")
print(f"Amplitude: {results['amplitude']:.4f}")
print(f"Phase shift: {results['phase_shift']:.2f} hours")
print(f"Coherence R̄: {results['coherence']:.4f}")
print(f"Status: {results['status']}")
```

**Output:**
```
Period: 23.45 hours
Amplitude: 0.6234
Phase shift: 2.15 hours
Coherence R̄: 0.7823
Status: GO
```

**Custom parameters:**

```python
from phaselab.circadian import simulate_sms_clock, SMSClockParams

# Custom clock parameters
params = SMSClockParams(
    tau_P=4.0,      # PER delay time constant (hours)
    alpha_P=2.0,    # PER suppression strength
    beta_R=0.5,     # RORα effect on BMAL1
    beta_V=0.5,     # REV-ERBα effect on BMAL1
    K_light=0.1,    # Light sensitivity
    omega_0=2*np.pi/24  # Base angular frequency
)

# Simulate with intervention (CRISPRa restoring RAI1)
results_baseline = simulate_sms_clock(rai1_level=0.5, params=params)
results_treated = simulate_sms_clock(rai1_level=0.85, params=params)

print(f"Baseline coherence: {results_baseline['coherence']:.4f}")
print(f"Treated coherence: {results_treated['coherence']:.4f}")
```

**Full time series output:**

```python
results = simulate_sms_clock(
    rai1_level=0.5,
    t_end=240.0,      # 10 days
    dt=0.1,           # Time step (hours)
    return_timeseries=True
)

# Access time series
t = results['time']
bmal1 = results['BMAL1']
per = results['PER']
rev_erb = results['REV_ERB']
ror = results['ROR']

# Plot
import matplotlib.pyplot as plt
plt.figure(figsize=(12, 6))
plt.plot(t, bmal1, label='BMAL1')
plt.plot(t, per, label='PER')
plt.xlabel('Time (hours)')
plt.ylabel('Expression level')
plt.legend()
plt.show()
```

### Kuramoto Oscillators

Model coupled oscillator networks:

```python
from phaselab.circadian import KuramotoNetwork

# Create network of 100 oscillators
network = KuramotoNetwork(
    n_oscillators=100,
    coupling_strength=0.5,
    natural_frequency_spread=0.1
)

# Simulate
results = network.simulate(t_end=100.0)

print(f"Order parameter R: {results['order_parameter']:.4f}")
print(f"Synchronization: {'Yes' if results['synchronized'] else 'No'}")
```

**Heterogeneous network:**

```python
from phaselab.circadian import KuramotoNetwork
import numpy as np

# Custom frequency distribution (bimodal)
frequencies = np.concatenate([
    np.random.normal(1.0, 0.1, 50),
    np.random.normal(1.1, 0.1, 50)
])

network = KuramotoNetwork(
    n_oscillators=100,
    coupling_strength=0.8,
    natural_frequencies=frequencies
)

# With external forcing (light-dark cycle)
def light_forcing(t):
    return 0.1 * np.sin(2 * np.pi * t / 24)

results = network.simulate(
    t_end=240.0,
    external_forcing=light_forcing
)
```

---

## Quantum Integration

PhaseLab integrates with IBM Quantum for hardware validation.

### Setting Up IBM Quantum

```python
import os
os.environ['IBM_QUANTUM_TOKEN'] = 'your-token-here'

# Or use .env file
# IBM_QUANTUM_TOKEN=your-token-here
```

### Running on Quantum Hardware

```python
from phaselab.quantum import QuantumCoherenceValidator

# Initialize validator
validator = QuantumCoherenceValidator(
    backend='ibm_torino',  # or 'ibm_brisbane', 'ibm_kyoto'
    shots=4096
)

# Validate a guide RNA
guide = "GAAGGAGAGCAAGAGCGCGA"
result = validator.validate(guide)

print(f"Hardware R̄: {result['coherence']:.4f}")
print(f"Status: {result['status']}")
print(f"Job ID: {result['job_id']}")
```

### Simulator Testing

```python
from phaselab.quantum import QuantumCoherenceValidator

# Use Aer simulator
validator = QuantumCoherenceValidator(backend='aer_simulator')

guides = [
    "GAAGGAGAGCAAGAGCGCGA",
    "AACTGCAAAGAAGTGGGCAC",
    "TACAGGAGCTTCCAGCGTCA"
]

for guide in guides:
    result = validator.validate(guide)
    print(f"{guide[:10]}... R̄={result['coherence']:.3f} [{result['status']}]")
```

---

## Complete Examples

### Example 1: SMS Gene Therapy Pipeline

```python
"""
Complete SMS gene therapy gRNA design pipeline.
"""
import phaselab as pl
from phaselab.crispr import design_guides, CRISPRConfig
from phaselab.circadian import simulate_sms_clock

# RAI1 promoter sequence (truncated for example)
rai1_promoter = """
GCGCGCTCGCGCGCTCGCGCGAAGGAGAGCAAGAGCGCGACGGCTAGCTAGCT
AGCTAGCTAGCTACAGGAGCTTCCAGCGTCAGGGCTAGCTAGCTAGCTAGCTA
GCTAGCTAGCTAACTGCAAAGAAGTGGGCACGCGCTAGCTAGCTAGCTAGCT
"""

# TSS position
tss_index = len(rai1_promoter) // 2

# Configure pipeline
config = CRISPRConfig(
    pam="NGG",
    window_upstream=400,
    window_downstream=50,
    min_gc=0.40,
    max_gc=0.70
)

# Design guides
guides = design_guides(rai1_promoter, tss_index, config=config)

# Filter GO candidates
go_candidates = guides[guides['go_no_go'] == 'GO']
print(f"Found {len(go_candidates)} GO candidates")

# Simulate therapeutic effect
for idx, row in go_candidates.head(3).iterrows():
    print(f"\n{row['sequence']}")
    print(f"  Position: {row['position']} bp from TSS")
    print(f"  Score: {row['combined_score']:.3f}")

    # Simulate with restored RAI1 (hypothetical 80% restoration)
    sim = simulate_sms_clock(rai1_level=0.8)
    print(f"  Predicted period: {sim['period']:.2f}h")
    print(f"  Predicted coherence: {sim['coherence']:.4f}")
```

### Example 2: Comparative Analysis

```python
"""
Compare coherence across multiple conditions.
"""
import phaselab as pl
from phaselab.circadian import simulate_sms_clock
import numpy as np

# RAI1 levels to test
rai1_levels = np.linspace(0.3, 1.0, 8)

results = []
for level in rai1_levels:
    sim = simulate_sms_clock(rai1_level=level)
    results.append({
        'rai1_level': level,
        'period': sim['period'],
        'coherence': sim['coherence'],
        'status': sim['status']
    })

# Find therapeutic threshold
import pandas as pd
df = pd.DataFrame(results)
go_threshold = df[df['status'] == 'GO']['rai1_level'].min()
print(f"Minimum RAI1 for GO status: {go_threshold:.1%}")

# Visualize
import matplotlib.pyplot as plt
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(df['rai1_level'], df['coherence'], 'bo-')
ax1.axhline(pl.core.constants.E_MINUS_2, color='r', linestyle='--', label='GO threshold')
ax1.set_xlabel('RAI1 Level')
ax1.set_ylabel('Coherence R̄')
ax1.legend()

ax2.plot(df['rai1_level'], df['period'], 'go-')
ax2.axhline(24, color='gray', linestyle='--', label='24h')
ax2.set_xlabel('RAI1 Level')
ax2.set_ylabel('Period (hours)')
ax2.legend()

plt.tight_layout()
plt.show()
```

### Example 3: Batch Processing

```python
"""
Batch process multiple sequences.
"""
from phaselab.crispr import design_guides
import pandas as pd

# Multiple target genes
targets = {
    'RAI1': 'ATGCGATCG...',
    'CLOCK': 'GCTAGCTAG...',
    'BMAL1': 'TAGCTAGCT...',
}

all_guides = []
for gene, sequence in targets.items():
    guides = design_guides(sequence, tss_index=len(sequence)//2)
    guides['target_gene'] = gene
    all_guides.append(guides)

combined = pd.concat(all_guides)
combined = combined.sort_values('combined_score', ascending=False)

# Top guide per gene
top_per_gene = combined.groupby('target_gene').first()
print(top_per_gene[['sequence', 'combined_score', 'go_no_go']])
```

---

## API Reference

### Core Functions

| Function | Description |
|----------|-------------|
| `coherence_score(data, mode='auto')` | Compute R̄ coherence score |
| `phase_variance(phases)` | Compute phase variance V_φ |
| `go_no_go(R_bar, threshold=E_MINUS_2)` | Binary GO/NO-GO classification |
| `classify_coherence(R_bar)` | Detailed category classification |

### CRISPR Functions

| Function | Description |
|----------|-------------|
| `design_guides(seq, tss, config, dnase_peaks)` | High-level guide design pipeline |
| `find_pam_sites(seq, pam, guide_length)` | Find PAM sites in sequence |
| `score_guide(guide, weights)` | Score individual guide |
| `gc_content(seq)` | Calculate GC fraction |
| `santalucia_delta_g(seq)` | Calculate thermodynamic ΔG |

### Circadian Functions

| Function | Description |
|----------|-------------|
| `simulate_sms_clock(rai1_level, params, t_end)` | Simulate SMS clock model |
| `KuramotoNetwork(n, coupling, frequencies)` | Create Kuramoto oscillator network |

### Constants

| Constant | Value | Description |
|----------|-------|-------------|
| `E_MINUS_2` | 0.1353... | e⁻² GO/NO-GO threshold |
| `FOUR_PI_SQUARED` | 39.478... | 4π² universal constant |

---

## Citation

If you use PhaseLab in your research, please cite:

```bibtex
@software{phaselab,
  author = {Vaca, Dylan},
  title = {PhaseLab: Phase-coherence analysis framework},
  year = {2024},
  url = {https://github.com/dylanvaca/phaselab}
}
```

---

## License

MIT License - see [LICENSE](LICENSE) for details.
