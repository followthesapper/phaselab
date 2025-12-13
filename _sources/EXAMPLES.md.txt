# PhaseLab Examples

Practical examples for common use cases.

## Basic Coherence Analysis

### Calculate Coherence from Experimental Data

```python
import phaselab as pl
import numpy as np

# Your experimental phase measurements
phase_data = np.array([0.12, 0.15, 0.11, 0.14, 0.13, 0.16, 0.12])

# Calculate coherence
R_bar = pl.coherence_score(phase_data)
V_phi = pl.phase_variance(phase_data)

print(f"Phase variance: {V_phi:.6f}")
print(f"Coherence R̄: {R_bar:.4f}")
print(f"Status: {pl.go_no_go(R_bar)}")
```

### From Quantum Expectation Values

```python
import phaselab as pl

# Expectation values from VQE or quantum measurement
expectations = [0.82, 0.85, 0.79, 0.88, 0.81]

R_bar = pl.coherence_score(expectations, mode='expectations')
print(f"Quantum coherence: {R_bar:.4f}")
```

---

## CRISPR Guide Design

### Simple Guide Design

```python
from phaselab.crispr import design_guides

# Your target promoter sequence
promoter = """
ATGCGATCGATCGATCGATCGAGGCGATCGATCGATCGATCGATCGATCG
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
GGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
"""

# TSS is at position 75
guides = design_guides(promoter, tss_index=75)

# Show top 5 guides
print(guides[['sequence', 'position', 'gc_content', 'combined_score', 'go_no_go']].head())
```

### Filter by Criteria

```python
from phaselab.crispr import design_guides

guides = design_guides(promoter, tss_index=75)

# Only GO candidates
go_guides = guides[guides['go_no_go'] == 'GO']

# GC content between 50-60%
optimal_gc = go_guides[(go_guides['gc_content'] >= 0.50) &
                        (go_guides['gc_content'] <= 0.60)]

# Within 200bp of TSS
proximal = optimal_gc[optimal_gc['position'].abs() <= 200]

print(f"Filtered from {len(guides)} to {len(proximal)} guides")
print(proximal.head())
```

### With Chromatin Accessibility

```python
from phaselab.crispr import design_guides

# DNase-seq peaks (accessible regions)
dnase_peaks = [
    (50, 120),   # Accessible region 1
    (180, 250),  # Accessible region 2
]

guides = design_guides(
    promoter,
    tss_index=75,
    dnase_peaks=dnase_peaks,
    verbose=True
)

# Guides in accessible regions get higher scores
print(guides.head(10))
```

---

## Circadian Clock Simulation

### Basic SMS Simulation

```python
from phaselab.circadian import simulate_sms_clock

# Normal individual (100% RAI1)
normal = simulate_sms_clock(rai1_level=1.0)

# SMS patient (50% RAI1)
sms = simulate_sms_clock(rai1_level=0.5)

print("Normal:")
print(f"  Period: {normal['period']:.2f}h")
print(f"  Coherence: {normal['coherence']:.4f}")
print(f"  Status: {normal['status']}")

print("\nSMS:")
print(f"  Period: {sms['period']:.2f}h")
print(f"  Coherence: {sms['coherence']:.4f}")
print(f"  Status: {sms['status']}")
```

### Therapeutic Dose-Response

```python
from phaselab.circadian import simulate_sms_clock
import numpy as np

# Test different RAI1 restoration levels
levels = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

print("RAI1 Level | Period | Coherence | Status")
print("-" * 45)
for level in levels:
    result = simulate_sms_clock(rai1_level=level)
    print(f"   {level:.0%}     | {result['period']:5.2f}h |  {result['coherence']:.4f}  | {result['status']}")
```

### Time Series Analysis

```python
from phaselab.circadian import simulate_sms_clock
import matplotlib.pyplot as plt

# Get full time series
result = simulate_sms_clock(
    rai1_level=0.5,
    t_end=168,  # 7 days
    return_timeseries=True
)

# Plot
fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

axes[0].plot(result['time'], result['BMAL1'], label='BMAL1', color='blue')
axes[0].plot(result['time'], result['PER'], label='PER', color='red')
axes[0].set_ylabel('Expression')
axes[0].legend()
axes[0].set_title('SMS Clock Genes (50% RAI1)')

axes[1].plot(result['time'], result['REV_ERB'], label='REV-ERBα', color='purple')
axes[1].plot(result['time'], result['ROR'], label='RORα', color='green')
axes[1].set_xlabel('Time (hours)')
axes[1].set_ylabel('Expression')
axes[1].legend()

plt.tight_layout()
plt.savefig('sms_clock_timeseries.png', dpi=150)
plt.show()
```

---

## Kuramoto Oscillator Networks

### Basic Synchronization

```python
from phaselab.circadian import KuramotoNetwork

# 50 coupled oscillators
network = KuramotoNetwork(
    n_oscillators=50,
    coupling_strength=0.5
)

result = network.simulate(t_end=100)

print(f"Order parameter: {result['order_parameter']:.4f}")
print(f"Synchronized: {result['synchronized']}")
```

### Critical Coupling Analysis

```python
from phaselab.circadian import KuramotoNetwork
import numpy as np

# Find critical coupling strength
couplings = np.linspace(0.1, 1.0, 10)
order_params = []

for K in couplings:
    network = KuramotoNetwork(n_oscillators=100, coupling_strength=K)
    result = network.simulate(t_end=50)
    order_params.append(result['order_parameter'])

# Find transition point
import matplotlib.pyplot as plt
plt.plot(couplings, order_params, 'bo-')
plt.xlabel('Coupling strength K')
plt.ylabel('Order parameter R')
plt.axhline(0.5, color='r', linestyle='--', alpha=0.5)
plt.title('Kuramoto Synchronization Transition')
plt.savefig('kuramoto_transition.png', dpi=150)
plt.show()
```

---

## IBM Quantum Hardware

### Simulator Test

```python
from phaselab.quantum import QuantumCoherenceValidator

# Use local simulator
validator = QuantumCoherenceValidator(backend='aer_simulator')

guide = "GAAGGAGAGCAAGAGCGCGA"
result = validator.validate(guide)

print(f"Guide: {guide}")
print(f"Coherence: {result['coherence']:.4f}")
print(f"Status: {result['status']}")
```

### Real Hardware Validation

```python
import os
os.environ['IBM_QUANTUM_TOKEN'] = 'your-token'

from phaselab.quantum import QuantumCoherenceValidator

# Connect to IBM Quantum
validator = QuantumCoherenceValidator(
    backend='ibm_torino',
    shots=4096
)

# Validate candidates
candidates = [
    "GAAGGAGAGCAAGAGCGCGA",
    "AACTGCAAAGAAGTGGGCAC",
    "TACAGGAGCTTCCAGCGTCA"
]

for guide in candidates:
    result = validator.validate(guide)
    print(f"{guide}: R̄={result['coherence']:.3f} [{result['status']}]")
```

---

## Complete Pipeline Example

### End-to-End Gene Therapy Design

```python
"""
Complete pipeline: sequence → guides → simulation → validation
"""
import phaselab as pl
from phaselab.crispr import design_guides
from phaselab.circadian import simulate_sms_clock

# 1. Design guides
print("=" * 60)
print("STEP 1: Guide RNA Design")
print("=" * 60)

promoter = """
GCGCGCTCGCGCGCTCGCGCGAAGGAGAGCAAGAGCGCGACGGCTAGCTAGCT
AGCTAGCTAGCTACAGGAGCTTCCAGCGTCAGGGCTAGCTAGCTAGCTAGCTA
GCTAGCTAGCTAACTGCAAAGAAGTGGGCACGCGCTAGCTAGCTAGCTAGCT
"""

guides = design_guides(promoter, tss_index=len(promoter)//2)
go_candidates = guides[guides['go_no_go'] == 'GO'].head(3)

print(f"Found {len(go_candidates)} top GO candidates:\n")
for i, (_, row) in enumerate(go_candidates.iterrows(), 1):
    print(f"  {i}. {row['sequence']}")
    print(f"     Position: {row['position']} bp | Score: {row['combined_score']:.3f}")

# 2. Simulate therapeutic effect
print("\n" + "=" * 60)
print("STEP 2: Therapeutic Simulation")
print("=" * 60)

# Baseline (untreated SMS)
baseline = simulate_sms_clock(rai1_level=0.5)
print(f"\nBaseline (50% RAI1):")
print(f"  Period: {baseline['period']:.2f}h")
print(f"  Coherence: {baseline['coherence']:.4f}")
print(f"  Status: {baseline['status']}")

# With CRISPRa treatment (estimated 80% restoration)
treated = simulate_sms_clock(rai1_level=0.8)
print(f"\nTreated (80% RAI1 restoration):")
print(f"  Period: {treated['period']:.2f}h")
print(f"  Coherence: {treated['coherence']:.4f}")
print(f"  Status: {treated['status']}")

# 3. Summary
print("\n" + "=" * 60)
print("STEP 3: Summary")
print("=" * 60)

improvement = treated['coherence'] - baseline['coherence']
print(f"\nCoherence improvement: +{improvement:.4f}")
print(f"Period correction: {abs(treated['period'] - 24):.2f}h → {abs(baseline['period'] - 24):.2f}h from 24h")

if treated['status'] == 'GO' and baseline['status'] != 'GO':
    print("\n✓ Treatment restores GO status!")
elif treated['status'] == 'GO':
    print("\n✓ Treatment maintains/improves GO status")
```

---

## Tips and Best Practices

### 1. Always Check GO Status

```python
# Before proceeding with any candidate
if pl.go_no_go(R_bar) == "NO-GO":
    print("Warning: Candidate below coherence threshold")
```

### 2. Use Batch Processing for Large Screens

```python
from concurrent.futures import ProcessPoolExecutor
from phaselab.crispr import design_guides

def process_sequence(seq_data):
    name, sequence = seq_data
    guides = design_guides(sequence, tss_index=len(sequence)//2)
    return name, guides[guides['go_no_go'] == 'GO']

sequences = [('gene1', seq1), ('gene2', seq2), ...]

with ProcessPoolExecutor(max_workers=4) as executor:
    results = list(executor.map(process_sequence, sequences))
```

### 3. Save Results for Reproducibility

```python
import json
import pandas as pd

# Save guide results
guides.to_csv('guide_candidates.csv', index=False)

# Save simulation results
results = simulate_sms_clock(rai1_level=0.5)
with open('simulation_results.json', 'w') as f:
    json.dump({k: v for k, v in results.items()
               if not isinstance(v, np.ndarray)}, f, indent=2)
```

### 4. Validate on Simulator Before Hardware

```python
from phaselab.quantum import QuantumCoherenceValidator

# Always test on simulator first
sim_validator = QuantumCoherenceValidator(backend='aer_simulator')
sim_result = sim_validator.validate(guide)

if sim_result['status'] == 'GO':
    # Then proceed to hardware
    hw_validator = QuantumCoherenceValidator(backend='ibm_torino')
    hw_result = hw_validator.validate(guide)
```
