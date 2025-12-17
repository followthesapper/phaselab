# Quantum Mode Guidance

## Overview

PhaseLab supports three quantum execution modes. This document specifies when each mode is appropriate.

---

## The Three Modes

| Mode | Meaning | Performance | Use Case |
|------|---------|-------------|----------|
| **OFF** | Classical only | Fastest | Standard analysis |
| **AUDIT** | Classical + quantum validation | Moderate | Research validation |
| **REQUIRED** | Quantum mandatory | Slowest | Fundamental research |

---

## The Decision Rule

> **If a classical experiment can falsify a claim, quantum mode should be OFF.**

This rule applies to 95%+ of PhaseLab use cases.

---

## Module-by-Module Guidance

### Biology Modules

| Module | Default | Recommended | Justification |
|--------|---------|-------------|---------------|
| `phaselab.spatial` | OFF | OFF | Classical coherence sufficient |
| `phaselab.surf` | OFF | OFF | Deconvolution is classical |
| `phaselab.crispr.pipeline` | OFF | OFF | Guide scoring is classical |
| `phaselab.crispr.scoring` | OFF | OFF | All metrics classical |

### Protein Modules

| Module | Default | Recommended | Justification |
|--------|---------|-------------|---------------|
| `phaselab.protein.mutscan` | OFF | OFF | DMS coherence is classical |
| `phaselab.protein.folding` | OFF | OFF | Ramachandran analysis is classical |
| `phaselab.protein.binding` | OFF | OFF | Affinity patterns are classical |

### Chemistry Modules

| Module | Default | Recommended | Justification |
|--------|---------|-------------|---------------|
| `phaselab.chem.binding` | OFF | OFF | Binding landscapes classical |
| `phaselab.chem.reaction` | OFF | OFF | Reaction optimization classical |
| `phaselab.chem.screening` | OFF | OFF | HTS analysis classical |

### Omics Modules

| Module | Default | Recommended | Justification |
|--------|---------|-------------|---------------|
| `phaselab.omics.expression` | OFF | OFF | RNA-seq analysis classical |
| `phaselab.omics.atac` | OFF | OFF | Accessibility classical |
| `phaselab.omics.chip` | OFF | OFF | Binding sites classical |

### Microbiology Modules

| Module | Default | Recommended | Justification |
|--------|---------|-------------|---------------|
| `phaselab.microbio.tnseq` | OFF | OFF | Essentiality classical |
| `phaselab.microbio.crispri` | OFF | OFF | Bacterial screens classical |
| `phaselab.microbio.drug` | OFF | OFF | Dose-response classical |

### Trial Modules

| Module | Default | Recommended | Justification |
|--------|---------|-------------|---------------|
| `phaselab.trials.sms.*` | OFF | OFF | All therapeutic trials classical |

### Circadian Modules

| Module | Default | Recommended | Justification |
|--------|---------|-------------|---------------|
| `phaselab.circadian.sms_model` | OFF | OFF | Kuramoto simulation classical |

---

## When to Use AUDIT Mode

AUDIT mode runs classical computation first, then validates a subset with quantum.

### Appropriate Use Cases

| Use Case | Justification |
|----------|---------------|
| Publishing coherence results | Independent validation |
| Cross-checking unusual findings | Catch classical edge cases |
| Methodology validation papers | Show classical ≈ quantum |
| Research-grade analysis | Extra rigor when needed |

### Configuration

```python
from phaselab.quantum import set_quantum_mode, QuantumMode

set_quantum_mode(QuantumMode.AUDIT)
# Now 10% of computations will be quantum-validated
```

---

## When to Use REQUIRED Mode

REQUIRED mode uses quantum for all coherence computation.

### Appropriate Use Cases

| Use Case | Justification |
|----------|---------------|
| Representation efficiency bounds | Quantum fundamentally different |
| Measurement limit experiments | Probing quantum-classical boundary |
| IR theoretical validation | Testing core theory |
| Fundamental physics research | Not applied biology |

### Configuration

```python
from phaselab.quantum import set_quantum_mode, QuantumMode

set_quantum_mode(QuantumMode.REQUIRED)
# All coherence via quantum simulation
# Raises error if quantum backend unavailable
```

---

## What Quantum Adds (and Doesn't)

### What Quantum Can Add

| Capability | When Relevant |
|------------|---------------|
| Representation efficiency bounds | Theoretical research |
| Phase transition classification | System characterization |
| Decoherence universality tests | Fundamental IR validation |
| Measurement limit probing | Physics experiments |

### What Quantum Does NOT Add

| Capability | Why Not |
|------------|---------|
| Better CRISPR guide selection | Classical coherence sufficient |
| Better DMS analysis | No quantum advantage |
| Better binding predictions | Classical physics governs |
| Better therapeutic decisions | Reliability is the question |

---

## Cost-Benefit Analysis

| Mode | Compute Cost | Value Added | When Worth It |
|------|--------------|-------------|---------------|
| OFF | 1x | Baseline | Always for applied work |
| AUDIT | 1.1x | Validation | Publications, unusual results |
| REQUIRED | 10-100x | Theoretical rigor | Fundamental research only |

---

## Common Mistakes

### Mistake 1: "Quantum is always better"

**Wrong**. Quantum is slower and adds no value when classical suffices.

### Mistake 2: "Use REQUIRED for therapeutic work"

**Wrong**. Therapeutic decisions need wet lab, not quantum.

### Mistake 3: "AUDIT mode for every analysis"

**Unnecessary**. Only use AUDIT when publishing or validating.

### Mistake 4: "Quantum makes results more accurate"

**Misleading**. Classical and quantum give same coherence values. Quantum validates the classical is correct.

---

## The Bottom Line

For 95%+ of PhaseLab use cases:

```python
from phaselab.quantum import set_quantum_mode

set_quantum_mode("off")  # This is the default. Leave it.
```

Quantum is for:
1. Methodology papers
2. Theoretical research
3. IR validation

Quantum is NOT for:
1. Guide selection
2. Protein analysis
3. Drug discovery
4. Any applied biology

---

*PhaseLab v1.0.0*
*December 2025*
