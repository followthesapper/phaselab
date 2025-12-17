# Two-Stage Perturbation Analysis Framework

## Overview

PhaseLab implements a two-stage perturbation analysis framework that enables spatial coherence analysis across multiple biological and chemical domains. This document describes the theoretical foundation, practical implementation, and cross-domain generalization of this approach.

## The Key Insight

> "The guide is the probe, not the structure."

Experiments E200-E211 showed that computing coherence from guide sequences (GC content, thermodynamic properties, structural features) has **no predictive value** for experimental outcomes (r ≈ 0).

Experiments E213-E216 validated that measuring **spatial coherence of response landscapes** predicts perturbation reliability:
- Correlation: r = -0.24 to -0.50 with outcome variance
- Variance reduction: 32-49% when selecting from stable regions
- Validation scale: 115,251 sgRNAs across 6 genes

## Two-Stage Framework

### Stage I: Landscape Feasibility (Pre-Tiling)

**Purpose**: Assess whether a perturbation landscape is likely to contain stable regions before expensive tiling experiments.

**Inputs**:
- Target sequence/structure
- Perturbation design parameters
- Prior knowledge (e.g., conserved domains, known functional regions)

**Outputs**:
- Feasibility assessment (GO/NO-GO)
- Predicted stable region candidates
- Recommended tiling experiment design

**What Stage I Does NOT Require**:
- Dense tiling screen data
- Experimental response measurements
- Large-scale screening results

**How It Works**:
1. **Structural Priors**: Uses evolutionary conservation, domain annotations, and sequence features to predict where coherent responses are likely
2. **Virtual Perturbation Lattice**: Models perturbation space geometry from design constraints
3. **Response Modeling**: Estimates expected response variance from prior screens or literature

### Stage II: Landscape Resolution (Tiling-Dependent)

**Purpose**: Empirically resolve the stability structure of a perturbation landscape.

**Inputs**:
- Actual tiling screen data
- Per-perturbation response measurements
- Replicate measurements (critical for variance estimation)

**Outputs**:
- Validated coherence profile
- Stable/Mixed/Amplifying region classification
- Quantitative variance reduction estimates

**Requirements**:
- Minimum 16-20 perturbations per gene/target
- At least 2-3 replicates per perturbation
- Replicate depth > perturbation density (prioritize quality over quantity)

## Stability Classification

PhaseLab classifies perturbation regions into four categories:

| Class | Coherence | Interpretation | Action |
|-------|-----------|----------------|--------|
| **STABLE** | > 0.7 | Low variance, reliable predictions | Select from this region |
| **MIXED** | 0.4-0.7 | Moderate variance, context-dependent | Use with caution |
| **AMPLIFYING** | < 0.4 | High variance, unreliable | Avoid or increase replicates |
| **IRRELEVANT** | - | Below response threshold | Ignore |

## Cross-Domain Generalization

The two-stage framework applies across multiple domains because perturbation-response relationships share common mathematical structure:

### Biology (CRISPR Screens)

| Component | Mapping |
|-----------|---------|
| Perturbation coordinate | Guide position (bp from TSS) |
| Response signal | Log2 fold-change, growth score |
| Coherence meaning | Neighboring guides yield consistent effects |
| Stable regions | Promoter elements with predictable modulation |

**Example**: CRISPRa RAI1 activation
- 8 guides identified in TSS-200 to TSS-50 window
- 1.36x fold change (68% expression target)
- Claim level: EXPLORATORY (pending tiling validation)

### Microbiology (TnSeq, CRISPRi)

| Component | Mapping |
|-----------|---------|
| Perturbation coordinate | Insertion position, guide position |
| Response signal | Fitness score, log2 fold-change |
| Coherence meaning | Essential domain structure |
| Stable regions | Core functional domains |

**Example**: CRISPRi dnaA essentiality
- Essential gene phenotype confirmed
- -2.19 log2FC (severe growth defect)
- Claim level: CONTEXT_DEPENDENT (known essential gene)

### Chemistry (Binding, SAR)

| Component | Mapping |
|-----------|---------|
| Perturbation coordinate | Residue position, modification site |
| Response signal | ΔΔG, IC50, Kd |
| Coherence meaning | Binding hot spot reliability |
| Stable regions | Consistent binding determinants |

**Example**: ABL1 binding landscape
- 11 positions analyzed
- 1.55 kcal/mol mean ΔΔG
- 2 hot spots identified in coherent regions

## Minimum Viable Tiling Experiments

One of the key practical insights: you don't need thousands of guides for Stage II validation.

### For CRISPRa/CRISPRi Targets

**Minimum Design**:
- 16-20 guides per gene
- Tiled at 20bp spacing across promoter
- 3 replicates per guide
- Total: ~50-60 data points per gene

**Why This Works**:
- Coherence is a local property (window = 5-10 guides)
- Edge effects minimal with proper windowing
- Replicate variance more important than guide density

### For Protein Binding

**Minimum Design**:
- 15-25 point mutations
- Covering predicted binding interface
- ΔΔG measurements with error bars
- Total: ~50-75 data points

### For Drug Dose-Response

**Minimum Design**:
- 8-12 concentration points per drug
- 3-4 replicates per concentration
- Include vehicle and positive controls
- Total: ~40-60 data points per drug

## API Usage

### Stage I: Feasibility Assessment

```python
from phaselab.landscapes import ResponseLandscape
from phaselab.spatial import analyze_tiling_coherence, TilingResult

# Create landscape from design (no experimental data yet)
landscape = ResponseLandscape(
    coords=np.array([50, 100, 150, 200, 250]),  # Planned positions
    responses=np.array([0.0, 0.0, 0.0, 0.0, 0.0]),  # Placeholder
)

# Stage I uses structural priors
# Returns feasibility assessment without requiring tiling data
```

### Stage II: Resolution with Data

```python
from phaselab.spatial import analyze_tiling_coherence, load_tiling_screen

# Load actual tiling screen data
landscape = load_tiling_screen(
    "tiling_screen_results.tsv",
    position_col="tss_distance",
    response_col="log2fc",
)

# Full coherence analysis
result = analyze_tiling_coherence(
    landscape,
    window=5,
    stable_threshold=0.7,
)

# Extract stable regions
for region in result.stable_regions:
    print(f"Stable: {region['start']}-{region['end']}")
    print(f"  Coherence: {region['coherence']:.3f}")
    print(f"  Mean response: {region['mean_response']:.3f}")
```

### Cross-Domain Examples

```python
# Biology: CRISPR tiling
from phaselab.spatial import analyze_tiling_coherence

# Microbiology: TnSeq essentiality
from phaselab.microbio import analyze_tnseq_coherence

# Chemistry: Binding landscapes
from phaselab.chem import analyze_binding_coherence

# Omics: Expression analysis
from phaselab.omics import analyze_expression_coherence
```

## Quantum Mode Integration

The two-stage framework integrates with PhaseLab's quantum computation modes:

```python
from phaselab.quantum import set_quantum_mode, QuantumMode

# Stage I: Classical only (fast feasibility)
set_quantum_mode(QuantumMode.OFF)

# Stage II: Quantum validation for production
set_quantum_mode(QuantumMode.AUDIT)

# Research-grade: Full quantum coherence
set_quantum_mode(QuantumMode.REQUIRED)
```

## Claim Level Propagation

Results carry claim levels that propagate through the pipeline:

| Level | Meaning | Requirements |
|-------|---------|--------------|
| **STRONG_COMPUTATIONAL** | High confidence | Validated coherence, multiple evidence sources |
| **CONTEXT_DEPENDENT** | Moderate confidence | Single evidence source, known target |
| **EXPLORATORY** | Low confidence | Structural priors only, no tiling data |
| **UNKNOWN** | Insufficient data | Cannot assess reliability |

## Experimental Validation Protocol

To move from EXPLORATORY to STRONG_COMPUTATIONAL:

1. **Design Minimum Viable Tiling**
   - 16-20 guides per gene
   - Cover predicted stable regions + flanking
   - Include positive/negative controls

2. **Execute Screen**
   - 3 replicates minimum
   - Consistent experimental conditions
   - Track batch effects

3. **Validate Coherence**
   - Compute spatial coherence profile
   - Compare to Stage I predictions
   - Assess variance reduction

4. **Update Claim Level**
   - If coherence validated: STRONG_COMPUTATIONAL
   - If partial: CONTEXT_DEPENDENT
   - If refuted: Return to design

## References

- E213-E216: Spatial coherence validation across 115,251 sgRNAs
- Chang et al. 2022: CRISPRa RAI1 activation in SMS models
- CRISPR-SURF: Deconvolution methods for tiling screens

---

*Framework version: 1.0.0*
*Last updated: December 2025*
