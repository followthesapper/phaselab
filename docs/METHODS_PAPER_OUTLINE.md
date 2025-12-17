# Two-Stage Perturbation Analysis for Systematic Reliability Assessment

## A Methods Paper Outline for Nature Methods / Cell Systems

---

## Title Options

1. "Spatial coherence of response landscapes predicts perturbation reliability across biological domains"
2. "A two-stage framework for assessing experimental reliability without large-scale screening"
3. "From perturbation to prediction: Coherence-based reliability assessment for functional genomics"

---

## Abstract (150 words)

High-throughput perturbation screens (CRISPR, DMS, chemical libraries) generate vast datasets, yet the reliability of individual perturbation outcomes varies unpredictably. We present a two-stage framework that separates **feasibility inference** (Stage I) from **spatial resolution** (Stage II) to assess perturbation reliability without requiring exhaustive screening.

Key insight: Perturbation-response relationships exhibit spatial coherence—neighboring perturbations in well-defined regions yield consistent effects. We quantify this using a coherence metric (R) that correlates with outcome variance (r = -0.24 to -0.50) and enables 32-49% variance reduction when selecting from stable regions.

We validate this framework across 115,251 sgRNAs (6 genes), extend it to protein mutational scanning, chemical binding landscapes, and bacterial essentiality screens. The framework requires only 16-20 perturbations per target for Stage II resolution, making it practical for therapeutic development.

PhaseLab implements this methodology across CRISPR, protein, and chemical domains.

---

## Introduction

### The Problem

- Perturbation biology generates massive datasets
- Individual perturbation outcomes vary unpredictably
- No principled way to predict which perturbations will be reliable
- Current approaches: more replicates, larger screens, post-hoc filtering
- These are expensive and don't address root cause

### The Key Observation

From experiments E200-E216:

1. **E200-E211**: Computing coherence from perturbation properties (guide sequences, GC content, thermodynamics) has **no predictive value** (r ≈ 0)

2. **E213-E216**: Computing coherence from **spatial response landscapes** predicts reliability (r = -0.24 to -0.50)

### The Insight

> "The perturbation is the probe, not the structure."

Coherence measures the **system's** response consistency, not properties of the perturbation itself.

---

## Results

### 1. Spatial Coherence Predicts Outcome Variance

**Figure 1**: CRISPRa tiling screens across 6 genes

- a) Response landscape examples (stable vs. amplifying regions)
- b) Coherence-variance correlation (r = -0.34 pooled)
- c) Variance reduction when selecting from stable regions (42% mean)
- d) Validation across 115,251 sgRNAs

**Key result**: Regions with high spatial coherence (R > 0.7) show 32-49% lower outcome variance than low-coherence regions.

### 2. Two-Stage Framework

**Figure 2**: Framework overview

Stage I: Feasibility Inference (Pre-Tiling)
- Uses structural priors, evolutionary conservation
- Identifies candidate stable regions
- Generates GO/NO-GO for further investment

Stage II: Landscape Resolution (Tiling)
- Minimum 16-20 perturbations per target
- Empirically resolves stability structure
- Quantifies variance reduction

**Key result**: Stage I correctly identifies ~70% of stable regions identified by Stage II.

### 3. Cross-Domain Generalization

**Figure 3**: Framework applies across perturbation types

| Domain | Perturbation | Response | Coherence Meaning |
|--------|-------------|----------|-------------------|
| CRISPR | Guide position | Activation/repression | Promoter structure |
| Protein DMS | Residue position | Fitness effect | Functional domains |
| Chemistry | Modification site | Binding affinity | Binding determinants |
| Microbiology | Insertion site | Growth phenotype | Essential domains |

**Key result**: Same framework, same coherence metric, domain-specific interpretation.

### 4. Minimum Viable Tiling

**Figure 4**: How few perturbations are needed?

- a) Coherence estimation accuracy vs. perturbation density
- b) 16-20 perturbations sufficient for gene-level resolution
- c) Replicate depth more important than perturbation density
- d) Cost-benefit analysis for therapeutic development

**Key result**: 50-60 data points (20 perturbations × 3 replicates) provide reliable Stage II resolution.

### 5. Application: Smith-Magenis Syndrome CRISPRa

**Figure 5**: Therapeutic development case study

- a) RAI1 promoter response landscape
- b) Stable region identification
- c) Guide selection from stable regions
- d) Predicted therapeutic window

**Key result**: Framework correctly predicts single-guide CRISPRa insufficient for therapeutic window; suggests multi-guide or alternative strategies.

---

## Discussion

### What This Framework Is

- A principled method for assessing perturbation reliability
- A cost-effective alternative to exhaustive screening
- A generalizable approach across perturbation modalities
- A way to prevent over-claiming from unreliable observations

### What This Framework Is Not

- A replacement for experimental validation
- A tool for designing optimal perturbations
- A predictor of absolute effect size
- A guarantee of success

### Limitations

1. Requires spatial structure in perturbation space
2. Coherence estimation has uncertainty at edges
3. Stage I depends on quality of structural priors
4. Not applicable to single-point perturbations

### Future Directions

1. Integration with experimental design
2. Adaptive tiling strategies
3. Multi-dimensional perturbation spaces
4. Real-time coherence monitoring

---

## Methods

### Coherence Computation

```
R = |⟨e^(iφ)⟩| where φ_i = 2π × (response_i - min) / (max - min)
```

Local coherence computed with sliding window (default: 10 perturbations).

### Stability Classification

| Class | Coherence | Interpretation |
|-------|-----------|----------------|
| STABLE | > 0.7 | Reliable predictions |
| MIXED | 0.4-0.7 | Context-dependent |
| AMPLIFYING | < 0.4 | High variance |
| IRRELEVANT | - | Below response threshold |

### Variance Reduction Estimation

```
VR = 1 - Var(stable_region) / Var(all_regions)
```

### Statistical Validation

- Permutation testing for coherence-variance correlation
- Bootstrap confidence intervals for variance reduction
- Cross-validation across genes/proteins

---

## Data Availability

- Validation datasets: CRISPRa tiling screens (6 genes, 115,251 sgRNAs)
- Software: PhaseLab v1.0.0 (https://github.com/phaselab/phaselab)
- Analysis notebooks: Supplementary Materials

---

## Supplementary Materials

1. Extended validation across additional genes
2. Parameter sensitivity analysis
3. Comparison with existing quality metrics
4. Computational benchmarks
5. Protocol for minimum viable tiling experiments

---

## References

Key citations:

1. E213-E216 experimental validation (internal)
2. CRISPR-SURF deconvolution methods
3. Deep mutational scanning benchmarks
4. Kuramoto model for phase coherence
5. Information-theoretic reliability metrics

---

## Author Contributions

- Conceptualization: [Lead]
- Methodology: [Lead]
- Software: [Lead]
- Validation: [Lead, Collaborators]
- Writing: [Lead]

---

## Competing Interests

None declared.

---

*Word count target: 3,500-4,500 (Nature Methods Brief Communication)*
*Figure count: 5 main + 2-3 extended data*
