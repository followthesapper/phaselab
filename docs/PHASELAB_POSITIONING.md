# Why PhaseLab Is Not a CRISPR Tool

## A Positioning Document

---

## The Misconception

When people first encounter PhaseLab, they often ask:

> "So it's like CRISPOR but better?"

This is a fundamental misunderstanding of what PhaseLab does and why it matters.

---

## What CRISPR Tools Do

Tools like CRISPOR, Benchling, CRISPRscan, and Cas-OFFinder are **perturbation design tools**. They:

1. Find PAM sites in a sequence
2. Score guide RNA properties (GC content, secondary structure)
3. Predict off-target binding sites
4. Rank guides by efficiency scores

These tools answer: **"What guide should I use?"**

They are essential, well-validated, and PhaseLab does not replace them.

---

## What PhaseLab Does

PhaseLab is a **perturbation reliability framework**. It:

1. Assesses whether ANY perturbation in a region will be reliable
2. Quantifies outcome variance before experiments are run
3. Identifies stable vs. amplifying regions
4. Propagates uncertainty through downstream analyses

PhaseLab answers: **"Can I trust the results I'll get?"**

---

## The Key Difference

### CRISPOR asks:
"Given this target, which guide is best?"

### PhaseLab asks:
"Given this target, is ANY guide going to give reliable results?"

These are orthogonal questions.

---

## Why This Matters

### Scenario 1: High-Scoring CRISPOR Guide Fails

You design a CRISPR knockout with a guide that scores:
- MIT: 98
- CFD: 97
- Doench: 0.72

You run the experiment. It fails. Why?

**CRISPOR can't tell you.** The guide was "good" by all metrics.

**PhaseLab can.** If the guide is in an amplifying region (low spatial coherence), outcome variance is high regardless of guide properties. The guide wasn't "bad"—the region is unreliable.

### Scenario 2: Medium-Scoring Guide Succeeds

A collaborator uses a "mediocre" guide:
- MIT: 65
- CFD: 70
- Doench: 0.45

It works beautifully, consistently across replicates.

**CRISPOR is surprised.** The guide "shouldn't" work well.

**PhaseLab isn't.** The guide is in a stable region. Spatial coherence is high. Individual guide properties matter less when the system responds consistently.

---

## The Experimental Evidence

### E200-E211: Guide-Sequence Coherence Fails

We computed coherence from guide properties:
- GC content
- Thermodynamic stability
- Secondary structure
- Sequence motifs

Result: **r ≈ 0** with experimental outcomes.

Guide properties do not predict reliability.

### E213-E216: Spatial Coherence Works

We computed coherence from response landscapes:
- Position along the target
- Effect at neighboring positions
- Local variance patterns

Result: **r = -0.24 to -0.50** with outcome variance.

Spatial structure predicts reliability.

---

## Domain Generalization

PhaseLab applies far beyond CRISPR:

| Domain | Perturbation | Tool Equivalent | PhaseLab Adds |
|--------|-------------|-----------------|---------------|
| CRISPR | Guide RNA | CRISPOR, CRISPick | Reliability assessment |
| Protein | Mutations | PolyPhen, EVE | Coherent domain identification |
| Chemistry | Modifications | Docking scores | Stable binding mode detection |
| Drug Screening | Compounds | IC50 ranking | Reliable hit identification |

CRISPOR is to PhaseLab as:
- PolyPhen is to PhaseLab (for proteins)
- Docking scores are to PhaseLab (for binding)
- Growth curves are to PhaseLab (for microbiology)

PhaseLab is the reliability layer that goes on top.

---

## Use Case Comparison

### When to Use CRISPOR

1. You need to design a guide for a specific target
2. You want to minimize off-target effects
3. You need standard efficiency predictions
4. You're doing routine CRISPR work

### When to Use PhaseLab

1. You need to know if your region is reliable
2. You're selecting from multiple candidate targets
3. You're designing a tiling experiment
4. You want to assess claim levels for therapeutic development
5. You're working outside CRISPR (protein, chemistry, etc.)

### When to Use Both

1. **First**: Run CRISPOR to get candidate guides
2. **Then**: Run PhaseLab to assess which regions are stable
3. **Select**: Guides from CRISPOR that fall in stable PhaseLab regions
4. **Claim**: Results with appropriate uncertainty levels

---

## The Two-Stage Framework

PhaseLab introduces a concept that has no equivalent in CRISPOR:

**Stage I: Feasibility** (Before tiling)
- Can we predict stable regions from structure alone?
- Should we invest in a full tiling experiment?
- GO/NO-GO decision

**Stage II: Resolution** (With tiling)
- Empirically resolve the stability landscape
- Identify stable vs. amplifying regions
- Quantify variance reduction

CRISPOR operates in neither stage—it's pre-Stage I (guide design).

---

## What PhaseLab Does NOT Do

1. **Not** a replacement for CRISPOR/CRISPick
2. **Not** an off-target predictor
3. **Not** an efficiency scorer
4. **Not** a guide ranking tool
5. **Not** specific to CRISPR

---

## What PhaseLab DOES Do

1. Assesses perturbation reliability across modalities
2. Quantifies outcome variance before experiments
3. Identifies stable/amplifying regions
4. Propagates uncertainty through analyses
5. Generates claim levels for honest reporting
6. Designs minimum viable validation experiments

---

## The Honest Position

If someone asks "Is PhaseLab better than CRISPOR?", the answer is:

> "They do different things. CRISPOR designs guides. PhaseLab assesses whether your results will be reliable. Use both."

If someone asks "Why do I need PhaseLab if I have CRISPOR?", the answer is:

> "CRISPOR tells you which guide is best. PhaseLab tells you whether 'best' will matter. A perfect guide in an amplifying region still gives unreliable results."

If someone asks "Does PhaseLab work outside CRISPR?", the answer is:

> "Yes. That's the point. PhaseLab is a perturbation reliability framework. CRISPR is one application. Protein engineering, drug discovery, and microbial screens are others."

---

## Summary Table

| Aspect | CRISPOR | PhaseLab |
|--------|---------|----------|
| Question answered | "Which guide?" | "Can I trust it?" |
| Input | Target sequence | Response landscape |
| Output | Ranked guides | Reliability assessment |
| Domain | CRISPR only | Any perturbation |
| Stage | Pre-experiment | Experiment design |
| Replaces | Manual PAM search | Nothing (new capability) |

---

## Conclusion

PhaseLab is not competing with CRISPOR. It's solving a different problem that CRISPOR cannot address: **perturbation reliability**.

The question "Is PhaseLab a CRISPR tool?" is like asking "Is a thermometer a cooking tool?"

A thermometer helps you cook better, but it doesn't replace your stove. PhaseLab helps you do perturbation biology better, but it doesn't replace your guide design tools.

Use both. They're complementary.

---

*PhaseLab v1.0.0*
*December 2025*
