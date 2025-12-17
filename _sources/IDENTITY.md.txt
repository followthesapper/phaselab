# PhaseLab Identity

## What PhaseLab Is

**PhaseLab is a reliability layer for perturbation science.**

That's it. Everything else is implementation.

---

## The One-Sentence Definition

> PhaseLab tells you whether your experimental results will be reproducible before you run the experiment.

---

## What PhaseLab Is NOT

| Not This | Why Not |
|----------|---------|
| A CRISPR tool | CRISPR is one application |
| A biology platform | Biology is one domain |
| A quantum project | Quantum is one mode |
| A theory engine | Theory is the foundation, not the product |
| A prediction system | It predicts reliability, not outcomes |

---

## The Core Abstraction

All perturbation experiments share this structure:

```
Perturbation → System → Response
```

PhaseLab asks: **Is the response reliable?**

| Domain | Perturbation | Response | Reliability Question |
|--------|-------------|----------|---------------------|
| CRISPR | Guide position | Activation/repression | Will this guide work consistently? |
| Protein | Mutation | Fitness effect | Is this residue functionally important? |
| Chemistry | Modification | Binding affinity | Is this interaction real? |
| Microbiology | Insertion | Growth phenotype | Is this gene essential? |

Same question. Same framework. Same math.

---

## The Value Proposition

### For Biologists

"Stop wasting time on experiments that won't replicate."

### For Chemists

"Find stable binding modes before synthesis."

### For Therapeutic Developers

"Know your claim levels before making promises."

### For Reviewers

"Reproducibility is quantified, not assumed."

---

## The Brand Position

### PhaseLab vs. Design Tools

| Tool | Question Answered |
|------|-------------------|
| CRISPOR | "Which guide is best?" |
| CRISPick | "Which guide is most efficient?" |
| **PhaseLab** | "Can I trust the result?" |

### PhaseLab vs. Analysis Tools

| Tool | Purpose |
|------|---------|
| DESeq2 | Differential expression statistics |
| MACS2 | Peak calling |
| **PhaseLab** | Reliability assessment |

### PhaseLab vs. Prediction Tools

| Tool | Predicts |
|------|----------|
| AlphaFold | Structure |
| EVE | Variant effects |
| **PhaseLab** | Prediction reliability |

---

## The Scientific Foundation

PhaseLab is built on one empirical observation:

> Spatial coherence of response landscapes correlates with outcome variance.

This observation was validated across 115,251 sgRNAs and extends to proteins, chemistry, and microbiology.

The theoretical foundation (Informational Relativity) explains *why* this works, but the empirical observation stands independently.

---

## What Makes PhaseLab Different

### 1. Domain-General

Most tools are domain-specific. PhaseLab works across:
- CRISPR screens
- Protein mutagenesis
- Chemical modification
- Microbial genetics

### 2. Pre-Experimental

Most tools analyze data after experiments. PhaseLab assesses reliability before or with minimal data.

### 3. Honest About Uncertainty

Most tools give point estimates. PhaseLab gives claim levels:
- UNKNOWN
- EXPLORATORY
- CONTEXT_DEPENDENT
- STRONG_COMPUTATIONAL

### 4. Falsifiable

Most tools are validated by success stories. PhaseLab specifies exactly what would prove it wrong.

---

## The Non-Negotiables

### Always True

1. PhaseLab assesses reliability, not outcomes
2. Claim levels never upgrade without evidence
3. Wet lab validation is always required for therapeutic claims
4. Quantum mode is OFF by default
5. The framework is falsifiable

### Never Claimed

1. PhaseLab never claims to predict biology
2. PhaseLab never claims to replace experiments
3. PhaseLab never claims quantum is necessary
4. PhaseLab never claims universal accuracy

---

## How to Talk About PhaseLab

### Good

- "PhaseLab helps you assess experimental reliability"
- "PhaseLab identifies stable vs. unreliable regions"
- "PhaseLab works across CRISPR, proteins, and chemistry"
- "PhaseLab quantifies uncertainty with claim levels"

### Bad

- "PhaseLab predicts which guides will work" (wrong focus)
- "PhaseLab uses quantum to improve accuracy" (misleading)
- "PhaseLab is a new CRISPR tool" (wrong category)
- "PhaseLab eliminates the need for experiments" (false)

---

## The Long-Term Vision

### Year 1: Establish Credibility

- Methods paper on two-stage framework
- Real dataset validation (MaveDB)
- First wet lab validation (RAI1 tiling)

### Year 2: Expand Domains

- Additional protein datasets
- Chemistry validation
- Microbiology validation

### Year 3: Enable Others

- Public API
- Integration with existing tools
- Community adoption

### Long-Term: Change the Standard

> Every perturbation experiment reports claim levels.
> Reproducibility is quantified, not assumed.
> Reliability assessment is standard practice.

---

## Summary

PhaseLab is simple:

1. **One abstraction**: Perturbation → Response → Reliability
2. **One metric**: Spatial coherence
3. **One claim**: Coherence predicts variance
4. **One purpose**: Help scientists trust their results

Everything else—the modules, the domains, the quantum modes—serves that purpose.

---

*PhaseLab v1.0.0*
*December 2025*
