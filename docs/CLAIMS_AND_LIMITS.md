# PhaseLab: Claims and Limits

## What PhaseLab Is

**PhaseLab is a reliability layer for perturbation science.**

It assesses whether experimental results from perturbation experiments (CRISPR, mutagenesis, chemical modification) will be reproducible before or with minimal data.

---

## The Core Scientific Claim

> **Spatial coherence of response landscapes predicts perturbation reliability across biological and chemical domains, enabling minimum viable experimental design.**

This claim is:
- **Narrow**: About reliability, not biology
- **Testable**: Coherence-variance correlation is measurable
- **Falsifiable**: If coherence doesn't correlate with variance, the claim fails
- **Domain-general**: Same math applies to CRISPR, proteins, chemistry

---

## What PhaseLab Guarantees

### Guaranteed (Mathematical)

| Guarantee | Basis |
|-----------|-------|
| Coherence computation is deterministic | Algorithm specification |
| StabilityClass assignment follows defined thresholds | Code implementation |
| Claim levels propagate conservatively | Worst-case aggregation |
| Cross-domain API consistency | Unified ResponseLandscape abstraction |

### Guaranteed (Empirical, E213-E216)

| Guarantee | Evidence |
|-----------|----------|
| Coherence correlates with outcome variance | r = -0.24 to -0.50 across 115,251 sgRNAs |
| Stable regions show lower variance | 32-49% variance reduction |
| Guide-sequence coherence does NOT work | r ≈ 0 in E200-E211 |

---

## What PhaseLab Estimates

### Estimates (With Uncertainty)

| Estimate | Uncertainty Source |
|----------|-------------------|
| Variance reduction magnitude | Depends on target, window size |
| Stage I stable region predictions | Structural prior quality |
| Therapeutic window boundaries | Model parameter uncertainty |
| Circadian rescue predictions | Simplified Kuramoto model |

### Estimates (Context-Dependent)

| Estimate | Context |
|----------|---------|
| Minimum viable tiling size | Gene-specific, 16-20 is typical |
| Coherence threshold for stability | 0.7 default, may vary by domain |
| Effect size thresholds | Application-specific |

---

## What Requires Wet Lab Validation

### Always Requires Validation

| Claim Type | Why |
|------------|-----|
| Absolute effect sizes | PhaseLab predicts reliability, not magnitude |
| Off-target activity | Requires experimental measurement |
| Cellular toxicity | Cannot be computed |
| In vivo efficacy | Delivery, biodistribution, immunogenicity |
| Therapeutic window | Patient-specific factors |

### Validation Upgrades Claim Level

| From | To | Requires |
|------|----|---------|
| EXPLORATORY | CONTEXT_DEPENDENT | Tiling data for target |
| CONTEXT_DEPENDENT | STRONG_COMPUTATIONAL | Cross-validated tiling |
| STRONG_COMPUTATIONAL | Therapeutic | Clinical validation |

---

## What Is Falsifiable

### Falsification Tests (Any of these failing invalidates the framework)

| Test | Failure Condition |
|------|-------------------|
| **A: Ranking validity** | PhaseLab-ranked guides do NOT outperform random controls |
| **B: Risk prediction** | High-risk flagged guides show no elevated failure rate |
| **C: Variance correlation** | Coherence does NOT correlate with outcome variance |
| **D: Stage I accuracy** | Stage I predictions match < 50% of Stage II results |

### What Would NOT Falsify PhaseLab

| Observation | Interpretation |
|-------------|----------------|
| A specific guide fails | Individual outcomes vary; framework predicts variance |
| A target has no stable regions | Correct behavior for amplifying targets |
| Effect sizes differ from prediction | PhaseLab predicts reliability, not magnitude |
| Quantum mode gives same result as classical | Classical often sufficient |

---

## What PhaseLab Does NOT Claim

### Explicit Non-Claims

| Non-Claim | Why Not |
|-----------|---------|
| "Predicts biology" | Predicts measurement reliability |
| "Replaces experiments" | Experiments remain essential |
| "Designs optimal guides" | Use CRISPOR/CRISPick for that |
| "Guarantees therapeutic success" | Reliability ≠ efficacy |
| "Works without any data" | Stage II requires tiling |
| "Quantum is always better" | Usually OFF is correct |

### What PhaseLab Is NOT

- Not a CRISPR design tool
- Not a biology predictor
- Not a quantum computing platform
- Not a replacement for experimental validation
- Not an accuracy guarantee

---

## Claim Level Semantics

| Level | Meaning | Requirements |
|-------|---------|--------------|
| **UNKNOWN** | Cannot assess reliability | Insufficient data |
| **EXPLORATORY** | Preliminary estimate only | Structural priors only |
| **CONTEXT_DEPENDENT** | Valid within tested context | Single tiling dataset |
| **STRONG_COMPUTATIONAL** | High confidence | Cross-validated, multiple sources |

### Claim Level Rules

1. Claims never upgrade without new evidence
2. Claims propagate at the minimum of inputs
3. UNKNOWN cannot be upgraded computationally
4. Only wet lab can reach therapeutic claims

---

## Quantum Mode Guidance

| Module | Default Mode | When to Change |
|--------|--------------|----------------|
| CRISPRa/i feasibility | OFF | Never for standard use |
| Protein DMS | OFF | Never for standard use |
| Binding landscapes | OFF | Never for standard use |
| Circadian simulation | OFF | Research exploration only |
| Representation bounds | AUDIT | Research validation |
| Measurement limits | REQUIRED | Fundamental research only |

**Rule**: If a classical experiment can falsify a claim, quantum mode should be OFF.

---

## Summary Statement

PhaseLab makes one claim, makes it carefully, and provides tools to falsify it.

The claim is:

> **Spatial coherence predicts reliability.**

Everything else—the modules, the domains, the applications—are implementations of that single insight.

If the correlation between coherence and variance disappears under proper testing, PhaseLab fails. That's the deal.

---

*PhaseLab v1.0.0*
*December 2025*
