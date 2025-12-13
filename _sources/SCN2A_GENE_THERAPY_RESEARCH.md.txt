# SCN2A Gene Therapy Research for Autism-Linked Haploinsufficiency

**CRISPRa Guide RNA Design Using Quantum-Validated Coherence Metrics**

*Dylan Vaca | December 2025*

---

## Abstract

We present a CRISPRa guide RNA design approach for SCN2A-related neurodevelopmental disorders, the most common monogenic cause of autism-linked haploinsufficiency. Using IBM Quantum hardware (Torino, 133 qubits) and genome-wide off-target analysis (CRISPOR), we identified and validated optimal guide RNAs for upregulating SCN2A expression. Our lead candidate (`GCTGACTGCTACATAGCCAA`) achieved excellent specificity scores (MIT 83, CFD 89) with zero off-targets at ≤2 mismatches, validated by quantum coherence R̄ = 0.970 on real hardware. This work extends the PhaseLab framework previously validated for RAI1/Smith-Magenis Syndrome to a second therapeutic target.

---

## 1. Background

### 1.1 SCN2A and Autism

SCN2A encodes the voltage-gated sodium channel NaV1.2, which is critical for excitatory neuron function during early brain development. Loss-of-function mutations leading to haploinsufficiency cause:

- **Autism Spectrum Disorder (ASD)**
- **Intellectual disability**
- **Epilepsy** (in a subset of patients)

SCN2A-related neurodevelopmental disorders (SCN2A-NDD) represent one of the most frequent monogenic causes of ASD.

**Genetic Basis:**
Most affected individuals have heterozygous loss-of-function variants, leaving approximately 50% of normal SCN2A expression. Unlike gain-of-function variants (which cause early-onset epilepsy), loss-of-function variants are ideal targets for CRISPRa upregulation.

### 1.2 Therapeutic Evidence

A landmark 2025 Nature paper demonstrated that CRISPRa upregulation of the remaining Scn2a allele in haploinsufficient mice:

- Restores Scn2a expression toward normal levels
- Rescues electrophysiological deficits in neurons
- Reduces seizures and improves behavioral phenotypes
- **Works even when delivered in adolescence** (critical for clinical translation)

Reference: [Nature 2025 - DOI: 10.1038/s41586-025-09522-w](https://www.nature.com/articles/s41586-025-09522-w)

### 1.3 The IR Coherence Framework

This research applies the **Informational Relativity (IR) framework** to assess simulation reliability:

| Component | Application |
|-----------|-------------|
| **R̄ = exp(-V_φ/2)** | Coherence metric for binding quality |
| **e⁻² threshold** | GO/NO-GO decision boundary (≈0.135) |
| **Pauli Hamiltonian** | DNA-gRNA binding encoding |
| **Hardware validation** | IBM Torino 133-qubit processor |

---

## 2. Methods

### 2.1 Genomic Data

**Source:** NCBI RefSeq NC_000002.12 (GRCh38.p14)

**Region Analyzed:**
```
Chromosome:  2
Cytoband:    2q24.3
Coordinates: 165,238,414 - 165,239,914 (1500 bp promoter)
TSS:         165,239,414
Gene:        SCN2A (sodium voltage-gated channel alpha subunit 2)
```

### 2.2 gRNA Design Pipeline

Our multi-layer validation pipeline (identical to RAI1/E200):

| Layer | Method | Purpose |
|-------|--------|---------|
| PAM Scanning | NGG (SpCas9) | Find Cas9 binding sites |
| GC Content | 40-70% filter | Optimal hybridization |
| Thermodynamics | SantaLucia ΔG | Binding energy prediction |
| Chromatin | DNase HS model (brain) | Accessibility assessment |
| MIT Score | Position-weighted | Off-target specificity |
| CFD Score | Mismatch penalty | Cutting frequency |
| **IR Coherence** | R̄ metric | Simulation reliability |

### 2.3 Quantum Simulation

**Hamiltonian Encoding:**
- DNA-RNA base pairing → Pauli ZZ terms
- Stacking interactions → XX/YY terms
- Nearest-neighbor energy model (SantaLucia parameters)

**Hardware:**
- IBM Torino (133 superconducting qubits)
- 7 circuits per guide (Pauli term measurements)
- 2,000 shots per circuit

### 2.4 CRISPRa Window

For transcriptional activation, guides must bind within:
- **Optimal window:** -400 to -50 bp from TSS
- **Mechanism:** dCas9-VP64 recruits transcriptional machinery
- **Chromatin:** Must be in accessible (open) region

---

## 3. Results

### 3.1 Candidate Identification

From 75 PAM sites in the SCN2A promoter, we identified 18 candidates in the optimal CRISPRa window.

**Initial Screening Results:**

| Metric | Value |
|--------|-------|
| Total PAM sites | 75 |
| CRISPRa window candidates | 18 |
| After quality filtering | 15 |
| GO classification rate | 100% (15/15) |
| Average coherence (R̄) | 0.944 |

### 3.2 IBM Quantum Hardware Validation

**Run Details:**
- Date: December 10, 2025
- Backend: IBM Torino (133 qubits)
- Cost: FREE (IBM Quantum Open Plan)

**Hardware Results:**

| Guide | Sequence | Position | Simulator R̄ | Hardware R̄ | Δ | Status |
|-------|----------|----------|--------------|-------------|-----|--------|
| 1 | `TTCCACTTTTGACCAGGAGA` | -64 bp | 0.944 | 0.971 | 2.9% | **GO** |
| 2 | `AGATGGTTCCACTTTTGACC` | -70 bp | 0.945 | 0.970 | 2.6% | **GO** |
| **3** | `GCTGACTGCTACATAGCCAA` | -104 bp | 0.945 | **0.970** | 2.7% | **GO** |

**Average simulator-hardware agreement: 2.7%**

**IBM Quantum Job IDs:**
- `d4t29p7t3pms7399q6ag` (Guide 1)
- `d4t29rkfitbs739jh28g` (Guide 2)
- `d4t29ts5fjns73d313ag` (Guide 3)

### 3.3 CRISPOR Genome-Wide Validation

**Analysis Details:**
- Tool: CRISPOR v5.2 (crispor.tefor.net)
- Genome: Human GRCh38/hg38
- PAM: NGG (SpCas9)

**CRISPOR Results:**

| Guide | Sequence | MIT | CFD | Off-targets (0-1-2-3-4 mm) | Total OT | Verdict |
|-------|----------|-----|-----|----------------------------|----------|---------|
| 1 | `TTCCACTTTTGACCAGGAGA` | 67 | 76 | 0-0-1-25-217 | 243 | Has 1 at 2mm |
| 2 | `AGATGGTTCCACTTTTGACC` | 79 | 89 | 0-0-3-9-110 | 122 | Has 3 at 2mm |
| **3** | `GCTGACTGCTACATAGCCAA` | **83** | **89** | **0-0-0-10-89** | 99 | **BEST** |

**Key Finding:** Guide 3 has zero off-targets at ≤2 mismatches, making it the safest candidate.

---

## 4. Lead Candidate

### 4.1 Primary Selection

Based on combined quantum and CRISPOR validation:

**Primary: Guide 3**
```
Sequence:    GCTGACTGCTACATAGCCAA
PAM:         AGG
Position:    -104 bp from TSS
GC Content:  50%
MIT Score:   83 (excellent)
CFD Score:   89 (excellent)
Off-targets: 0 at ≤2 mismatches
Hardware R̄: 0.970
Status:      GO
```

### 4.2 Selection Rationale

1. **Highest specificity:** MIT 83 = excellent genome-wide specificity
2. **Zero critical off-targets:** No off-targets at ≤2 mismatches (safest)
3. **Optimal GC:** 50% GC content in ideal range
4. **Open chromatin:** Located in DNase hypersensitive region
5. **Validated coherence:** R̄ = 0.970 on IBM Torino

---

## 5. Comparison with RAI1/SMS (E200)

| Metric | E200 (RAI1/SMS) | E201 (SCN2A/Autism) |
|--------|-----------------|---------------------|
| Disease | Smith-Magenis Syndrome | Autism-linked NDD |
| Haploinsufficiency | Yes | Yes |
| CRISPRa mouse model | Yes (PMC9762195) | Yes (Nature 2025) |
| Total PAM sites | 274 | 75 |
| CRISPRa candidates | 91 | 18 |
| **Lead guide MIT** | 83 | 83 |
| **Lead guide CFD** | 93 | 89 |
| **Off-targets ≤2mm** | 0 | 0 |
| **Hardware R̄** | 0.839 | 0.970 |
| Hardware platform | IBM Torino | IBM Torino |

**Key observation:** Both targets yielded lead candidates with:
- MIT specificity scores of 83 (identical)
- Zero off-targets at ≤2 mismatches
- Confirmed GO classification on IBM hardware

---

## 6. Discussion

### 6.1 Framework Generalization

This study demonstrates that the PhaseLab/IR framework generalizes across:

1. **Different genes** (RAI1 → SCN2A)
2. **Different chromosomes** (chr17 → chr2)
3. **Different diseases** (SMS → Autism-linked NDD)

The consistent performance (100% GO rate, similar MIT scores, excellent hardware agreement) suggests the framework captures fundamental properties of gRNA-DNA binding.

### 6.2 Hardware Coherence

The SCN2A guides achieved higher hardware coherence (0.970) than RAI1 guides (0.839). This may reflect:
- Different sequence composition
- Position within promoter architecture
- Thermodynamic stability differences

Importantly, both are well above the GO threshold (e⁻² ≈ 0.135).

### 6.3 Clinical Relevance

The Nature 2025 paper establishing CRISPRa efficacy for Scn2a used AAV delivery to mouse brain. Our lead candidate provides a human-optimized guide for:

- iPSC-derived neuron testing
- Human clinical trial design
- Companion diagnostic development

---

## 7. Recommended Next Steps

### 7.1 Wet Lab Validation

1. **Order Synthesis**
   - Primary: `GCTGACTGCTACATAGCCAA` (Guide 3)
   - Source: IDT, Synthego, or similar

2. **Clone into Vector**
   - dCas9-VP64 system (e.g., Addgene #135338)
   - Alternative: dCas9-VPR for stronger activation

3. **Cell Culture Testing**
   - Initial: HEK293T (high transfection efficiency)
   - Disease-relevant: iPSC-derived neurons from SCN2A patients

4. **Readouts**
   - qPCR for SCN2A mRNA (48-72h post-transfection)
   - Western blot for NaV1.2 protein
   - Patch clamp electrophysiology (functional validation)

### 7.2 Research Partnerships

We recommend contacting:
- **SFARI** (Simons Foundation Autism Research Initiative)
- The lab that published the Nature 2025 paper
- AAV/gene therapy companies with CNS programs

---

## 8. Limitations

1. **Computational only** - No wet lab validation yet
2. **Chromatin model** - Brain-specific but not patient-derived
3. **Off-target prediction** - Algorithmic, not empirical GUIDE-seq
4. **Binding ≠ Activation** - Strong binding doesn't guarantee efficient activation
5. **Delivery** - In vivo CNS delivery remains challenging

---

## 9. Conclusion

We have identified and validated a lead CRISPRa guide RNA candidate for SCN2A upregulation in autism-linked neurodevelopmental disorders:

**`GCTGACTGCTACATAGCCAA AGG`**

This candidate was validated through:
- Quantum coherence R̄ = 0.970 on IBM Torino
- CRISPOR specificity (MIT 83, CFD 89)
- Zero off-targets with ≤2 mismatches
- Optimal position (-104 bp) and GC content (50%)

The PhaseLab/IR coherence framework successfully generalized from RAI1/SMS to SCN2A/Autism, demonstrating its utility as a cross-target gRNA reliability metric. This is the second haploinsufficiency disorder for which we have provided quantum-validated CRISPRa guide designs.

---

## References

1. Nature 2025 - CRISPR activation for SCN2A-related neurodevelopmental disorders
   - DOI: 10.1038/s41586-025-09522-w

2. SFARI - Development of CRISPR activation therapeutics to rescue SCN2A function
   - https://www.sfari.org/funded-project/development-of-crispr-activation-therapeutics-to-rescue-scn2a-function/

3. Hsu et al. (2013) - MIT specificity score
4. Doench et al. (2016) - CFD score methodology
5. SantaLucia (1998) - Nearest-neighbor thermodynamic parameters
6. CRISPOR - https://crispor.tefor.net/
7. IBM Quantum - https://quantum.ibm.com/

---

## Appendix A: Data Availability

| Resource | Location |
|----------|----------|
| Quantum results | IBM Quantum Jobs: `d4t29p7t3pms7399q6ag`, `d4t29rkfitbs739jh28g`, `d4t29ts5fjns73d313ag` |
| CRISPOR analysis | http://crispor.tefor.net/ (GRCh38, SCN2A promoter) |
| PhaseLab package | https://github.com/followthesapper/phaselab |
| Experiment data | `experiments/E201_SCN2A_CRISPRa/Data/` |

## Appendix B: Code Availability

All code is available in the PhaseLab Python package:

```bash
pip install phaselab[quantum]
```

Example usage:

```python
from phaselab.targets import load_target_config
from phaselab.crispr import design_guides

# Load SCN2A configuration
scn2a = load_target_config("SCN2A")
print(f"Gene: {scn2a.gene_symbol}")
print(f"TSS: chr{scn2a.chrom}:{scn2a.tss_genomic}")

# Design guides (with promoter sequence)
guides = design_guides(scn2a_promoter, tss_index=scn2a.tss_index)
top_candidates = guides[guides['go_no_go'] == 'GO'].head(5)
print(top_candidates[['sequence', 'position', 'coherence_R', 'go_no_go']])
```

---

*This research was conducted as part of the Informational Relativity experimental series.*

*Hardware validation: IBM Torino, December 2025*

*Extends E200 (RAI1/SMS) methodology to a second therapeutic target.*
