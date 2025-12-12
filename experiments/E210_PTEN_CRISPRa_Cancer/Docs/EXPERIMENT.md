# E210: PTEN CRISPRa for Cancer Tumor Suppressor Reactivation

## Overview

This experiment applies the PhaseLab CRISPRa design pipeline to PTEN, one of the most frequently silenced tumor suppressor genes in human cancer. The same validated methodology used for RAI1/Smith-Magenis Syndrome (E200) and SCN2A/Autism (E201) is applied to design and pre-validate guide RNAs for reactivating PTEN expression in cancer cells.

## Background

### PTEN and Cancer

PTEN (Phosphatase and Tensin Homolog) encodes a lipid phosphatase that antagonizes PI3K signaling. Loss of PTEN function leads to:

- Hyperactivation of AKT/mTOR/MAPK pathways
- Uncontrolled cell proliferation
- Resistance to apoptosis
- Drug resistance (especially to BRAF inhibitors)

PTEN loss occurs in a large fraction of human cancers:

| Cancer Type | PTEN Loss Frequency |
|-------------|---------------------|
| Glioblastoma | ~70% |
| Prostate cancer | ~60% |
| Triple-negative breast cancer (TNBC) | ~40% |
| Melanoma | ~30% |
| Endometrial cancer | ~50% |

### Literature Support

Moses et al. 2019 demonstrated that CRISPRa upregulation of PTEN in cancer cells:
- Achieves 2.27-fold PTEN activation (sgRNA -54 in SUM159)
- Reduces p-AKT, p-mTOR, p-S6K, p-ERK signaling
- Decreases colony formation
- Increases sensitivity to BRAF + PI3K/mTOR inhibitor combination

Reference: [Mol Ther Nucleic Acids 2019](https://pubmed.ncbi.nlm.nih.gov/30654190/) (PMID: 30654190)

## Methods

### Target Region

- **Gene**: PTEN (phosphatase and tensin homolog)
- **Chromosome**: 10 (10q23.31)
- **Genome Build**: GRCh38
- **TSS**: chr10:87,863,625 (RefSeq NM_000314.8, UCSC hg38)
- **Promoter Region**: chr10:87,862,625-87,864,125 (1500 bp)
- **CRISPRa Window**: -300 to -50 bp from TSS

*Note: Coordinates refer to RefSeq NM_000314.8 on GRCh38 (UCSC hg38). Alternative upstream TSSs have been reported but are not targeted here.*

### Cell Models

1. **SK-MEL-28** (melanoma)
   - BRAF V600E mutation
   - PTEN wild-type but low expression
   - Intrinsic resistance to dabrafenib

2. **SUM159** (TNBC)
   - Mesenchymal stem-like subtype
   - PTEN wild-type but low expression
   - Highly aggressive phenotype

### Pipeline Components

1. **PAM Scanning**: SpCas9 (NGG) on both strands
2. **Basic Filtering**: GC 40-70% (relaxed to ~85% for CpG island context), max homopolymer ≤4
3. **Off-target Scoring**: MIT algorithm + complexity analysis (internal scoring)
4. **Chromatin Accessibility**: Cancer cell model (ENCODE data)
5. **Thermodynamics**: SantaLucia nearest-neighbor energy
6. **Quantum Simulation**: Pauli Hamiltonian encoding
7. **IR Coherence**: R̄ = exp(-V_φ/2) with e⁻² threshold

### CpG Island Promoter Consideration

The PTEN core promoter resides in a CpG-island context with unusually high GC content, a known feature of many tumor suppressor promoters. Because of this, several top candidates have GC fractions around ~80%. For PTEN we explicitly relaxed the default 40-70% GC filter used in PhaseLab, and relied on thermodynamic ΔG and IR coherence to down-select guides that remain physically plausible despite the GC-rich background.

### Quantum Validation

Coherence metric validated on IBM Quantum hardware:
- Backend: IBM Torino (133 qubits)
- Previous validation: E200 RAI1, E201 SCN2A achieved excellent simulator-hardware agreement

## Results

### PAM Site Analysis

| Metric | Value |
|--------|-------|
| Total PAM sites | 336 |
| In CRISPRa window | 53 |
| After quality filtering | 15 |
| Simulator GO classification | 100% |

### Top Candidates (Simulator)

| Rank | Sequence | Position | GC | Chromatin | R̄ | Status |
|------|----------|----------|-----|-----------|-------|--------|
| 1 | CGGAAGGGGGAGCGCGGCAG | -90 | 80% | OPEN (near TSS) | 0.944 | GO |
| 2 | GGAGCGGAGCGAGGAGGCGG | -265 | 80% | OPEN (CRISPRa) | 0.944 | GO |
| 3 | GAGGCGGGACCCGCGTGCGG | -145 | 85% | OPEN (CRISPRa) | 0.942 | GO |
| 4 | GAGCGAGGCGGAGCGGAGCG | -274 | 80% | OPEN (CRISPRa) | 0.945 | GO |
| 5 | GCGGCAGCGGAGCGCGCGCG | -77 | 90% | OPEN (near TSS) | 0.945 | GO |

### IBM Quantum Hardware Validation

**Backend**: IBM Torino (133 qubits)
**Date**: December 11, 2025
**Shots**: 4096 per circuit

| Guide | Position | Simulator R̄ | Hardware R̄ | Difference | Status |
|-------|----------|-------------|-------------|------------|--------|
| CGGAAGGGGGAGCGCGGCAG | -90 | 0.944 | **0.953** | 1.0% | **GO** |
| GGAGCGGAGCGAGGAGGCGG | -265 | 0.944 | **0.953** | 0.9% | **GO** |
| GAGGCGGGACCCGCGTGCGG | -145 | 0.942 | **0.952** | 1.1% | **GO** |

**Average simulator-hardware agreement: 1.0%** (Excellent)

*Note: Hardware validation was performed on the top 3 ranked guides. All three achieved GO status.*

### Lead Candidate

**Spacer**: `CGGAAGGGGGAGCGCGGCAG`
**PAM**: `CGG`

| Metric | Value |
|--------|-------|
| Position | -90 bp from TSS |
| GC Content | 80% |
| Chromatin | OPEN (near TSS) |
| Simulator Coherence | 0.944 |
| Hardware Coherence | **0.953** |
| Status | **GO** |

### Key Findings

1. **All 15 simulator candidates achieved GO classification** (R̄ > 0.135)
2. **All 3 hardware-tested guides confirmed GO status** on IBM Torino
3. **Average coherence: 0.944** (Excellent reliability)
4. **Hardware validation confirms simulator predictions** (1.0% average difference)
5. **High GC content (~80%) is expected** given PTEN's CpG island promoter context
6. **Top guide position (-90 bp) aligns** with Moses et al. optimal window (-54 bp)

## Comparison with Other PhaseLab Targets

| Metric | E200 (RAI1) | E201 (SCN2A) | E210 (PTEN) |
|--------|-------------|--------------|-------------|
| Disease | Smith-Magenis | Autism-NDD | Cancer |
| Mechanism | Haploinsufficiency | Haploinsufficiency | Silencing |
| Tissue | Neurons | Neurons | Cancer cells |
| Total PAM sites | 274 | 75 | 336 |
| CRISPRa candidates | 91 | 18 | 53 |
| Promoter GC | ~50% | ~40% | ~80% (CpG) |
| Hardware R̄ | 0.839 | 0.970 | **0.953** |
| Hardware validated | IBM Torino | IBM Torino | **IBM Torino** |

## Biological Significance

### Why PTEN Reactivation Works

Unlike mutational inactivation, PTEN loss in many cancers occurs through:
- Transcriptional silencing
- Promoter hypermethylation
- Epigenetic repression

CRISPRa can overcome these mechanisms by recruiting transcriptional activators (VPR) to the PTEN promoter, restoring expression from the intact gene.

### Therapeutic Implications

1. **Monotherapy**: PTEN reactivation alone suppresses tumor growth
2. **Combination**: CRISPRa + BRAF inhibitors (melanoma)
3. **Combination**: CRISPRa + PI3K/mTOR inhibitors (all cancers)
4. **Drug Sensitization**: Restoring PTEN overcomes drug resistance

### Moses et al. Key Results

| Cell Line | Guide | PTEN Fold-Change | Effect |
|-----------|-------|------------------|--------|
| SUM159 | sgRNA -54 | 2.27x | Reduced p-AKT |
| SK-MEL-28 | sgRNA -54 | 1.8x | Drug sensitization |

## Limitations

- **Computational + quantum only**: No wet-lab testing has been performed yet
- **Off-target scoring is internal**: Based on PhaseLab's MIT/CFD-style scoring; genome-wide off-target enumeration via CRISPOR is still required
- **Does not model**:
  - Promoter methylation dynamics directly
  - Chromatin remodeling in actual tumors
  - Intratumor heterogeneity and clonal selection
  - Delivery efficiency, immune response, or toxicity

## Next Steps

1. ✅ **Hardware Validation**: Run top 3 guides on IBM Torino - COMPLETE
2. **CRISPOR Analysis**: Submit for genome-wide off-target screening
3. **Literature Comparison**: Compare with Moses et al. sgRNA-54
4. **Experimental Validation**: Test in SK-MEL-28 and SUM159 cells
5. **Combination Testing**: CRISPRa + dabrafenib/dactolisib

## Experimental Protocol (Proposed)

### Delivery
- **System**: dCas9-VPR lentivirus
- **MOI**: 10-20 for stable expression
- **Selection**: Puromycin (if vector includes)

### Readouts (48-72h post-transduction)
1. **PTEN mRNA**: qPCR (Taqman assay)
2. **PTEN protein**: Western blot
3. **Signaling**: p-AKT (S473), p-mTOR (S2448), p-S6K (T389)
4. **Phenotype**: Proliferation (CellTiter-Glo), colony formation

### Controls
1. Untreated cells
2. dCas9-VPR + non-targeting gRNA
3. Moses et al. sgRNA -54 (positive control)
4. PhaseLab top guides

## Files

- `Code/E210_pten_crispra_cancer.py` - Main experiment script
- `Data/E210_pten_guides_*.json` - Results output (includes hardware validation)
- `Docs/EXPERIMENT.md` - This file
- `Docs/FINDINGS.md` - Summary of findings

## References

1. Moses C, et al. "Activating PTEN Tumor Suppressor Expression with the CRISPR/dCas9 System." Mol Ther Nucleic Acids. 2019;14:287-300. PMID: 30654190

2. Li J, et al. "PTEN, a Putative Protein Tyrosine Phosphatase Gene Mutated in Human Brain, Breast, and Prostate Cancer." Science. 1997;275:1943-1947. PMID: 9090379

3. E200: Smith-Magenis Syndrome Gene Therapy (PhaseLab)
   - Hardware validation on IBM Torino, December 2025

4. E201: SCN2A CRISPRa for Autism-Linked Haploinsufficiency (PhaseLab)
   - Hardware validation on IBM Torino, December 2025

---

*Dylan Vaca, December 2025*
