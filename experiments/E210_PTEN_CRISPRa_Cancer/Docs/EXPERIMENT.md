# E210: PTEN CRISPRa for Cancer Tumor Suppressor Reactivation

## Overview

This experiment applies PhaseLab's CRISPRa guide design pipeline to reactivate the
PTEN tumor suppressor gene in cancer cells. Based on Moses et al. 2019 (Mol Ther
Nucleic Acids), which demonstrated that dCas9-VPR CRISPRa can activate PTEN
expression and suppress oncogenic signaling in melanoma and TNBC cell lines.

## Scientific Background

### PTEN Tumor Suppressor

PTEN (Phosphatase and Tensin Homolog) is one of the most frequently mutated or
silenced tumor suppressor genes in human cancer:

- **Function**: Antagonizes PI3K signaling by dephosphorylating PIP3 to PIP2
- **Loss consequences**: Hyperactivation of AKT/mTOR/MAPK pathways
- **Cancer prevalence**:
  - ~40% of triple-negative breast cancer (TNBC)
  - ~30% of melanoma
  - ~70% of glioblastoma
  - ~60% of prostate cancer

### Therapeutic Rationale

Unlike traditional gene therapy (which adds exogenous genes), CRISPRa reactivates
the endogenous PTEN locus. This approach:

1. Uses the native gene with proper regulation
2. Restores physiological expression levels
3. Avoids insertional mutagenesis risks
4. Can potentially sensitize tumors to targeted therapies

### Literature Foundation

**Moses et al. 2019** (PMID: 30654190):
- Designed sgRNAs targeting PTEN proximal promoter (-300 to -50 bp from TSS)
- Best guide (sgRNA -54) achieved 2.27-fold activation in SUM159 cells
- Activation reduced p-AKT, p-mTOR, p-S6K, p-ERK signaling
- Reduced colony formation and increased sensitivity to BRAF + PI3K/mTOR inhibitors

## Experimental Design

### Target Gene
- **Gene**: PTEN
- **Coordinates**: chr10:87,862,625-87,864,125 (GRCh38)
- **TSS**: 87,863,625
- **CRISPRa window**: -300 to -50 bp from TSS

### Cell Models
1. **SK-MEL-28** (melanoma)
   - BRAF V600E mutation
   - PTEN wild-type but low expression
   - Intrinsic resistance to dabrafenib

2. **SUM159** (TNBC)
   - Mesenchymal stem-like subtype
   - PTEN wild-type but low expression
   - Highly aggressive phenotype

### PhaseLab Pipeline

1. **PAM site scanning**: NGG sites in CRISPRa window
2. **Sequence filtering**: GC content, homopolymers, complexity
3. **Off-target scoring**: MIT algorithm, CFD score
4. **Chromatin modeling**: ENCODE cancer cell accessibility data
5. **Thermodynamic scoring**: SantaLucia nearest-neighbor model
6. **Quantum simulation**: VQE-based binding energy
7. **IR coherence**: R̄ metric with e⁻² threshold
8. **Hardware validation**: IBM Quantum backend

### Output Metrics

For each guide candidate:
- Position relative to TSS
- GC content (overall and seed region)
- Off-target specificity score
- Chromatin accessibility score
- Thermodynamic binding energy (ΔG)
- Quantum binding energy
- IR coherence (R̄) with GO/NO-GO classification

## Expected Results

### Guide Selection Criteria
- GC content: 40-70%
- Position: Preferably -50 to -100 bp from TSS (per Moses et al.)
- Chromatin: Open regions (DNase HS)
- Coherence: R̄ > 0.135 (e⁻² threshold)
- No homopolymer runs ≥5 bp

### Validation Plan
1. **Computational**: CRISPOR genome-wide off-target analysis
2. **In silico**: Compare with Moses et al. sgRNA-54 performance
3. **Experimental**: PTEN activation in SK-MEL-28 and SUM159 cells

## Clinical Implications

### Potential Applications
1. **Monotherapy**: PTEN reactivation alone
2. **Combination therapy**: CRISPRa + BRAF inhibitors (melanoma)
3. **Combination therapy**: CRISPRa + PI3K/mTOR inhibitors
4. **Biomarker**: PTEN expression as response predictor

### Delivery Considerations
- **Vector**: AAV or lentivirus
- **Serotype**: AAV9 (broad tropism) or tissue-specific
- **Packaging**: dCas9-VPR + sgRNA within 4.7kb limit
- **Immunogenicity**: Monitor for anti-Cas9 responses

## References

1. Moses C, et al. "Activating PTEN Tumor Suppressor Expression with the
   CRISPR/dCas9 System." Mol Ther Nucleic Acids. 2019;14:287-300.
   PMID: 30654190

2. Li J, et al. "PTEN, a Putative Protein Tyrosine Phosphatase Gene Mutated
   in Human Brain, Breast, and Prostate Cancer." Science. 1997;275:1943-1947.
   PMID: 9090379

3. Song MS, et al. "The functions and regulation of the PTEN tumour suppressor."
   Nat Rev Mol Cell Biol. 2012;13:283-296. PMID: 22473468

## Files

- `Code/E210_pten_crispra_cancer.py` - Main experiment script
- `Data/` - Output JSON files with guide rankings
- `Figures/` - Visualization outputs
- `Docs/EXPERIMENT.md` - This documentation
- `Docs/FINDINGS.md` - Results and conclusions (post-run)
