# E210 FINDINGS: PTEN CRISPRa for Cancer

## Experiment Status: ✅ COMPLETE (Hardware Validated)

**Date**: December 11, 2025
**Backend**: IBM Torino (133 qubits)
**Shots**: 4096 per circuit

---

## Executive Summary

PhaseLab successfully designed and pre-validated CRISPRa guides for PTEN tumor suppressor reactivation in cancer cells. All three top-ranked guides achieved GO status on IBM Quantum hardware, providing strong candidates for experimental validation.

**Bottom Line**: PhaseLab can design and pre-validate CRISPRa guides for tumor suppressor reactivation in cancer models, providing strong candidates for gene-therapy-oriented experiments.

---

## Key Results

### Hardware Validation (IBM Torino)

| Guide | Position | Sim R̄ | Hardware R̄ | Agreement |
|-------|----------|--------|-------------|-----------|
| CGGAAGGGGGAGCGCGGCAG | -90 bp | 0.944 | **0.953** | 1.0% |
| GGAGCGGAGCGAGGAGGCGG | -265 bp | 0.944 | **0.953** | 0.9% |
| GAGGCGGGACCCGCGTGCGG | -145 bp | 0.942 | **0.952** | 1.1% |

**All 3 top-ranked guides: GO status confirmed on real quantum hardware**

### Pipeline Statistics

| Metric | Value |
|--------|-------|
| Total PAM sites scanned | 336 |
| CRISPRa window candidates | 53 |
| Evaluated candidates | 15 |
| Simulator GO classification | **100%** |
| Hardware-tested guides | 3 |
| Hardware GO classification | **100%** |
| Average coherence (R̄) | 0.944 |
| Hardware-simulator agreement | **1.0%** |

---

## Top Recommendation

### Lead Guide

**Sequence**: `CGGAAGGGGGAGCGCGGCAG`
**PAM**: CGG
**Position**: -90 bp from TSS

| Metric | Value |
|--------|-------|
| GC Content | 80% |
| Chromatin State | OPEN (near TSS) |
| Simulator R̄ | 0.944 |
| **Hardware R̄** | **0.953** |
| Status | **GO** |

### Why This Guide?

1. **Position**: -90 bp is in the optimal CRISPRa window (Moses et al. found -54 bp best)
2. **Chromatin**: Open region near TSS - accessible for dCas9-VPR binding
3. **Coherence**: 0.953 on hardware - exceeds e⁻² threshold by 7x
4. **Validated**: Confirmed on real IBM quantum hardware

---

## Comparison with Literature

### Moses et al. 2019 (PMID: 30654190)

| Parameter | Moses et al. | PhaseLab E210 |
|-----------|--------------|---------------|
| Best position | -54 bp | -90 bp |
| Activation | 2.27-fold | TBD (experimental) |
| CRISPRa window | -300 to -50 bp | -300 to -50 bp |
| Cell lines | SK-MEL-28, SUM159 | Same recommended |

### Comparison with Other PhaseLab Targets

| Experiment | Target | Disease | Hardware R̄ | Status |
|------------|--------|---------|-------------|--------|
| E200 | RAI1 | Smith-Magenis | 0.839 | GO |
| E201 | SCN2A | Autism | 0.970 | GO |
| **E210** | **PTEN** | **Cancer** | **0.953** | **GO** |

---

## Biological Context

### What PTEN Does

```
Normal: PTEN active → PI3K suppressed → Controlled growth
Cancer: PTEN silenced → PI3K hyperactive → Uncontrolled growth
CRISPRa: PTEN reactivated → PI3K suppressed → Growth control partially restored
```

### Expected Downstream Effects

When PTEN is reactivated by CRISPRa:

1. **Signaling**: ↓ p-AKT, ↓ p-mTOR, ↓ p-S6K, ↓ p-ERK
2. **Phenotype**: ↓ Proliferation, ↓ Colony formation, ↑ Apoptosis
3. **Drug Response**: ↑ Sensitivity to BRAF/PI3K/mTOR inhibitors

### Cancer Types That Could Benefit

| Cancer | PTEN Loss | CRISPRa Potential |
|--------|-----------|-------------------|
| Glioblastoma | ~70% | High |
| Prostate | ~60% | High |
| TNBC | ~40% | Medium-High |
| Melanoma | ~30% | Medium |
| Endometrial | ~50% | High |

---

## Technical Notes

### High GC Content and CpG Island Context

The PTEN core promoter resides in a CpG-island context with unusually high GC content, a known feature of many tumor suppressor promoters. Because of this, several top candidates have GC fractions around ~80%. For PTEN we explicitly relaxed the default 40-70% GC filter used in PhaseLab, and relied on thermodynamic ΔG and IR coherence to down-select guides that remain physically plausible despite the GC-rich background.

- CpG islands are common in tumor suppressors
- High GC = strong secondary structure, but still accessible
- All hardware-tested guides achieved GO status

### Hardware Validation Details

- **Backend**: ibm_torino
- **Qubits used**: 12
- **Circuit depth**: ~20 gates
- **Shots**: 4096
- **Job completion**: Successful

---

## Limitations

- **Computational + quantum only**: These results are purely computational plus quantum-hardware validation; no wet-lab testing has been performed yet.
- **Off-target scoring is internal**: PhaseLab's scoring does not replace genome-wide off-target enumeration via tools like CRISPOR, which is still required for any experimental work.
- **Does not model**:
  - Promoter methylation dynamics directly
  - Chromatin remodeling in actual tumors
  - Intratumor heterogeneity and clonal selection
  - Delivery efficiency, immune response, or toxicity

---

## Recommended Next Steps

### Immediate (Computational)

1. **CRISPOR validation**: Submit top 5 guides for genome-wide off-target analysis
2. **Literature comparison**: Order Moses et al. sgRNA-54 as positive control
3. **Chromatin profiling**: Check ENCODE data for SK-MEL-28/SUM159 specifically

### Experimental (Wet Lab)

1. **Clone guides** into dCas9-VPR lentiviral vector
2. **Transduce** SK-MEL-28 and SUM159 cells
3. **Readouts** at 48-72h:
   - qPCR: PTEN mRNA
   - Western: PTEN protein, p-AKT, p-mTOR
   - Phenotype: Proliferation, colony formation

### Combination Therapy

1. **Test**: CRISPRa + dabrafenib (BRAF inhibitor) in SK-MEL-28
2. **Test**: CRISPRa + dactolisib (PI3K/mTOR inhibitor) in both lines
3. **Expected**: Synergistic tumor suppression

---

## Files Generated

| File | Description |
|------|-------------|
| `Data/E210_pten_guides_20251211_173147.json` | Full results with hardware validation |
| `Docs/EXPERIMENT.md` | Complete experimental documentation |
| `Docs/FINDINGS.md` | This summary document |
| `Code/E210_pten_crispra_cancer.py` | Experiment script |

---

## Conclusion

E210 demonstrates that PhaseLab can design and pre-validate CRISPRa guides for tumor suppressor reactivation in cancer models. The PTEN guides achieved excellent coherence scores on IBM Quantum hardware, matching or exceeding previous results for RAI1 (Smith-Magenis) and SCN2A (Autism).

This opens the door for:
- Systematic CRISPRa design for other tumor suppressors (CDKN2A, RB1, TP53)
- Combination approaches with targeted cancer therapies
- Experimental validation toward potential clinical translation for PTEN-deficient cancers

---

*Dylan Vaca, December 2025*
*Hardware validated on IBM Torino*
