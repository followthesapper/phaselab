# E210 FINDINGS: PTEN CRISPRa for Cancer

## Experiment Status: COMPLETE

**Date**: December 11, 2025
**Validation**: Simulator (AerSimulator) - IBM hardware fallback

## Key Results

### Guide Rankings

All 15 evaluated guides achieved **GO** status (R̄ > e⁻² = 0.135).

| Rank | Sequence | Position | GC% | R̄ | Status |
|------|----------|----------|-----|------|--------|
| 1 | CGGAAGGGGGAGCGCGGCAG | -90 | 80% | 0.944 | GO |
| 2 | GGAGCGGAGCGAGGAGGCGG | -265 | 80% | 0.944 | GO |
| 3 | GAGGCGGGACCCGCGTGCGG | -145 | 85% | 0.944 | GO |
| 4 | ACCCGCGTGCGGCGGAGGAG | -137 | 80% | 0.945 | GO |
| 5 | GAGCGAGGCGGAGCGGAGCG | -274 | 80% | 0.944 | GO |

### Top Recommendation

**Guide**: `CGGAAGGGGGAGCGCGGCAG`
- **Position**: -90 bp from TSS (close to Moses et al. -54)
- **Coherence**: R̄ = 0.944 [GO]
- **Chromatin**: OPEN (near TSS)
- **GC Content**: 80% (high due to CpG island)

### Comparison with Literature

Moses et al. 2019 found optimal activation at position **-54 bp** from TSS.
Our top guide at **-90 bp** is in the same high-activity region.

| Parameter | Moses et al. | PhaseLab |
|-----------|--------------|----------|
| Best position | -54 bp | -90 bp |
| Activation | 2.27-fold | TBD |
| Window | -300 to -50 | -300 to -50 |

### Chromatin Accessibility

The PTEN promoter is a **CpG island** with consistently high GC content (~80-90%).
All guides in the -300 to -50 bp window showed:
- Open chromatin (DNase HS signal)
- High accessibility scores (0.70-0.85)

### Validation Notes

**Warning flags observed**:
- High GC content (>80%) - inherent to PTEN CpG island promoter
- Some guides with Gx5 homopolymer runs - acceptable for CpG islands

**All guides passed**:
- IR coherence threshold (R̄ > 0.135)
- Sequence complexity requirements
- Chromatin accessibility criteria

## Biological Implications

### Therapeutic Potential

1. **Monotherapy**: PTEN reactivation alone may reduce tumor growth
2. **Combination**: CRISPRa + BRAF inhibitors (melanoma)
3. **Combination**: CRISPRa + PI3K/mTOR inhibitors (TNBC)

### Expected Phenotypes (per Moses et al.)

- Reduced p-AKT levels
- Reduced p-mTOR levels
- Reduced p-S6K levels
- Reduced p-ERK levels
- Decreased colony formation
- Increased drug sensitivity

## Recommended Next Steps

### Immediate
1. Submit top 5 guides to CRISPOR for genome-wide off-target analysis
2. Compare CRISPOR scores with PhaseLab rankings
3. Design experimental validation protocol

### Experimental Validation
1. **Cell lines**: SK-MEL-28, SUM159
2. **Vector**: dCas9-VPR lentivirus
3. **Readouts**:
   - qPCR (PTEN mRNA, 48-72h)
   - Western blot (PTEN protein, p-AKT, p-mTOR)
   - Proliferation assay (CellTiter-Glo)
   - Colony formation (2-3 weeks)

### Comparative Analysis
1. Order Moses et al. sgRNA-54 as positive control
2. Test PhaseLab guides vs literature guides
3. Assess IR coherence correlation with activation efficiency

## Files Generated

- `Data/E210_pten_guides_20251211_173147.json` - Full results
- `Docs/EXPERIMENT.md` - Experimental design
- `Docs/FINDINGS.md` - This document

## Conclusion

PhaseLab successfully designed CRISPRa guides for PTEN tumor suppressor
reactivation in cancer cells. All evaluated guides achieved GO status with
high coherence scores. The top recommendation targets the -90 bp position,
within the optimal CRISPRa activation window identified by Moses et al.

The high GC content (~80%) is inherent to the PTEN promoter's CpG island
nature and does not indicate poor guide quality. Experimental validation
in SK-MEL-28 and SUM159 cells is recommended to confirm activation efficiency.
