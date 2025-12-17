# Minimal RAI1 Tiling Experiment Specification

## Purpose

Design the minimum viable tiling experiment to validate PhaseLab's Stage I predictions for RAI1 CRISPRa activation in Smith-Magenis Syndrome.

---

## Experimental Design

### Target Information

| Parameter | Value |
|-----------|-------|
| Gene | RAI1 (Retinoic Acid Induced 1) |
| Chromosome | 17p11.2 |
| TSS | chr17:17,584,792 (GRCh38) |
| CRISPRa Window | TSS-400 to TSS+100 |
| Therapeutic Goal | 70-110% of normal RAI1 expression |
| SMS Baseline | ~50% (haploinsufficiency) |

### Guide Design

#### Total Guides: 20

| Category | Count | Spacing | Purpose |
|----------|-------|---------|---------|
| Core tiling | 12 | 40bp | Cover TSS-350 to TSS-50 |
| Flanking | 4 | 100bp | Cover TSS-400 to TSS-350, TSS-50 to TSS+100 |
| Positive controls | 2 | - | Chang et al. validated guides |
| Negative controls | 2 | - | Non-targeting, off-target guides |

#### Guide Positions (TSS-relative)

```
Core tiling (40bp spacing):
  -350, -310, -270, -230, -190, -150, -110, -70, -30, +10, +50, +90

Flanking (100bp spacing):
  -400, -380 (upstream)
  +100, +120 (downstream)

Controls:
  Chang sg2: TSS-80 (positive control)
  Chang sg4: TSS-120 (positive control)
  NT1: Non-targeting
  OT1: Off-TSS control
```

### Replicates

| Replicate Type | Count | Purpose |
|----------------|-------|---------|
| Biological | 3 | Independent transfections |
| Technical | 2 | Same transfection, separate wells |

**Total data points**: 20 guides × 3 biological × 2 technical = **120 measurements**

---

## Cell Model

### Primary Model: iPSC-Derived Neurons

| Parameter | Specification |
|-----------|--------------|
| Cell type | iPSC-derived cortical neurons |
| Genotype | SMS patient (17p11.2 deletion) or isogenic |
| Passage | P15-P25 |
| Culture | Neuronal maintenance medium |
| Density | 50,000 cells/well (96-well) |

### Alternative Model: HEK293T + Reporter

| Parameter | Specification |
|-----------|--------------|
| Cell type | HEK293T-RAI1-GFP |
| Reporter | RAI1 promoter driving destabilized GFP |
| Density | 20,000 cells/well (96-well) |
| Advantage | Faster, cheaper screening |

---

## CRISPRa System

### Recommended: dCas9-VP64-p65-Rta (VPR)

| Component | Vector | Dose |
|-----------|--------|------|
| dCas9-VPR | Lenti (stable line) | Integrated |
| sgRNA | Lenti or plasmid | MOI 0.5 or 200ng/well |

### Alternative: SAM System

| Component | Vector | Dose |
|-----------|--------|------|
| dCas9-VP64 | Lenti | Integrated |
| MS2-p65-HSF1 | Lenti | Integrated |
| sgRNA-MS2 | Plasmid | 200ng/well |

---

## Readouts

### Primary: RAI1 Expression

| Timepoint | Method | Metric |
|-----------|--------|--------|
| 48h | RT-qPCR | RAI1 mRNA (normalized to GAPDH) |
| 72h | RT-qPCR | RAI1 mRNA (sustained expression) |

**qPCR Primers**:
```
RAI1-F: 5'-AGCCTGAGCAAGAAGAAGCA-3'
RAI1-R: 5'-TCTGCTGCTGCTGCTGCTG-3'
GAPDH-F: 5'-GAAGGTGAAGGTCGGAGTC-3'
GAPDH-R: 5'-GAAGATGGTGATGGGATTTC-3'
```

### Secondary: Protein Level

| Timepoint | Method | Metric |
|-----------|--------|--------|
| 72h | Western blot | RAI1 protein (normalized to β-actin) |
| 72h | Immunofluorescence | RAI1 nuclear localization |

### Tertiary: Functional (if available)

| Timepoint | Method | Metric |
|-----------|--------|--------|
| 96h | PER2::Luc reporter | Circadian period/amplitude |

---

## Analysis Plan

### Stage 1: Quality Control

1. Exclude wells with transfection efficiency < 70%
2. Exclude outliers > 3 SD from replicate mean
3. Confirm positive controls show expected activation

### Stage 2: Spatial Coherence Analysis

```python
from phaselab.spatial import analyze_tiling_coherence

# Load data
landscape = ResponseLandscape(
    coords=guide_positions,
    responses=log2_fold_changes,
)

# Analyze
result = analyze_tiling_coherence(
    landscape,
    window=5,  # 5-guide window
    stable_threshold=0.7,
)

# Extract stable regions
stable_regions = result.stable_regions
```

### Stage 3: Validation Metrics

| Metric | Target | Pass Criteria |
|--------|--------|---------------|
| Coherence-variance correlation | r < -0.2 | Validates framework |
| Stage I prediction accuracy | > 60% | Stage I is useful |
| Variance reduction in stable regions | > 25% | Practical benefit |
| Top guide fold-change | > 1.5x | Biological effect |

---

## Success Criteria

### GO Criteria

1. **Primary**: At least 3 guides achieve > 1.5x RAI1 activation
2. **Secondary**: Stable regions identified with coherence > 0.7
3. **Tertiary**: Stage I predictions match > 60% of empirical stable regions

### NO-GO Criteria

1. No guides achieve > 1.3x activation
2. All regions show amplifying behavior (coherence < 0.4)
3. Positive controls fail to replicate published effects

---

## Timeline

| Week | Activity |
|------|----------|
| 1 | Guide cloning, cell prep |
| 2 | Transfection, 48h timepoint |
| 3 | 72h timepoint, qPCR, Western |
| 4 | Data analysis, coherence computation |
| 5 | Validation, report generation |

---

## Budget Estimate

| Item | Cost |
|------|------|
| Guide synthesis (20 guides) | $400 |
| Cloning/sequencing | $300 |
| Cell culture (96-well × 6 plates) | $200 |
| qPCR reagents | $500 |
| Western blot | $300 |
| Personnel (40h @ $50/h) | $2,000 |
| **Total** | **$3,700** |

---

## Expected Outcomes

### If Successful

1. Empirical validation of PhaseLab Stage I predictions
2. Identification of stable regions in RAI1 promoter
3. Top guide candidates for therapeutic development
4. Publication-ready validation data

### If Partial Success

1. Coherence framework validated but Stage I accuracy limited
2. Need larger tiling for better resolution
3. Alternative CRISPRa strategies indicated

### If Failure

1. RAI1 promoter may be uniformly amplifying
2. CRISPRa may be insufficient for therapeutic dosage
3. Consider multi-guide strategies or alternative modalities

---

## Contingency Plans

### If Coherence Not Validated

- Increase replicate depth (n=5 biological)
- Reduce window size (window=3)
- Check for batch effects

### If No Guides Achieve Threshold

- Try SAM system instead of VPR
- Design guides for enhancer regions
- Consider CRISPRa + small molecule combination

---

## Contact

For questions about this experiment specification, contact the PhaseLab development team.

---

*Version: 1.0*
*Last updated: December 2025*
