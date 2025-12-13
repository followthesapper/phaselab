# PhaseLab Cancer CRISPRa Workflow

## Overview

This document describes how to use PhaseLab for CRISPRa-based reactivation of
tumor suppressor genes in cancer cells. The workflow combines computational
guide design with IR coherence validation to select reliable guides for
cancer gene therapy applications.

## Scientific Background

### CRISPRa in Cancer Therapy

CRISPRa (CRISPR activation) uses catalytically dead Cas9 (dCas9) fused to
transcriptional activators (e.g., VPR, VP64) to upregulate gene expression.
In cancer, this approach can:

1. **Reactivate tumor suppressors** (PTEN, CDKN2A, RB1, TP53)
2. **Restore checkpoint genes** (p21, p16)
3. **Enhance immune recognition** (HLA, PD-L1 downregulators)

### Key Advantages Over Gene Addition

- Uses endogenous gene with native regulation
- Restores physiological expression levels
- Avoids insertional mutagenesis
- Can target epigenetically silenced genes

## Workflow Steps

### Step 1: Target Selection

Choose cancer-relevant tumor suppressors based on:
- Frequency of loss in target cancer type
- Mechanism of silencing (transcriptional vs. mutational)
- Therapeutic relevance (drug sensitization potential)

**Example targets for different cancers**:

| Cancer Type | Target Gene | Loss Frequency |
|-------------|-------------|----------------|
| Melanoma | PTEN | ~30% |
| TNBC | PTEN | ~40% |
| Glioblastoma | PTEN | ~70% |
| Prostate | PTEN | ~60% |
| Multiple | CDKN2A (p16) | ~50% |
| Multiple | RB1 | ~30% |

### Step 2: Configure Target in PhaseLab

Create a YAML configuration file in `src/phaselab/targets/`:

```yaml
# GENE_NAME.yaml
gene_symbol: "GENE_NAME"
genome_build: "GRCh38"
gene_id: "ENSGXXXXXXXXXXX"
refseq_id: "NM_XXXXXXXX"

promoter:
  chrom: "chrX"
  start: XXXXXXXX      # TSS - 1000 bp
  end:   XXXXXXXX      # TSS + 500 bp
  tss:   XXXXXXXX      # TSS position

crispr:
  system: "SpCas9"
  pam: "NGG"
  activator: "VPR"
  window:
    start_offset: -400  # CRISPRa optimal
    end_offset: -50

filters:
  gc_min: 0.40
  gc_max: 0.70
  max_homopolymer: 4

scoring:
  coherence_threshold: 0.135  # e^-2

notes:
  disease: "Cancer type"
  therapeutic_goal: "Description"
```

### Step 3: Run Guide Design Pipeline

```python
from phaselab.targets import load_target_config
from phaselab.crispr import design_guides, GuideDesignConfig

# Load target
target = load_target_config("PTEN")

# Design guides
guides = design_guides(
    sequence=promoter_sequence,
    tss_index=target.tss_index,
    config=GuideDesignConfig(
        cas_system="SpCas9",
        pam="NGG",
        mode="CRISPRa",
        window_start=-400,
        window_end=-50,
        gc_min=0.4,
        gc_max=0.7,
        require_open_chromatin=True,
    ),
)

# Rank by coherence and specificity
top_guides = guides.sort_values(
    ["go_no_go", "coherence_R", "mit_score"],
    ascending=[False, False, False]
).head(10)
```

### Step 4: Validate with IR Coherence

The pipeline automatically computes:
- **R̄ (coherence)**: Phase coherence metric
- **GO/NO-GO**: Binary classification (R̄ > e⁻² = 0.135)

Guides with GO status are recommended for experimental validation.

### Step 5: Off-Target Analysis

Export to CRISPOR for genome-wide off-target validation:

```python
from phaselab.io import export_crispor_batch

export_crispor_batch(
    sequence=promoter_sequence,
    outfile="guides_crispor.fa",
    guides=top_guides
)
```

Upload to https://crispor.tefor.net/ and verify:
- 0 off-targets at 0-1 mismatches
- Few off-targets at 2 mismatches
- No exonic off-targets

### Step 6: Experimental Validation

**Recommended protocol**:

1. **Cell lines**: Select cancer lines with:
   - Wild-type but silenced target gene
   - Confirmed low expression (qPCR/Western)
   - Established culture protocols

2. **Delivery system**: dCas9-VPR via:
   - Lentivirus (stable expression)
   - Transfection (transient screening)

3. **Controls**:
   - Untreated cells
   - dCas9-VPR + non-targeting gRNA
   - Literature-validated gRNA (if available)

4. **Readouts** (48-72h post-transduction):
   - Target mRNA (qPCR)
   - Target protein (Western blot)
   - Downstream signaling (pathway-specific)
   - Phenotype (proliferation, apoptosis, colony formation)

## Example: PTEN Reactivation

### Cell Models
- **SK-MEL-28**: BRAF V600E melanoma
- **SUM159**: Triple-negative breast cancer

### PhaseLab Results (E210)

Top guide: `CGGAAGGGGGAGCGCGGCAG`
- Position: -90 bp from TSS
- Coherence: R̄ = 0.944 [GO]

### Expected Outcomes (per Moses et al. 2019)
- 2-3 fold PTEN upregulation
- Reduced p-AKT, p-mTOR signaling
- Decreased colony formation
- Increased drug sensitivity

## Combination Therapy Applications

CRISPRa tumor suppressor reactivation synergizes with:

### PTEN + Targeted Therapy
- BRAF inhibitors (dabrafenib, vemurafenib)
- PI3K/mTOR inhibitors (dactolisib, BEZ235)
- AKT inhibitors (MK-2206)

### CDKN2A + CDK Inhibitors
- CDK4/6 inhibitors (palbociclib, ribociclib)

### RB1 + E2F Targeting
- E2F inhibitors
- HDAC inhibitors

## Delivery Considerations

### AAV Selection (using PhaseLab)

```python
from phaselab.delivery.aav import select_optimal_serotype

# For brain tumors (glioblastoma)
serotype = select_optimal_serotype(
    target_tissue="brain",
    species="human",
    require_bbb_penetration=True
)
# Returns: AAV9 or AAVrh10

# For liver metastases
serotype = select_optimal_serotype(
    target_tissue="liver",
    species="human"
)
# Returns: AAV8
```

### Packaging Constraints

dCas9-VPR cassette: ~4.2 kb
- Fits within AAV packaging limit (~4.7 kb)
- Use compact promoters (CMV mini, EF1α short)
- Consider dual-AAV split-intein systems for larger constructs

### Immunogenicity Assessment

```python
from phaselab.delivery.immunogenicity import assess_immunogenic_risk

risk = assess_immunogenic_risk(
    guide_sequence="CGGAAGGGGGAGCGCGGCAG",
    cas_protein="SpCas9",
    synthesis_method="chemical",
    population="general"
)
# Returns risk assessment with recommendations
```

## Reporting Results

Use PhaseLab's report generator:

```python
from phaselab.io.report import generate_therapeutic_report

report = generate_therapeutic_report(
    guides=top_guides,
    target="PTEN",
    disease="melanoma",
    include_coherence=True,
    include_delivery=True
)

# Export as HTML
export_report(report, "pten_crispra_report.html", format="html")
```

## References

1. Moses C, et al. "Activating PTEN Tumor Suppressor Expression with the
   CRISPR/dCas9 System." Mol Ther Nucleic Acids. 2019;14:287-300.

2. Gilbert LA, et al. "Genome-Scale CRISPR-Mediated Control of Gene
   Repression and Activation." Cell. 2014;159:647-661.

3. Konermann S, et al. "Genome-scale transcriptional activation by an
   engineered CRISPR-Cas9 complex." Nature. 2015;517:583-588.

4. Matharu N, et al. "CRISPR-mediated activation of a promoter or enhancer
   rescues obesity caused by haploinsufficiency." Science. 2019;363:eaau0629.
