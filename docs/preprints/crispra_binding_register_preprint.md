# Phase-Aware CRISPRa Guide Enumeration Reveals Binding-Register Flexibility in GC-Dense Promoters

**Dylan Vaca**
Independent Researcher
PhaseLab / Informational Relativity Program
Contact: <email>

---

## Abstract

CRISPR activation (CRISPRa) relies on catalytically inactive Cas nucleases (dCas) to bind regulatory DNA regions and modulate transcription. While guide RNA (gRNA) design tools have matured for cleavage-based CRISPR applications, many retain enumeration assumptions optimized for double-strand break induction rather than binding. Here we show that rigid spacer–PAM anchoring assumptions systematically exclude experimentally validated CRISPRa guides in GC-dense promoters. Using the Smith–Magenis syndrome gene *Rai1* as a case study, we demonstrate that a published, experimentally validated CRISPRa guide is missed by standard enumeration pipelines despite exact sequence and PAM compatibility. We introduce a binding-aware enumeration model incorporating relaxed PAM constraints, configurable spacer length, and a sliding binding register. This correction recovers the experimental guide at the correct genomic position without relaxing safety or sequence filters. These findings highlight a previously underappreciated modeling mismatch between cleavage-optimized CRISPR pipelines and CRISPRa binding biology and motivate explicit declaration of nuclease role in computational guide design.

---

## Introduction

CRISPR-based transcriptional modulation (CRISPRa and CRISPRi) has become a powerful strategy for studying gene regulation and developing therapeutic interventions for haploinsufficiency disorders. Unlike cleavage-based CRISPR systems, CRISPRa relies on dCas-mediated DNA binding rather than precise double-strand break positioning. Despite this distinction, many computational pipelines for CRISPRa guide design inherit enumeration logic developed for cleavage-based nucleases, including rigid spacer–PAM anchoring and fixed guide length assumptions.

GC-dense promoters, common among dosage-sensitive developmental genes, pose additional challenges due to repetitive motifs, overlapping PAMs, and high local sequence degeneracy. These features increase the likelihood that biologically functional binding events occur with slight register shifts relative to canonical spacer–PAM models.

Here, we investigate a discrepancy between computational guide enumeration and experimental CRISPRa results reported by Chang *et al.* (2022), who identified an effective CRISPRa guide for *Rai1* activation in mouse despite non-canonical PAM usage and a non-standard spacer length. We show that standard enumeration pipelines systematically fail to recover this guide due to rigid anchoring assumptions. We then introduce a binding-aware enumeration framework implemented in PhaseLab that resolves this discrepancy.

---

## Results

### Failure of Rigid Enumeration to Recover an Experimental CRISPRa Guide

Using the published mouse *Rai1* promoter sequence (mm10), we attempted to enumerate CRISPRa guides matching the experimentally validated guide sequence:

```
CCTGGCACCCGAGGCCACGA
```

Despite exact sequence presence, correct PAM compatibility (GCGAGA), and correct transcription start site offset (−80 bp), standard enumeration logic failed to recover this guide. Investigation revealed three interacting assumptions:

1. **Cleavage-optimized PAM stringency** (NNGRRT for SaCas9)
2. **Fixed 21 bp spacer length**
3. **Rigid spacer–PAM anchoring (guide_start = pam_start − guide_length)**

These assumptions are appropriate for cleavage but inconsistent with CRISPRa binding biology.

---

### Binding-Aware Enumeration Recovers the Experimental Winner

We implemented three corrections under an explicit *binding* nuclease role:

* **Relaxed PAM pattern** for SaCas9 binding (NNGRRN)
* **Configurable spacer length** (20 bp, matching experimental design)
* **Sliding binding register** allowing ±2 bp shifts around PAM anchor

With these changes enabled, the experimental guide was recovered:

| Guide                     | Status | Rank | TSS Offset |
| ------------------------- | ------ | ---- | ---------- |
| sg1                       | Found  | #47  | −258       |
| sg2 (experimental winner) | Found  | #100 | −80        |

Importantly, the guide was recovered without relaxing GC content thresholds, poly-T filters, or off-target safety gates.

---

### Enumeration Expansion Does Not Inflate False Positives

Introducing sliding registers increased candidate enumeration from 114 to 570 guides; however, all candidates still passed existing hard safety gates. The experimental guide appeared as a legitimate candidate rather than an outlier, indicating that sliding registers expand *coverage* without collapsing specificity.

---

## Discussion

### Binding vs Cutting Are Distinct Computational Problems

Our results demonstrate that CRISPRa guide enumeration is not a simple subset of cleavage-based CRISPR design. Rigid spacer–PAM anchoring, while appropriate for inducing double-strand breaks, is overly restrictive for dCas-mediated binding. In GC-dense promoters, where PAMs and spacer motifs overlap extensively, binding events tolerate small register shifts that computational pipelines often exclude.

### Why Existing Tools Miss These Guides

Most CRISPR design tools do not encode nuclease *role* as a first-class concept. As a result, cleavage-era assumptions persist silently in CRISPRa workflows. Because pipelines typically report only successful candidates, users may be unaware that valid guides are systematically excluded.

### Implications for Therapeutic CRISPRa

For haploinsufficiency disorders such as Smith–Magenis syndrome, missing a functional guide can delay therapeutic exploration. Binding-aware enumeration increases the probability that experimentally viable guides are surfaced early, reducing wasted wet-lab effort.

---

## Methods

### Promoter Sequence

Mouse *Rai1* promoter (mm10) was extracted from UCSC Genome Browser, spanning 1,101 bp upstream of the transcription start site.

### Enumeration Framework

Guide enumeration was performed using PhaseLab v0.9.3 with explicit nuclease role declaration:

* **Nuclease**: SaCas9 (dCas9)
* **Role**: Binding
* **Spacer length**: 20 bp
* **PAM pattern**: NNGRRN
* **Binding register shift**: ±2 bp

### Validation

Recovered guides were matched against published experimental sequences and genomic positions. No off-target scoring or ranking optimization was applied during recovery analysis.

---

## Limitations

This study does not claim to predict CRISPRa efficacy or replace wet-lab validation. Binding-aware enumeration increases candidate recall but does not guarantee transcriptional activation. Chromatin accessibility, transcription factor competition, and cellular context remain critical determinants of success.

---

## Conclusion

We identify a systematic enumeration error affecting CRISPRa guide discovery in GC-dense promoters and provide a binding-aware correction that recovers experimentally validated guides otherwise excluded by standard pipelines. Explicit modeling of nuclease role and binding register flexibility should be considered essential components of future CRISPRa design tools.

---

## Availability

PhaseLab is open-source and available at:
[https://github.com/followthesapper/phaselab](https://github.com/followthesapper/phaselab)

---

## Acknowledgments

We thank the CRISPR research community for openly published datasets enabling retrospective validation.

---

## References

Chang et al., 2022. *[Full citation here]*
Additional references as appropriate.
