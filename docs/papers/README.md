# PhaseLab Research Papers

Three publishable academic papers establishing PhaseLab and demonstrating its applications.

---

## Paper Overview

| # | Title | Focus | Target Journals |
|---|-------|-------|-----------------|
| **1** | [PhaseLab Framework](PAPER_1_PHASELAB_FRAMEWORK.md) | Mathematical foundation, cross-domain validation | *Nature Computational Science*, *NPJ Quantum Information*, *PLOS Computational Biology* |
| **2** | [CRISPRa gRNA Design](PAPER_2_CRISPRA_GRNA_DESIGN.md) | Quantum-informed guide selection for RAI1 | *Nature Biotechnology*, *Nucleic Acids Research*, *CRISPR Journal* |
| **3** | [Circadian Modeling](PAPER_3_CIRCADIAN_MODELING.md) | Phase dynamics of SMS clock dysregulation | *Cell Systems*, *eLife*, *Journal of Biological Rhythms* |

---

## Key Contributions

### Paper 1: Framework Foundation
- Introduces IR coherence metric: $\bar{R} = e^{-V_\phi/2}$
- Validates $e^{-2}$ threshold across quantum and biological systems
- Demonstrates 1-2% hardware-simulator agreement
- Establishes PhaseLab as domain-general tool

### Paper 2: CRISPRa Application
- First quantum-informed CRISPR guide design
- Identifies optimal RAI1 guide: `TACAGGAGCTTCCAGCGTCA`
- MIT 83, CFD 93, zero off-targets ≤2mm
- Shows quantum energy correlates with specificity

### Paper 3: Systems Biology
- Extended Kuramoto model with PER delay, REV/ROR modulation
- Predicts non-monotonic RAI1 response (80% optimal, not 100%)
- Identifies therapeutic window: +20-60% boost
- Explains inverted melatonin in SMS

---

## Figures to Generate

### Paper 1
1. PhaseLab architecture diagram
2. IR coherence function curve
3. gRNA Hamiltonian encoding
4. Hardware vs simulator scatter plot
5. Kuramoto network visualization
6. Cross-domain coherence comparison

### Paper 2
1. RAI1 promoter schematic
2. PAM site distribution
3. Scoring heatmaps (ΔG, MIT/CFD)
4. Coherence vs specificity correlation
5. Hardware vs simulator agreement
6. CRISPOR off-target profile
7. Pipeline flowchart

### Paper 3
1. Circadian network diagram
2. PER/ROR/REV dynamics
3. R̄ vs RAI1 dosage curve
4. Phase synchrony heatmap
5. SMS vs normal phase trajectories
6. Therapeutic window analysis
7. Melatonin rhythm prediction
8. Parameter sensitivity analysis

---

## Hardware Validation Data

| System | Backend | Job ID | R̄ |
|--------|---------|--------|-----|
| H₂ VQE | IBM Brisbane | (prior work) | 0.891 |
| gRNA_1 | IBM Torino | d4suvd7t3pms7399mq8g | 0.854 |
| gRNA_2 | IBM Torino | d4suvd7t3pms7399mq8g | 0.840 |
| gRNA_3 | IBM Torino | d4suvd7t3pms7399mq8g | 0.839 |

---

## CRISPOR Validation

- **Batch ID:** QkJQAaLNZB1vU53T54MY
- **URL:** https://crispor.gi.ucsc.edu/crispor.py?batchId=QkJQAaLNZB1vU53T54MY
- **Best candidate:** TACAGGAGCTTCCAGCGTCA (MIT 83, CFD 93)

---

## Citation

If using these papers or PhaseLab in research:

```bibtex
@software{phaselab2025,
  author = {Vaca, Dylan},
  title = {PhaseLab: Phase-coherence analysis for quantum and biological systems},
  year = {2025},
  url = {https://github.com/followthesapper/phaselab}
}
```

---

## Status

| Paper | Draft | Figures | Review | Submission |
|-------|-------|---------|--------|------------|
| 1 | ✅ Complete | 🔲 Pending | 🔲 | 🔲 |
| 2 | ✅ Complete | 🔲 Pending | 🔲 | 🔲 |
| 3 | ✅ Complete | 🔲 Pending | 🔲 | 🔲 |

---

*Dylan Vaca, December 2025*
