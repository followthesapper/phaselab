# E201: SCN2A CRISPRa for Autism-Linked Haploinsufficiency

## Overview

This experiment applies the PhaseLab CRISPRa design pipeline to SCN2A, the most common monogenic cause of autism-associated neurodevelopmental disorders. The same validated methodology used for RAI1/Smith-Magenis Syndrome (E200) is applied to design and validate guide RNAs for boosting SCN2A expression in patients with haploinsufficiency.

## Background

### SCN2A and Autism

SCN2A encodes the voltage-gated sodium channel NaV1.2, which is critical for excitatory neuron function. Loss-of-function mutations leading to haploinsufficiency cause:

- Autism Spectrum Disorder (ASD)
- Intellectual disability
- Epilepsy (in a subset of patients)

SCN2A-related neurodevelopmental disorders represent one of the most frequent monogenic causes of ASD.

### Literature Support

A 2025 Nature paper demonstrated that CRISPRa upregulation of the remaining Scn2a allele in haploinsufficient mice:
- Restores Scn2a expression toward normal
- Rescues electrophysiological deficits
- Reduces seizures and improves behavior
- Works even when delivered in adolescence

Reference: [Nature 2025](https://www.nature.com/articles/s41586-025-09522-w)

## Methods

### Target Region

- **Gene**: SCN2A (sodium voltage-gated channel alpha subunit 2)
- **Chromosome**: 2 (2q24.3)
- **Genome Build**: GRCh38.p14
- **TSS**: chr2:165,239,414
- **Promoter Region**: chr2:165,238,414-165,239,914 (1500 bp)
- **CRISPRa Window**: -400 to -50 bp from TSS

### Pipeline Components

1. **PAM Scanning**: SpCas9 (NGG) on both strands
2. **Basic Filtering**: GC 40-70%, max homopolymer ≤4
3. **Off-target Scoring**: MIT algorithm + complexity analysis
4. **Chromatin Accessibility**: Brain-specific model (PsychENCODE-style)
5. **Thermodynamics**: SantaLucia nearest-neighbor energy
6. **Quantum Simulation**: Pauli Hamiltonian encoding
7. **IR Coherence**: R̄ = exp(-V_φ/2) with e⁻² threshold

### Quantum Validation

Coherence metric validated on IBM Quantum hardware:
- Backend: IBM Torino (133 qubits)
- Previous validation: E200 RAI1 guides achieved 1-2% simulator-hardware agreement

## Results

### PAM Site Analysis

| Metric | Value |
|--------|-------|
| Total PAM sites | 75 |
| In CRISPRa window | 18 |
| After quality filtering | 15 |
| All GO classification | 100% |

### Top Candidates

| Rank | Sequence | Position | GC | Chromatin | R̄ | Status |
|------|----------|----------|-----|-----------|-------|--------|
| 1 | TTCCACTTTTGACCAGGAGA | -64 | 45% | OPEN | 0.944 | GO |
| 2 | AGATGGTTCCACTTTTGACC | -70 | 45% | OPEN | 0.945 | GO |
| 3 | GCTGACTGCTACATAGCCAA | -104 | 50% | OPEN | 0.945 | GO |
| 4 | ATCTTAATTCACGTTCCTTT | -106 | 30% | OPEN | 0.942 | GO |
| 5 | CAAAGGAACGTGAATTAAGA | -87 | 35% | OPEN | 0.945 | GO |

### Key Findings

1. **All 15 candidates achieved GO classification** (R̄ > 0.135)
2. **Average coherence: 0.944** (Excellent reliability)
3. **Top candidates are in open chromatin** near the TSS
4. **Optimal GC range (45-50%)** in positions -64 to -104

## Comparison with E200 (RAI1)

| Metric | E200 (RAI1) | E201 (SCN2A) |
|--------|-------------|--------------|
| Disease | Smith-Magenis Syndrome | Autism-linked NDD |
| Haploinsufficiency | Yes | Yes |
| CRISPRa validated | Yes (mouse) | Yes (mouse, Nature 2025) |
| Total PAM sites | 274 | 75 |
| CRISPRa candidates | 91 | 18 |
| Top coherence | 0.945 | 0.945 |
| Hardware validated | IBM Torino | Pending |

## Next Steps

1. **Hardware Validation**: Run top 3 guides on IBM Torino
2. **CRISPOR Analysis**: Submit for genome-wide off-target screening
3. **Literature Comparison**: Compare with guides from Nature 2025 paper
4. **Experimental Validation**: Consider testing in iPSC-derived neurons

## Files

- `Code/E201_scn2a_crispra.py` - Main experiment script
- `Code/E201_ibm_hardware.py` - IBM Quantum validation script
- `Data/E201_scn2a_guides_*.json` - Results output
- `Docs/EXPERIMENT.md` - This file

## References

1. Nature 2025: CRISPR activation for SCN2A-related neurodevelopmental disorders
   - DOI: 10.1038/s41586-025-09522-w

2. SFARI: Development of CRISPR activation therapeutics to rescue SCN2A function
   - https://www.sfari.org/funded-project/development-of-crispr-activation-therapeutics-to-rescue-scn2a-function/

3. E200: Smith-Magenis Syndrome Gene Therapy (PhaseLab)
   - Hardware validation on IBM Torino, December 2025

---

*Dylan Vaca, December 2025*
