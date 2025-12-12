Changelog
=========

Version 0.6.1 (December 2025)
-----------------------------

**Coherence Mode Parameter**

- Added ``mode="heuristic"`` (fast) vs ``mode="quantum"`` (VQE) parameter
- Heuristic mode is now default for speed
- Quantum mode provides research-grade accuracy

**Honest Coherence Weighting**

- Heuristic coherence demoted to tie-breaker weight (0.05 vs 0.30)
- Prevents over-reliance on proxy metrics
- Two-stage scoring: hard safety gates + soft ranking

**Risk Mass Metrics**

- Added ``risk_mass_close``: Off-targets within 100bp of TSS
- Added ``risk_mass_exonic``: Off-targets in exonic regions
- Added ``tail_risk_score``: Aggregate tail risk metric

**Evidence Levels**

- Level A: Hardware-validated (IBM Quantum)
- Level B: VQE-simulated (quantum mode)
- Level C: Heuristic only (capped influence)

**Score Capping**

- Unvalidated guides capped to prevent misleading rankings
- Evidence level affects maximum achievable score

Version 0.6.0 (November 2025)
-----------------------------

**ATLAS-Q Integration**

- Full integration with ATLAS-Q tensor network simulator
- IR measurement grouping (5x variance reduction)
- Real circular statistics coherence (replaces heuristic)

**Coherence-Aware VQE**

- VQE optimization with real-time coherence tracking
- GO/NO-GO classification during optimization
- Optional GPU acceleration via Triton kernels

**Rust Backend Support**

- Optional Rust backend for 30-77x faster simulation
- Automatic fallback to Python if Rust unavailable

**CRISPOR Integration**

- IR-enhanced off-target analysis
- Off-target entropy metrics
- Coherence contrast scoring

**Unified ATLAS-Q Coherence**

- Single coherence computation path for all CRISPR modalities
- Consistent API across CRISPRa, CRISPRi, knockout, editing

Version 0.5.0 (October 2025)
----------------------------

**Real ATAC-seq Integration**

- BigWig file support for tissue-specific accessibility
- CpG methylation modeling for CRISPRa efficiency

**Nucleosome Occupancy**

- NuPoP-like algorithm for nucleosome prediction
- Integration with guide scoring

**Multi-Guide Synergy**

- Combinatorial CRISPR design
- Pairwise synergy prediction
- Guide set optimization

**Enhancer Targeting**

- Enhancer identification and scoring
- CRISPRa enhancer guide design
- Promoter vs enhancer comparison

**AAV Delivery Modeling**

- Serotype selection
- Delivery efficiency prediction
- Immunogenicity assessment

**Validation Reports**

- Comprehensive validation report generation
- Evidence summary and confidence scoring

Version 0.4.0 (September 2025)
------------------------------

**Complete CRISPR Toolkit**

- CRISPR knockout (Cas9 cutting)
- CRISPRi (transcriptional repression)
- All modalities hardware-validated on IBM Torino

**Therapeutic Dosage Optimization**

- Haploinsufficiency models
- Dose-response prediction
- Therapeutic window estimation

Version 0.3.0 (August 2025)
---------------------------

**Multi-Tissue Circadian Models**

- Inter-tissue coupling
- Tissue-specific parameters
- SCN-peripheral synchronization

**Drug Response Modeling**

- Chronotherapy optimization
- Dosing schedule prediction
- Response curve modeling

**Expanded CRISPR Editors**

- Base editing (ABE/CBE)
- Prime editing (pegRNA design)
- Bystander prediction

Version 0.2.0 (July 2025)
-------------------------

**Protein Folding Coherence**

- Folding reliability assessment
- Coherence-structure correlation
- Engineering applications

**Tissue-Specific Chromatin**

- ENCODE integration
- Cell-type specific accessibility
- Tissue-aware guide scoring

Version 0.1.0 (June 2025)
-------------------------

**Initial Release**

- Core coherence metrics (R, V_phi)
- GO/NO-GO classification
- Basic CRISPRa guide design
- SMS circadian clock model
- IBM Quantum integration
