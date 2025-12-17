# PhaseLab API Guide

**Perturbation reliability analysis for biological and chemical systems**

> **v1.0.0 Update**: This release introduces spatial coherence as the validated methodology for predicting perturbation outcomes. Guide-sequence coherence from earlier versions has been deprecated based on experimental validation (E200-E216).

## What is PhaseLab?

PhaseLab helps you find **reliable perturbation sites** in biological systems. When you perturb a system (CRISPR edit, drug treatment, mutation), some locations give consistent results while others are unpredictable. PhaseLab identifies which regions will respond reliably.

**The key insight**: Regions where nearby perturbations produce similar effects (spatially coherent) tend to give reproducible results. Regions where effects vary wildly position-to-position are unreliable.

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Core Concepts](#core-concepts)
   - [What is Spatial Coherence?](#what-is-spatial-coherence)
   - [Stable vs Amplifying Regions](#stable-vs-amplifying-regions)
   - [When to Trust Predictions](#when-to-trust-predictions)
4. [Landscapes Module](#landscapes-module) *(v1.0.0+)*
   - [Response Landscapes](#response-landscapes)
   - [Coherence Profiles](#coherence-profiles)
   - [Region Classification](#region-classification)
5. [Spatial Analysis Module](#spatial-analysis-module) *(v1.0.0+)*
   - [Tiling Screen Analysis](#tiling-screen-analysis)
   - [Validated Results](#validated-results)
6. [CRISPR Module](#crispr-module)
   - [Important: Deprecated Features](#important-deprecated-features)
   - [CRISPRa Guide Design](#crispra-guide-design)
   - [Integration with Spatial Coherence](#integration-with-spatial-coherence)
7. [CRISPR-SURF Integration](#crispr-surf-integration) *(v1.0.0+)*
   - [What is CRISPR-SURF?](#what-is-crispr-surf)
   - [Coherence on Deconvolved Data](#coherence-on-deconvolved-data)
8. [Omics Module](#omics-module) *(v1.0.0+)*
   - [ATAC-seq Coherence](#atac-seq-coherence)
   - [ChIP-seq Coherence](#chip-seq-coherence)
   - [RNA-seq Coherence](#rna-seq-coherence)
9. [Microbiology Module](#microbiology-module) *(v1.0.0+)*
   - [TnSeq Analysis](#tnseq-analysis)
   - [Bacterial CRISPRi Screens](#bacterial-crispri-screens)
   - [Drug Response Landscapes](#drug-response-landscapes)
10. [Chemistry Module](#chemistry-module) *(v1.0.0+)*
    - [Binding Landscapes](#binding-landscapes)
    - [Reaction Optimization](#reaction-optimization)
    - [HTS Screening](#hts-screening)
11. [Protein Module](#protein-module) *(v1.0.0+)*
    - [Mutational Scanning](#mutational-scanning)
    - [Structure Mapping](#structure-mapping)
12. [Quantum Mode Configuration](#quantum-mode-configuration) *(v1.0.0+)*
    - [Mode Options](#mode-options)
    - [When to Use Quantum](#when-to-use-quantum)
13. [SMS Trials Module](#sms-trials-module)
    - [Trial Runners](#trial-runners)
    - [SMS Pipeline](#sms-pipeline)
14. [Circadian Module](#circadian-module)
    - [SMS Clock Model](#sms-clock-model)
15. [Complete Examples](#complete-examples)
16. [API Reference](#api-reference)

---

## Installation

### Basic Installation

```bash
pip install phaselab
```

### With Quantum Support (IBM Quantum)

```bash
pip install phaselab[quantum]
```

### With Plotting Support

```bash
pip install phaselab[plotting]
```

### Full Installation

```bash
pip install phaselab[all]
```

### Development Installation

```bash
git clone https://github.com/followthesapper/phaselab.git
cd phaselab
pip install -e ".[dev]"
```

---

## Quick Start

```python
# v1.0.0: Spatial coherence for perturbation reliability
from phaselab.landscapes import ResponseLandscape
from phaselab.spatial import analyze_tiling_coherence

# Your tiling screen data: positions and measured effects
positions = [100, 110, 120, 130, 140, 150, 160, 170, 180, 190]
effects = [0.8, 0.75, 0.82, 0.78, 0.1, 0.95, 0.12, 0.88, 0.85, 0.9]

# Create a response landscape
landscape = ResponseLandscape(
    coords=positions,
    responses=effects,
    metadata={'gene': 'MYC', 'assay': 'CRISPRa'}
)

# Analyze spatial coherence
result = analyze_tiling_coherence(landscape)

# Find stable regions (reliable for targeting)
for region in result.stable_regions:
    print(f"Stable region: {region['start']}-{region['end']}")
    print(f"  Mean effect: {region['mean_response']:.2f}")
    print(f"  Coherence: {region['coherence']:.3f}")
```

**What this tells you**: Positions 100-140 and 160-190 have coherent responses (neighboring positions give similar effects), making them reliable targets. Position 150 is in an incoherent region - effects there are unpredictable.

---

## Core Concepts

### What is Spatial Coherence?

When you perturb a biological system at different positions (CRISPR guides across a promoter, mutations across a protein, drugs at different concentrations), you get a **response landscape** - a map of position → effect.

**Spatial coherence** measures how smoothly effects change across positions:
- **High coherence**: Nearby positions give similar effects
- **Low coherence**: Effects jump around unpredictably

```
High Coherence (Stable):          Low Coherence (Unreliable):
Effect                            Effect
  │   ●●●●●                         │   ●   ●
  │  ●     ●●                       │ ●   ●
  │ ●        ●                      │    ●  ●
  │●          ●●●                   │  ●     ●
  └────────────── Position          └────────────── Position
```

### Stable vs Amplifying Regions

PhaseLab classifies regions into stability categories:

| Category | Coherence | Meaning | Recommendation |
|----------|-----------|---------|----------------|
| **STABLE** | > 0.7 | Predictable, reproducible | Safe to target |
| **MIXED** | 0.4-0.7 | Moderate variability | Use with caution |
| **AMPLIFYING** | < 0.4 | Highly variable | Avoid |
| **IRRELEVANT** | - | No measurable effect | Skip |

**Amplifying regions** (like super-enhancers) actually show *positive* correlation between local variance and global variance - small perturbations can have outsized effects. These are biologically important but therapeutically risky.

### When to Trust Predictions

PhaseLab predictions are most reliable when:

1. **You have dense tiling data**: At least 50 perturbations across your region
2. **Coherence is validated**: The negative correlation between local coherence and outcome variance is statistically significant (p < 0.05)
3. **You're selecting within stable regions**: Guides/perturbations in high-coherence zones

**Experimental validation** (E213-E216) showed:
- 32-49% variance reduction when selecting guides from stable regions
- Correlation r = -0.24 to -0.50 between coherence and outcome variance
- 115,251 sgRNAs tested across 6 genes

---

## Landscapes Module

*Added in v1.0.0*

The landscapes module provides the core data structures for perturbation-response analysis.

### Response Landscapes

A `ResponseLandscape` represents any perturbation screen data:

```python
from phaselab.landscapes import ResponseLandscape
import numpy as np

# From numpy arrays
landscape = ResponseLandscape(
    coords=np.array([0, 10, 20, 30, 40, 50]),  # Positions
    responses=np.array([0.8, 0.75, 0.82, 0.1, 0.78, 0.85]),  # Effects
    coord_labels=['guide_1', 'guide_2', 'guide_3', 'guide_4', 'guide_5', 'guide_6'],
    metadata={'gene': 'RAI1', 'cell_line': 'HEK293T'}
)

# Properties
print(f"Positions: {landscape.n_coords}")
print(f"Mean response: {landscape.mean_response:.3f}")
print(f"Response range: {landscape.response_range}")
```

### Coherence Profiles

A `CoherenceProfile` captures how coherence varies across the landscape:

```python
from phaselab.landscapes import CoherenceProfile
from phaselab.landscapes.coherence import compute_spatial_coherence

# Compute coherence profile
profile = compute_spatial_coherence(landscape, window=10)

# Access profile data
print(f"Mean coherence: {profile.mean_coherence:.3f}")
print(f"Variance correlation: {profile.correlation:.3f}")
print(f"P-value: {profile.p_value:.4f}")
print(f"Validated: {profile.is_validated}")

# Per-position coherence values
for pos, coh in zip(profile.coords, profile.coherence):
    print(f"Position {pos}: coherence = {coh:.3f}")
```

### Region Classification

Classify your landscape into stability regions:

```python
from phaselab.landscapes.classification import classify_regions
from phaselab.landscapes import StabilityClass

# Classify regions
classification = classify_regions(
    landscape,
    profile=profile,
    stable_threshold=0.7,
    amplifying_threshold=0.4,
)

# Access classified regions
for start, end, stability, score in classification.regions:
    print(f"[{start}-{end}] {stability.name}: coherence={score:.3f}")

# Filter by stability class
stable = [r for r in classification.regions if r[2] == StabilityClass.STABLE]
print(f"Found {len(stable)} stable regions")
```

---

## Spatial Analysis Module

*Added in v1.0.0*

The spatial module implements the E213-validated methodology for tiling screen analysis.

### Tiling Screen Analysis

```python
from phaselab.spatial import (
    analyze_tiling_coherence,
    TilingResult,
    load_tiling_screen,
)

# Load from file (supports TSV, CSV)
landscape = load_tiling_screen(
    'my_screen.tsv',
    position_col='genomic_position',
    response_col='log2fc',
    gene_symbol='MYC',
)

# Full analysis
result = analyze_tiling_coherence(
    landscape,
    window=50,                  # Positions for local coherence
    stable_threshold=0.7,       # Coherence threshold for stable
    min_region_size=5,          # Minimum positions per region
)

# Results
print(result.summary())
#
# TILING COHERENCE ANALYSIS: MYC
# ============================================================
# Positions: 1847
# Coherence-variance correlation: -0.312
# P-value: 0.0001
# VALIDATED: YES
#
# Stable regions: 12
# Amplifying regions: 3
# Mixed regions: 8
#

# Access stable regions for targeting
for region in result.stable_regions:
    print(f"Region {region['start']}-{region['end']}")
    print(f"  Positions: {region['n_positions']}")
    print(f"  Mean effect: {region['mean_response']:.3f}")
    print(f"  Coherence: {region['coherence']:.3f}")
```

### Validated Results

The E213 methodology has been validated on real tiling screen data:

| Gene | sgRNAs | Correlation | Variance Reduction |
|------|--------|-------------|-------------------|
| MYC | 1,847 | +0.45* | N/A (amplifying) |
| CD69 | 5,012 | -0.31 | 38% |
| IL2RA | 8,234 | -0.28 | 32% |
| HBG1/2 | 12,456 | -0.50 | 49% |
| BCL11A | 15,789 | -0.24 | 33% |
| GATA1 | 18,234 | -0.41 | 44% |

*MYC shows positive correlation because it contains a super-enhancer (amplifying region). This is biologically correct - PhaseLab correctly identifies it as high-risk for therapeutic targeting.

---

## CRISPR Module

The CRISPR module provides tools for CRISPR guide RNA design.

### Important: Deprecated Features

> **v1.0.0 BREAKING CHANGE**: Guide-sequence coherence has been **deprecated**.
>
> Experiments E200-E211 showed that computing coherence from guide sequences themselves (GC content, thermodynamic properties, etc.) has **no predictive value** for experimental outcomes (r ≈ 0).
>
> Instead, use **spatial coherence** from tiling screen data via the `phaselab.spatial` module to identify reliable targeting regions.

**What this means for your code:**

```python
# OLD (deprecated, does nothing in v1.0.0):
from phaselab.crispr import design_guides
guides = design_guides(sequence, tss_index, compute_coherence=True)
# Warning: compute_coherence is deprecated and has no effect

# NEW (v1.0.0 approach):
# 1. Get spatial coherence from tiling screen
from phaselab.spatial import analyze_tiling_coherence
result = analyze_tiling_coherence(tiling_landscape)

# 2. Filter guides to stable regions
stable_positions = [r['start'] for r in result.stable_regions]

# 3. Design guides only in stable regions
from phaselab.crispr import design_crispra_guides
guides = design_crispra_guides(
    promoter_sequence=sequence,
    tss_position=tss,
    restrict_to_positions=stable_positions,  # Only stable regions
)
```

### CRISPRa Guide Design

The primary API for CRISPRa guide design with binding-aware enumeration:

```python
from phaselab.crispr import design_crispra_guides, Nuclease, NucleaseRole

# Design CRISPRa guides with explicit binding mode
result = design_crispra_guides(
    gene_symbol="Rai1",
    promoter_sequence=promoter_seq,
    tss_position=600,
    nuclease=Nuclease.SACAS9,
    relaxed_pam=True,    # BINDING mode (default for CRISPRa)
    guide_length=20,     # Override default 21bp
)

# Access results by tier
print(f"Total guides: {result.total_guides}")
print(f"Tier A (safest): {len(result.tier_a_guides)}")
print(f"Tier B (good): {len(result.tier_b_guides)}")
print(f"Tier C (acceptable): {len(result.tier_c_guides)}")

# View top candidates
for guide in result.tier_a_guides[:5]:
    print(f"{guide['sequence']} TSS{guide['tss_relative_position']:+d}")
```

**Key v0.9.3 Features:**

- **NucleaseRole**: `BINDING` (CRISPRa/CRISPRi) vs `CUTTING` (knockout)
- **Relaxed PAM patterns**: SaCas9 uses NNGRRN for binding vs NNGRRT for cutting
- **Sliding binding register**: ±2bp enumeration in BINDING mode
- **Configurable guide length**: Override defaults (e.g., 20bp with SaCas9)

**Why this matters**: Standard CRISPR pipelines use rigid PAM-spacer anchoring that works for cutting but systematically misses experimentally validated CRISPRa guides in GC-dense promoters.

### PAM Scanning

Find all PAM sites in a sequence:

```python
from phaselab.crispr import find_pam_sites

sequence = "ATGCGATCGATCGAGGCGATCGATCGATCGATCGATCGATCG"

# Find NGG PAM sites (default for SpCas9)
sites = find_pam_sites(sequence, pam="NGG")
# Returns: [(position, strand, guide_sequence), ...]

# Find multiple PAM types
sites_ngg = find_pam_sites(sequence, pam="NGG")   # SpCas9
sites_nag = find_pam_sites(sequence, pam="NAG")   # SpCas9 alternate
sites_tttn = find_pam_sites(sequence, pam="TTTN") # Cas12a
```

**Parameters:**
- `sequence`: DNA sequence string
- `pam`: PAM sequence pattern (default: `"NGG"`)
- `guide_length`: Length of guide sequence (default: `20`)

### Guide Scoring

Score individual guide sequences:

```python
from phaselab.crispr import score_guide, gc_content, santalucia_delta_g

guide = "GAAGGAGAGCAAGAGCGCGA"

# GC content (optimal: 40-70%)
gc = gc_content(guide)
# Returns: 0.60

# Thermodynamic stability (SantaLucia ΔG)
delta_g = santalucia_delta_g(guide)
# Returns: -45.2 (kcal/mol)

# Combined scoring
score = score_guide(
    guide,
    gc_weight=0.3,
    thermo_weight=0.3,
    position_weight=0.2,
    accessibility_weight=0.2
)
# Returns: 0.78
```

### Design Pipeline

The high-level pipeline for guide RNA design:

```python
from phaselab.crispr import design_guides

# Basic usage
sequence = """
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
AGGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
"""
tss_index = 50  # Transcription start site position

guides = design_guides(sequence, tss_index)
print(guides.head())
```

**Output DataFrame columns:**
- `sequence`: 20bp guide sequence
- `position`: Distance from TSS (negative = upstream)
- `strand`: `+` or `-`
- `gc_content`: GC fraction (0-1)
- `delta_g`: Thermodynamic stability
- `combined_score`: Weighted composite score
- `go_no_go`: Classification status

**Advanced usage with custom configuration:**

```python
from phaselab.crispr import design_guides, CRISPRConfig

# Custom configuration
config = CRISPRConfig(
    pam="NGG",
    guide_length=20,
    window_upstream=500,    # Search 500bp upstream of TSS
    window_downstream=100,  # Search 100bp downstream
    min_gc=0.40,
    max_gc=0.70,
    weights={
        'gc': 0.25,
        'thermo': 0.25,
        'position': 0.25,
        'accessibility': 0.25
    }
)

# With DNase accessibility peaks
dnase_peaks = [(45, 85), (120, 180)]  # Accessible regions

guides = design_guides(
    sequence,
    tss_index,
    config=config,
    dnase_peaks=dnase_peaks,
    verbose=True
)

# Filter top candidates
top_guides = guides[guides['go_no_go'] == 'GO'].head(10)
```

### Integration with Spatial Coherence

The recommended workflow combines spatial coherence analysis with guide design:

```python
from phaselab.spatial import analyze_tiling_coherence, load_tiling_screen
from phaselab.crispr import design_crispra_guides

# 1. Load your tiling screen data
tiling_data = load_tiling_screen(
    'rai1_tiling_screen.tsv',
    position_col='position',
    response_col='log2fc',
    gene_symbol='RAI1',
)

# 2. Analyze spatial coherence
coherence_result = analyze_tiling_coherence(tiling_data)

# 3. Get stable region boundaries
stable_regions = coherence_result.stable_regions
print(f"Found {len(stable_regions)} stable regions")

# 4. Design guides restricted to stable regions
guides = design_crispra_guides(
    promoter_sequence=rai1_promoter,
    tss_position=600,
    # Only consider positions in stable regions
    position_filter=lambda pos: any(
        r['start'] <= pos <= r['end']
        for r in stable_regions
    ),
)

# 5. Annotate guides with local coherence
for guide in guides.candidates:
    pos = guide['position']
    # Find coherence at this position
    local_coh = coherence_result.get_coherence_at(pos)
    guide['local_coherence'] = local_coh
    guide['region_status'] = 'STABLE' if local_coh > 0.7 else 'MIXED'
```

---

## CRISPR-SURF Integration

*Added in v1.0.0*

### What is CRISPR-SURF?

[CRISPR-SURF](https://github.com/pinellolab/CRISPR-SURF) is a deconvolution algorithm that separates true regulatory signal from guide-specific noise in tiling screens. PhaseLab integrates with SURF output to provide coherence analysis on cleaner data.

**Why use SURF + PhaseLab together?**
- SURF removes guide-specific artifacts (efficiency, accessibility)
- PhaseLab then measures spatial coherence of the deconvolved signal
- Result: More reliable identification of regulatory regions

### Coherence on Deconvolved Data

```python
from phaselab.surf import (
    parse_surf_output,
    compute_surf_coherence,
    SURFPipeline,
)

# Parse SURF output files
surf_output = parse_surf_output(
    beta_file='surf_results/beta_profile.tsv',
    regions_file='surf_results/significant_regions.bed',
    gene_symbol='RAI1',
)

# Compute coherence on deconvolved data
result = compute_surf_coherence(
    surf_output,
    window=50,
    stable_threshold=0.7,
)

# Compare raw vs deconvolved coherence
print(f"Raw data coherence: {result.raw_coherence:.3f}")
print(f"SURF deconvolved coherence: {result.deconvolved_coherence:.3f}")
print(f"Improvement: {result.coherence_improvement:.1%}")

# High-confidence targets (stable in deconvolved data)
for target in result.high_confidence_targets:
    print(f"  {target['position']}: β={target['beta']:.3f}, coh={target['coherence']:.3f}")
```

### Full SURF Pipeline

```python
from phaselab.surf import SURFPipeline, SURFPipelineConfig

# Configure pipeline
config = SURFPipelineConfig(
    surf_path='/path/to/CRISPR-SURF',  # SURF installation
    window=50,
    stable_threshold=0.7,
    run_surf=True,  # Run SURF or use existing output
)

# Create and run pipeline
pipeline = SURFPipeline(config)
result = pipeline.run(
    input_file='tiling_screen.tsv',
    output_dir='surf_analysis/',
)

# Results include both SURF and coherence analysis
print(result.summary())
```

---

## Omics Module

*Added in v1.0.0*

The omics module applies spatial coherence to common genomics assays.

### ATAC-seq Coherence

Identify stably accessible chromatin regions:

```python
from phaselab.omics import (
    load_atac_data,
    analyze_atac_coherence,
    ATACLandscape,
)

# Load ATAC-seq signal
landscape = load_atac_data(
    'atac_signal.tsv',
    position_col='position',
    signal_col='signal',
    gene_symbol='RAI1',
    cell_type='neurons',
)

# Analyze coherence
result = analyze_atac_coherence(
    landscape,
    window=50,
    stable_threshold=0.7,
    peak_threshold=2.0,  # std above mean for peaks
)

# Find stable accessible regions (good for CRISPR targeting)
for region in result.accessible_stable:
    print(f"Stable peak: {region['start']}-{region['end']}")
    print(f"  Mean signal: {region['mean_signal']:.2f}")
    print(f"  Coherence: {region['coherence']:.3f}")
```

### ChIP-seq Coherence

Identify stable protein-DNA binding sites:

```python
from phaselab.omics import load_chip_data, analyze_chip_coherence

# Load ChIP-seq signal
landscape = load_chip_data(
    'h3k27ac_chip.tsv',
    target='H3K27ac',
    gene_symbol='MYC',
    cell_type='K562',
)

# Analyze
result = analyze_chip_coherence(landscape)

# Stable binding sites
print(f"Total peaks: {len(result.peak_regions)}")
print(f"Stable binding sites: {result.n_stable_peaks}")
```

### RNA-seq Coherence

Identify genes with reliable expression changes:

```python
from phaselab.omics import load_expression_data, analyze_expression_coherence

# Load expression data (genes as "positions")
landscape = load_expression_data(
    'differential_expression.tsv',
    gene_col='gene_id',
    expression_col='log2fc',
    condition='treatment_vs_control',
)

# Analyze
result = analyze_expression_coherence(landscape)

# Reliable changes (high coherence among similar genes)
for change in result.reliable_changes[:10]:
    print(f"{change['gene_id']}: {change['log2fc']:.2f} (coh={change['coherence']:.3f})")
```

---

## Microbiology Module

*Added in v1.0.0*

Apply spatial coherence to microbial screens and fitness assays.

### TnSeq Analysis

Identify essential genes/domains from transposon screens:

```python
from phaselab.microbio import (
    load_tnseq_data,
    analyze_tnseq_coherence,
    identify_essential_domains,
)

# Load TnSeq fitness data
landscape = load_tnseq_data(
    'tnseq_fitness.tsv',
    position_col='insertion_site',
    fitness_col='fitness_score',
    gene_symbol='essential_gene',
    organism='E. coli',
)

# Analyze coherence
result = analyze_tnseq_coherence(landscape)

# Find essential domains (low fitness, high coherence)
essential = identify_essential_domains(
    result,
    fitness_threshold=-2.0,  # Log2 fitness
    coherence_threshold=0.7,
)

for domain in essential:
    print(f"Essential domain: {domain['start']}-{domain['end']}")
    print(f"  Mean fitness: {domain['mean_fitness']:.2f}")
    print(f"  Coherence: {domain['coherence']:.3f}")
```

### Bacterial CRISPRi Screens

```python
from phaselab.microbio import load_crispri_screen, analyze_crispri_coherence

# Load bacterial CRISPRi screen
landscape = load_crispri_screen(
    'crispri_screen.tsv',
    gene='dnaA',
    organism='E. coli',
)

# Analyze
result = analyze_crispri_coherence(landscape)

# Rank guides by coherence-weighted score
for guide in result.ranked_guides[:10]:
    print(f"{guide['position']}: score={guide['score']:.3f}, coh={guide['coherence']:.3f}")
```

### Drug Response Landscapes

Identify stable dosing windows:

```python
from phaselab.microbio import (
    load_dose_response,
    analyze_drug_coherence,
    identify_stable_dosing_window,
)

# Load dose-response data
landscape = load_dose_response(
    'dose_response.tsv',
    concentration_col='concentration',
    response_col='viability',
    drug_name='Compound_X',
)

# Analyze
result = analyze_drug_coherence(landscape)

# Find stable dosing window (reliable response)
window = identify_stable_dosing_window(result)
print(f"Stable dosing window: {window['min_dose']:.2f} - {window['max_dose']:.2f}")
print(f"Expected response: {window['mean_response']:.2f} ± {window['response_std']:.2f}")
```

---

## Chemistry Module

*Added in v1.0.0*

Apply spatial coherence to chemical/biochemical systems.

### Binding Landscapes

Analyze protein-ligand or protein-protein binding:

```python
from phaselab.chem import (
    load_binding_data,
    analyze_binding_coherence,
    hot_spot_analysis,
)

# Load binding data (position = residue, response = ΔΔG)
landscape = load_binding_data(
    'alanine_scan.tsv',
    position_col='residue',
    binding_col='ddG',
    protein='antibody_CDR',
)

# Analyze coherence
result = analyze_binding_coherence(landscape)

# Find stable binding hot spots
hot_spots = hot_spot_analysis(result)
for spot in hot_spots:
    print(f"Hot spot: residues {spot['start']}-{spot['end']}")
    print(f"  Mean ΔΔG: {spot['mean_ddG']:.2f} kcal/mol")
    print(f"  Coherence: {spot['coherence']:.3f}")
```

### Reaction Optimization

Find stable reaction conditions:

```python
from phaselab.chem import (
    load_reaction_data,
    analyze_reaction_coherence,
    identify_stable_conditions,
)

# Load reaction optimization data
landscape = load_reaction_data(
    'reaction_screen.tsv',
    condition_col='temperature',  # or pH, concentration, etc.
    yield_col='yield',
    reaction='Suzuki_coupling',
)

# Analyze
result = analyze_reaction_coherence(landscape)

# Find stable operating window
stable = identify_stable_conditions(result)
print(f"Optimal temperature range: {stable['min']:.0f}°C - {stable['max']:.0f}°C")
print(f"Expected yield: {stable['mean_yield']:.0f}% ± {stable['yield_std']:.0f}%")
```

### HTS Screening

Identify reliable hits from high-throughput screens:

```python
from phaselab.chem import (
    load_screening_data,
    analyze_screening_coherence,
)

# Load HTS data
landscape = load_screening_data(
    'hts_plate.tsv',
    position_col='well',
    activity_col='inhibition',
    assay_name='kinase_inhibition',
)

# Analyze
result = analyze_screening_coherence(landscape)

# Reliable hits (active + in coherent region)
print(f"Total hits: {len(result.reliable_hits)}")
print(f"Hit rate: {result.hit_rate:.1%}")

for hit in result.reliable_hits[:10]:
    compound = hit.get('compound_id', f"well_{hit['position']}")
    print(f"{compound}: activity={hit['activity']:.1f}%, coh={hit['coherence']:.3f}")
```

---

## Protein Module

*Added in v1.0.0*

Analyze mutational scanning and protein fitness landscapes.

### Mutational Scanning

```python
from phaselab.protein import (
    load_mutscan_data,
    analyze_mutscan_coherence,
    local_coherence_profile,
)

# Load deep mutational scan data
landscape = load_mutscan_data(
    'dms_data.tsv',
    position_col='position',
    fitness_col='fitness',
    protein='GFP',
)

# Analyze coherence
result = analyze_mutscan_coherence(landscape)

# Find functional domains (conserved, high coherence)
for domain in result.functional_domains:
    print(f"Functional domain: {domain['start']}-{domain['end']}")
    print(f"  Mean fitness effect: {domain['mean_fitness']:.3f}")
    print(f"  Conservation: {domain['conservation']:.2f}")
    print(f"  Coherence: {domain['coherence']:.3f}")

# Local coherence profile for structure mapping
profile = local_coherence_profile(landscape, window=5)
```

### Structure Mapping

Map coherence to 3D structure:

```python
from phaselab.protein import map_coherence_to_structure

# Map coherence values to B-factor column in PDB
map_coherence_to_structure(
    result,
    pdb_file='protein.pdb',
    output_file='protein_coherence.pdb',
    chain='A',
)

# Visualize in PyMOL:
# spectrum b, blue_white_red
# Low coherence (variable) = blue
# High coherence (stable) = red
```

---

## Quantum Mode Configuration

*Added in v1.0.0*

Configure how PhaseLab uses quantum computation for coherence validation.

### Mode Options

```python
from phaselab.quantum import (
    QuantumMode,
    set_quantum_mode,
    get_quantum_mode,
    quantum_status,
)

# Available modes
# QuantumMode.OFF      - Classical only (fastest, default)
# QuantumMode.AUDIT    - Classical + quantum validation on subset
# QuantumMode.REQUIRED - Quantum mandatory for all coherence

# Set mode
set_quantum_mode(QuantumMode.OFF)      # Default, fastest
set_quantum_mode(QuantumMode.AUDIT)    # Validate with quantum
set_quantum_mode(QuantumMode.REQUIRED) # Full quantum (slowest)

# Check current mode
mode = get_quantum_mode()
print(f"Current mode: {mode.name}")

# Get detailed status
status = quantum_status()
print(status)
# Quantum Mode: AUDIT
# ATLAS-Q Available: Yes
# IBM Quantum Token: Configured
# Default Backend: ibm_torino
```

### Full Configuration

```python
from phaselab.quantum import configure_quantum, QuantumConfig

# Full configuration
config = QuantumConfig(
    mode=QuantumMode.AUDIT,
    audit_fraction=0.1,        # Validate 10% of calculations
    backend='ibm_torino',      # IBM Quantum backend
    shots=1000,                # Measurements per circuit
    use_gpu=True,              # GPU acceleration if available
    cache_results=True,        # Cache quantum results
)

configure_quantum(config)
```

### When to Use Quantum

| Mode | Use Case | Speed |
|------|----------|-------|
| **OFF** | Development, exploration, large screens | Fastest |
| **AUDIT** | Publication, validation, spot-checking | Medium |
| **REQUIRED** | High-stakes decisions, regulatory | Slowest |

Most users should use **OFF** for exploration and **AUDIT** for final analysis. **REQUIRED** mode is only needed when you need quantum-validated coherence for every calculation.

---

## Circadian Module

The circadian module models biological clock systems for Smith-Magenis Syndrome research.

### SMS Clock Model

Simulate the SMS circadian clock with RAI1 haploinsufficiency:

```python
from phaselab.circadian import simulate_sms_clock, SMSClockParams

# Basic simulation with 50% RAI1 (SMS condition)
results = simulate_sms_clock(rai1_level=0.5)

print(f"Period: {results['period']:.2f} hours")
print(f"Amplitude: {results['amplitude']:.4f}")
print(f"Phase shift: {results['phase_shift']:.2f} hours")
print(f"Coherence R̄: {results['coherence']:.4f}")
print(f"Status: {results['status']}")
```

**Output:**
```
Period: 23.45 hours
Amplitude: 0.6234
Phase shift: 2.15 hours
Coherence R̄: 0.7823
Status: GO
```

**Custom parameters:**

```python
from phaselab.circadian import simulate_sms_clock, SMSClockParams

# Custom clock parameters
params = SMSClockParams(
    tau_P=4.0,      # PER delay time constant (hours)
    alpha_P=2.0,    # PER suppression strength
    beta_R=0.5,     # RORα effect on BMAL1
    beta_V=0.5,     # REV-ERBα effect on BMAL1
    K_light=0.1,    # Light sensitivity
    omega_0=2*np.pi/24  # Base angular frequency
)

# Simulate with intervention (CRISPRa restoring RAI1)
results_baseline = simulate_sms_clock(rai1_level=0.5, params=params)
results_treated = simulate_sms_clock(rai1_level=0.85, params=params)

print(f"Baseline coherence: {results_baseline['coherence']:.4f}")
print(f"Treated coherence: {results_treated['coherence']:.4f}")
```

**Full time series output:**

```python
results = simulate_sms_clock(
    rai1_level=0.5,
    t_end=240.0,      # 10 days
    dt=0.1,           # Time step (hours)
    return_timeseries=True
)

# Access time series
t = results['time']
bmal1 = results['BMAL1']
per = results['PER']
rev_erb = results['REV_ERB']
ror = results['ROR']

# Plot
import matplotlib.pyplot as plt
plt.figure(figsize=(12, 6))
plt.plot(t, bmal1, label='BMAL1')
plt.plot(t, per, label='PER')
plt.xlabel('Time (hours)')
plt.ylabel('Expression level')
plt.legend()
plt.show()
```

### Kuramoto Oscillators

Model coupled oscillator networks:

```python
from phaselab.circadian import KuramotoNetwork

# Create network of 100 oscillators
network = KuramotoNetwork(
    n_oscillators=100,
    coupling_strength=0.5,
    natural_frequency_spread=0.1
)

# Simulate
results = network.simulate(t_end=100.0)

print(f"Order parameter R: {results['order_parameter']:.4f}")
print(f"Synchronization: {'Yes' if results['synchronized'] else 'No'}")
```

**Heterogeneous network:**

```python
from phaselab.circadian import KuramotoNetwork
import numpy as np

# Custom frequency distribution (bimodal)
frequencies = np.concatenate([
    np.random.normal(1.0, 0.1, 50),
    np.random.normal(1.1, 0.1, 50)
])

network = KuramotoNetwork(
    n_oscillators=100,
    coupling_strength=0.8,
    natural_frequencies=frequencies
)

# With external forcing (light-dark cycle)
def light_forcing(t):
    return 0.1 * np.sin(2 * np.pi * t / 24)

results = network.simulate(
    t_end=240.0,
    external_forcing=light_forcing
)
```

---

## SMS Trials Module

*Added in v0.9.0*

The SMS Trials module provides a complete therapeutic trial framework for Smith-Magenis Syndrome gene therapy development.

### Trial Runners

Individual trial runners for each therapeutic modality:

```python
from phaselab.trials.sms import (
    run_sms_crispra_trial,
    run_sms_crispri_trial,
    run_sms_knockout_trial,
    run_sms_base_editing_trial,
    run_sms_prime_editing_trial,
    run_circadian_rescue_simulation,
    run_delivery_assessment,
    SMSTrialConfig,
)

# Configure trials
config = SMSTrialConfig(
    therapeutic_window=(0.70, 1.10),  # 70-110% of normal RAI1
    optimal_expression=0.80,          # Target 80%
    baseline_expression=0.50,         # SMS baseline is ~50%
    use_virtual_assay=True,           # Use enhanced pipeline
    coherence_mode="heuristic",       # Fast mode
    verbose=True,
)

# CRISPRa trial for RAI1 activation
crispra_result = run_sms_crispra_trial(
    promoter_sequence=rai1_promoter,  # Your promoter sequence
    config=config,
)

print(f"Status: {crispra_result.status}")
print(f"Candidates: {crispra_result.n_candidates}")
print(f"Claim level: {crispra_result.claim_level}")
if crispra_result.best_candidate:
    print(f"Best guide: {crispra_result.best_candidate['sequence']}")
```

#### CRISPRi Modifier Suppression

```python
# Target circadian modifier genes (PER1, PER2, CRY1)
crispri_result = run_sms_crispri_trial(
    target_gene="PER1",
    config=config,
)

# Check suppression levels (target: 30-60% suppression)
if crispri_result.best_candidate:
    supp = crispri_result.best_candidate['expected_suppression']
    print(f"Expected suppression: {supp:.0%}")
```

#### Circadian Rescue Simulation

```python
# Predict circadian rescue from RAI1 boost
circadian_result = run_circadian_rescue_simulation(
    predicted_rai1_expression=0.80,  # Expected RAI1 level after treatment
    config=config,
)

print(f"Rescue status: {circadian_result.metrics['rescue_status']}")
print(f"Final R̄: {circadian_result.metrics['final_R_bar']:.3f}")
print(f"Sleep quality: {circadian_result.metrics['sleep_quality_prediction']}")
```

#### Delivery Assessment

```python
# Assess AAV delivery feasibility for CNS
delivery_result = run_delivery_assessment(
    modality="CRISPRa_VP64",
    target_tissue="brain",
    config=config,
)

print(f"Feasibility: {delivery_result.metrics['delivery_feasibility']}")
print(f"Payload size: {delivery_result.metrics['payload_size']}bp")
print(f"Recommended serotype: {delivery_result.best_candidate['recommended_serotype']}")
print(f"Capacity utilization: {delivery_result.metrics['packaging_utilization']:.0%}")
```

### SMS Pipeline

The SMSPipeline orchestrates all trials and provides integrated GO/NO-GO decisions:

```python
from phaselab.trials.sms import SMSPipeline, SMSTrialConfig

# Configure pipeline
config = SMSTrialConfig(
    therapeutic_window=(0.70, 1.10),
    use_virtual_assay=True,
    simulation_hours=120.0,
    n_circadian_trials=5,
)

# Create and run pipeline
pipeline = SMSPipeline(config=config)
result = pipeline.run_full_pipeline(
    include_modifiers=True,   # Include CRISPRi trials
    include_editing=True,     # Include base/prime editing
)

# Overall decision
print(f"Overall GO/NO-GO: {result.overall_go_nogo}")
print(f"Claim level: {result.overall_claim_level}")

# Individual results
print(f"CRISPRa candidates: {result.crispra_result.n_candidates}")
print(f"Circadian rescue: {result.circadian_result.metrics['rescue_status']}")
print(f"Delivery feasible: {result.delivery_result.metrics['delivery_feasibility']}")

# Wet lab recommendations
for rec in result.wet_lab_recommendations:
    print(f"- {rec}")
```

### Falsification Tests

The pipeline automatically generates falsification tests for wet-lab validation:

```python
# Get falsification tests
for test in result.falsification_tests:
    print(f"\nTest {test['id']}: {test['name']}")
    print(f"  Description: {test['description']}")
    print(f"  Failure condition: {test['failure_condition']}")
    print(f"  Required for: {test['required_for']}")
```

**Available Tests:**

| Test | Name | Purpose |
|------|------|---------|
| **A** | Ranking Validity | PhaseLab-ranked guides must outperform random CRISPOR-acceptable controls |
| **B** | Risk Prediction | High-risk (CAUTION) guides should fail safety screen more often |
| **C** | Dosage Prediction | Predicted expression levels should correlate (r > 0.6) with observed |
| **D** | UNKNOWN Bucket | UNKNOWN-labeled guides should fail at approximately random rate |

### Trial Result Structure

All trial runners return `SMSTrialResult` with consistent structure:

```python
@dataclass
class SMSTrialResult:
    trial_type: TrialType          # Type of trial (CRISPRA_RAI1, CRISPRI_MODIFIER, etc.)
    status: TrialStatus            # COMPLETED, FAILED, PENDING
    summary: str                   # Human-readable summary
    candidates: List[Dict]         # Ranked candidates
    best_candidate: Optional[Dict] # Top candidate
    metrics: Dict[str, Any]        # Trial-specific metrics
    claim_level: str               # STRONG_COMPUTATIONAL, CONTEXT_DEPENDENT, etc.
    claim_description: str         # Explanation of claim level
    warnings: List[str]            # Important warnings
    errors: List[str]              # Any errors encountered

    # Helper properties
    @property
    def is_successful(self) -> bool: ...
    @property
    def has_viable_candidates(self) -> bool: ...
    @property
    def n_candidates(self) -> int: ...
```

---

## Quantum Discriminator

*Added in v0.9.5*

The Quantum Discriminator is a late-stage guide selection tool that uses quantum chemistry on IBM Quantum hardware to discriminate between elite CRISPRa guides whose predicted effectiveness is classically indistinguishable.

**Scientific Claim (defensible):**
> "IR-enhanced quantum VQE on current IBM hardware can resolve binding energy differences between CRISPRa guides that are indistinguishable under classical scoring, providing a physically grounded late-stage discriminator for therapeutic guide selection."

**What this IS:**
- Quantum resolves energetic degeneracy
- Quantum increases hit rate per wet-lab experiment
- Quantum reduces wasted biological trials

**What this is NOT:**
- Not "quantum finds cures"
- Not "quantum replaces biology"
- Not "quantum predicts expression directly"

### When to Use Quantum Discrimination

Use the quantum discriminator when:

1. **Classical methods saturate**: Top guides have combined scores within 5% of each other
2. **Elite tier only**: Only for guides that pass all classical gates
3. **Late-stage**: After multi-evidence pipeline (binding, phase, geometry) has ranked candidates

```python
from phaselab.crispr import DISCRIMINATOR_GATES

# Pre-quantum gates that must be passed:
print(DISCRIMINATOR_GATES)
# {
#     'min_mit_score': 50,
#     'max_exonic_ot': 0,
#     'min_delta_r': 0.30,
#     'min_phase_coherence': 0.90,
#     'min_guides_for_quantum': 2,
# }
```

### Effective Binding Hamiltonian

The discriminator constructs an effective Hamiltonian for RNA-DNA binding:

$$H = H_{HB} + H_{stack} + H_{charge} + H_{constraint}$$

Where:
- **H_HB**: Watson-Crick hydrogen bonding (G-C: -0.18 eV, A-T: -0.12 eV)
- **H_stack**: π-π stacking stabilization between adjacent bases
- **H_charge**: Backbone electrostatic interactions with 0.7 screening
- **H_constraint**: Prevents unphysical strand separation

The Hamiltonian uses the seed region (12bp PAM-proximal) encoded as Pauli operators:

```python
from phaselab.crispr.quantum_discriminator import _build_effective_binding_hamiltonian

guide = "GCGCGCGCGCGCGCGCGCGC"
target = "CGCGCGCGCGCGCGCGCGCG"

coefficients, paulis, n_qubits = _build_effective_binding_hamiltonian(guide, target)

print(f"Qubits: {n_qubits}")        # 12
print(f"Pauli terms: {len(paulis)}") # ~44
print(f"First term: {coefficients[0]:.4f} * {paulis[0]}")  # -0.1800 * ZIIIIIIIIIII
```

### Running the Discriminator

#### Basic Usage (Simulation)

```python
from phaselab.crispr import run_quantum_discriminator, DiscriminatorStatus

# Elite guides from classical pipeline (scores are degenerate)
guides = [
    {'sequence': 'GCGCGCGCGCGCGCGCGCGC', 'combined_score': 0.92},
    {'sequence': 'ATCGATCGATCGATCGATCG', 'combined_score': 0.91},
    {'sequence': 'CCGGCCGGCCGGCCGGCCGG', 'combined_score': 0.90},
]

# DNA target context
dna_context = "CGCGCGCGCGCGCGCGCGCG..."  # Promoter region

# Run discriminator (simulation mode)
result = run_quantum_discriminator(
    guides=guides,
    dna_context=dna_context,
    backend_name="ibm_torino",
    use_hardware=False,  # Simulation
    shots=1000,
    max_iterations=30,
    degeneracy_threshold=0.05,
)

print(result.summary())
# ======================================================================
# QUANTUM DISCRIMINATOR RESULT
# ======================================================================
# Status: quantum_success
# Guides evaluated: 3
# Quantum advantage: YES
#
# Ranking by quantum binding energy:
#   #1: CCGGCCGGCCGGCCG... E=-1.993470 Ha [GO]
#   #2: ATCGATCGATCGATC... E=-1.681406 Ha [GO]
#   #3: GCGCGCGCGCGCGCG... E=-0.556050 Ha [GO]
```

#### Running on IBM Quantum Hardware

```python
import os
os.environ['IBM_QUANTUM_TOKEN'] = 'your-token-here'

result = run_quantum_discriminator(
    guides=guides,
    dna_context=dna_context,
    backend_name="ibm_torino",
    use_hardware=True,  # Real hardware
    shots=1000,
    max_iterations=30,
)

# Check individual guide results
for g in result.ranked_guides:
    print(f"{g.guide_sequence[:15]}... E={g.binding_energy:.6f} Ha, "
          f"R̄={g.coherence:.3f}, {'GO' if g.is_go else 'NO-GO'}")
```

#### Understanding Results

```python
# Result structure
print(f"Status: {result.status}")              # DiscriminatorStatus.QUANTUM_SUCCESS
print(f"Quantum advantage: {result.quantum_advantage}")  # True if ordering differs

# Energy separations between guide pairs
for pair, delta_e in result.energy_separations.items():
    is_sig = result.significant_separations[pair]
    print(f"{pair}: ΔE = {delta_e:.6f} Ha {'*' if is_sig else ''}")

# Classical vs quantum comparison
print(f"Classical scores: {result.classical_scores}")
```

### PhaseLab API Integration

The high-level API integrates quantum discrimination into the full design pipeline:

```python
from phaselab.crispr import design_guides_with_quantum_discriminator

# After running classical pipeline...
result = design_guides_with_quantum_discriminator(
    gene="RAI1",
    guides=top_classical_guides,  # From multi-evidence pipeline
    dna_context=rai1_promoter,
    quantum_stage="late",          # Only use at end
    quantum_backend="ibm_torino",
    use_hardware=False,            # True for real hardware
    max_quantum_guides=3,          # Limit quantum evaluations
)

print(f"Gene: {result['gene']}")
print(f"Quantum advantage: {result['quantum_advantage']}")
print(f"Recommendation: {result['recommendation']}")

# Final ranking (quantum-ordered if advantage detected)
for i, g in enumerate(result['final_ranking'], 1):
    if g.get('ranking_source') == 'quantum':
        print(f"#{i}: {g['sequence'][:15]}... E={g['quantum_energy']:.6f} Ha")
    else:
        print(f"#{i}: {g['sequence'][:15]}... score={g['classical_score']:.3f}")
```

### Discriminator Status Codes

| Status | Meaning |
|--------|---------|
| `QUANTUM_SUCCESS` | VQE completed, results available |
| `NO_DEGENERACY` | Classical scores differ enough, quantum not needed |
| `INSUFFICIENT_GUIDES` | Need ≥2 guides to discriminate |
| `QUANTUM_FAILED` | VQE failed, fell back to classical |

### GO/NO-GO Threshold

Quantum execution quality is validated using IR coherence:

$$\text{GO if } \bar{R} > e^{-2} \approx 0.135$$

Guides with coherence below this threshold during VQE execution are marked NO-GO and excluded from final ranking.

### API Reference

| Function | Description |
|----------|-------------|
| `run_quantum_discriminator(guides, dna_context, ...)` | Main discriminator entry point |
| `design_guides_with_quantum_discriminator(gene, guides, ...)` | High-level PhaseLab API |
| `DiscriminatorStatus` | Enum for result status |
| `QuantumGuideResult` | Single guide quantum result |
| `DiscriminatorResult` | Complete discriminator result |
| `DISCRIMINATOR_GATES` | Pre-quantum gate thresholds |

---

## Quantum Integration

PhaseLab integrates with IBM Quantum for hardware validation.

### Setting Up IBM Quantum

```python
import os
os.environ['IBM_QUANTUM_TOKEN'] = 'your-token-here'

# Or use .env file
# IBM_QUANTUM_TOKEN=your-token-here
```

### Running on Quantum Hardware

```python
from phaselab.quantum import QuantumCoherenceValidator

# Initialize validator
validator = QuantumCoherenceValidator(
    backend='ibm_torino',  # or 'ibm_brisbane', 'ibm_kyoto'
    shots=4096
)

# Validate a guide RNA
guide = "GAAGGAGAGCAAGAGCGCGA"
result = validator.validate(guide)

print(f"Hardware R̄: {result['coherence']:.4f}")
print(f"Status: {result['status']}")
print(f"Job ID: {result['job_id']}")
```

### Simulator Testing

```python
from phaselab.quantum import QuantumCoherenceValidator

# Use Aer simulator
validator = QuantumCoherenceValidator(backend='aer_simulator')

guides = [
    "GAAGGAGAGCAAGAGCGCGA",
    "AACTGCAAAGAAGTGGGCAC",
    "TACAGGAGCTTCCAGCGTCA"
]

for guide in guides:
    result = validator.validate(guide)
    print(f"{guide[:10]}... R̄={result['coherence']:.3f} [{result['status']}]")
```

---

## Complete Examples

### Example 1: SMS Gene Therapy Pipeline

```python
"""
Complete SMS gene therapy gRNA design pipeline.
"""
import phaselab as pl
from phaselab.crispr import design_guides, CRISPRConfig
from phaselab.circadian import simulate_sms_clock

# RAI1 promoter sequence (truncated for example)
rai1_promoter = """
GCGCGCTCGCGCGCTCGCGCGAAGGAGAGCAAGAGCGCGACGGCTAGCTAGCT
AGCTAGCTAGCTACAGGAGCTTCCAGCGTCAGGGCTAGCTAGCTAGCTAGCTA
GCTAGCTAGCTAACTGCAAAGAAGTGGGCACGCGCTAGCTAGCTAGCTAGCT
"""

# TSS position
tss_index = len(rai1_promoter) // 2

# Configure pipeline
config = CRISPRConfig(
    pam="NGG",
    window_upstream=400,
    window_downstream=50,
    min_gc=0.40,
    max_gc=0.70
)

# Design guides
guides = design_guides(rai1_promoter, tss_index, config=config)

# Filter GO candidates
go_candidates = guides[guides['go_no_go'] == 'GO']
print(f"Found {len(go_candidates)} GO candidates")

# Simulate therapeutic effect
for idx, row in go_candidates.head(3).iterrows():
    print(f"\n{row['sequence']}")
    print(f"  Position: {row['position']} bp from TSS")
    print(f"  Score: {row['combined_score']:.3f}")

    # Simulate with restored RAI1 (hypothetical 80% restoration)
    sim = simulate_sms_clock(rai1_level=0.8)
    print(f"  Predicted period: {sim['period']:.2f}h")
    print(f"  Predicted coherence: {sim['coherence']:.4f}")
```

### Example 2: Comparative Analysis

```python
"""
Compare coherence across multiple conditions.
"""
import phaselab as pl
from phaselab.circadian import simulate_sms_clock
import numpy as np

# RAI1 levels to test
rai1_levels = np.linspace(0.3, 1.0, 8)

results = []
for level in rai1_levels:
    sim = simulate_sms_clock(rai1_level=level)
    results.append({
        'rai1_level': level,
        'period': sim['period'],
        'coherence': sim['coherence'],
        'status': sim['status']
    })

# Find therapeutic threshold
import pandas as pd
df = pd.DataFrame(results)
go_threshold = df[df['status'] == 'GO']['rai1_level'].min()
print(f"Minimum RAI1 for GO status: {go_threshold:.1%}")

# Visualize
import matplotlib.pyplot as plt
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(df['rai1_level'], df['coherence'], 'bo-')
ax1.axhline(pl.core.constants.E_MINUS_2, color='r', linestyle='--', label='GO threshold')
ax1.set_xlabel('RAI1 Level')
ax1.set_ylabel('Coherence R̄')
ax1.legend()

ax2.plot(df['rai1_level'], df['period'], 'go-')
ax2.axhline(24, color='gray', linestyle='--', label='24h')
ax2.set_xlabel('RAI1 Level')
ax2.set_ylabel('Period (hours)')
ax2.legend()

plt.tight_layout()
plt.show()
```

### Example 3: Batch Processing

```python
"""
Batch process multiple sequences.
"""
from phaselab.crispr import design_guides
import pandas as pd

# Multiple target genes
targets = {
    'RAI1': 'ATGCGATCG...',
    'CLOCK': 'GCTAGCTAG...',
    'BMAL1': 'TAGCTAGCT...',
}

all_guides = []
for gene, sequence in targets.items():
    guides = design_guides(sequence, tss_index=len(sequence)//2)
    guides['target_gene'] = gene
    all_guides.append(guides)

combined = pd.concat(all_guides)
combined = combined.sort_values('combined_score', ascending=False)

# Top guide per gene
top_per_gene = combined.groupby('target_gene').first()
print(top_per_gene[['sequence', 'combined_score', 'go_no_go']])
```

---

## API Reference

### Core Functions

| Function | Description |
|----------|-------------|
| `coherence_score(data, mode='auto')` | Compute R̄ coherence score |
| `phase_variance(phases)` | Compute phase variance V_φ |
| `go_no_go(R_bar, threshold=E_MINUS_2)` | Binary GO/NO-GO classification |
| `classify_coherence(R_bar)` | Detailed category classification |

### CRISPR Functions

| Function | Description |
|----------|-------------|
| `design_guides(seq, tss, config, dnase_peaks)` | High-level guide design pipeline |
| `find_pam_sites(seq, pam, guide_length)` | Find PAM sites in sequence |
| `score_guide(guide, weights)` | Score individual guide |
| `gc_content(seq)` | Calculate GC fraction |
| `santalucia_delta_g(seq)` | Calculate thermodynamic ΔG |

### Circadian Functions

| Function | Description |
|----------|-------------|
| `simulate_sms_clock(rai1_level, params, t_end)` | Simulate SMS clock model |
| `KuramotoNetwork(n, coupling, frequencies)` | Create Kuramoto oscillator network |

### Constants

| Constant | Value | Description |
|----------|-------|-------------|
| `E_MINUS_2` | 0.1353... | e⁻² GO/NO-GO threshold |
| `FOUR_PI_SQUARED` | 39.478... | 4π² universal constant |

---

## Citation

If you use PhaseLab in your research, please cite:

```bibtex
@software{phaselab,
  author = {Vaca, Dylan},
  title = {PhaseLab: Phase-coherence analysis framework},
  year = {2024},
  url = {https://github.com/followthesapper/phaselab}
}
```

---

## License

MIT License - see [LICENSE](LICENSE) for details.
