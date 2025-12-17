# PhaseLab

**A reliability layer for perturbation science.**

[![PyPI version](https://badge.fury.io/py/phaselab.svg)](https://badge.fury.io/py/phaselab)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://github.com/followthesapper/phaselab/actions/workflows/publish.yml/badge.svg)](https://github.com/followthesapper/phaselab/actions)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://followthesapper.github.io/phaselab/)

PhaseLab assesses whether experimental results from perturbation experiments will be reproducible—before you run the experiment.

## The Core Insight

> **Spatial coherence of response landscapes predicts perturbation reliability.**

Validated across 115,251 sgRNAs (6 genes):
- Correlation with outcome variance: r = -0.24 to -0.50
- Variance reduction in stable regions: 32-49%

## Quick Start

```bash
pip install phaselab
```

### Assess CRISPR Tiling Reliability

```python
from phaselab.spatial import analyze_tiling_coherence
from phaselab.landscapes import ResponseLandscape

# Your tiling screen data
landscape = ResponseLandscape(
    coords=[100, 150, 200, 250, 300],  # Guide positions
    responses=[0.8, 0.75, 0.1, 0.82, 0.78],  # Log2 fold-changes
)

# Analyze spatial coherence
result = analyze_tiling_coherence(landscape)

# Find stable regions
for region in result.stable_regions:
    print(f"Stable: {region['start']}-{region['end']}")
    print(f"  Coherence: {region['coherence']:.3f}")
    print(f"  Recommendation: Select guides from this region")
```

### Analyze Protein Mutational Scanning

```python
from phaselab.protein.mutscan import MutScanLandscape, analyze_mutscan_coherence

# Your DMS data
landscape = MutScanLandscape(
    positions=residue_positions,
    effects=fitness_effects,
    protein_id="TEM1",
    protein_name="TEM-1 β-Lactamase",
)

# Find functional domains
result = analyze_mutscan_coherence(landscape)

print(f"Essential domains: {len(result.essential_regions)}")
print(f"Variable regions: {len(result.variable_regions)}")
```

### Assess Binding Landscape Reliability

```python
from phaselab.chem import analyze_binding_coherence, BindingLandscape

# Your binding data
landscape = BindingLandscape(
    positions=mutation_positions,
    affinities=delta_delta_g_values,
    target="ABL1",
    ligand="Imatinib",
)

# Find reliable binding determinants
result = analyze_binding_coherence(landscape)

for hotspot in result.hot_spots:
    print(f"Position {hotspot['position']}: ΔΔG={hotspot['effect']:.2f}")
```

## What PhaseLab Is

| PhaseLab Is | PhaseLab Is NOT |
|-------------|-----------------|
| A reliability assessment framework | A CRISPR design tool |
| Domain-general (CRISPR, protein, chemistry) | A predictor of biological outcomes |
| Pre-experimental guidance | A replacement for validation |
| Honest about uncertainty | A guarantee of success |

**PhaseLab answers:** "Can I trust the results I'll get?"

**CRISPOR answers:** "Which guide should I use?"

Use both. They're complementary.

## Two-Stage Framework

### Stage I: Feasibility (Before Tiling)

- Uses structural priors and response modeling
- Identifies candidate stable regions
- GO/NO-GO decision for experimental investment

### Stage II: Resolution (With Tiling)

- Minimum 16-20 perturbations per target
- Empirically resolves stability structure
- Quantifies variance reduction

## Modules

| Module | Domain | Purpose |
|--------|--------|---------|
| `phaselab.spatial` | CRISPR | Tiling screen coherence |
| `phaselab.protein.mutscan` | Protein | DMS functional domains |
| `phaselab.protein.folding` | Protein | Structure prediction QC |
| `phaselab.chem.binding` | Chemistry | Binding hot spots |
| `phaselab.omics` | Genomics | ATAC/ChIP/RNA-seq reliability |
| `phaselab.microbio` | Microbiology | TnSeq/CRISPRi essentiality |
| `phaselab.trials.sms` | Therapeutic | SMS gene therapy pipeline |

## Claim Levels

PhaseLab reports uncertainty honestly:

| Level | Meaning | Requirements |
|-------|---------|--------------|
| `UNKNOWN` | Cannot assess | Insufficient data |
| `EXPLORATORY` | Preliminary | Structural priors only |
| `CONTEXT_DEPENDENT` | Valid in context | Single tiling dataset |
| `STRONG_COMPUTATIONAL` | High confidence | Cross-validated |

## Quantum Mode

For most use cases, quantum mode should be OFF (default):

```python
from phaselab.quantum import set_quantum_mode, QuantumMode

set_quantum_mode(QuantumMode.OFF)  # Classical only (fastest, default)
set_quantum_mode(QuantumMode.AUDIT)  # Validation subset
set_quantum_mode(QuantumMode.REQUIRED)  # Research only
```

**Rule:** If a classical experiment can falsify a claim, quantum is optional.

## Installation

```bash
# Basic installation
pip install phaselab

# With quantum support
pip install phaselab[quantum]

# Full installation
pip install phaselab[all]
```

## Documentation

- **[API Guide](docs/API_GUIDE.md)** - Complete API reference
- **[Two-Stage Framework](docs/TWO_STAGE_PERTURBATION_FRAMEWORK.md)** - Methodology
- **[Claims and Limits](docs/CLAIMS_AND_LIMITS.md)** - What PhaseLab guarantees
- **[Quantum Mode Guidance](docs/QUANTUM_MODE_GUIDANCE.md)** - When to use quantum
- **[Positioning](docs/PHASELAB_POSITIONING.md)** - Why PhaseLab is not a CRISPR tool

### Research Documentation

- **[SMS Gene Therapy](docs/SMS_GENE_THERAPY_RESEARCH.md)** - RAI1 CRISPRa trials
- **[Methods Paper Outline](docs/METHODS_PAPER_OUTLINE.md)** - Publication pathway
- **[RAI1 Tiling Spec](docs/RAI1_TILING_EXPERIMENT_SPEC.md)** - Wet lab protocol

## Validation

| Dataset | Scale | Correlation | Variance Reduction |
|---------|-------|-------------|-------------------|
| CRISPRa tiling (6 genes) | 115,251 sgRNAs | r = -0.34 | 42% mean |
| E213-E216 experiments | Multiple screens | r = -0.24 to -0.50 | 32-49% |

## Example: SMS Therapeutic Pipeline

```python
from phaselab.trials.sms import SMSPipeline, SMSTrialConfig

# Configure pipeline
config = SMSTrialConfig(
    therapeutic_window=(0.70, 1.10),
    verbose=True,
)

# Run full assessment
pipeline = SMSPipeline(config=config)
result = pipeline.run_full_pipeline()

print(f"GO/NO-GO: {result.overall_go_nogo}")
print(f"Claim level: {result.overall_claim_level}")

# Get falsification tests for wet lab
for test in result.falsification_tests:
    print(f"Test {test['id']}: {test['name']}")
```

## The Science

**E200-E211**: Guide-sequence coherence does NOT work (r ≈ 0)

**E213-E216**: Spatial coherence of response landscapes DOES work

The key insight:
> "The guide is the probe, not the structure. Coherence measures the system's response consistency, not properties of the perturbation itself."

## Citation

```bibtex
@software{phaselab2025,
  author = {Vaca, Dylan},
  title = {PhaseLab: A reliability layer for perturbation science},
  year = {2025},
  version = {1.0.0},
  url = {https://github.com/followthesapper/phaselab}
}
```

## Contributing

Contributions welcome! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

MIT License - see [LICENSE](LICENSE) for details.

---

*PhaseLab v1.0.0 - Spatial Coherence Paradigm*
*315 tests passing | Production ready*
