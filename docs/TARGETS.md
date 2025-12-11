# PhaseLab Target Library

PhaseLab supports gene target configurations for CRISPRa/CRISPRi experiments. Each target is defined by a YAML configuration file containing genomic coordinates, CRISPRa parameters, and metadata.

## Available Targets

| Target | Disease | Mode | Status |
|--------|---------|------|--------|
| **RAI1** | Smith-Magenis Syndrome | Haploinsufficiency | Hardware validated (IBM Torino) |
| **SCN2A** | Autism-linked NDD | Haploinsufficiency | Hardware validated (IBM Torino) |

### Validated Guide Candidates

| Target | Lead Candidate | MIT | CFD | Off-targets ≤2mm | Hardware R̄ |
|--------|----------------|-----|-----|------------------|-------------|
| RAI1 | `TACAGGAGCTTCCAGCGTCA` | 83 | 93 | 0 | 0.839 |
| SCN2A | `GCTGACTGCTACATAGCCAA` | 83 | 89 | 0 | 0.970 |

## Adding a New Target

### Step 1: Find Genomic Coordinates

Use Ensembl, NCBI, or UCSC Genome Browser to identify:
- Canonical transcript and TSS
- Chromosome and coordinates
- Promoter region (~1000bp upstream, ~500bp downstream of TSS)

### Step 2: Create YAML Configuration

Create a file in `src/phaselab/targets/GENE_NAME.yaml`:

```yaml
gene_symbol: "GENE_NAME"
genome_build: "GRCh38"

# Gene identifiers
gene_id: "ENSG..."
refseq_id: "NM_..."

promoter:
  chrom: "chrN"
  start: 123456789     # TSS - 1000 bp
  end:   123458289     # TSS + 500 bp
  tss:   123457789     # Transcription start site

crispr:
  system: "SpCas9"
  pam: "NGG"
  window:
    start_offset: -400
    end_offset: -50

filters:
  gc_min: 0.40
  gc_max: 0.70
  max_homopolymer: 4

scoring:
  coherence_threshold: 0.135

chromatin:
  source: "PsychENCODE"  # or "ENCODE", "None"
  tissue: "human_cortex"
  cell_type: "neurons"

notes:
  disease: "Disease name"
  mode: "haploinsufficiency"
  therapeutic_goal: "Description of therapeutic goal"
```

### Step 3: Load and Use

```python
from phaselab.targets import load_target_config

# Load configuration
target = load_target_config("GENE_NAME")

# Access properties
print(f"Gene: {target.gene_symbol}")
print(f"Chromosome: {target.chrom}")
print(f"TSS: {target.tss_genomic}")
print(f"TSS index in promoter: {target.tss_index}")
print(f"CRISPRa window: {target.crispr_window}")
```

## Target Configuration Schema

### Required Fields

| Field | Type | Description |
|-------|------|-------------|
| `gene_symbol` | string | Gene symbol (e.g., "SCN2A") |
| `genome_build` | string | Genome assembly (e.g., "GRCh38") |
| `promoter.chrom` | string | Chromosome (e.g., "chr2") |
| `promoter.start` | int | Promoter start coordinate |
| `promoter.end` | int | Promoter end coordinate |
| `promoter.tss` | int | TSS genomic coordinate |

### Optional Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `gene_id` | string | None | Ensembl gene ID |
| `refseq_id` | string | None | RefSeq transcript ID |
| `crispr.system` | string | "SpCas9" | CRISPR system |
| `crispr.pam` | string | "NGG" | PAM sequence |
| `crispr.window.start_offset` | int | -400 | Window start (relative to TSS) |
| `crispr.window.end_offset` | int | -50 | Window end (relative to TSS) |
| `filters.gc_min` | float | 0.40 | Minimum GC content |
| `filters.gc_max` | float | 0.70 | Maximum GC content |
| `filters.max_homopolymer` | int | 4 | Max homopolymer run |
| `scoring.coherence_threshold` | float | 0.135 | IR coherence threshold |
| `chromatin.source` | string | None | Chromatin data source |

## API Reference

### `load_target_config(name: str) -> TargetConfig`

Load a target configuration from YAML.

```python
from phaselab.targets import load_target_config

scn2a = load_target_config("SCN2A")
rai1 = load_target_config("RAI1")
```

### `list_available_targets() -> List[str]`

List all available target configurations.

```python
from phaselab.targets import list_available_targets

targets = list_available_targets()
print(targets)  # ['RAI1', 'SCN2A']
```

### `TargetConfig`

Dataclass containing all target configuration.

```python
@dataclass
class TargetConfig:
    gene_symbol: str
    genome_build: str
    chrom: str
    promoter_start: int
    promoter_end: int
    tss_genomic: int
    crispr_pam: str
    crispr_window: Tuple[int, int]
    # ... and more

    @property
    def tss_index(self) -> int:
        """TSS position as 0-based index within promoter."""
        return self.tss_genomic - self.promoter_start
```

## Example: Full Workflow

```python
from phaselab.targets import load_target_config
from phaselab.crispr import design_guides, GuideDesignConfig

# 1. Load target
target = load_target_config("SCN2A")

# 2. Fetch promoter sequence (your implementation)
from your_genome_lib import fetch_sequence
promoter_seq = fetch_sequence(
    target.genome_build,
    target.chrom,
    target.promoter_start,
    target.promoter_end
)

# 3. Configure pipeline
config = GuideDesignConfig(
    pam=target.crispr_pam,
    crispr_window=target.crispr_window,
    min_gc=target.gc_min,
    max_gc=target.gc_max,
)

# 4. Design guides
guides = design_guides(
    sequence=promoter_seq,
    tss_index=target.tss_index,
    config=config
)

# 5. View results
print(guides[['sequence', 'position', 'coherence_R', 'go_no_go']])
```

## Supported Diseases

### Haploinsufficiency Disorders

PhaseLab is particularly suited for disorders where:
- One gene copy is non-functional
- Boosting expression from the remaining allele is therapeutic
- CRISPRa can provide 50-200% increase in expression

Current targets:
- **RAI1** → Smith-Magenis Syndrome
- **SCN2A** → Autism-linked NDD, epilepsy

### Future Targets

Candidates for future PhaseLab targets:
- **CHD8** → Autism (SINEUP upregulation validated)
- **ASH1L** → Autism-like phenotypes
- **SHANK3** → Phelan-McDermid Syndrome
- **NRXN1** → Autism, schizophrenia

---

*See individual target YAML files in `src/phaselab/targets/` for complete configurations.*
