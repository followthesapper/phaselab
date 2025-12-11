"""
Target configuration loader for PhaseLab.

Loads gene target configurations from YAML files and provides
a unified interface for accessing genomic coordinates and parameters.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import yaml


@dataclass
class TargetConfig:
    """
    Configuration for a gene target.

    Contains all information needed to run a PhaseLab CRISPRa
    experiment on a specific gene.
    """

    # Gene identifiers
    gene_symbol: str
    genome_build: str
    gene_id: Optional[str] = None
    refseq_id: Optional[str] = None
    uniprot_id: Optional[str] = None

    # Genomic coordinates
    chrom: str = ""
    promoter_start: int = 0
    promoter_end: int = 0
    tss_genomic: int = 0

    # CRISPR settings
    crispr_system: str = "SpCas9"
    crispr_pam: str = "NGG"
    crispr_window: Tuple[int, int] = (-400, -50)

    # Filtering thresholds
    gc_min: float = 0.40
    gc_max: float = 0.70
    max_homopolymer: int = 4
    min_complexity: float = 0.5

    # Scoring settings
    coherence_threshold: float = 0.135  # e^-2

    # Chromatin metadata
    chromatin_source: Optional[str] = None
    chromatin_tissue: Optional[str] = None
    chromatin_cell_type: Optional[str] = None
    chromatin_tracks: List[str] = field(default_factory=list)

    # Notes and metadata
    disease: Optional[str] = None
    omim: Optional[str] = None
    mode: Optional[str] = None
    therapeutic_goal: Optional[str] = None
    references: List[Dict[str, Any]] = field(default_factory=list)
    hardware_validation: Optional[Dict[str, Any]] = None

    @property
    def promoter_length(self) -> int:
        """Length of the promoter region in base pairs."""
        return self.promoter_end - self.promoter_start

    @property
    def tss_index(self) -> int:
        """TSS position as 0-based index within promoter sequence."""
        return self.tss_genomic - self.promoter_start

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "gene_symbol": self.gene_symbol,
            "genome_build": self.genome_build,
            "gene_id": self.gene_id,
            "refseq_id": self.refseq_id,
            "chrom": self.chrom,
            "promoter_start": self.promoter_start,
            "promoter_end": self.promoter_end,
            "tss_genomic": self.tss_genomic,
            "tss_index": self.tss_index,
            "crispr_system": self.crispr_system,
            "crispr_pam": self.crispr_pam,
            "crispr_window": self.crispr_window,
            "gc_min": self.gc_min,
            "gc_max": self.gc_max,
            "max_homopolymer": self.max_homopolymer,
            "coherence_threshold": self.coherence_threshold,
            "chromatin_source": self.chromatin_source,
            "disease": self.disease,
            "therapeutic_goal": self.therapeutic_goal,
        }


def get_target_path(name: str) -> Path:
    """
    Get the path to a target YAML file.

    Args:
        name: Target name (e.g., 'SCN2A', 'RAI1')

    Returns:
        Path to the YAML file.

    Raises:
        FileNotFoundError: If target file doesn't exist.
    """
    # Look in the package targets directory
    package_dir = Path(__file__).parent
    yaml_path = package_dir / f"{name}.yaml"

    if yaml_path.exists():
        return yaml_path

    # Also check without extension
    if not name.endswith(".yaml"):
        yaml_path = package_dir / f"{name}.yaml"
        if yaml_path.exists():
            return yaml_path

    raise FileNotFoundError(f"Target configuration not found: {name}")


def list_available_targets() -> List[str]:
    """
    List all available target configurations.

    Returns:
        List of target names (without .yaml extension).
    """
    package_dir = Path(__file__).parent
    yaml_files = list(package_dir.glob("*.yaml"))
    return sorted([f.stem for f in yaml_files])


def load_target_config(name: str) -> TargetConfig:
    """
    Load a target configuration from YAML.

    Args:
        name: Target name (e.g., 'SCN2A', 'RAI1')

    Returns:
        TargetConfig with all settings loaded.

    Raises:
        FileNotFoundError: If target file doesn't exist.
        ValueError: If YAML is malformed.

    Example:
        >>> from phaselab.targets import load_target_config
        >>> scn2a = load_target_config("SCN2A")
        >>> print(scn2a.gene_symbol, scn2a.chrom, scn2a.tss_genomic)
        SCN2A chr2 165239414
    """
    yaml_path = get_target_path(name)

    with open(yaml_path, "r") as f:
        cfg = yaml.safe_load(f)

    if cfg is None:
        raise ValueError(f"Empty or invalid YAML file: {yaml_path}")

    # Extract nested values
    promoter = cfg.get("promoter", {})
    crispr = cfg.get("crispr", {})
    filters = cfg.get("filters", {})
    scoring = cfg.get("scoring", {})
    chromatin = cfg.get("chromatin", {})
    notes = cfg.get("notes", {})

    # Build window tuple
    window_cfg = crispr.get("window", {})
    crispr_window = (
        window_cfg.get("start_offset", -400),
        window_cfg.get("end_offset", -50),
    )

    return TargetConfig(
        # Gene identifiers
        gene_symbol=cfg.get("gene_symbol", name),
        genome_build=cfg.get("genome_build", "GRCh38"),
        gene_id=cfg.get("gene_id"),
        refseq_id=cfg.get("refseq_id"),
        uniprot_id=cfg.get("uniprot_id"),
        # Genomic coordinates
        chrom=promoter.get("chrom", ""),
        promoter_start=promoter.get("start", 0),
        promoter_end=promoter.get("end", 0),
        tss_genomic=promoter.get("tss", 0),
        # CRISPR settings
        crispr_system=crispr.get("system", "SpCas9"),
        crispr_pam=crispr.get("pam", "NGG"),
        crispr_window=crispr_window,
        # Filtering thresholds
        gc_min=filters.get("gc_min", 0.40),
        gc_max=filters.get("gc_max", 0.70),
        max_homopolymer=filters.get("max_homopolymer", 4),
        min_complexity=filters.get("min_complexity", 0.5),
        # Scoring
        coherence_threshold=scoring.get("coherence_threshold", 0.135),
        # Chromatin
        chromatin_source=chromatin.get("source"),
        chromatin_tissue=chromatin.get("tissue"),
        chromatin_cell_type=chromatin.get("cell_type"),
        chromatin_tracks=chromatin.get("tracks", []),
        # Notes
        disease=notes.get("disease"),
        omim=notes.get("omim"),
        mode=notes.get("mode"),
        therapeutic_goal=notes.get("therapeutic_goal"),
        references=notes.get("references", []),
        hardware_validation=notes.get("hardware_validation"),
    )
