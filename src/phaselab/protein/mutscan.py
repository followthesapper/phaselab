"""
PhaseLab Protein: Mutational scanning coherence analysis.

Applies the E213-validated spatial coherence methodology to:
- Deep mutational scanning (DMS) data
- Saturation mutagenesis
- Alanine scanning
- Evolutionary coupling analysis

Key insight:
- Residue position = perturbation coordinate
- Mutational effect (fitness, function) = response signal
- Spatial coherence identifies structurally/functionally coherent domains

This is the protein equivalent of E213 CRISPRa tiling analysis.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Union, Tuple
from pathlib import Path

from ..landscapes.core import (
    ResponseLandscape,
    CoherenceProfile,
    StabilityClass,
    RegionClassification,
)
from ..landscapes.coherence import compute_spatial_coherence
from ..landscapes.classification import classify_regions


@dataclass
class MutScanLandscape:
    """
    Mutational scanning landscape for a protein.

    Represents how mutation effects vary across protein sequence.

    Attributes:
        positions: Residue positions (1-indexed)
        effects: Mutational effects (fitness, function score)
        protein_id: Protein identifier
        protein_name: Protein name
        wildtype_sequence: Wild-type amino acid sequence
        mutations: Optional list of specific mutations
        effect_type: Type of effect measured
        metadata: Additional metadata
    """
    positions: np.ndarray
    effects: np.ndarray
    protein_id: str
    protein_name: Optional[str] = None
    wildtype_sequence: Optional[str] = None
    mutations: Optional[List[str]] = None
    effect_type: str = "fitness"
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_positions(self) -> int:
        """Number of positions."""
        return len(self.positions)

    @property
    def protein_length(self) -> Optional[int]:
        """Protein length from sequence."""
        if self.wildtype_sequence:
            return len(self.wildtype_sequence)
        return int(self.positions.max()) if len(self.positions) > 0 else None

    @property
    def mean_effect(self) -> float:
        """Mean mutational effect."""
        return float(np.nanmean(self.effects))

    @property
    def effect_variance(self) -> float:
        """Effect variance."""
        return float(np.nanvar(self.effects))

    def to_response_landscape(self) -> ResponseLandscape:
        """Convert to generic ResponseLandscape."""
        return ResponseLandscape(
            coords=self.positions,
            responses=self.effects,
            coord_labels=self.mutations,
            metadata={
                'type': 'mutscan',
                'protein_id': self.protein_id,
                'protein_name': self.protein_name,
                'effect_type': self.effect_type,
                **self.metadata,
            },
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'protein_id': self.protein_id,
            'protein_name': self.protein_name,
            'n_positions': self.n_positions,
            'protein_length': self.protein_length,
            'mean_effect': self.mean_effect,
            'effect_type': self.effect_type,
        }


@dataclass
class FunctionalDomain:
    """
    A functional domain identified from mutational scanning.

    Attributes:
        start: Start residue
        end: End residue
        stability: Stability classification
        coherence: Coherence score
        mean_effect: Mean mutational effect
        functional_class: Classification of function
        n_positions: Number of positions
    """
    start: int
    end: int
    stability: StabilityClass
    coherence: float
    mean_effect: float
    functional_class: str  # "essential", "neutral", "beneficial"
    n_positions: int

    @property
    def is_essential(self) -> bool:
        """Whether domain is functionally essential."""
        return self.functional_class == "essential"

    @property
    def is_stable(self) -> bool:
        """Whether domain shows stable behavior."""
        return self.stability == StabilityClass.STABLE

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'start': self.start,
            'end': self.end,
            'stability': self.stability.value,
            'coherence': self.coherence,
            'mean_effect': self.mean_effect,
            'functional_class': self.functional_class,
            'n_positions': self.n_positions,
        }


@dataclass
class MutScanResult:
    """
    Result of mutational scanning coherence analysis.

    Attributes:
        landscape: Original mutational scanning landscape
        profile: Coherence profile
        functional_domains: Identified functional domains
        essential_regions: High-confidence essential regions
        variable_regions: Regions tolerant to mutation
        validation: Validation statistics
    """
    landscape: MutScanLandscape
    profile: CoherenceProfile
    functional_domains: List[FunctionalDomain]
    essential_regions: List[FunctionalDomain]
    variable_regions: List[FunctionalDomain]
    validation: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_validated(self) -> bool:
        """Whether coherence model is validated."""
        return self.profile.is_validated

    @property
    def correlation(self) -> float:
        """Coherence-variance correlation."""
        return self.profile.correlation

    def summary(self) -> str:
        """Generate summary string."""
        lines = [
            "=" * 60,
            f"MUTATIONAL SCANNING: {self.landscape.protein_name or self.landscape.protein_id}",
            "=" * 60,
            "",
            f"Positions analyzed: {self.landscape.n_positions}",
            f"Protein length: {self.landscape.protein_length or 'unknown'}",
            f"Effect type: {self.landscape.effect_type}",
            "",
            "COHERENCE:",
            f"  Correlation: {self.correlation:.3f}",
            f"  Validated: {'YES' if self.is_validated else 'NO'}",
            f"  Variance reduction: {self.profile.variance_reduction_estimate:.1%}",
            "",
            "DOMAINS:",
            f"  Functional domains: {len(self.functional_domains)}",
            f"  Essential regions: {len(self.essential_regions)}",
            f"  Variable regions: {len(self.variable_regions)}",
        ]

        if self.essential_regions:
            lines.extend(["", "ESSENTIAL REGIONS:"])
            for region in self.essential_regions[:5]:
                lines.append(
                    f"  [{region.start}-{region.end}] "
                    f"coh={region.coherence:.3f} "
                    f"effect={region.mean_effect:.3f}"
                )

        lines.append("=" * 60)
        return "\n".join(lines)


def load_mutscan_data(
    filepath: Union[str, Path],
    protein_id: str,
    position_col: str = "position",
    effect_col: str = "effect",
    mutation_col: Optional[str] = "mutation",
    delimiter: str = "\t",
    protein_name: Optional[str] = None,
    effect_type: str = "fitness",
) -> MutScanLandscape:
    """
    Load mutational scanning data from a file.

    Supports common DMS output formats.

    Args:
        filepath: Path to data file.
        protein_id: Protein identifier.
        position_col: Position column name.
        effect_col: Effect column name.
        mutation_col: Mutation label column name.
        delimiter: Column delimiter.
        protein_name: Protein name.
        effect_type: Type of effect measured.

    Returns:
        MutScanLandscape for analysis.

    Example:
        >>> landscape = load_mutscan_data(
        ...     'dms_results.tsv',
        ...     protein_id='P53',
        ...     protein_name='p53',
        ...     effect_type='fitness',
        ... )
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"MutScan file not found: {filepath}")

    try:
        import pandas as pd
        df = pd.read_csv(filepath, delimiter=delimiter)

        positions = df[position_col].values
        effects = df[effect_col].values

        mutations = None
        if mutation_col and mutation_col in df.columns:
            mutations = df[mutation_col].tolist()

        return MutScanLandscape(
            positions=positions.astype(float),
            effects=effects.astype(float),
            protein_id=protein_id,
            protein_name=protein_name,
            mutations=mutations,
            effect_type=effect_type,
            metadata={'source_file': str(filepath)},
        )

    except ImportError:
        raise ImportError("pandas required for mutscan loading")


def analyze_mutscan_coherence(
    landscape: MutScanLandscape,
    window: int = 10,
    stable_threshold: float = 0.7,
    essential_threshold: float = -1.0,
    beneficial_threshold: float = 0.5,
) -> MutScanResult:
    """
    Analyze spatial coherence in mutational scanning data.

    Identifies:
    - Functional domains with coherent mutational effects
    - Essential regions (stable + strong negative effect)
    - Variable regions (stable + neutral/positive effect)

    Args:
        landscape: MutScan landscape to analyze.
        window: Window size for coherence.
        stable_threshold: Threshold for stable regions.
        essential_threshold: Effect threshold for essential.
        beneficial_threshold: Effect threshold for beneficial.

    Returns:
        MutScanResult with complete analysis.

    Example:
        >>> result = analyze_mutscan_coherence(landscape)
        >>> essential = result.essential_regions
        >>> print(f"Found {len(essential)} essential domains")
    """
    response_landscape = landscape.to_response_landscape()

    # Compute coherence
    profile = compute_spatial_coherence(
        response_landscape,
        window=window,
    )

    # Classify regions
    classification = classify_regions(
        response_landscape,
        profile=profile,
        stable_threshold=stable_threshold,
    )

    # Build functional domains
    functional_domains = []
    essential_regions = []
    variable_regions = []

    for region in classification.regions:
        start, end, stability, score = region

        # Get effects in this region
        mask = (landscape.positions >= start) & (landscape.positions <= end)
        region_effects = landscape.effects[mask]

        if len(region_effects) == 0:
            continue

        mean_effect = float(np.nanmean(region_effects))

        # Classify functional impact
        if mean_effect < essential_threshold:
            functional_class = "essential"
        elif mean_effect > beneficial_threshold:
            functional_class = "beneficial"
        else:
            functional_class = "neutral"

        domain = FunctionalDomain(
            start=int(start),
            end=int(end),
            stability=stability,
            coherence=float(score),
            mean_effect=mean_effect,
            functional_class=functional_class,
            n_positions=int(np.sum(mask)),
        )

        functional_domains.append(domain)

        # Categorize based on stability and function
        if stability == StabilityClass.STABLE:
            if functional_class == "essential":
                essential_regions.append(domain)
            elif functional_class in ("neutral", "beneficial"):
                variable_regions.append(domain)

    # Sort by coherence
    functional_domains.sort(key=lambda d: d.coherence, reverse=True)
    essential_regions.sort(key=lambda d: d.mean_effect)  # Most essential first
    variable_regions.sort(key=lambda d: d.coherence, reverse=True)

    validation = {
        'coherence_variance_correlation': profile.correlation,
        'p_value': profile.p_value,
        'is_validated': profile.is_validated,
        'variance_reduction': profile.variance_reduction_estimate,
    }

    return MutScanResult(
        landscape=landscape,
        profile=profile,
        functional_domains=functional_domains,
        essential_regions=essential_regions,
        variable_regions=variable_regions,
        validation=validation,
    )


def local_coherence_profile(
    landscape: MutScanLandscape,
    window: int = 5,
) -> Dict[str, np.ndarray]:
    """
    Compute per-residue local coherence.

    Returns a coherence value for each position, useful for:
    - Mapping coherence onto 3D structure
    - Identifying position-specific reliability
    - Comparing to structural features

    Args:
        landscape: MutScan landscape.
        window: Window size for local coherence.

    Returns:
        Dictionary with position and coherence arrays.
    """
    response_landscape = landscape.to_response_landscape()
    profile = compute_spatial_coherence(response_landscape, window=window)

    return {
        'positions': profile.coords,
        'coherence': profile.coherence,
        'local_variance': profile.local_variance,
    }


def map_coherence_to_structure(
    result: MutScanResult,
    pdb_path: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Map coherence values to protein structure.

    Creates a mapping suitable for visualization in molecular viewers.

    Args:
        result: MutScanResult from analyze_mutscan_coherence.
        pdb_path: Optional PDB file for B-factor replacement.

    Returns:
        Dictionary with structural mapping information.
    """
    # Create per-residue coherence mapping
    positions = result.landscape.positions
    coherence_values = result.profile.coherence
    coh_positions = result.profile.coords

    # Interpolate coherence to all positions
    residue_coherence = {}
    for pos in positions:
        idx = np.argmin(np.abs(coh_positions - pos))
        residue_coherence[int(pos)] = float(coherence_values[idx])

    # Domain annotations
    domain_annotations = []
    for domain in result.functional_domains:
        domain_annotations.append({
            'start': domain.start,
            'end': domain.end,
            'type': domain.functional_class,
            'stability': domain.stability.value,
            'color_by_stability': (
                'green' if domain.stability == StabilityClass.STABLE
                else 'red' if domain.stability == StabilityClass.AMPLIFYING
                else 'yellow'
            ),
        })

    output = {
        'residue_coherence': residue_coherence,
        'domain_annotations': domain_annotations,
        'n_residues': len(residue_coherence),
    }

    # Optionally write B-factor replaced PDB
    if pdb_path:
        try:
            bfactor_pdb = _write_bfactor_pdb(pdb_path, residue_coherence)
            output['bfactor_pdb'] = bfactor_pdb
        except Exception as e:
            output['bfactor_pdb_error'] = str(e)

    return output


def _write_bfactor_pdb(pdb_path: str, residue_coherence: Dict[int, float]) -> str:
    """Write PDB with coherence as B-factors."""
    try:
        from Bio.PDB import PDBParser, PDBIO
    except ImportError:
        raise ImportError("BioPython required for PDB manipulation")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    # Replace B-factors with coherence values
    for model in structure:
        for chain in model:
            for residue in chain:
                res_num = residue.get_id()[1]
                coherence = residue_coherence.get(res_num, 0.5)
                # Scale to 0-100 for B-factor range
                scaled_coherence = coherence * 100

                for atom in residue:
                    atom.bfactor = scaled_coherence

    # Write modified PDB
    output_path = pdb_path.replace('.pdb', '_coherence.pdb')
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path)

    return output_path


def compare_mutscan_datasets(
    landscape_a: MutScanLandscape,
    landscape_b: MutScanLandscape,
    window: int = 10,
) -> Dict[str, Any]:
    """
    Compare two mutational scanning datasets.

    Useful for:
    - Comparing conditions (e.g., +/- drug)
    - Comparing experimental repeats
    - Cross-validating DMS studies

    Args:
        landscape_a: First MutScan landscape.
        landscape_b: Second MutScan landscape.
        window: Coherence window size.

    Returns:
        Comparison results dictionary.
    """
    result_a = analyze_mutscan_coherence(landscape_a, window=window)
    result_b = analyze_mutscan_coherence(landscape_b, window=window)

    # Find common positions
    common_positions = np.intersect1d(landscape_a.positions, landscape_b.positions)

    mask_a = np.isin(landscape_a.positions, common_positions)
    mask_b = np.isin(landscape_b.positions, common_positions)

    if len(common_positions) > 2:
        effect_correlation = np.corrcoef(
            landscape_a.effects[mask_a],
            landscape_b.effects[mask_b],
        )[0, 1]
    else:
        effect_correlation = np.nan

    return {
        'comparison': {
            'n_common_positions': len(common_positions),
            'effect_correlation': float(effect_correlation),
            'coherence_a': result_a.correlation,
            'coherence_b': result_b.correlation,
            'essential_domains_a': len(result_a.essential_regions),
            'essential_domains_b': len(result_b.essential_regions),
        },
        'result_a': result_a,
        'result_b': result_b,
    }
