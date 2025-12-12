"""
PhaseLab Methylation: CpG methylation modeling for CRISPR efficiency.

Implements:
- CpG island prediction and scoring
- Methylation-dependent efficiency suppression
- Promoter methylation state modeling
- Integration with CRISPRa/CRISPRi predictions

References:
- C-RNNCrispr: DNA methylation affects sgRNA efficiency
- ENCODE methylation data (RRBS, WGBS)
- CpGIMethPred: CpG island methylation prediction
- Research shows CRISPRa efficiency drops in methylated promoters

Version: 0.5.0
"""

import numpy as np
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass, field
import re
import logging

logger = logging.getLogger(__name__)


# CpG density thresholds based on Gardiner-Garden & Frommer criteria
CPG_ISLAND_CRITERIA = {
    'min_length': 200,          # Minimum 200bp
    'min_gc_content': 0.50,     # ≥50% GC
    'min_cpg_ratio': 0.60,      # Observed/Expected CpG ≥0.6
}

# Methylation impact on CRISPR efficiency
# Based on C-RNNCrispr and experimental data
METHYLATION_EFFICIENCY_FACTORS = {
    'unmethylated': 1.0,        # Full efficiency
    'low': 0.85,                # Minor impact
    'moderate': 0.60,           # Significant reduction
    'high': 0.35,               # Major reduction
    'hypermethylated': 0.15,    # Severe suppression
}


@dataclass
class CpGIsland:
    """CpG island annotation."""
    chrom: str
    start: int
    end: int
    gc_content: float
    cpg_count: int
    cpg_ratio: float  # Observed/Expected
    length: int

    @property
    def cpg_density(self) -> float:
        """CpGs per 100bp."""
        return (self.cpg_count / self.length) * 100 if self.length > 0 else 0


@dataclass
class MethylationScore:
    """Methylation assessment for a genomic region."""
    position: int
    chrom: str

    # Methylation state
    methylation_level: float  # 0-1 (0=unmethylated, 1=fully methylated)
    methylation_class: str    # unmethylated, low, moderate, high, hypermethylated

    # CpG context
    in_cpg_island: bool
    cpg_density: float        # CpGs per 100bp in region
    cpg_count: int            # Number of CpGs in region

    # Impact on CRISPR
    efficiency_factor: float  # Multiplier for CRISPR efficiency
    suppression_risk: str     # LOW, MODERATE, HIGH

    # Data source
    data_source: str = "predicted"  # predicted, ENCODE, RRBS, WGBS


@dataclass
class MethylationConfig:
    """Configuration for methylation analysis."""

    # Analysis window
    window_size: int = 500

    # CpG island detection
    detect_islands: bool = True
    island_criteria: Dict = field(default_factory=lambda: CPG_ISLAND_CRITERIA.copy())

    # Methylation estimation (when no experimental data)
    use_sequence_prediction: bool = True
    promoter_default_methylation: float = 0.2  # Promoters often unmethylated
    gene_body_default_methylation: float = 0.6  # Gene bodies often methylated

    # Efficiency impact
    apply_efficiency_penalty: bool = True


def count_cpg(sequence: str) -> int:
    """
    Count CpG dinucleotides in a sequence.

    Args:
        sequence: DNA sequence

    Returns:
        Number of CpG sites
    """
    sequence = sequence.upper()
    return len(re.findall(r'CG', sequence))


def calculate_cpg_ratio(sequence: str) -> float:
    """
    Calculate observed/expected CpG ratio.

    Expected CpG = (C count × G count) / length
    Observed CpG = actual CpG count

    Args:
        sequence: DNA sequence

    Returns:
        CpG ratio (O/E)
    """
    sequence = sequence.upper()
    length = len(sequence)

    if length == 0:
        return 0.0

    c_count = sequence.count('C')
    g_count = sequence.count('G')
    cpg_count = count_cpg(sequence)

    # Expected CpG under random model
    expected = (c_count * g_count) / length if length > 0 else 0

    if expected == 0:
        return 0.0

    return cpg_count / expected


def gc_content(sequence: str) -> float:
    """Calculate GC content of sequence."""
    sequence = sequence.upper()
    if len(sequence) == 0:
        return 0.0
    gc = sequence.count('G') + sequence.count('C')
    return gc / len(sequence)


def is_cpg_island(
    sequence: str,
    criteria: Optional[Dict] = None
) -> Tuple[bool, Dict]:
    """
    Determine if sequence qualifies as a CpG island.

    Uses Gardiner-Garden & Frommer criteria:
    - Length ≥200bp
    - GC content ≥50%
    - Observed/Expected CpG ≥0.6

    Args:
        sequence: DNA sequence to analyze
        criteria: Optional custom criteria

    Returns:
        (is_island, metrics_dict)
    """
    if criteria is None:
        criteria = CPG_ISLAND_CRITERIA

    sequence = sequence.upper()
    length = len(sequence)

    metrics = {
        'length': length,
        'gc_content': gc_content(sequence),
        'cpg_count': count_cpg(sequence),
        'cpg_ratio': calculate_cpg_ratio(sequence),
    }

    is_island = (
        length >= criteria['min_length'] and
        metrics['gc_content'] >= criteria['min_gc_content'] and
        metrics['cpg_ratio'] >= criteria['min_cpg_ratio']
    )

    return is_island, metrics


def find_cpg_islands(
    sequence: str,
    min_length: int = 200,
    step: int = 50,
    chrom: str = "unknown"
) -> List[CpGIsland]:
    """
    Scan sequence for CpG islands using sliding window.

    Args:
        sequence: DNA sequence
        min_length: Minimum island length
        step: Sliding window step size
        chrom: Chromosome name for annotation

    Returns:
        List of CpGIsland objects
    """
    sequence = sequence.upper()
    islands = []
    length = len(sequence)

    if length < min_length:
        return islands

    # Sliding window scan
    i = 0
    while i < length - min_length:
        window = sequence[i:i + min_length]
        is_island, metrics = is_cpg_island(window)

        if is_island:
            # Extend island as far as possible
            end = i + min_length
            while end < length:
                extended = sequence[i:end + step]
                still_island, _ = is_cpg_island(extended)
                if still_island:
                    end += step
                else:
                    break

            # Create island object
            island_seq = sequence[i:end]
            island = CpGIsland(
                chrom=chrom,
                start=i,
                end=end,
                gc_content=gc_content(island_seq),
                cpg_count=count_cpg(island_seq),
                cpg_ratio=calculate_cpg_ratio(island_seq),
                length=end - i,
            )
            islands.append(island)

            # Skip past this island
            i = end
        else:
            i += step

    return islands


def predict_methylation_from_sequence(
    sequence: str,
    position_type: str = "unknown"
) -> float:
    """
    Predict methylation level from sequence features.

    This is a simplified model for when experimental data isn't available.
    Uses CpG density and context as proxies.

    Args:
        sequence: DNA sequence
        position_type: "promoter", "gene_body", "intergenic", "unknown"

    Returns:
        Predicted methylation level 0-1
    """
    sequence = sequence.upper()

    # CpG island regions are typically unmethylated
    is_island, metrics = is_cpg_island(sequence)

    if is_island:
        # CpG islands at promoters are usually unmethylated
        return 0.15

    cpg_density = (count_cpg(sequence) / len(sequence)) * 100 if sequence else 0

    # Base prediction on position type
    if position_type == "promoter":
        # Promoters tend to be unmethylated (active genes)
        base = 0.20
    elif position_type == "gene_body":
        # Gene bodies tend to be methylated
        base = 0.65
    elif position_type == "intergenic":
        # Intergenic regions variable
        base = 0.50
    else:
        base = 0.40

    # Adjust based on CpG density
    # High CpG density often correlates with lower methylation
    if cpg_density > 5:  # High density
        base *= 0.5
    elif cpg_density > 2:  # Moderate density
        base *= 0.8

    return min(1.0, max(0.0, base))


def classify_methylation(level: float) -> str:
    """
    Classify methylation level into category.

    Args:
        level: Methylation level 0-1

    Returns:
        Category string
    """
    if level < 0.1:
        return "unmethylated"
    elif level < 0.3:
        return "low"
    elif level < 0.6:
        return "moderate"
    elif level < 0.85:
        return "high"
    else:
        return "hypermethylated"


def methylation_efficiency_factor(level: float) -> float:
    """
    Calculate CRISPR efficiency factor from methylation level.

    Based on experimental observations that CRISPRa/CRISPRi
    efficiency decreases in methylated regions.

    Args:
        level: Methylation level 0-1

    Returns:
        Efficiency multiplier 0-1
    """
    # Logistic decay model
    # Efficiency drops significantly above ~40% methylation
    return 1.0 / (1.0 + np.exp(5 * (level - 0.5)))


def score_methylation(
    sequence: str,
    chrom: str = "unknown",
    position: int = 0,
    position_type: str = "unknown",
    experimental_methylation: Optional[float] = None,
    config: Optional[MethylationConfig] = None,
) -> MethylationScore:
    """
    Compute methylation score for a genomic region.

    Args:
        sequence: DNA sequence of the region
        chrom: Chromosome
        position: Genomic position
        position_type: "promoter", "gene_body", "intergenic"
        experimental_methylation: Optional measured methylation level
        config: MethylationConfig

    Returns:
        MethylationScore with assessment

    Example:
        >>> score = score_methylation(
        ...     sequence=promoter_seq,
        ...     position_type="promoter"
        ... )
        >>> print(f"Efficiency factor: {score.efficiency_factor:.2f}")
    """
    if config is None:
        config = MethylationConfig()

    sequence = sequence.upper()

    # Determine methylation level
    if experimental_methylation is not None:
        meth_level = experimental_methylation
        data_source = "experimental"
    elif config.use_sequence_prediction:
        meth_level = predict_methylation_from_sequence(sequence, position_type)
        data_source = "predicted"
    else:
        meth_level = 0.3  # Default moderate
        data_source = "default"

    # CpG analysis
    is_island, metrics = is_cpg_island(sequence)
    cpg_count = count_cpg(sequence)
    cpg_density = (cpg_count / len(sequence)) * 100 if sequence else 0

    # Classification
    meth_class = classify_methylation(meth_level)
    eff_factor = methylation_efficiency_factor(meth_level)

    # Suppression risk
    if eff_factor >= 0.8:
        risk = "LOW"
    elif eff_factor >= 0.5:
        risk = "MODERATE"
    else:
        risk = "HIGH"

    return MethylationScore(
        position=position,
        chrom=chrom,
        methylation_level=meth_level,
        methylation_class=meth_class,
        in_cpg_island=is_island,
        cpg_density=cpg_density,
        cpg_count=cpg_count,
        efficiency_factor=eff_factor,
        suppression_risk=risk,
        data_source=data_source,
    )


def adjust_crispra_efficiency(
    base_efficiency: float,
    methylation_score: MethylationScore
) -> Dict[str, float]:
    """
    Adjust CRISPRa efficiency prediction based on methylation.

    Args:
        base_efficiency: Predicted efficiency without methylation
        methylation_score: MethylationScore for the target region

    Returns:
        Dictionary with adjusted efficiency and factors

    Example:
        >>> meth = score_methylation(promoter_seq, position_type="promoter")
        >>> result = adjust_crispra_efficiency(0.85, meth)
        >>> print(f"Adjusted efficiency: {result['adjusted_efficiency']:.2f}")
    """
    adjusted = base_efficiency * methylation_score.efficiency_factor

    return {
        'base_efficiency': base_efficiency,
        'adjusted_efficiency': adjusted,
        'methylation_factor': methylation_score.efficiency_factor,
        'methylation_level': methylation_score.methylation_level,
        'methylation_class': methylation_score.methylation_class,
        'suppression_risk': methylation_score.suppression_risk,
        'recommendation': _get_methylation_recommendation(methylation_score),
    }


def _get_methylation_recommendation(score: MethylationScore) -> str:
    """Generate recommendation based on methylation score."""
    if score.suppression_risk == "LOW":
        return "PROCEED: Low methylation, CRISPRa should work well."
    elif score.suppression_risk == "MODERATE":
        return (
            "CAUTION: Moderate methylation may reduce efficiency. "
            "Consider demethylation or alternative targets."
        )
    else:
        return (
            "WARNING: High methylation likely to suppress CRISPRa. "
            "Recommend: (1) Use DNMT inhibitors, (2) Target enhancer, "
            "or (3) Choose different guide location."
        )


@dataclass
class MethylationProfile:
    """
    Methylation profile for a genomic region or gene.

    Can be loaded from ENCODE/Roadmap bisulfite sequencing data.
    """

    name: str
    tissue: str
    cell_type: Optional[str] = None

    # Methylation data: position -> beta value (0-1)
    _methylation_data: Dict[str, Dict[int, float]] = field(
        default_factory=dict, repr=False
    )

    # CpG islands in this profile
    cpg_islands: List[CpGIsland] = field(default_factory=list)

    # Data source
    data_source: str = "unknown"  # ENCODE, Roadmap, RRBS, WGBS

    def add_methylation_data(
        self,
        chrom: str,
        position: int,
        beta_value: float
    ):
        """Add methylation measurement at a position."""
        if chrom not in self._methylation_data:
            self._methylation_data[chrom] = {}
        self._methylation_data[chrom][position] = beta_value

    def get_methylation(
        self,
        chrom: str,
        position: int,
        window: int = 100
    ) -> Optional[float]:
        """
        Get methylation level at a position.

        Args:
            chrom: Chromosome
            position: Genomic position
            window: Window to average over

        Returns:
            Methylation level or None if no data
        """
        if chrom not in self._methylation_data:
            return None

        data = self._methylation_data[chrom]

        # Find values in window
        values = []
        for pos, beta in data.items():
            if abs(pos - position) <= window // 2:
                values.append(beta)

        if not values:
            return None

        return float(np.mean(values))

    def load_bedmethyl(self, path: str) -> int:
        """
        Load methylation data from bedMethyl format.

        bedMethyl format (ENCODE standard):
        chrom, start, end, name, score, strand, thickStart, thickEnd,
        itemRgb, coverage, percentage

        Args:
            path: Path to bedMethyl file

        Returns:
            Number of CpG sites loaded
        """
        count = 0
        try:
            with open(path, 'r') as f:
                for line in f:
                    if line.startswith('#') or line.startswith('track'):
                        continue
                    fields = line.strip().split('\t')
                    if len(fields) >= 11:
                        chrom = fields[0]
                        position = int(fields[1])
                        percentage = float(fields[10])
                        beta = percentage / 100.0
                        self.add_methylation_data(chrom, position, beta)
                        count += 1
        except Exception as e:
            logger.error(f"Error loading bedMethyl: {e}")

        logger.info(f"Loaded {count} methylation sites from {path}")
        return count


# Pre-computed methylation profiles for common genes
GENE_METHYLATION_DEFAULTS = {
    # Haploinsufficiency genes - typically unmethylated promoters
    'RAI1': {'promoter': 0.15, 'gene_body': 0.60},
    'SCN2A': {'promoter': 0.12, 'gene_body': 0.55},
    'SHANK3': {'promoter': 0.18, 'gene_body': 0.58},
    'CHD8': {'promoter': 0.10, 'gene_body': 0.52},
    'MECP2': {'promoter': 0.08, 'gene_body': 0.65},

    # Tumor suppressors - can be hypermethylated in cancer
    'TP53': {'promoter': 0.10, 'gene_body': 0.50},
    'BRCA1': {'promoter': 0.12, 'gene_body': 0.55},

    # Default for unknown genes
    'DEFAULT': {'promoter': 0.20, 'gene_body': 0.55},
}


def get_gene_methylation(
    gene: str,
    region: str = "promoter"
) -> float:
    """
    Get default methylation level for a gene region.

    Args:
        gene: Gene symbol
        region: "promoter" or "gene_body"

    Returns:
        Default methylation level
    """
    gene = gene.upper()
    if gene in GENE_METHYLATION_DEFAULTS:
        return GENE_METHYLATION_DEFAULTS[gene].get(region, 0.3)
    return GENE_METHYLATION_DEFAULTS['DEFAULT'].get(region, 0.3)


def score_guide_methylation(
    guide_sequence: str,
    target_chrom: str,
    target_position: int,
    gene: Optional[str] = None,
    position_type: str = "promoter",
    methylation_profile: Optional[MethylationProfile] = None,
) -> Dict[str, any]:
    """
    Score CRISPR guide considering methylation context.

    This is the main interface for integrating methylation
    into guide scoring pipelines.

    Args:
        guide_sequence: Guide RNA sequence
        target_chrom: Chromosome
        target_position: Target genomic position
        gene: Optional gene symbol for defaults
        position_type: "promoter", "gene_body", etc.
        methylation_profile: Optional experimental data

    Returns:
        Dictionary with methylation assessment

    Example:
        >>> result = score_guide_methylation(
        ...     "GCGACTGCTACATAGCCAGG",
        ...     "chr17", 17681458,
        ...     gene="RAI1",
        ...     position_type="promoter"
        ... )
        >>> print(f"Risk: {result['suppression_risk']}")
    """
    # Get methylation level
    if methylation_profile:
        meth_level = methylation_profile.get_methylation(
            target_chrom, target_position
        )
        data_source = methylation_profile.data_source
    else:
        meth_level = None
        data_source = "predicted"

    # Fall back to gene defaults or sequence prediction
    if meth_level is None:
        if gene:
            meth_level = get_gene_methylation(gene, position_type)
        else:
            meth_level = predict_methylation_from_sequence(
                guide_sequence, position_type
            )

    # Score the methylation
    meth_score = score_methylation(
        sequence=guide_sequence,
        chrom=target_chrom,
        position=target_position,
        position_type=position_type,
        experimental_methylation=meth_level,
    )

    return {
        'methylation_level': meth_score.methylation_level,
        'methylation_class': meth_score.methylation_class,
        'efficiency_factor': meth_score.efficiency_factor,
        'suppression_risk': meth_score.suppression_risk,
        'in_cpg_island': meth_score.in_cpg_island,
        'cpg_density': meth_score.cpg_density,
        'data_source': data_source,
        'gene': gene,
        'recommendation': _get_methylation_recommendation(meth_score),
    }
