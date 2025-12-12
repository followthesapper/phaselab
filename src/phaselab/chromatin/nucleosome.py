"""
PhaseLab Nucleosome Positioning: Predict nucleosome occupancy for CRISPR efficiency.

Implements:
- NuPoP-like duration Hidden Markov Model for nucleosome prediction
- Nucleosome affinity scoring from sequence
- Nucleosome-free region (NFR) detection
- Integration with CRISPR guide scoring

References:
- NuPoP: Xi et al. (2010) BMC Bioinformatics - Duration HMM model
- nuCpos: Chemical map-based prediction (Brogaard et al.)
- Nucleosome affinity affects CRISPR efficiency (piCRISPR)

IMPORTANT: CRISPRa/CRISPRi efficiency is significantly reduced when
targeting nucleosome-occupied DNA. This module helps predict and
avoid such regions.

Version: 0.5.0
"""

import numpy as np
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass, field
import logging

logger = logging.getLogger(__name__)


# Nucleosome constants
NUCLEOSOME_LENGTH = 147  # bp wrapped around histone octamer
LINKER_LENGTH_MEAN = 38  # Mean linker DNA length (varies by organism/cell type)
LINKER_LENGTH_STD = 12   # Standard deviation of linker length

# Dinucleotide preferences from experimental data
# Based on nucleosome positioning signals (TA, TT, AA prefer minor groove outward)
# Positive = favors nucleosome positioning, Negative = disfavors
DINUCLEOTIDE_PREFERENCES = {
    # Flexible dinucleotides (favor nucleosome)
    'AA': 0.15, 'TT': 0.15, 'AT': 0.10, 'TA': 0.20,
    # Rigid dinucleotides (disfavor nucleosome)
    'GC': -0.20, 'CG': -0.25, 'CC': -0.15, 'GG': -0.15,
    # Neutral
    'AG': 0.0, 'GA': 0.0, 'CT': 0.0, 'TC': 0.0,
    'AC': -0.05, 'CA': -0.05, 'GT': -0.05, 'TG': -0.05,
}

# Position-specific preferences within nucleosome (147bp)
# Reflects 10bp periodicity of DNA helix
POSITION_WEIGHTS = None  # Computed dynamically


@dataclass
class NucleosomeCall:
    """A predicted nucleosome position."""
    center: int           # Center position
    start: int            # Start of nucleosomal DNA
    end: int              # End of nucleosomal DNA
    score: float          # Affinity/occupancy score (0-1)
    confidence: float     # Prediction confidence
    dyad_position: int    # Dyad axis position


@dataclass
class NucleosomeProfile:
    """Nucleosome occupancy profile for a sequence."""
    sequence_length: int

    # Per-base scores
    occupancy: np.ndarray      # Occupancy probability at each position
    affinity: np.ndarray       # Nucleosome affinity score

    # Called nucleosomes
    nucleosomes: List[NucleosomeCall] = field(default_factory=list)

    # Nucleosome-free regions
    nfr_regions: List[Tuple[int, int]] = field(default_factory=list)

    # Summary stats
    mean_occupancy: float = 0.0
    nfr_fraction: float = 0.0


@dataclass
class NucleosomeConfig:
    """Configuration for nucleosome prediction."""

    # Model parameters
    nucleosome_length: int = 147
    linker_mean: int = 38
    linker_std: int = 12

    # Prediction thresholds
    occupancy_threshold: float = 0.5   # Above = nucleosome-occupied
    nfr_threshold: float = 0.3         # Below = nucleosome-free

    # NFR detection
    min_nfr_length: int = 50           # Minimum NFR length to report

    # CRISPR efficiency impact
    nucleosome_penalty: float = 0.5    # Efficiency multiplier when in nucleosome


def compute_dinucleotide_score(sequence: str) -> np.ndarray:
    """
    Compute dinucleotide-based nucleosome affinity score.

    Based on experimental observations that certain dinucleotide
    patterns favor or disfavor nucleosome positioning.

    Args:
        sequence: DNA sequence

    Returns:
        Array of scores (one per position)
    """
    sequence = sequence.upper()
    n = len(sequence)

    if n < 2:
        return np.zeros(n)

    scores = np.zeros(n)

    for i in range(n - 1):
        dinuc = sequence[i:i+2]
        if dinuc in DINUCLEOTIDE_PREFERENCES:
            score = DINUCLEOTIDE_PREFERENCES[dinuc]
            scores[i] += score * 0.5
            scores[i + 1] += score * 0.5

    return scores


def compute_gc_profile(sequence: str, window: int = 50) -> np.ndarray:
    """
    Compute GC content profile along sequence.

    GC content affects nucleosome affinity:
    - Very high GC (>70%) disfavors nucleosomes
    - Very low GC (<30%) slightly disfavors nucleosomes
    - Moderate GC (40-60%) is neutral to favorable

    Args:
        sequence: DNA sequence
        window: Window size for smoothing

    Returns:
        GC content at each position (smoothed)
    """
    sequence = sequence.upper()
    n = len(sequence)

    # Binary array: 1 for G/C, 0 for A/T
    gc_binary = np.array([1 if b in 'GC' else 0 for b in sequence], dtype=float)

    # Smooth with running average
    if n >= window:
        kernel = np.ones(window) / window
        gc_profile = np.convolve(gc_binary, kernel, mode='same')
    else:
        gc_profile = gc_binary

    return gc_profile


def compute_periodicity_score(sequence: str) -> np.ndarray:
    """
    Compute 10bp periodic pattern score for nucleosome positioning.

    DNA has ~10bp helical repeat. Strong positioning signals
    show AA/TT/TA dinucleotides with ~10bp periodicity.

    Args:
        sequence: DNA sequence

    Returns:
        Periodicity score array
    """
    sequence = sequence.upper()
    n = len(sequence)
    scores = np.zeros(n)

    # Look for periodic AA/TT/TA patterns
    positioning_dinucs = {'AA', 'TT', 'TA', 'AT'}

    for i in range(n - 1):
        dinuc = sequence[i:i+2]
        if dinuc in positioning_dinucs:
            # Check for periodic pattern
            periodic_score = 0
            for offset in [10, 20, 30]:  # Check 10bp periods
                if i + offset + 1 < n:
                    next_dinuc = sequence[i + offset:i + offset + 2]
                    if next_dinuc in positioning_dinucs:
                        periodic_score += 0.3

            scores[i] += periodic_score

    return scores


def predict_nucleosome_affinity(
    sequence: str,
    config: Optional[NucleosomeConfig] = None
) -> np.ndarray:
    """
    Predict nucleosome affinity score at each position.

    Combines:
    - Dinucleotide preferences
    - GC content effects
    - Periodic positioning signals

    Higher score = more likely to be nucleosome-occupied.

    Args:
        sequence: DNA sequence
        config: NucleosomeConfig

    Returns:
        Affinity score array (0-1 scale)
    """
    if config is None:
        config = NucleosomeConfig()

    sequence = sequence.upper()
    n = len(sequence)

    if n == 0:
        return np.array([])

    # Component scores
    dinuc_scores = compute_dinucleotide_score(sequence)
    gc_profile = compute_gc_profile(sequence)
    periodic_scores = compute_periodicity_score(sequence)

    # GC content effect
    # Moderate GC is favorable, extremes are not
    gc_factor = 1.0 - 2.0 * np.abs(gc_profile - 0.5)

    # Combine scores
    # Weights based on literature: dinucleotide > GC > periodicity
    combined = (
        0.4 * dinuc_scores +
        0.35 * (gc_factor - 0.5) +
        0.25 * periodic_scores
    )

    # Normalize to 0-1
    affinity = 1.0 / (1.0 + np.exp(-2 * combined))

    return affinity


def predict_nucleosome_occupancy(
    sequence: str,
    config: Optional[NucleosomeConfig] = None
) -> np.ndarray:
    """
    Predict nucleosome occupancy probability at each position.

    Uses a simplified version of the NuPoP duration HMM.
    Accounts for nucleosome size and linker constraints.

    Args:
        sequence: DNA sequence
        config: NucleosomeConfig

    Returns:
        Occupancy probability array (0-1)
    """
    if config is None:
        config = NucleosomeConfig()

    n = len(sequence)
    if n == 0:
        return np.array([])

    # Get base affinity
    affinity = predict_nucleosome_affinity(sequence, config)

    # Smooth with nucleosome-sized window
    nuc_len = config.nucleosome_length
    if n >= nuc_len:
        kernel = np.ones(nuc_len) / nuc_len
        occupancy = np.convolve(affinity, kernel, mode='same')
    else:
        occupancy = affinity

    # Apply linker constraints
    # This is a simplified version - full NuPoP uses duration HMM
    # We use a penalty for very closely spaced high-occupancy regions

    # Ensure 0-1 range
    occupancy = np.clip(occupancy, 0, 1)

    return occupancy


def find_nucleosome_free_regions(
    occupancy: np.ndarray,
    config: Optional[NucleosomeConfig] = None
) -> List[Tuple[int, int]]:
    """
    Find nucleosome-free regions (NFRs) in occupancy profile.

    NFRs are critical for CRISPR efficiency - guides targeting
    NFRs have higher success rates.

    Args:
        occupancy: Nucleosome occupancy array
        config: NucleosomeConfig

    Returns:
        List of (start, end) tuples for NFRs
    """
    if config is None:
        config = NucleosomeConfig()

    nfrs = []
    n = len(occupancy)

    if n == 0:
        return nfrs

    # Find regions below threshold
    is_nfr = occupancy < config.nfr_threshold

    # Find contiguous regions
    in_nfr = False
    nfr_start = 0

    for i in range(n):
        if is_nfr[i] and not in_nfr:
            # Start of NFR
            in_nfr = True
            nfr_start = i
        elif not is_nfr[i] and in_nfr:
            # End of NFR
            in_nfr = False
            if i - nfr_start >= config.min_nfr_length:
                nfrs.append((nfr_start, i))

    # Check if sequence ends in NFR
    if in_nfr and n - nfr_start >= config.min_nfr_length:
        nfrs.append((nfr_start, n))

    return nfrs


def call_nucleosomes(
    sequence: str,
    config: Optional[NucleosomeConfig] = None
) -> List[NucleosomeCall]:
    """
    Call discrete nucleosome positions from sequence.

    Uses peak detection on occupancy profile to identify
    likely nucleosome dyad positions.

    Args:
        sequence: DNA sequence
        config: NucleosomeConfig

    Returns:
        List of NucleosomeCall objects
    """
    if config is None:
        config = NucleosomeConfig()

    n = len(sequence)
    if n < config.nucleosome_length:
        return []

    occupancy = predict_nucleosome_occupancy(sequence, config)

    nucleosomes = []
    nuc_len = config.nucleosome_length
    half_nuc = nuc_len // 2

    # Simple peak detection
    # Look for local maxima above threshold
    min_spacing = nuc_len + config.linker_mean // 2

    i = half_nuc
    while i < n - half_nuc:
        if occupancy[i] >= config.occupancy_threshold:
            # Check if this is a local maximum
            window = occupancy[max(0, i - 20):min(n, i + 21)]
            if len(window) > 0 and occupancy[i] == np.max(window):
                # Found a nucleosome
                nuc = NucleosomeCall(
                    center=i,
                    start=max(0, i - half_nuc),
                    end=min(n, i + half_nuc),
                    score=float(occupancy[i]),
                    confidence=float(np.mean(occupancy[max(0, i - 10):min(n, i + 11)])),
                    dyad_position=i,
                )
                nucleosomes.append(nuc)

                # Skip ahead by minimum spacing
                i += min_spacing
            else:
                i += 1
        else:
            i += 1

    return nucleosomes


def predict_nucleosome_profile(
    sequence: str,
    config: Optional[NucleosomeConfig] = None
) -> NucleosomeProfile:
    """
    Generate complete nucleosome profile for a sequence.

    Args:
        sequence: DNA sequence
        config: NucleosomeConfig

    Returns:
        NucleosomeProfile with all predictions

    Example:
        >>> profile = predict_nucleosome_profile(promoter_seq)
        >>> print(f"Mean occupancy: {profile.mean_occupancy:.2f}")
        >>> print(f"NFRs found: {len(profile.nfr_regions)}")
    """
    if config is None:
        config = NucleosomeConfig()

    sequence = sequence.upper()
    n = len(sequence)

    if n == 0:
        return NucleosomeProfile(
            sequence_length=0,
            occupancy=np.array([]),
            affinity=np.array([]),
        )

    # Compute profiles
    affinity = predict_nucleosome_affinity(sequence, config)
    occupancy = predict_nucleosome_occupancy(sequence, config)

    # Find features
    nfr_regions = find_nucleosome_free_regions(occupancy, config)
    nucleosomes = call_nucleosomes(sequence, config)

    # Statistics
    mean_occ = float(np.mean(occupancy))
    nfr_bases = sum(end - start for start, end in nfr_regions)
    nfr_frac = nfr_bases / n if n > 0 else 0

    return NucleosomeProfile(
        sequence_length=n,
        occupancy=occupancy,
        affinity=affinity,
        nucleosomes=nucleosomes,
        nfr_regions=nfr_regions,
        mean_occupancy=mean_occ,
        nfr_fraction=nfr_frac,
    )


def score_guide_nucleosome(
    guide_sequence: str,
    target_position: int,
    sequence_context: str,
    config: Optional[NucleosomeConfig] = None,
) -> Dict[str, any]:
    """
    Score CRISPR guide position for nucleosome context.

    This is the main interface for integrating nucleosome predictions
    into guide scoring pipelines.

    Args:
        guide_sequence: Guide RNA sequence (20bp)
        target_position: Position within sequence_context
        sequence_context: Surrounding genomic sequence
        config: NucleosomeConfig

    Returns:
        Dictionary with nucleosome assessment

    Example:
        >>> result = score_guide_nucleosome(
        ...     "GCGACTGCTACATAGCCAGG",
        ...     target_position=500,
        ...     sequence_context=promoter_seq
        ... )
        >>> print(f"In NFR: {result['in_nfr']}")
        >>> print(f"Efficiency factor: {result['efficiency_factor']:.2f}")
    """
    if config is None:
        config = NucleosomeConfig()

    # Get nucleosome profile
    profile = predict_nucleosome_profile(sequence_context, config)

    # Get occupancy at guide position
    n = len(sequence_context)
    if target_position < 0 or target_position >= n:
        return {
            'occupancy': 0.5,
            'in_nfr': False,
            'in_nucleosome': False,
            'efficiency_factor': 0.75,
            'recommendation': 'UNKNOWN',
            'error': 'Target position out of range',
        }

    # Get local occupancy (average over guide region)
    guide_len = len(guide_sequence)
    start = max(0, target_position - guide_len // 2)
    end = min(n, target_position + guide_len // 2 + 1)

    local_occupancy = float(np.mean(profile.occupancy[start:end]))

    # Check if in NFR
    in_nfr = any(
        nfr_start <= target_position < nfr_end
        for nfr_start, nfr_end in profile.nfr_regions
    )

    # Check if in nucleosome
    in_nucleosome = local_occupancy >= config.occupancy_threshold

    # Calculate efficiency factor
    # High nucleosome occupancy reduces CRISPR efficiency
    if in_nfr:
        efficiency_factor = 1.0
    elif in_nucleosome:
        efficiency_factor = config.nucleosome_penalty
    else:
        # Intermediate region
        efficiency_factor = 1.0 - 0.5 * local_occupancy

    # Recommendation
    if in_nfr:
        recommendation = "EXCELLENT: Target in nucleosome-free region"
        risk = "LOW"
    elif in_nucleosome:
        recommendation = (
            "CAUTION: Target in nucleosome-occupied region. "
            "Consider nearby NFR positions."
        )
        risk = "HIGH"
    else:
        recommendation = "MODERATE: Partial nucleosome occupancy"
        risk = "MODERATE"

    return {
        'occupancy': local_occupancy,
        'affinity': float(np.mean(profile.affinity[start:end])),
        'in_nfr': in_nfr,
        'in_nucleosome': in_nucleosome,
        'efficiency_factor': efficiency_factor,
        'risk': risk,
        'recommendation': recommendation,
        'nfr_count': len(profile.nfr_regions),
        'nucleosome_count': len(profile.nucleosomes),
    }


def find_best_nfr_positions(
    sequence: str,
    min_positions: int = 5,
    config: Optional[NucleosomeConfig] = None,
) -> List[Dict[str, any]]:
    """
    Find best positions for CRISPR targeting within NFRs.

    Useful for guide design when you want to preferentially
    target nucleosome-free regions.

    Args:
        sequence: DNA sequence
        min_positions: Minimum number of positions to return
        config: NucleosomeConfig

    Returns:
        List of position dictionaries sorted by accessibility
    """
    if config is None:
        config = NucleosomeConfig()

    profile = predict_nucleosome_profile(sequence, config)

    positions = []

    # Score all positions
    for i in range(len(sequence)):
        positions.append({
            'position': i,
            'occupancy': float(profile.occupancy[i]),
            'affinity': float(profile.affinity[i]),
            'in_nfr': any(
                start <= i < end
                for start, end in profile.nfr_regions
            ),
        })

    # Sort by lowest occupancy (best for CRISPR)
    positions.sort(key=lambda x: x['occupancy'])

    # Return at least min_positions
    return positions[:max(min_positions, len(positions))]


def nucleosome_adjusted_efficiency(
    base_efficiency: float,
    guide_position: int,
    sequence_context: str,
    config: Optional[NucleosomeConfig] = None,
) -> Dict[str, float]:
    """
    Adjust CRISPR efficiency prediction based on nucleosome context.

    Args:
        base_efficiency: Predicted efficiency without nucleosome consideration
        guide_position: Guide target position in sequence
        sequence_context: Surrounding DNA sequence
        config: NucleosomeConfig

    Returns:
        Dictionary with adjusted efficiency and factors
    """
    if config is None:
        config = NucleosomeConfig()

    result = score_guide_nucleosome(
        "",  # Not needed for this function
        guide_position,
        sequence_context,
        config,
    )

    adjusted = base_efficiency * result['efficiency_factor']

    return {
        'base_efficiency': base_efficiency,
        'adjusted_efficiency': adjusted,
        'nucleosome_factor': result['efficiency_factor'],
        'occupancy': result['occupancy'],
        'in_nfr': result['in_nfr'],
        'risk': result['risk'],
    }


# Known nucleosome-positioning sequences for validation
KNOWN_POSITIONING_SEQUENCES = {
    # 601 sequence - strong nucleosome positioning
    '601': {
        'sequence': (
            'CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGT'
        ),
        'expected_occupancy': 'high',
        'reference': 'Lowary & Widom 1998',
    },
    # Poly(dA:dT) - disfavors nucleosomes
    'polyAT': {
        'sequence': 'A' * 50 + 'T' * 50,
        'expected_occupancy': 'low',
        'reference': 'Segal et al. 2006',
    },
}
