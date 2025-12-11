"""
PhaseLab CRISPR Scoring: Guide RNA quality metrics.

Implements:
- GC content calculation
- Homopolymer run detection
- SantaLucia thermodynamic ΔG
- MIT specificity algorithm
- CFD (Cutting Frequency Determination) score
- Chromatin accessibility modeling
"""

import numpy as np
from typing import Tuple, Optional, Dict
import re


# SantaLucia nearest-neighbor parameters (kcal/mol)
# From: SantaLucia & Hicks (2004) Annu. Rev. Biophys. Biomol. Struct.
NN_PARAMS = {
    'AA': (-1.00, -0.0027),  # (ΔH, ΔS)
    'TT': (-1.00, -0.0027),
    'AT': (-0.88, -0.0024),
    'TA': (-0.58, -0.0015),
    'CA': (-1.45, -0.0039),
    'TG': (-1.45, -0.0039),
    'GT': (-1.44, -0.0037),
    'AC': (-1.44, -0.0037),
    'CT': (-1.28, -0.0033),
    'AG': (-1.28, -0.0033),
    'GA': (-1.30, -0.0032),
    'TC': (-1.30, -0.0032),
    'CG': (-2.17, -0.0055),
    'GC': (-2.24, -0.0056),
    'GG': (-1.84, -0.0046),
    'CC': (-1.84, -0.0046),
}

# Initiation parameters
NN_INIT = {
    'G': (0.98, 0.0024),
    'C': (0.98, 0.0024),
    'A': (1.03, 0.0027),
    'T': (1.03, 0.0027),
}

# MIT position weights for off-target scoring
# Higher weight = more important for specificity
MIT_POSITION_WEIGHTS = [
    0, 0, 0.014, 0, 0,       # positions 1-5 (PAM-distal)
    0.395, 0.317, 0, 0.389, 0.079,  # positions 6-10
    0.445, 0.508, 0.613, 0.851, 0.732,  # positions 11-15
    0.828, 0.615, 0.804, 0.685, 0.583,  # positions 16-20 (PAM-proximal)
]


def gc_content(sequence: str) -> float:
    """
    Calculate GC content of a sequence.

    Args:
        sequence: DNA/RNA sequence.

    Returns:
        GC fraction (0.0 to 1.0).
    """
    sequence = sequence.upper()
    gc = sum(1 for b in sequence if b in 'GC')
    return gc / len(sequence) if sequence else 0.0


def max_homopolymer_run(sequence: str) -> int:
    """
    Find the longest homopolymer run in a sequence.

    Args:
        sequence: DNA/RNA sequence.

    Returns:
        Length of longest single-nucleotide repeat.
    """
    if not sequence:
        return 0

    sequence = sequence.upper()
    max_run = 1
    current_run = 1

    for i in range(1, len(sequence)):
        if sequence[i] == sequence[i - 1]:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            current_run = 1

    return max_run


def delta_g_santalucia(
    sequence: str,
    temperature: float = 37.0,
    na_conc: float = 0.1,
) -> float:
    """
    Calculate ΔG of hybridization using SantaLucia nearest-neighbor model.

    Args:
        sequence: DNA/RNA sequence (assumes binding to perfect complement).
        temperature: Temperature in Celsius.
        na_conc: Na+ concentration in M.

    Returns:
        ΔG in kcal/mol (negative = favorable binding).
    """
    sequence = sequence.upper().replace('U', 'T')
    T_kelvin = temperature + 273.15

    if len(sequence) < 2:
        return 0.0

    # Sum nearest-neighbor contributions
    delta_H = 0.0
    delta_S = 0.0

    for i in range(len(sequence) - 1):
        dinuc = sequence[i:i+2]
        if dinuc in NN_PARAMS:
            dH, dS = NN_PARAMS[dinuc]
            delta_H += dH
            delta_S += dS

    # Initiation
    first_base = sequence[0]
    last_base = sequence[-1]
    if first_base in NN_INIT:
        dH, dS = NN_INIT[first_base]
        delta_H += dH
        delta_S += dS
    if last_base in NN_INIT:
        dH, dS = NN_INIT[last_base]
        delta_H += dH
        delta_S += dS

    # Salt correction (simplified)
    delta_S_corrected = delta_S + 0.368 * len(sequence) * np.log(na_conc) / 1000

    # ΔG = ΔH - TΔS
    delta_G = delta_H - T_kelvin * delta_S_corrected

    return delta_G


def sequence_complexity(sequence: str) -> float:
    """
    Calculate sequence complexity (0 = repetitive, 1 = complex).

    Uses linguistic complexity based on unique k-mers.

    Args:
        sequence: DNA sequence.

    Returns:
        Complexity score (0.0 to 1.0).
    """
    sequence = sequence.upper()
    n = len(sequence)

    if n < 3:
        return 1.0

    # Count unique k-mers for k=1,2,3
    total_possible = 0
    total_unique = 0

    for k in [1, 2, 3]:
        kmers = set()
        for i in range(n - k + 1):
            kmers.add(sequence[i:i+k])
        possible = min(4**k, n - k + 1)
        total_possible += possible
        total_unique += len(kmers)

    return total_unique / total_possible if total_possible > 0 else 1.0


def mit_specificity_score(
    guide_seq: str,
    off_target_count: int = 0,
    avg_mismatches: float = 4.0,
) -> float:
    """
    Calculate MIT specificity score (simplified).

    The full MIT algorithm requires genome-wide alignment.
    This provides an estimate based on guide sequence properties.

    Higher score = more specific (fewer predicted off-targets).

    Args:
        guide_seq: 20bp guide sequence.
        off_target_count: Number of known off-targets (if available).
        avg_mismatches: Average mismatches to off-targets.

    Returns:
        MIT specificity score (0-100).
    """
    guide_seq = guide_seq.upper()

    # Base score from sequence complexity
    complexity = sequence_complexity(guide_seq)

    # Penalize low GC or very high GC
    gc = gc_content(guide_seq)
    gc_penalty = 0.0
    if gc < 0.4 or gc > 0.7:
        gc_penalty = 10 * abs(gc - 0.55)

    # Penalize homopolymers
    max_homo = max_homopolymer_run(guide_seq)
    homo_penalty = max(0, (max_homo - 3) * 5)

    # Seed region (PAM-proximal 12nt) importance
    seed = guide_seq[-12:]
    seed_complexity = sequence_complexity(seed)

    # Estimate score
    base_score = 100 * complexity * seed_complexity
    score = base_score - gc_penalty - homo_penalty

    # Adjust for known off-targets
    if off_target_count > 0:
        score -= min(30, off_target_count * 2)

    return max(0, min(100, score))


def cfd_score(
    guide_seq: str,
    target_seq: Optional[str] = None,
) -> float:
    """
    Calculate CFD (Cutting Frequency Determination) score.

    CFD predicts how likely an off-target site will be cut.
    For on-target (no mismatches), returns 100.

    Args:
        guide_seq: 20bp guide sequence.
        target_seq: Target sequence (if different from perfect match).

    Returns:
        CFD score (0-100, higher = more cutting).
    """
    guide_seq = guide_seq.upper()

    if target_seq is None:
        # On-target: perfect match
        return 100.0

    target_seq = target_seq.upper()

    if len(guide_seq) != len(target_seq):
        return 0.0

    # Count mismatches and their positions
    mismatches = []
    for i, (g, t) in enumerate(zip(guide_seq, target_seq)):
        if g != t:
            mismatches.append(i)

    if not mismatches:
        return 100.0

    # CFD penalty increases with mismatches and seed region location
    cfd = 100.0
    for pos in mismatches:
        # Higher penalty for PAM-proximal (seed) mismatches
        if pos >= 8:  # Seed region
            cfd *= 0.5
        else:
            cfd *= 0.7

    return max(0, cfd)


def chromatin_accessibility_score(
    position: int,
    tss_position: int,
    dnase_peaks: Optional[list] = None,
) -> Tuple[str, float]:
    """
    Estimate chromatin accessibility at a genomic position.

    Without experimental data, uses heuristic based on TSS proximity.
    Near TSS = more likely to be open chromatin.

    Args:
        position: Genomic position.
        tss_position: Transcription start site position.
        dnase_peaks: Optional list of (start, end) DNase HS peaks.

    Returns:
        (state, accessibility_score)
        state: "OPEN", "MODERATE", or "CLOSED"
        accessibility_score: 0.0 to 1.0
    """
    rel_pos = position - tss_position

    # Check if in provided DNase peaks
    if dnase_peaks:
        for start, end in dnase_peaks:
            if start <= position <= end:
                return ("OPEN", 0.9)

    # Heuristic based on TSS proximity
    abs_dist = abs(rel_pos)

    if abs_dist < 200:
        # Very close to TSS - likely open
        score = 0.8 - 0.001 * abs_dist
        state = "OPEN"
    elif abs_dist < 500:
        score = 0.6 - 0.0005 * (abs_dist - 200)
        state = "MODERATE"
    else:
        score = max(0.2, 0.5 - 0.0002 * (abs_dist - 500))
        state = "MODERATE" if score > 0.35 else "CLOSED"

    return (state, score)
