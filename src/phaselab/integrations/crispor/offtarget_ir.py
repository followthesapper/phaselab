"""
Off-Target IR Analysis Module.

Applies Informational Relativity metrics to CRISPOR off-target data
to provide structure-aware risk assessment beyond simple counting.

Key Metrics:
1. Off-target entropy - Is risk concentrated or diffuse?
2. Coherence contrast (ΔR̄) - Are off-targets "too good" physically?
3. Phase clustering - Do off-targets form dangerous families?
4. Energy spectrum - Tail-risk from favorable binding energies

These metrics complement CRISPOR's MIT/CFD scores by asking
"how dangerous is the off-target *structure*?" not just "how many?"
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class OffTargetSite:
    """Represents a single off-target site with IR metrics."""
    sequence: str
    mismatches: int
    cfd_score: float
    gene: Optional[str] = None
    region: Optional[str] = None  # exon, intron, intergenic

    # IR metrics (computed)
    r_bar: float = 0.0
    binding_energy: float = 0.0
    phase_mean: float = 0.0

    @property
    def is_exonic(self) -> bool:
        return self.region and 'exon' in self.region.lower()

    @property
    def is_close(self) -> bool:
        """0-2 mismatches are considered 'close' and dangerous."""
        return self.mismatches <= 2


@dataclass
class OffTargetIRAnalysis:
    """
    Complete IR analysis of a guide's off-target landscape.

    This provides structure-aware metrics that CRISPOR alone doesn't compute.
    """
    guide_sequence: str

    # Core IR metrics
    r_bar_on_target: float = 0.0
    r_bar_off_max: float = 0.0  # Highest R̄ among off-targets
    delta_r_bar: float = 0.0    # R̄_on - R̄_off_max (coherence contrast)

    # Off-target entropy (risk distribution)
    entropy: float = 0.0  # Shannon entropy of risk weights
    entropy_normalized: float = 0.0  # 0-1 scale

    # Energy spectrum
    energy_on_target: float = 0.0
    energy_off_mean: float = 0.0
    energy_off_max: float = 0.0  # Most favorable off-target energy
    energy_tail_risk: float = 0.0  # 95th percentile

    # Clustering (if computed)
    n_clusters: int = 0
    largest_cluster_size: int = 0
    cluster_concentration: float = 0.0

    # Summary
    risk_concentrated: bool = False  # Low entropy = few dangerous sites
    off_target_competitive: bool = False  # ΔR̄ < threshold
    has_exonic_risk: bool = False
    has_exonic_close_ot: bool = False  # Exonic off-targets with ≤2 mismatches

    # IR-enhanced score (pre-computed for convenience)
    ir_enhanced_score: float = 0.0  # Score penalty/bonus from IR analysis

    # Raw data
    n_offtargets_analyzed: int = 0
    close_offtargets: List[OffTargetSite] = field(default_factory=list)

    def to_dict(self) -> Dict:
        return {
            'guide_sequence': self.guide_sequence,
            'r_bar_on_target': float(self.r_bar_on_target),
            'r_bar_off_max': float(self.r_bar_off_max),
            'delta_r_bar': float(self.delta_r_bar),
            'entropy': float(self.entropy),
            'entropy_normalized': float(self.entropy_normalized),
            'energy_on_target': float(self.energy_on_target),
            'energy_off_mean': float(self.energy_off_mean),
            'energy_tail_risk': float(self.energy_tail_risk),
            'n_clusters': self.n_clusters,
            'risk_concentrated': self.risk_concentrated,
            'off_target_competitive': self.off_target_competitive,
            'has_exonic_risk': self.has_exonic_risk,
            'has_exonic_close_ot': self.has_exonic_close_ot,
            'ir_enhanced_score': float(self.ir_enhanced_score),
            'n_offtargets_analyzed': self.n_offtargets_analyzed,
        }


# Thermodynamic parameters for binding energy estimation
NN_PARAMS = {
    'AA': -1.00, 'AT': -0.88, 'AG': -1.28, 'AC': -1.44,
    'TA': -0.58, 'TT': -1.00, 'TG': -1.45, 'TC': -1.30,
    'GA': -1.30, 'GT': -1.44, 'GG': -1.84, 'GC': -2.24,
    'CA': -1.45, 'CT': -1.28, 'CG': -2.17, 'CC': -1.84,
}


def compute_binding_energy(guide: str, target: str) -> float:
    """
    Estimate binding energy using nearest-neighbor thermodynamics.

    More negative = more stable binding.
    """
    if len(guide) != len(target):
        return 0.0

    energy = 0.0
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    for i in range(len(guide) - 1):
        dinuc = guide[i:i+2]
        base_energy = NN_PARAMS.get(dinuc, -1.0)

        # Perfect match contribution
        match_i = guide[i] == complement.get(target[i], 'N')
        match_j = guide[i+1] == complement.get(target[i+1], 'N')

        if match_i and match_j:
            energy += base_energy
        elif match_i or match_j:
            energy += base_energy * 0.3  # Partial match
        else:
            energy += abs(base_energy) * 0.5  # Mismatch penalty

    return energy


def compute_coherence_simple(expectations: np.ndarray) -> float:
    """
    Compute phase coherence R̄ from expectation values.

    R̄ = |⟨e^(iφ)⟩| where φ = arccos(⟨O⟩)
    """
    phases = np.arccos(np.clip(expectations, -1, 1))
    return float(np.abs(np.mean(np.exp(1j * phases))))


def compute_offtarget_entropy(
    offtargets: List[OffTargetSite],
    use_cfd_weights: bool = True
) -> Tuple[float, float]:
    """
    Compute Shannon entropy of off-target risk distribution.

    Low entropy = risk concentrated in few sites (dangerous)
    High entropy = risk diffuse across many weak sites (safer)

    Args:
        offtargets: List of off-target sites
        use_cfd_weights: If True, weight by CFD score; otherwise by mismatches

    Returns:
        (entropy, normalized_entropy) where normalized is 0-1 scale
    """
    if not offtargets:
        return 0.0, 1.0  # No off-targets = maximum "safety"

    # Compute weights for each off-target
    weights = []
    for ot in offtargets:
        if use_cfd_weights:
            # CFD is 0-1, higher = more likely to cut
            w = ot.cfd_score if ot.cfd_score > 0 else 0.01
        else:
            # Weight inversely by mismatches (0mm = highest weight)
            w = 1.0 / (1.0 + ot.mismatches)

        # Boost weight for exonic off-targets
        if ot.is_exonic:
            w *= 2.0

        weights.append(w)

    weights = np.array(weights)

    # Normalize to probability distribution
    if weights.sum() == 0:
        return 0.0, 1.0

    p = weights / weights.sum()

    # Shannon entropy
    # H = -Σ p_i log(p_i)
    p_nonzero = p[p > 0]
    entropy = -np.sum(p_nonzero * np.log(p_nonzero))

    # Normalize to 0-1 (max entropy = log(n))
    max_entropy = np.log(len(offtargets)) if len(offtargets) > 1 else 1.0
    normalized = entropy / max_entropy if max_entropy > 0 else 1.0

    return float(entropy), float(normalized)


def compute_coherence_contrast(
    guide_sequence: str,
    on_target_r_bar: float,
    offtargets: List[OffTargetSite],
    top_k: int = 20
) -> Tuple[float, float]:
    """
    Compute coherence contrast: ΔR̄ = R̄_on - max(R̄_off).

    This measures whether off-targets are "too good" physically.

    Args:
        guide_sequence: The guide RNA sequence
        on_target_r_bar: Pre-computed on-target coherence
        offtargets: List of off-target sites
        top_k: Only analyze top K most dangerous off-targets

    Returns:
        (r_bar_off_max, delta_r_bar)
    """
    if not offtargets:
        return 0.0, on_target_r_bar

    # Sort by danger (low mismatches + high CFD)
    sorted_ots = sorted(
        offtargets,
        key=lambda x: (x.mismatches, -x.cfd_score)
    )[:top_k]

    r_bar_values = []
    for ot in sorted_ots:
        # Simple coherence estimate based on sequence similarity
        # In full implementation, would run quantum simulation
        mismatch_penalty = 1.0 - (ot.mismatches * 0.1)  # ~10% drop per mismatch
        estimated_r_bar = on_target_r_bar * max(0.5, mismatch_penalty)

        # Adjust by CFD (higher CFD = closer to on-target behavior)
        estimated_r_bar *= (0.7 + 0.3 * ot.cfd_score)

        r_bar_values.append(estimated_r_bar)
        ot.r_bar = estimated_r_bar

    r_bar_off_max = max(r_bar_values) if r_bar_values else 0.0
    delta_r_bar = on_target_r_bar - r_bar_off_max

    return r_bar_off_max, delta_r_bar


def compute_energy_spectrum(
    guide_sequence: str,
    offtargets: List[OffTargetSite]
) -> Tuple[float, float, float]:
    """
    Compute energy spectrum of off-targets.

    Returns:
        (mean_energy, max_energy, tail_risk_95)
    """
    if not offtargets:
        return 0.0, 0.0, 0.0

    energies = []
    for ot in offtargets:
        energy = compute_binding_energy(guide_sequence, ot.sequence)
        ot.binding_energy = energy
        energies.append(energy)

    energies = np.array(energies)

    mean_energy = float(np.mean(energies))
    max_energy = float(np.min(energies))  # Most negative = most stable
    tail_risk = float(np.percentile(energies, 5))  # 5th percentile (most stable)

    return mean_energy, max_energy, tail_risk


def compute_offtarget_clustering(
    offtargets: List[OffTargetSite],
    similarity_threshold: float = 0.8,
) -> Tuple[int, int, float]:
    """
    Cluster off-targets by sequence similarity to detect dangerous "families".

    Uses a simple hierarchical approach based on sequence identity.
    Off-targets that are similar to each other may indicate a systematic
    vulnerability (e.g., pseudogenes, gene families).

    Args:
        offtargets: List of off-target sites
        similarity_threshold: Minimum similarity to be in same cluster (0-1)

    Returns:
        (n_clusters, largest_cluster_size, cluster_concentration)
        - cluster_concentration: fraction of off-targets in largest cluster
    """
    if len(offtargets) < 2:
        return len(offtargets), len(offtargets), 1.0 if offtargets else 0.0

    # Simple single-linkage clustering based on sequence similarity
    n = len(offtargets)
    cluster_id = list(range(n))  # Each off-target starts in its own cluster

    def sequence_similarity(seq1: str, seq2: str) -> float:
        """Compute fraction of matching positions."""
        if len(seq1) != len(seq2):
            min_len = min(len(seq1), len(seq2))
            seq1 = seq1[:min_len]
            seq2 = seq2[:min_len]
        if not seq1:
            return 0.0
        matches = sum(1 for a, b in zip(seq1.upper(), seq2.upper()) if a == b)
        return matches / len(seq1)

    def find_root(i: int) -> int:
        """Union-find: find cluster root."""
        while cluster_id[i] != i:
            cluster_id[i] = cluster_id[cluster_id[i]]  # Path compression
            i = cluster_id[i]
        return i

    def union(i: int, j: int):
        """Union-find: merge clusters."""
        root_i = find_root(i)
        root_j = find_root(j)
        if root_i != root_j:
            cluster_id[root_i] = root_j

    # Cluster similar sequences
    for i in range(n):
        for j in range(i + 1, n):
            sim = sequence_similarity(offtargets[i].sequence, offtargets[j].sequence)
            if sim >= similarity_threshold:
                union(i, j)

    # Count cluster sizes
    from collections import Counter
    roots = [find_root(i) for i in range(n)]
    cluster_sizes = Counter(roots)

    n_clusters = len(cluster_sizes)
    largest_cluster = max(cluster_sizes.values())
    concentration = largest_cluster / n

    return n_clusters, largest_cluster, concentration


def compute_region_difficulty(
    sequence: str,
    k: int = 3,
) -> Dict[str, float]:
    """
    Compute region difficulty metrics ("soup index").

    GC-rich and repetitive regions are inherently harder to design
    specific guides for. This normalizes expectations.

    Args:
        sequence: DNA sequence to analyze
        k: k-mer size for entropy calculation

    Returns:
        Dict with difficulty metrics:
        - gc_content: GC fraction
        - kmer_entropy: k-mer Shannon entropy (higher = more complex)
        - repeat_fraction: fraction of sequence in homopolymers
        - difficulty_index: composite score (higher = harder)
    """
    sequence = sequence.upper().replace('N', '')

    if len(sequence) < k:
        return {
            'gc_content': 0.5,
            'kmer_entropy': 0.0,
            'repeat_fraction': 0.0,
            'difficulty_index': 0.5,
        }

    # GC content
    gc_count = sequence.count('G') + sequence.count('C')
    gc_content = gc_count / len(sequence)

    # k-mer entropy
    from collections import Counter
    kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
    kmer_counts = Counter(kmers)
    total_kmers = sum(kmer_counts.values())

    kmer_entropy = 0.0
    for count in kmer_counts.values():
        p = count / total_kmers
        if p > 0:
            kmer_entropy -= p * np.log2(p)

    # Normalize by max possible entropy (all unique k-mers)
    max_entropy = np.log2(min(4**k, total_kmers))
    kmer_entropy_norm = kmer_entropy / max_entropy if max_entropy > 0 else 0

    # Repeat fraction (homopolymers of length >= 4)
    import re
    homopolymers = re.findall(r'(.)\1{3,}', sequence)
    repeat_length = sum(len(h) + 3 for h in homopolymers)  # +3 because findall returns just the char
    repeat_fraction = repeat_length / len(sequence)

    # Composite difficulty index
    # High GC, low k-mer diversity, and repeats all make design harder
    gc_difficulty = abs(gc_content - 0.5) * 2  # 0 at 50% GC, 1 at 0% or 100%
    diversity_difficulty = 1 - kmer_entropy_norm  # Low entropy = repetitive = hard
    repeat_difficulty = repeat_fraction

    # Weighted combination
    difficulty_index = (
        0.4 * gc_difficulty +
        0.4 * diversity_difficulty +
        0.2 * repeat_difficulty
    )

    return {
        'gc_content': gc_content,
        'kmer_entropy': kmer_entropy_norm,
        'repeat_fraction': repeat_fraction,
        'difficulty_index': difficulty_index,
    }


def analyze_offtarget_landscape(
    guide_sequence: str,
    on_target_r_bar: float,
    offtargets: List[OffTargetSite],
    on_target_energy: Optional[float] = None,
    delta_r_bar_threshold: float = 0.05
) -> OffTargetIRAnalysis:
    """
    Comprehensive IR analysis of a guide's off-target landscape.

    This is the main entry point that computes all IR-enhanced metrics.

    Args:
        guide_sequence: 20bp guide sequence
        on_target_r_bar: Pre-computed on-target coherence
        offtargets: List of OffTargetSite objects from CRISPOR
        on_target_energy: Pre-computed on-target binding energy (optional)
        delta_r_bar_threshold: Threshold below which off-targets are "competitive"

    Returns:
        OffTargetIRAnalysis with all computed metrics
    """
    analysis = OffTargetIRAnalysis(
        guide_sequence=guide_sequence,
        r_bar_on_target=on_target_r_bar,
        n_offtargets_analyzed=len(offtargets)
    )

    if not offtargets:
        analysis.entropy_normalized = 1.0  # No off-targets = safe
        analysis.delta_r_bar = on_target_r_bar
        return analysis

    # 1. Compute off-target entropy
    entropy, entropy_norm = compute_offtarget_entropy(offtargets)
    analysis.entropy = entropy
    analysis.entropy_normalized = entropy_norm
    analysis.risk_concentrated = entropy_norm < 0.5  # Low entropy = concentrated risk

    # 2. Compute coherence contrast
    r_bar_off_max, delta_r_bar = compute_coherence_contrast(
        guide_sequence, on_target_r_bar, offtargets
    )
    analysis.r_bar_off_max = r_bar_off_max
    analysis.delta_r_bar = delta_r_bar
    analysis.off_target_competitive = delta_r_bar < delta_r_bar_threshold

    # 3. Compute energy spectrum
    mean_e, max_e, tail_e = compute_energy_spectrum(guide_sequence, offtargets)
    analysis.energy_off_mean = mean_e
    analysis.energy_off_max = max_e
    analysis.energy_tail_risk = tail_e

    if on_target_energy is not None:
        analysis.energy_on_target = on_target_energy

    # 4. Check for exonic risk
    close_exonic = [ot for ot in offtargets if ot.is_close and ot.is_exonic]
    analysis.has_exonic_risk = len(close_exonic) > 0
    analysis.has_exonic_close_ot = len(close_exonic) > 0  # Alias for clarity

    # 5. Store close off-targets for inspection
    analysis.close_offtargets = [ot for ot in offtargets if ot.is_close]

    # 6. Compute off-target clustering (phase families)
    n_clusters, largest_cluster, concentration = compute_offtarget_clustering(offtargets)
    analysis.n_clusters = n_clusters
    analysis.largest_cluster_size = largest_cluster
    analysis.cluster_concentration = concentration

    # 7. Pre-compute IR enhanced score for convenience
    # This gives the penalty/bonus that would be applied to a base score
    ir_score, _ = compute_ir_enhanced_score(0.0, analysis)
    analysis.ir_enhanced_score = ir_score

    return analysis


def compute_risk_mass(
    offtargets: List[OffTargetSite],
) -> Dict[str, float]:
    """
    Compute risk mass metrics: sum(CFD) by mismatch bucket and annotation.

    Risk mass is more informative than off-target counts because it weights
    by actual cutting probability (CFD score).

    Args:
        offtargets: List of off-target sites

    Returns:
        Dict with:
        - risk_mass_close: sum(CFD) for 0-2mm off-targets
        - risk_mass_distant: sum(CFD) for 3-4mm off-targets
        - risk_mass_exonic: sum(CFD) for exonic off-targets
        - risk_mass_by_mm: {0: sum, 1: sum, 2: sum, 3: sum, 4: sum}
    """
    risk_mass_by_mm = {i: 0.0 for i in range(5)}
    risk_mass_exonic = 0.0

    for ot in offtargets:
        mm = min(ot.mismatches, 4)  # Cap at 4
        risk_mass_by_mm[mm] += ot.cfd_score

        if ot.is_exonic:
            risk_mass_exonic += ot.cfd_score

    return {
        'risk_mass_close': sum(risk_mass_by_mm[i] for i in range(3)),  # 0-2mm
        'risk_mass_distant': sum(risk_mass_by_mm[i] for i in range(3, 5)),  # 3-4mm
        'risk_mass_exonic': risk_mass_exonic,
        'risk_mass_by_mm': risk_mass_by_mm,
    }


def compute_tail_risk(
    offtargets: List[OffTargetSite],
    focus_exonic: bool = True,
) -> Tuple[float, Optional[OffTargetSite]]:
    """
    Compute tail-risk score: the worst single off-target.

    Tail risk dominates safety perception more than total OT count.
    A single high-CFD exonic off-target is more dangerous than
    100 low-CFD intergenic off-targets.

    Args:
        offtargets: List of off-target sites
        focus_exonic: If True, only consider exonic off-targets for max

    Returns:
        (max_cfd_score, worst_offtarget)
    """
    if not offtargets:
        return 0.0, None

    if focus_exonic:
        exonic = [ot for ot in offtargets if ot.is_exonic]
        if exonic:
            worst = max(exonic, key=lambda x: x.cfd_score)
            return worst.cfd_score, worst

    # Fallback to all off-targets
    worst = max(offtargets, key=lambda x: x.cfd_score)
    return worst.cfd_score, worst


def compute_concentration_measures(
    offtargets: List[OffTargetSite],
) -> Dict[str, float]:
    """
    Compute concentration measures for off-target risk.

    These complement entropy by capturing different aspects of risk distribution:
    - Gini coefficient: 0 = all equal, 1 = all in one site
    - Herfindahl-Hirschman Index (HHI): 0-1, higher = more concentrated

    Args:
        offtargets: List of off-target sites

    Returns:
        Dict with gini_coefficient, herfindahl_index
    """
    if not offtargets:
        return {'gini_coefficient': 0.0, 'herfindahl_index': 0.0}

    # Get CFD scores as "market shares"
    cfd_scores = np.array([ot.cfd_score for ot in offtargets])

    if cfd_scores.sum() == 0:
        return {'gini_coefficient': 0.0, 'herfindahl_index': 0.0}

    # Normalize to proportions
    proportions = cfd_scores / cfd_scores.sum()

    # Gini coefficient
    # G = (2 * sum_i(i * x_i) / (n * sum(x_i))) - (n + 1) / n
    n = len(cfd_scores)
    sorted_scores = np.sort(cfd_scores)
    cumsum = np.cumsum(sorted_scores)
    gini = (2 * np.sum((np.arange(1, n + 1) * sorted_scores))) / (n * cfd_scores.sum()) - (n + 1) / n

    # Herfindahl-Hirschman Index
    # HHI = sum(p_i^2) where p_i is market share
    hhi = float(np.sum(proportions ** 2))

    return {
        'gini_coefficient': float(max(0, min(1, gini))),
        'herfindahl_index': hhi,
    }


def compute_ir_enhanced_score(
    base_score: float,
    ir_analysis: OffTargetIRAnalysis,
    weights: Optional[Dict[str, float]] = None
) -> Tuple[float, List[str]]:
    """
    Adjust a guide's score based on IR analysis of off-targets.

    This adds/subtracts from the base score based on:
    - Entropy bonus/penalty
    - Coherence contrast
    - Exonic risk
    - Energy tail risk

    Args:
        base_score: Initial score from standard pipeline
        ir_analysis: Results from analyze_offtarget_landscape
        weights: Optional custom weights

    Returns:
        (adjusted_score, list of adjustment reasons)
    """
    if weights is None:
        weights = {
            'entropy_bonus': 0.05,      # Bonus for diffuse risk
            'entropy_penalty': 0.10,    # Penalty for concentrated risk
            'delta_r_bar_penalty': 0.08,  # If off-targets competitive
            'exonic_risk_penalty': 0.15,  # Exonic close off-targets
            'energy_tail_penalty': 0.05,  # Very stable off-targets
        }

    adjustments = []
    score = base_score

    # 1. Entropy adjustment
    if ir_analysis.entropy_normalized > 0.7:
        # Diffuse risk = bonus
        bonus = weights['entropy_bonus']
        score += bonus
        adjustments.append(f"+{bonus:.2f} diffuse risk (entropy={ir_analysis.entropy_normalized:.2f})")
    elif ir_analysis.entropy_normalized < 0.3:
        # Concentrated risk = penalty
        penalty = weights['entropy_penalty']
        score -= penalty
        adjustments.append(f"-{penalty:.2f} concentrated risk (entropy={ir_analysis.entropy_normalized:.2f})")

    # 2. Coherence contrast
    if ir_analysis.off_target_competitive:
        penalty = weights['delta_r_bar_penalty']
        score -= penalty
        adjustments.append(f"-{penalty:.2f} competitive off-targets (ΔR̄={ir_analysis.delta_r_bar:.3f})")

    # 3. Exonic risk
    if ir_analysis.has_exonic_risk:
        penalty = weights['exonic_risk_penalty']
        score -= penalty
        n_close_exonic = len([ot for ot in ir_analysis.close_offtargets if ot.is_exonic])
        adjustments.append(f"-{penalty:.2f} exonic close off-targets ({n_close_exonic})")

    # 4. Energy tail risk
    # Penalize if off-targets have very favorable binding energy
    if ir_analysis.energy_tail_risk < -15:  # Very stable
        penalty = weights['energy_tail_penalty']
        score -= penalty
        adjustments.append(f"-{penalty:.2f} high-affinity off-targets (E_tail={ir_analysis.energy_tail_risk:.1f})")

    # Bound score
    score = max(0.0, min(1.0, score))

    return score, adjustments
