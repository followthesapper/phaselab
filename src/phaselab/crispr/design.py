"""
PhaseLab Guide Design Pipeline (v0.9.2+).

The complete pipeline for guide RNA design:
1. Region Declaration - Define where to look (promoter hypothesis)
2. Enumeration - Find all candidate protospacers
3. Evaluation - Compute risk structure (off-targets, constraints)
4. Policy Ranking - Tier and rank under explicit policy

KEY PRINCIPLES:
- Every assumption is declared, not guessed
- Policies explain outcomes
- Results are reproducible (manifest)
- Biology is user-declared, not inferred

Example:
    >>> from phaselab.crispr.design import design_crispra_guides
    >>>
    >>> result = design_crispra_guides(
    ...     gene_symbol="RAI1",
    ...     promoter_sequence=promoter_seq,
    ...     tss_position=500,
    ...     policy=RankingPolicy.BINDING_STRICT,
    ... )
    >>>
    >>> print(result.summary())
    >>> print(result.top_tier_a_guides())
"""

from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple
from datetime import datetime
import json
import hashlib
import logging

from .region import (
    Region, RegionSet, RegionBuilder, Window,
    GenomeBuild, Modality, TSSSource, TSSAnnotation,
    build_regions_for_gene, get_tss_for_gene,
)
from .enumerate import (
    Nuclease, NucleaseRole, NUCLEASE_CONFIGS, CandidateGuide,
    EnumerationResult, PAMScanner,
    enumerate_guides, enumerate_from_sequence,
)
from .scoring import (
    RankingPolicy, POLICY_CONFIGS, PolicyConfig,
    GuideTier, GateResult,
    rank_guides, apply_hard_gates, emit_manifest, print_ranking_report,
    gc_content, max_homopolymer_run, sequence_complexity,
    poly_t_penalty, is_repeat_region, u6_compatibility_check,
    delta_g_santalucia, mit_specificity_score,
)

logger = logging.getLogger(__name__)


# =============================================================================
# DESIGN RESULT
# =============================================================================

@dataclass
class DesignResult:
    """
    Complete result of guide design pipeline.

    Contains all intermediate outputs for transparency.
    """
    # Core results
    ranked_guides: List[Dict[str, Any]]
    tiers: Dict[str, List[Dict[str, Any]]]

    # Configuration
    gene_symbol: str
    genome_build: GenomeBuild
    modality: Modality
    nuclease: Nuclease
    policy: RankingPolicy

    # Intermediate data
    region_set: Optional[RegionSet] = None
    enumeration: Optional[EnumerationResult] = None
    total_enumerated: int = 0
    total_after_gates: int = 0

    # Manifest
    manifest: Dict[str, Any] = field(default_factory=dict)

    @property
    def tier_a_guides(self) -> List[Dict[str, Any]]:
        """Get Tier A guides (0/0/0 off-targets)."""
        return self.tiers.get('A', [])

    @property
    def tier_b_guides(self) -> List[Dict[str, Any]]:
        """Get Tier B guides."""
        return self.tiers.get('B', [])

    @property
    def top_guide(self) -> Optional[Dict[str, Any]]:
        """Get the top-ranked guide."""
        return self.ranked_guides[0] if self.ranked_guides else None

    def summary(self) -> str:
        """Generate human-readable summary."""
        lines = [
            "=" * 70,
            f"PHASELAB DESIGN RESULT: {self.gene_symbol}",
            "=" * 70,
            f"Genome: {self.genome_build.value}",
            f"Modality: {self.modality.value}",
            f"Nuclease: {self.nuclease.value}",
            f"Policy: {self.policy.value}",
            "",
            f"Candidates enumerated: {self.total_enumerated}",
            f"Passed hard gates: {self.total_after_gates}",
            "",
            f"Tier A (0/0/0): {len(self.tier_a_guides)} guides",
            f"Tier B (0/0/1-2): {len(self.tier_b_guides)} guides",
            f"Tier C (other): {len(self.tiers.get('C', []))} guides",
        ]

        if self.top_guide:
            lines.extend([
                "",
                "TOP RECOMMENDATION:",
                f"  Sequence: {self.top_guide.get('sequence')}",
                f"  Tier: {self.top_guide.get('tier')}",
                f"  Position: {self.top_guide.get('tss_relative_position')} from TSS",
            ])

        lines.append("=" * 70)
        return '\n'.join(lines)

    def to_dataframe(self):
        """Convert ranked guides to pandas DataFrame (if pandas available)."""
        try:
            import pandas as pd
            return pd.DataFrame(self.ranked_guides)
        except ImportError:
            raise ImportError("pandas required for to_dataframe()")


# =============================================================================
# EVALUATION (sequence quality + basic scoring)
# =============================================================================

def evaluate_candidate(
    candidate: CandidateGuide,
    check_u6: bool = True,
) -> Dict[str, Any]:
    """
    Evaluate a candidate guide for sequence quality.

    This applies local heuristics. Full off-target evaluation
    requires CRISPOR integration.

    Args:
        candidate: CandidateGuide from enumeration.
        check_u6: Check U6/Pol III compatibility.

    Returns:
        Dict with evaluation metrics.
    """
    seq = candidate.sequence

    # Basic metrics
    gc = gc_content(seq)
    homo = max_homopolymer_run(seq)
    complexity = sequence_complexity(seq)
    delta_g = delta_g_santalucia(seq)

    # Heuristic specificity (not real off-target count)
    mit_estimate = mit_specificity_score(seq)

    # U6 compatibility
    u6_compatible = True
    u6_warnings = []
    if check_u6:
        u6_compatible, u6_warnings = u6_compatibility_check(seq)

    # Repeat region check
    is_repeat, repeat_reason = is_repeat_region(seq)

    # Build result
    result = {
        **candidate.to_dict(),
        'gc_content': gc,
        'homopolymer_run': homo,
        'complexity': complexity,
        'delta_g': delta_g,
        'mit_score_estimate': mit_estimate,
        'u6_compatible': u6_compatible,
        'u6_warnings': u6_warnings,
        'is_repeat': is_repeat,
        'repeat_reason': repeat_reason if is_repeat else None,

        # Placeholder for CRISPOR data (to be filled by integration)
        'off_targets': {},
        'crispor_validated': False,
        'is_unscorable': False,
    }

    return result


def evaluate_candidates(
    candidates: List[CandidateGuide],
    check_u6: bool = True,
) -> List[Dict[str, Any]]:
    """Evaluate a list of candidates."""
    return [evaluate_candidate(c, check_u6) for c in candidates]


# =============================================================================
# MAIN DESIGN FUNCTIONS
# =============================================================================

def design_crispra_guides(
    gene_symbol: str,
    promoter_sequence: Optional[str] = None,
    tss_position: Optional[int] = None,
    genome_build: GenomeBuild = GenomeBuild.HG38,
    nuclease: Nuclease = Nuclease.SPCAS9,
    policy: RankingPolicy = RankingPolicy.BINDING_STRICT,
    window: Optional[Window] = None,
    check_u6: bool = True,
    relaxed_pam: bool = True,  # v0.9.3: Use relaxed PAM for CRISPRa by default
    guide_length: Optional[int] = None,  # v0.9.3: Override guide length
) -> DesignResult:
    """
    Design CRISPRa guides for a target gene.

    This is the main entry point for CRISPRa guide design.

    v0.9.3: CRISPRa uses BINDING mode by default (relaxed PAM constraints).
            This is experimentally validated - dCas9 binding tolerates
            non-canonical PAMs that would never support cutting.
            Guide length can be overridden (some papers use 20bp with SaCas9).

    Args:
        gene_symbol: Target gene name (e.g., "RAI1").
        promoter_sequence: Promoter DNA sequence (required).
        tss_position: Position of TSS within promoter_sequence.
        genome_build: Genome build for annotation lookup.
        nuclease: CRISPR nuclease to use.
        policy: Ranking policy (default: BINDING_STRICT for CRISPRa).
        window: Custom window (uses default CRISPRa window if None).
        check_u6: Check U6/Pol III compatibility.
        relaxed_pam: Use binding-mode PAM (default True for CRISPRa).
        guide_length: Override default guide length (e.g., 20bp with SaCas9).

    Returns:
        DesignResult with ranked guides.

    Example:
        >>> # Match Chang et al. 2022 parameters
        >>> result = design_crispra_guides(
        ...     gene_symbol="Rai1",
        ...     promoter_sequence=rai1_promoter,
        ...     tss_position=600,
        ...     nuclease=Nuclease.SACAS9,
        ...     guide_length=20,  # They used 20bp guides
        ... )
        >>> for guide in result.tier_a_guides[:5]:
        ...     print(f"{guide['sequence']} (Tier A)")
    """
    if promoter_sequence is None:
        raise ValueError(
            "promoter_sequence is required. "
            "PhaseLab does not fetch sequences - provide the promoter DNA."
        )

    if tss_position is None:
        raise ValueError(
            "tss_position is required. "
            "Specify the TSS position within your promoter_sequence."
        )

    # v0.9.3: CRISPRa uses BINDING mode (relaxed PAM)
    role = NucleaseRole.BINDING if relaxed_pam else NucleaseRole.CUTTING

    logger.info(f"Designing CRISPRa guides for {gene_symbol} (PAM mode: {role.value}, guide_length: {guide_length or 'default'})")

    # Step 1: Enumerate candidates with appropriate PAM mode
    enumeration = enumerate_from_sequence(
        sequence=promoter_sequence,
        gene_symbol=gene_symbol,
        tss_position=tss_position,
        nuclease=nuclease,
        role=role,  # v0.9.3: Pass role for PAM selection
        guide_length=guide_length,  # v0.9.3: Allow guide length override
    )

    logger.info(f"Enumerated {enumeration.total_candidates} candidates")

    # Step 2: Evaluate candidates
    evaluated = evaluate_candidates(enumeration.candidates, check_u6=check_u6)

    logger.info(f"Evaluated {len(evaluated)} candidates")

    # Step 3: Apply policy ranking
    # Note: Without CRISPOR, off_targets is empty, so all guides will pass OT gates
    # This is a LIMITATION - for real ranking, integrate CRISPOR
    ranking_result = rank_guides(evaluated, policy, return_excluded=True)

    logger.info(f"Ranked {len(ranking_result['ranked'])} guides")

    # Build manifest
    manifest = {
        'version': '0.9.3',
        'timestamp': datetime.now().isoformat(),
        'gene_symbol': gene_symbol,
        'genome_build': genome_build.value,
        'modality': Modality.CRISPRA.value,
        'nuclease': nuclease.value,
        'nuclease_role': role.value,  # v0.9.3: Track PAM mode
        'policy': policy.value,
        'tss_position': tss_position,
        'sequence_length': len(promoter_sequence),
        'total_enumerated': enumeration.total_candidates,
        'total_ranked': len(ranking_result['ranked']),
        'tier_counts': {
            'A': len(ranking_result['tiers']['A']),
            'B': len(ranking_result['tiers']['B']),
            'C': len(ranking_result['tiers']['C']),
        },
        'crispor_validated': False,  # Mark that this lacks real OT data
        'relaxed_pam': relaxed_pam,
        'guide_length': guide_length,  # v0.9.3: Track guide length override
    }

    return DesignResult(
        ranked_guides=ranking_result['ranked'],
        tiers=ranking_result['tiers'],
        gene_symbol=gene_symbol,
        genome_build=genome_build,
        modality=Modality.CRISPRA,
        nuclease=nuclease,
        policy=policy,
        region_set=enumeration.region_set,
        enumeration=enumeration,
        total_enumerated=enumeration.total_candidates,
        total_after_gates=len(ranking_result['ranked']),
        manifest=manifest,
    )


def design_knockout_guides(
    gene_symbol: str,
    exon_sequence: str,
    nuclease: Nuclease = Nuclease.SPCAS9,
    policy: RankingPolicy = RankingPolicy.CUTTING_STRICT,
    check_u6: bool = True,
) -> DesignResult:
    """
    Design knockout guides for a target gene.

    For knockout, we typically target early exons.

    Args:
        gene_symbol: Target gene name.
        exon_sequence: Exon DNA sequence to target.
        nuclease: CRISPR nuclease.
        policy: Ranking policy (default: CUTTING_STRICT for safety).
        check_u6: Check U6/Pol III compatibility.

    Returns:
        DesignResult with ranked guides.
    """
    logger.info(f"Designing knockout guides for {gene_symbol}")

    # Enumerate
    enumeration = enumerate_from_sequence(
        sequence=exon_sequence,
        gene_symbol=gene_symbol,
        tss_position=0,  # No TSS concept for exon targeting
        nuclease=nuclease,
    )

    # Evaluate
    evaluated = evaluate_candidates(enumeration.candidates, check_u6=check_u6)

    # Rank with strict policy
    ranking_result = rank_guides(evaluated, policy, return_excluded=True)

    manifest = {
        'version': '0.9.2',
        'timestamp': datetime.now().isoformat(),
        'gene_symbol': gene_symbol,
        'modality': Modality.KNOCKOUT.value,
        'nuclease': nuclease.value,
        'policy': policy.value,
        'sequence_length': len(exon_sequence),
        'total_enumerated': enumeration.total_candidates,
        'total_ranked': len(ranking_result['ranked']),
    }

    return DesignResult(
        ranked_guides=ranking_result['ranked'],
        tiers=ranking_result['tiers'],
        gene_symbol=gene_symbol,
        genome_build=GenomeBuild.HG38,
        modality=Modality.KNOCKOUT,
        nuclease=nuclease,
        policy=policy,
        enumeration=enumeration,
        total_enumerated=enumeration.total_candidates,
        total_after_gates=len(ranking_result['ranked']),
        manifest=manifest,
    )


# =============================================================================
# BENCHMARK MODE
# =============================================================================

@dataclass
class BenchmarkResult:
    """Result of benchmarking PhaseLab against published guides."""
    gene_symbol: str
    published_guides: List[Dict[str, Any]]
    phaselab_rankings: List[Dict[str, Any]]

    # Pass/fail
    passed: bool
    pass_rate: float
    details: List[str]

    def summary(self) -> str:
        lines = [
            "=" * 70,
            f"BENCHMARK RESULT: {self.gene_symbol}",
            "=" * 70,
            f"Published guides tested: {len(self.published_guides)}",
            f"Pass rate: {self.pass_rate:.1%}",
            f"Overall: {'PASS' if self.passed else 'FAIL'}",
            "",
            "Details:",
        ]
        for d in self.details:
            lines.append(f"  {d}")
        lines.append("=" * 70)
        return '\n'.join(lines)


def benchmark_against_published(
    published_guides: List[Dict[str, Any]],
    design_result: DesignResult,
    require_tier_a: bool = False,
    require_top_n: int = 10,
) -> BenchmarkResult:
    """
    Benchmark PhaseLab rankings against published experimental results.

    Args:
        published_guides: List of dicts with 'sequence' and 'experimental_winner'.
        design_result: DesignResult from design_crispra_guides or similar.
        require_tier_a: Require winners to be Tier A.
        require_top_n: Require winners to be in top N.

    Returns:
        BenchmarkResult with pass/fail analysis.

    Example:
        >>> published = [
        ...     {'sequence': 'CCTGGCACCCGAGGCCACGA', 'experimental_winner': True},
        ...     {'sequence': 'GTCTAAGTCCAAAATCCTCA', 'experimental_winner': False},
        ... ]
        >>> result = benchmark_against_published(published, design_result)
        >>> print(result.summary())
    """
    details = []
    passes = 0
    total_winners = 0

    # Index PhaseLab results by sequence
    phaselab_index = {
        g['sequence'].upper(): g
        for g in design_result.ranked_guides
    }

    for pub in published_guides:
        seq = pub['sequence'].upper()
        is_winner = pub.get('experimental_winner', False)

        if not is_winner:
            continue

        total_winners += 1

        if seq not in phaselab_index:
            details.append(f"FAIL: {seq[:15]}... not found in PhaseLab output")
            continue

        pl_guide = phaselab_index[seq]
        rank = pl_guide.get('rank', 999)
        tier = pl_guide.get('tier', 'X')

        # Check criteria
        tier_ok = (tier == 'A') if require_tier_a else (tier in ['A', 'B', 'C'])
        rank_ok = rank <= require_top_n

        if tier_ok and rank_ok:
            passes += 1
            details.append(f"PASS: {seq[:15]}... ranked #{rank}, Tier {tier}")
        elif not tier_ok:
            details.append(f"FAIL: {seq[:15]}... Tier {tier} (expected A)")
        else:
            details.append(f"WARN: {seq[:15]}... ranked #{rank} > top {require_top_n}")

    pass_rate = passes / total_winners if total_winners > 0 else 0
    passed = pass_rate >= 0.8  # 80% threshold

    return BenchmarkResult(
        gene_symbol=design_result.gene_symbol,
        published_guides=published_guides,
        phaselab_rankings=design_result.ranked_guides,
        passed=passed,
        pass_rate=pass_rate,
        details=details,
    )


# =============================================================================
# CONVENIENCE EXPORTS
# =============================================================================

__all__ = [
    # Main design functions
    'design_crispra_guides',
    'design_knockout_guides',

    # Results
    'DesignResult',
    'BenchmarkResult',

    # Evaluation
    'evaluate_candidate',
    'evaluate_candidates',

    # Benchmarking
    'benchmark_against_published',

    # Re-exports from region
    'GenomeBuild',
    'Modality',
    'Window',
    'TSSAnnotation',
    'RegionBuilder',

    # Re-exports from enumerate
    'Nuclease',
    'NucleaseRole',  # v0.9.3: BINDING vs CUTTING
    'CandidateGuide',

    # Re-exports from scoring
    'RankingPolicy',
    'GuideTier',
]
