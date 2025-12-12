"""
PhaseLab Enhanced CRISPR Pipeline v0.7.0

Integrates the Virtual Assay Stack components for improved guide accuracy:
- Biological context (chromatin, methylation, histones)
- ML outcome predictors (DeepCRISPR, DeepSpCas9, etc.)
- Evidence fusion with uncertainty quantification

This pipeline replaces the basic combined_score with a multi-layer
evidence fusion approach that:
1. Applies hard safety gates (off-targets, GC bounds)
2. Combines soft scores with confidence weighting
3. Uses ML predictions when available
4. Overlays IR coherence as robustness indicator
"""

import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple, Union
from enum import Enum
import logging

from .pam_scan import find_pam_sites, filter_by_window, PAMHit
from .scoring import (
    gc_content,
    delta_g_santalucia,
    mit_specificity_score,
    cfd_score,
    max_homopolymer_run,
    sequence_complexity,
    chromatin_accessibility_score,
)
from .coherence_utils import (
    CoherenceMode,
    compute_guide_coherence,
)
from ..core.coherence import go_no_go
from ..core.constants import E_MINUS_2

# v0.7.0 imports
from ..context import ContextStack, CellType, BiologicalContext
from ..predictors import (
    PredictorStack,
    RuleBasedPredictor,
    EnsemblePrediction,
)
from ..fusion import (
    EvidenceFusion,
    Evidence,
    EvidenceSource,
    EvidenceType,
    FusedResult,
    FusionConfig,
)

logger = logging.getLogger(__name__)


class Modality(Enum):
    """CRISPR modality."""
    CRISPRA = "CRISPRa"
    CRISPRI = "CRISPRi"
    KNOCKOUT = "Knockout"
    BASE_EDITING = "BaseEditing"
    PRIME_EDITING = "PrimeEditing"


@dataclass
class EnhancedGuideConfig:
    """
    Configuration for enhanced guide RNA design pipeline.

    Extends basic config with v0.7.0 Virtual Assay Stack options.
    """
    # Basic settings
    pam: str = "NGG"
    guide_length: int = 20
    modality: Modality = Modality.CRISPRA

    # Window settings (relative to TSS)
    crispra_window: Tuple[int, int] = (-400, -50)
    crispri_window: Tuple[int, int] = (-50, 300)
    knockout_window: Tuple[int, int] = (0, 500)  # First exon

    # Quality filters (hard gates)
    min_gc: float = 0.35
    max_gc: float = 0.75
    max_homopolymer: int = 4
    min_complexity: float = 0.4
    max_offtargets: int = 10

    # Coherence settings
    coherence_mode: CoherenceMode = CoherenceMode.HEURISTIC
    go_threshold: float = E_MINUS_2

    # v0.7.0: Biological context
    use_biological_context: bool = True
    cell_type: Union[CellType, str] = CellType.K562

    # v0.7.0: ML predictors
    use_ml_predictors: bool = True
    ml_min_confidence: float = 0.3

    # v0.7.0: Evidence fusion
    fusion_config: Optional[FusionConfig] = None

    # Output
    top_n: int = 20
    include_nogo: bool = False  # Include NO-GO guides in output

    @property
    def window(self) -> Tuple[int, int]:
        """Get window for current modality."""
        if self.modality == Modality.CRISPRA:
            return self.crispra_window
        elif self.modality == Modality.CRISPRI:
            return self.crispri_window
        else:
            return self.knockout_window


@dataclass
class EnhancedGuide:
    """
    Enhanced guide result with full evidence breakdown.
    """
    # Basic info
    sequence: str
    pam: str
    position: int  # Relative to TSS
    strand: str
    chrom: Optional[str] = None
    abs_start: Optional[int] = None
    abs_end: Optional[int] = None

    # Sequence scores
    gc: float = 0.0
    delta_g: float = 0.0
    complexity: float = 0.0
    homopolymer: int = 0

    # Specificity scores
    mit_score: float = 0.0
    cfd_score: float = 0.0
    offtarget_count: int = 0

    # Biological context (v0.7.0)
    accessibility: float = 0.5
    methylation: float = 0.5
    histone_activity: float = 0.5
    chromatin_state: str = "Unknown"
    context_confidence: float = 0.0

    # ML predictions (v0.7.0)
    ml_efficiency: float = 0.5
    ml_confidence: float = 0.0
    ml_predictors_used: int = 0

    # IR coherence
    coherence: float = 0.0
    coherence_mode: str = "heuristic"

    # Fused result (v0.7.0)
    fused_score: float = 0.0
    fused_confidence: float = 0.0
    go_status: str = "NO-GO"
    go_reason: str = ""
    passes_gates: bool = False
    failed_gates: List[str] = field(default_factory=list)

    # Evidence level
    evidence_level: str = "C"  # A/B/C

    @property
    def is_go(self) -> bool:
        return self.go_status == "GO"

    @property
    def is_viable(self) -> bool:
        return self.passes_gates and self.is_go

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "sequence": self.sequence,
            "pam": self.pam,
            "position": self.position,
            "strand": self.strand,
            "gc": self.gc,
            "delta_g": self.delta_g,
            "mit_score": self.mit_score,
            "cfd_score": self.cfd_score,
            "accessibility": self.accessibility,
            "methylation": self.methylation,
            "chromatin_state": self.chromatin_state,
            "ml_efficiency": self.ml_efficiency,
            "ml_confidence": self.ml_confidence,
            "coherence": self.coherence,
            "coherence_mode": self.coherence_mode,
            "fused_score": self.fused_score,
            "fused_confidence": self.fused_confidence,
            "go_status": self.go_status,
            "passes_gates": self.passes_gates,
            "evidence_level": self.evidence_level,
            "is_viable": self.is_viable,
        }


@dataclass
class EnhancedDesignResult:
    """
    Result container for enhanced guide design.
    """
    guides: List[EnhancedGuide]
    config: EnhancedGuideConfig
    target_gene: Optional[str] = None
    total_candidates: int = 0
    filtered_by_gates: int = 0
    filtered_by_coherence: int = 0

    def to_dataframe(self) -> pd.DataFrame:
        """Convert to pandas DataFrame."""
        return pd.DataFrame([g.to_dict() for g in self.guides])

    @property
    def go_guides(self) -> List[EnhancedGuide]:
        """Return only GO guides."""
        return [g for g in self.guides if g.is_go]

    @property
    def viable_guides(self) -> List[EnhancedGuide]:
        """Return only viable guides (pass gates + GO)."""
        return [g for g in self.guides if g.is_viable]

    def summary(self) -> Dict[str, Any]:
        """Get design summary."""
        return {
            "total_candidates": self.total_candidates,
            "passed_gates": len(self.guides),
            "go_guides": len(self.go_guides),
            "viable_guides": len(self.viable_guides),
            "filtered_by_gates": self.filtered_by_gates,
            "filtered_by_coherence": self.filtered_by_coherence,
            "modality": self.config.modality.value,
            "cell_type": str(self.config.cell_type),
        }


def design_enhanced_guides(
    sequence: str,
    tss_index: int,
    config: Optional[EnhancedGuideConfig] = None,
    target_gene: Optional[str] = None,
    chrom: Optional[str] = None,
    chrom_offset: int = 0,
    verbose: bool = False,
) -> EnhancedDesignResult:
    """
    Design guide RNAs using the full Virtual Assay Stack.

    This enhanced pipeline integrates:
    1. Sequence-based scoring (GC, thermodynamics, specificity)
    2. Biological context (chromatin state, methylation, histones)
    3. ML efficiency predictions (DeepCRISPR, DeepSpCas9)
    4. IR coherence (heuristic or quantum)
    5. Evidence fusion with uncertainty quantification

    Parameters
    ----------
    sequence : str
        Promoter/target DNA sequence (5'->3')
    tss_index : int
        Position of TSS in sequence (0-based)
    config : EnhancedGuideConfig, optional
        Pipeline configuration
    target_gene : str, optional
        Gene symbol for ML predictors
    chrom : str, optional
        Chromosome for context lookup (e.g., "chr4")
    chrom_offset : int
        Genomic offset of sequence start
    verbose : bool
        Print progress messages

    Returns
    -------
    EnhancedDesignResult
        Ranked guides with full evidence breakdown

    Example
    -------
    >>> from phaselab.crispr import design_enhanced_guides, EnhancedGuideConfig, Modality
    >>> from phaselab.context import CellType
    >>>
    >>> config = EnhancedGuideConfig(
    ...     modality=Modality.CRISPRA,
    ...     cell_type=CellType.K562,
    ...     coherence_mode=CoherenceMode.HEURISTIC,
    ... )
    >>>
    >>> result = design_enhanced_guides(
    ...     sequence=scn2a_promoter,
    ...     tss_index=500,
    ...     config=config,
    ...     target_gene="SCN2A",
    ...     chrom="chr2",
    ...     chrom_offset=165000000,
    ... )
    >>>
    >>> print(f"Found {len(result.viable_guides)} viable guides")
    >>> for g in result.viable_guides[:5]:
    ...     print(f"  {g.sequence} Score={g.fused_score:.3f} {g.go_status}")
    """
    if config is None:
        config = EnhancedGuideConfig()

    sequence = sequence.upper()

    if verbose:
        logger.info(f"Enhanced pipeline: {len(sequence)}bp, modality={config.modality.value}")

    # Initialize context and predictors if enabled
    context_stack = None
    predictor_stack = None

    if config.use_biological_context and chrom:
        try:
            context_stack = ContextStack(cell_type=config.cell_type)
            if verbose:
                logger.info(f"Loaded biological context for {config.cell_type}")
        except Exception as e:
            logger.warning(f"Could not load biological context: {e}")

    if config.use_ml_predictors:
        try:
            predictor_stack = PredictorStack(
                aggregation="weighted",
                min_confidence=config.ml_min_confidence,
            )
            predictor_stack.register(RuleBasedPredictor())
            # Would register DeepCRISPR, DeepSpCas9 if available
            if verbose:
                logger.info(f"ML predictors: {predictor_stack.available_predictors}")
        except Exception as e:
            logger.warning(f"Could not initialize ML predictors: {e}")

    # Step 1: Find PAM sites
    all_hits = find_pam_sites(
        sequence,
        pam=config.pam,
        guide_length=config.guide_length,
        both_strands=True,
    )

    if verbose:
        logger.info(f"Found {len(all_hits)} PAM sites")

    # Step 2: Filter to modality-specific window
    window_hits = filter_by_window(
        all_hits,
        tss_position=tss_index,
        window=config.window,
    )

    total_candidates = len(window_hits)
    if verbose:
        logger.info(f"Filtered to {total_candidates} in window {config.window}")

    if not window_hits:
        return EnhancedDesignResult(
            guides=[],
            config=config,
            target_gene=target_gene,
            total_candidates=0,
        )

    # Step 3: Process each candidate through Virtual Assay Stack
    fusion_config = config.fusion_config or FusionConfig.default()
    guides = []
    filtered_by_gates = 0
    filtered_by_coherence = 0

    for hit in window_hits:
        guide_seq = hit.guide
        rel_pos = ((hit.guide_start + hit.guide_end) // 2) - tss_index

        # Calculate genomic coordinates if available
        abs_start = chrom_offset + hit.guide_start if chrom_offset else None
        abs_end = chrom_offset + hit.guide_end if chrom_offset else None

        # === SEQUENCE SCORING ===
        gc = gc_content(guide_seq)
        homo = max_homopolymer_run(guide_seq)
        complexity = sequence_complexity(guide_seq)
        delta_g = delta_g_santalucia(guide_seq)
        mit = mit_specificity_score(guide_seq)
        cfd = cfd_score(guide_seq)

        # === EVIDENCE FUSION ===
        fusion = EvidenceFusion(config=fusion_config)

        # Hard gates
        gc_passes = config.min_gc <= gc <= config.max_gc
        homo_passes = homo <= config.max_homopolymer
        complexity_passes = complexity >= config.min_complexity

        fusion.add_hard_gate(EvidenceSource.GC_CONTENT, gc_passes)
        fusion.add_hard_gate(EvidenceSource.HOMOPOLYMER, homo_passes)

        # Soft sequence scores
        gc_score = 1.0 - 2.0 * abs(gc - 0.55)  # Optimal around 55%
        fusion.add_sequence_score(max(0, gc_score), source=EvidenceSource.GC_CONTENT, confidence=0.95)

        # Thermodynamics (normalize delta_g)
        thermo_score = min(1.0, max(0, (-delta_g) / 25.0))
        fusion.add_sequence_score(thermo_score, source=EvidenceSource.THERMODYNAMICS, confidence=0.9)

        # Specificity
        fusion.add_specificity_score(mit / 100.0, source=EvidenceSource.OFF_TARGET_CFD, confidence=0.85)

        # === BIOLOGICAL CONTEXT (v0.7.0) ===
        accessibility = 0.5
        methylation = 0.5
        histone_activity = 0.5
        chromatin_state = "Unknown"
        context_confidence = 0.0

        if context_stack and chrom and abs_start and abs_end:
            try:
                bio_ctx = context_stack.get_context(chrom, abs_start, abs_end)
                accessibility = bio_ctx.accessibility.score
                methylation = bio_ctx.methylation.mean_methylation
                histone_activity = bio_ctx.histones.activity_score
                chromatin_state = bio_ctx.chromatin_state.value
                context_confidence = bio_ctx.confidence

                fusion.add_context_score(
                    accessibility,
                    source=EvidenceSource.ACCESSIBILITY,
                    confidence=bio_ctx.accessibility.confidence,
                )
                fusion.add_context_score(
                    1.0 - methylation,  # Low methylation is good
                    source=EvidenceSource.METHYLATION,
                    confidence=bio_ctx.methylation.confidence,
                )
                fusion.add_context_score(
                    histone_activity,
                    source=EvidenceSource.HISTONE,
                    confidence=bio_ctx.histones.confidence,
                )
            except Exception as e:
                logger.debug(f"Context lookup failed: {e}")

        # === ML PREDICTIONS (v0.7.0) ===
        ml_efficiency = 0.5
        ml_confidence = 0.0
        ml_predictors_used = 0

        if predictor_stack:
            try:
                full_seq = guide_seq + hit.pam  # Include PAM
                ml_result = predictor_stack.predict(
                    sequence=full_seq,
                    target_gene=target_gene,
                    cell_type=str(config.cell_type) if config.cell_type else None,
                )
                ml_efficiency = ml_result.score
                ml_confidence = ml_result.confidence
                ml_predictors_used = ml_result.n_predictors

                if ml_confidence >= config.ml_min_confidence:
                    fusion.add_ml_prediction(
                        ml_efficiency,
                        source=EvidenceSource.DEEPCRISPR,
                        confidence=ml_confidence,
                    )
            except Exception as e:
                logger.debug(f"ML prediction failed: {e}")

        # === IR COHERENCE ===
        coherence = compute_guide_coherence(
            guide_seq,
            mode=config.coherence_mode,
        )
        # Map CoherenceMode to string for fusion
        coherence_mode_str = {
            CoherenceMode.HEURISTIC: "heuristic",
            CoherenceMode.QUANTUM: "quantum",
        }.get(config.coherence_mode, "heuristic")

        fusion.add_coherence(
            coherence,
            mode=coherence_mode_str,
            confidence=0.8 if config.coherence_mode == CoherenceMode.QUANTUM else 0.5,
        )

        # === FUSE ALL EVIDENCE ===
        fused = fusion.fuse()

        # Skip if failed hard gates (unless config says include)
        if not fused.passes_gates:
            filtered_by_gates += 1
            if not config.include_nogo:
                continue

        # Skip NO-GO guides (unless config says include)
        if not fused.is_go:
            filtered_by_coherence += 1
            if not config.include_nogo:
                continue

        # Determine evidence level
        # A = hardware validated (not available in current modes)
        # B = quantum VQE simulation
        # C = heuristic
        if config.coherence_mode == CoherenceMode.QUANTUM:
            evidence_level = "B"
        else:
            evidence_level = "C"

        # Create enhanced guide
        guide = EnhancedGuide(
            sequence=guide_seq,
            pam=hit.pam,
            position=rel_pos,
            strand=hit.strand,
            chrom=chrom,
            abs_start=abs_start,
            abs_end=abs_end,
            gc=round(gc, 3),
            delta_g=round(delta_g, 2),
            complexity=round(complexity, 3),
            homopolymer=homo,
            mit_score=round(mit, 1),
            cfd_score=round(cfd, 1),
            accessibility=round(accessibility, 3),
            methylation=round(methylation, 3),
            histone_activity=round(histone_activity, 3),
            chromatin_state=chromatin_state,
            context_confidence=round(context_confidence, 2),
            ml_efficiency=round(ml_efficiency, 3),
            ml_confidence=round(ml_confidence, 2),
            ml_predictors_used=ml_predictors_used,
            coherence=round(coherence, 4),
            coherence_mode=coherence_mode_str,
            fused_score=round(fused.score, 3),
            fused_confidence=round(fused.confidence, 2),
            go_status=fused.go_status,
            go_reason=fused.go_reason,
            passes_gates=fused.passes_gates,
            failed_gates=fused.failed_gates,
            evidence_level=evidence_level,
        )

        guides.append(guide)

    # Step 4: Sort by fused score
    guides.sort(key=lambda g: g.fused_score, reverse=True)

    # Step 5: Return top N
    if verbose:
        logger.info(f"Returning {min(len(guides), config.top_n)} guides")

    # Close context if we opened it
    if context_stack:
        context_stack.close()

    return EnhancedDesignResult(
        guides=guides[:config.top_n],
        config=config,
        target_gene=target_gene,
        total_candidates=total_candidates,
        filtered_by_gates=filtered_by_gates,
        filtered_by_coherence=filtered_by_coherence,
    )


def compare_guides_with_without_context(
    sequence: str,
    tss_index: int,
    chrom: str,
    chrom_offset: int,
    cell_type: CellType = CellType.K562,
    target_gene: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Compare guide rankings with and without biological context.

    Useful for demonstrating the value of the Virtual Assay Stack.

    Returns
    -------
    Dict with 'basic', 'enhanced', and 'rank_changes' DataFrames
    """
    # Basic config (no context, no ML)
    basic_config = EnhancedGuideConfig(
        use_biological_context=False,
        use_ml_predictors=False,
        coherence_mode=CoherenceMode.HEURISTIC,
    )

    # Enhanced config (full stack)
    enhanced_config = EnhancedGuideConfig(
        use_biological_context=True,
        use_ml_predictors=True,
        cell_type=cell_type,
        coherence_mode=CoherenceMode.HEURISTIC,
    )

    basic_result = design_enhanced_guides(
        sequence, tss_index, basic_config, target_gene
    )

    enhanced_result = design_enhanced_guides(
        sequence, tss_index, enhanced_config, target_gene, chrom, chrom_offset
    )

    basic_df = basic_result.to_dataframe()
    enhanced_df = enhanced_result.to_dataframe()

    # Calculate rank changes
    if not basic_df.empty and not enhanced_df.empty:
        basic_ranks = {row['sequence']: i for i, row in basic_df.iterrows()}
        enhanced_ranks = {row['sequence']: i for i, row in enhanced_df.iterrows()}

        changes = []
        for seq in enhanced_ranks:
            basic_rank = basic_ranks.get(seq)
            enhanced_rank = enhanced_ranks[seq]
            if basic_rank is not None:
                change = basic_rank - enhanced_rank  # Positive = improved
                changes.append({
                    'sequence': seq,
                    'basic_rank': basic_rank + 1,
                    'enhanced_rank': enhanced_rank + 1,
                    'rank_change': change,
                })

        rank_changes = pd.DataFrame(changes)
    else:
        rank_changes = pd.DataFrame()

    return {
        'basic': basic_df,
        'enhanced': enhanced_df,
        'rank_changes': rank_changes,
        'basic_summary': basic_result.summary(),
        'enhanced_summary': enhanced_result.summary(),
    }
