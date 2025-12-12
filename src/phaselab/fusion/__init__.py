"""
Evidence Fusion System for PhaseLab Virtual Assay Stack.

This module reconciles multiple evidence sources into a unified
guide score with proper uncertainty quantification.

Evidence Sources
----------------
1. Sequence-based scores (GC, thermodynamics, position)
2. Biological context (chromatin, methylation, histones)
3. ML predictions (DeepCRISPR, DeepSpCas9, etc.)
4. IR coherence (heuristic or quantum)
5. Off-target analysis

Fusion Strategy
---------------
The fusion system uses a hierarchical approach:

1. **Hard gates**: Binary pass/fail (off-target safety, GC bounds)
2. **Soft scores**: Weighted combination with confidence
3. **Coherence overlay**: IR coherence as reliability indicator

Example
-------
>>> from phaselab.fusion import EvidenceFusion, EvidenceSource
>>>
>>> fusion = EvidenceFusion()
>>>
>>> # Add evidence from different sources
>>> fusion.add_sequence_score(0.75, confidence=0.9)
>>> fusion.add_context_score(0.80, confidence=0.7)
>>> fusion.add_ml_prediction(0.72, confidence=0.85)
>>> fusion.add_coherence(0.45, mode="heuristic")
>>>
>>> # Get fused result
>>> result = fusion.fuse()
>>> print(f"Score: {result.score:.2f}")
>>> print(f"GO status: {result.go_status}")
"""

from phaselab.fusion.evidence import (
    Evidence,
    EvidenceSource,
    EvidenceType,
)
from phaselab.fusion.fuser import (
    EvidenceFusion,
    FusedResult,
    FusionConfig,
)
from phaselab.fusion.calibration import (
    Calibrator,
    CalibrationCurve,
)

__all__ = [
    # Evidence
    "Evidence",
    "EvidenceSource",
    "EvidenceType",
    # Fusion
    "EvidenceFusion",
    "FusedResult",
    "FusionConfig",
    # Calibration
    "Calibrator",
    "CalibrationCurve",
]
