"""
ML Outcome Predictor Protocol for PhaseLab Virtual Assay Stack.

This module provides a protocol-based interface for integrating external
ML prediction models (DeepCRISPR, DeepSpCas9, Enformer, etc.) into
the PhaseLab scoring pipeline.

Design Principles
-----------------
1. **Protocol-based**: Predictors implement a common interface
2. **Lazy loading**: Models load only when called
3. **Graceful fallback**: Missing dependencies return neutral scores
4. **Confidence tracking**: Each prediction includes confidence bounds

Components
----------
- OutcomePredictor: Abstract protocol for ML predictors
- PredictionResult: Standardized prediction output
- DeepCRISPRAdapter: Adapter for DeepCRISPR models
- DeepSpCas9Adapter: Adapter for DeepSpCas9
- EnformerAdapter: Adapter for Enformer expression predictions

Example
-------
>>> from phaselab.predictors import PredictorStack, DeepCRISPRAdapter
>>>
>>> # Create predictor stack
>>> stack = PredictorStack()
>>> stack.register(DeepCRISPRAdapter())
>>>
>>> # Get predictions
>>> result = stack.predict(sequence="ATGCGATCGATCGATCGATCNGG", target_gene="CLOCK")
>>> print(f"Efficiency: {result.efficiency:.2f}")
>>> print(f"Confidence: {result.confidence:.2f}")
"""

from phaselab.predictors.protocol import (
    OutcomePredictor,
    PredictionResult,
    PredictorMetadata,
)
from phaselab.predictors.adapters import (
    DeepCRISPRAdapter,
    DeepSpCas9Adapter,
    EnformerAdapter,
    RuleBasedPredictor,
)
from phaselab.predictors.stack import (
    PredictorStack,
    EnsemblePrediction,
)

__all__ = [
    # Protocol
    "OutcomePredictor",
    "PredictionResult",
    "PredictorMetadata",
    # Adapters
    "DeepCRISPRAdapter",
    "DeepSpCas9Adapter",
    "EnformerAdapter",
    "RuleBasedPredictor",
    # Stack
    "PredictorStack",
    "EnsemblePrediction",
]
