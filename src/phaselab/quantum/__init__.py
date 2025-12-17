"""
PhaseLab Quantum: ATLAS-Q integration and quantum computation modes.

This module provides:
- Quantum mode configuration (off/audit/required)
- Integration with ATLAS-Q for IR measurement grouping
- VQE optimization for coherence validation
- IBM Quantum hardware support
- GPU acceleration (when available)

Quantum Modes (v1.0.0):
- OFF: Classical-only (default, fastest)
- AUDIT: Classical + quantum validation on subset
- REQUIRED: Quantum mandatory for all coherence

All features are optional and gracefully degrade if atlas-quantum
is not installed.
"""

from typing import TYPE_CHECKING

# Lazy imports for optional atlas-quantum dependency
_ATLAS_Q_AVAILABLE = None


def is_atlas_q_available() -> bool:
    """Check if atlas-quantum is installed."""
    global _ATLAS_Q_AVAILABLE
    if _ATLAS_Q_AVAILABLE is None:
        try:
            import atlas_q
            _ATLAS_Q_AVAILABLE = True
        except ImportError:
            _ATLAS_Q_AVAILABLE = False
    return _ATLAS_Q_AVAILABLE


def get_atlas_q_version() -> str:
    """Get atlas-quantum version if installed."""
    if is_atlas_q_available():
        import atlas_q
        return atlas_q.__version__
    return "not installed"


# Configuration exports (eagerly loaded for convenience)
from .config import (
    QuantumMode,
    QuantumConfig,
    get_quantum_config,
    set_quantum_config,
    set_quantum_mode,
    get_quantum_mode,
    configure_quantum,
    quantum_status,
)


# Lazy imports for submodules
def __getattr__(name: str):
    """Lazy import submodules."""
    if name == "grouping":
        from . import grouping
        return grouping
    elif name == "coherence":
        from . import coherence
        return coherence
    elif name == "vqe":
        from . import vqe
        return vqe
    elif name == "backend":
        from . import backend
        return backend
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    # Availability checks
    "is_atlas_q_available",
    "get_atlas_q_version",
    # Quantum mode configuration
    "QuantumMode",
    "QuantumConfig",
    "get_quantum_config",
    "set_quantum_config",
    "set_quantum_mode",
    "get_quantum_mode",
    "configure_quantum",
    "quantum_status",
    # Submodules
    "grouping",
    "coherence",
    "vqe",
    "backend",
]
