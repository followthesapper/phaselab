"""
PhaseLab Quantum: ATLAS-Q integration for advanced quantum simulation.

This module provides integration with ATLAS-Q for:
- IR measurement grouping (5× variance reduction)
- Real circular statistics coherence (vs heuristic)
- VQE optimization for guide validation
- GPU acceleration (when available)
- Rust backend for fast simulation

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
    "is_atlas_q_available",
    "get_atlas_q_version",
    "grouping",
    "coherence",
    "vqe",
    "backend",
]
