"""
PhaseLab Quantum: Configuration for quantum computation modes.

Provides three modes for quantum/classical computation:

1. OFF (default): Classical-only computation
   - No quantum simulation
   - Fastest, no special requirements
   - Use when quantum features not needed

2. AUDIT: Classical with quantum validation
   - Primary computation is classical
   - Quantum used to audit/validate subset of results
   - Good for development and validation

3. REQUIRED: Quantum computation mandatory
   - All coherence computed via quantum simulation
   - Requires atlas-quantum or IBM Quantum access
   - Use for research-grade analysis

The mode can be set globally or per-analysis.
"""

from enum import Enum
from typing import Optional, Dict, Any
from dataclasses import dataclass
import os
import logging

logger = logging.getLogger(__name__)


class QuantumMode(Enum):
    """
    Quantum computation mode.

    OFF: Classical-only (default, fastest)
    AUDIT: Classical + quantum validation subset
    REQUIRED: Quantum mandatory for all coherence
    """
    OFF = "off"
    AUDIT = "audit"
    REQUIRED = "required"


@dataclass
class QuantumConfig:
    """
    Configuration for quantum computation.

    Attributes:
        mode: Quantum computation mode (off/audit/required)
        backend: Quantum backend to use (simulator/ibm_torino/etc.)
        shots: Number of measurement shots
        audit_fraction: Fraction of results to audit (for AUDIT mode)
        error_on_unavailable: Raise error if required quantum unavailable
        ibm_token: IBM Quantum API token (optional)
    """
    mode: QuantumMode = QuantumMode.OFF
    backend: str = "simulator"
    shots: int = 4096
    audit_fraction: float = 0.1
    error_on_unavailable: bool = True
    ibm_token: Optional[str] = None

    def __post_init__(self):
        """Validate configuration."""
        if self.mode == QuantumMode.REQUIRED:
            if not self.is_quantum_available():
                msg = (
                    "Quantum mode REQUIRED but quantum backend unavailable. "
                    "Install atlas-quantum or set IBMQ_TOKEN environment variable."
                )
                if self.error_on_unavailable:
                    raise RuntimeError(msg)
                else:
                    logger.warning(msg)

    def is_quantum_available(self) -> bool:
        """Check if quantum computation is available."""
        # Check for atlas-quantum
        try:
            import atlas_q
            return True
        except ImportError:
            pass

        # Check for IBM Quantum token
        if self.ibm_token or os.getenv("IBMQ_TOKEN"):
            return True

        return False

    def should_use_quantum(self) -> bool:
        """Determine if quantum should be used for this config."""
        if self.mode == QuantumMode.OFF:
            return False
        if self.mode == QuantumMode.REQUIRED:
            return True
        # AUDIT mode - quantum available and configured
        return self.is_quantum_available()

    def get_backend_info(self) -> Dict[str, Any]:
        """Get information about configured backend."""
        return {
            'mode': self.mode.value,
            'backend': self.backend,
            'shots': self.shots,
            'quantum_available': self.is_quantum_available(),
            'will_use_quantum': self.should_use_quantum(),
        }


# Global configuration
_global_config: Optional[QuantumConfig] = None


def get_quantum_config() -> QuantumConfig:
    """Get the global quantum configuration."""
    global _global_config
    if _global_config is None:
        _global_config = QuantumConfig()
    return _global_config


def set_quantum_config(config: QuantumConfig) -> None:
    """Set the global quantum configuration."""
    global _global_config
    _global_config = config


def set_quantum_mode(mode) -> None:
    """
    Set the global quantum mode.

    Args:
        mode: QuantumMode enum or string ("off", "audit", or "required")

    Example:
        >>> from phaselab.quantum import set_quantum_mode, QuantumMode
        >>> set_quantum_mode("audit")  # String form
        >>> set_quantum_mode(QuantumMode.AUDIT)  # Enum form
    """
    # Accept both QuantumMode enum and string
    if isinstance(mode, QuantumMode):
        mode_enum = mode
    elif isinstance(mode, str):
        mode_enum = QuantumMode(mode.lower())
    else:
        raise TypeError(f"mode must be QuantumMode or str, got {type(mode).__name__}")

    global _global_config
    if _global_config is None:
        _global_config = QuantumConfig(mode=mode_enum)
    else:
        _global_config.mode = mode_enum


def get_quantum_mode() -> QuantumMode:
    """Get the current quantum mode as QuantumMode enum."""
    return get_quantum_config().mode


def configure_quantum(
    mode: str = "off",
    backend: str = "simulator",
    shots: int = 4096,
    ibm_token: Optional[str] = None,
) -> QuantumConfig:
    """
    Configure quantum computation settings.

    Args:
        mode: "off", "audit", or "required"
        backend: Backend name (simulator, ibm_torino, etc.)
        shots: Number of measurement shots
        ibm_token: IBM Quantum API token

    Returns:
        QuantumConfig object

    Example:
        >>> config = configure_quantum(
        ...     mode="required",
        ...     backend="ibm_torino",
        ...     shots=8192,
        ... )
    """
    config = QuantumConfig(
        mode=QuantumMode(mode.lower()),
        backend=backend,
        shots=shots,
        ibm_token=ibm_token or os.getenv("IBMQ_TOKEN"),
    )
    set_quantum_config(config)
    return config


def quantum_status() -> Dict[str, Any]:
    """
    Get comprehensive quantum status.

    Returns:
        Dictionary with quantum configuration and availability.
    """
    config = get_quantum_config()

    # Check atlas-quantum
    atlas_available = False
    atlas_version = None
    try:
        import atlas_q
        atlas_available = True
        atlas_version = atlas_q.__version__
    except ImportError:
        pass

    # Check IBM Quantum
    ibm_available = bool(config.ibm_token or os.getenv("IBMQ_TOKEN"))

    return {
        'config': config.get_backend_info(),
        'backends': {
            'atlas_q': {
                'available': atlas_available,
                'version': atlas_version,
            },
            'ibm_quantum': {
                'available': ibm_available,
                'token_configured': ibm_available,
            },
        },
        'recommendation': _get_recommendation(config, atlas_available, ibm_available),
    }


def _get_recommendation(
    config: QuantumConfig,
    atlas_available: bool,
    ibm_available: bool,
) -> str:
    """Generate recommendation based on current status."""
    if config.mode == QuantumMode.OFF:
        return "Quantum disabled. Use set_quantum_mode('audit') to enable validation."

    if config.mode == QuantumMode.AUDIT:
        if atlas_available:
            return "AUDIT mode active with ATLAS-Q. A subset of results will be quantum-validated."
        elif ibm_available:
            return "AUDIT mode active with IBM Quantum. A subset of results will be quantum-validated."
        else:
            return "AUDIT mode but no quantum backend. Results will be classical-only."

    if config.mode == QuantumMode.REQUIRED:
        if atlas_available:
            return "REQUIRED mode with ATLAS-Q. All coherence computed via quantum simulation."
        elif ibm_available:
            return "REQUIRED mode with IBM Quantum. All coherence computed via quantum hardware."
        else:
            return "WARNING: REQUIRED mode but no quantum backend available!"

    return "Unknown configuration state."
