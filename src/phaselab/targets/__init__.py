"""
PhaseLab Target Configurations.

Load and manage gene target configurations for CRISPRa experiments.
Each target defines genomic coordinates, CRISPRa parameters, and metadata.
"""

from .config import (
    TargetConfig,
    load_target_config,
    list_available_targets,
    get_target_path,
)

__all__ = [
    "TargetConfig",
    "load_target_config",
    "list_available_targets",
    "get_target_path",
]
