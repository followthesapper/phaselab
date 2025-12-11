"""
PhaseLab I/O: File formats, export utilities, and report generation.

Version: 0.5.0
"""

from .export import (
    export_guides_json,
    export_guides_fasta,
    export_crispor_batch,
)

# v0.5.0: Report generation
from .report import (
    generate_guide_report,
    generate_therapeutic_report,
    generate_coherence_report,
    export_report,
    generate_crispor_link,
    ReportMetadata,
)

__all__ = [
    # Export
    "export_guides_json",
    "export_guides_fasta",
    "export_crispor_batch",
    # Reports (v0.5.0)
    "generate_guide_report",
    "generate_therapeutic_report",
    "generate_coherence_report",
    "export_report",
    "generate_crispor_link",
    "ReportMetadata",
]
