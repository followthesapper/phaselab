"""
CRISPOR Output Parsers.

Parse CRISPOR TSV output files into structured data.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List, Dict, Any


@dataclass
class CrisporGuideRow:
    """Parsed guide from CRISPOR output."""

    guide_id: str
    sequence: str  # 20bp guide sequence
    pam: str
    strand: str
    chrom: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None

    # On-target scores
    doench_2016: Optional[float] = None  # Recommended for U6 promoter
    moreno_mateos: Optional[float] = None  # For T7 in vitro
    out_of_frame: Optional[float] = None

    # Specificity scores
    mit_specificity: Optional[float] = None  # 0-100, higher = more specific
    cfd_specificity: Optional[float] = None

    # Off-target counts by mismatch
    ot_0mm: int = 0
    ot_1mm: int = 0
    ot_2mm: int = 0
    ot_3mm: int = 0
    ot_4mm: int = 0

    # Exonic off-targets (critical for safety)
    ot_exonic_0mm: int = 0
    ot_exonic_1mm: int = 0
    ot_exonic_2mm: int = 0
    ot_exonic_3mm: int = 0

    # Total off-targets
    total_offtargets: int = 0

    # Raw row for debugging
    raw: Dict[str, Any] = field(default_factory=dict)

    @property
    def total_close_offtargets(self) -> int:
        """Off-targets with 0-2 mismatches (most dangerous)."""
        return self.ot_0mm + self.ot_1mm + self.ot_2mm

    @property
    def has_exonic_close_offtargets(self) -> bool:
        """Whether there are exonic off-targets with 0-2 mismatches."""
        return (self.ot_exonic_0mm + self.ot_exonic_1mm + self.ot_exonic_2mm) > 0


@dataclass
class CrisporOffTarget:
    """A single off-target site from CRISPOR."""

    guide_id: str
    chrom: str
    start: int
    end: int
    strand: str
    sequence: str
    mismatches: int
    mismatch_positions: List[int]
    cfd_score: float
    gene: Optional[str] = None
    gene_region: Optional[str] = None  # exon, intron, intergenic


def _get_float(row: Dict, *keys: str) -> Optional[float]:
    """Get float value from row, trying multiple column names."""
    for k in keys:
        if k in row and row[k] not in ("", None, "NA", "na", "-"):
            try:
                return float(row[k])
            except (ValueError, TypeError):
                continue
    return None


def _get_int(row: Dict, *keys: str) -> int:
    """Get int value from row, trying multiple column names."""
    for k in keys:
        if k in row and row[k] not in ("", None, "NA", "na", "-"):
            try:
                return int(float(row[k]))
            except (ValueError, TypeError):
                continue
    return 0


def _get_str(row: Dict, *keys: str) -> str:
    """Get string value from row, trying multiple column names."""
    for k in keys:
        if k in row and row[k] not in (None,):
            return str(row[k])
    return ""


def parse_guides_tsv(path: Path) -> List[CrisporGuideRow]:
    """
    Parse CRISPOR guides.tsv output file.

    CRISPOR column names vary between versions, so we try multiple names.
    """
    if not path.exists():
        return []

    rows: List[CrisporGuideRow] = []

    with path.open(newline="", encoding="utf-8") as f:
        # Skip comment lines
        lines = [l for l in f if not l.startswith("#")]

    if not lines:
        return []

    # Parse as TSV
    reader = csv.DictReader(lines, delimiter="\t")

    for r in reader:
        # Guide identification
        guide_id = _get_str(r, "guideId", "id", "name", "Guide")
        sequence = _get_str(r, "guideSeq", "guide", "seq", "targetSeq", "Sequence")
        pam = _get_str(r, "pam", "PAM")
        strand = _get_str(r, "strand", "Strand", "guideStrand")

        # Location
        chrom = _get_str(r, "chrom", "chr", "chromosome") or None
        start = _get_int(r, "start", "pos", "position", "chromStart") or None
        end = _get_int(r, "end", "chromEnd") or None

        # On-target scores
        doench = _get_float(r, "doenchScore", "Doench2016", "doench2016", "doench",
                           "Doench '16", "fusi", "Fusi")
        moreno = _get_float(r, "morenoScore", "Moreno-Mateos", "moreno", "crisprScan")
        oof = _get_float(r, "outOfFrame", "OOF", "oof")

        # Specificity scores
        mit = _get_float(r, "mitSpecScore", "MIT", "mitSpec", "specificity",
                        "MIT Spec.", "mitSpecificity")
        cfd = _get_float(r, "cfdSpecScore", "CFD", "cfdSpec", "CFD Spec.")

        # Off-target counts
        ot0 = _get_int(r, "numOfftarget0Mm", "ot0", "OT0", "0mm")
        ot1 = _get_int(r, "numOfftarget1Mm", "ot1", "OT1", "1mm")
        ot2 = _get_int(r, "numOfftarget2Mm", "ot2", "OT2", "2mm")
        ot3 = _get_int(r, "numOfftarget3Mm", "ot3", "OT3", "3mm")
        ot4 = _get_int(r, "numOfftarget4Mm", "ot4", "OT4", "4mm")

        # Parse the summary string if available (format: "0 - 1 - 4 - 27 - 203")
        ot_summary = _get_str(r, "offtargetCount", "offTargets", "OT summary")
        if ot_summary and " - " in ot_summary:
            parts = ot_summary.split(" - ")
            if len(parts) >= 5:
                try:
                    ot0 = int(parts[0])
                    ot1 = int(parts[1])
                    ot2 = int(parts[2])
                    ot3 = int(parts[3])
                    ot4 = int(parts[4])
                except ValueError:
                    pass

        total_ot = ot0 + ot1 + ot2 + ot3 + ot4

        rows.append(CrisporGuideRow(
            guide_id=guide_id,
            sequence=sequence.upper().replace(" ", ""),
            pam=pam.upper(),
            strand=strand,
            chrom=chrom,
            start=start,
            end=end,
            doench_2016=doench,
            moreno_mateos=moreno,
            out_of_frame=oof,
            mit_specificity=mit,
            cfd_specificity=cfd,
            ot_0mm=ot0,
            ot_1mm=ot1,
            ot_2mm=ot2,
            ot_3mm=ot3,
            ot_4mm=ot4,
            total_offtargets=total_ot,
            raw=dict(r),
        ))

    return rows


def parse_offtargets_tsv(path: Path) -> List[CrisporOffTarget]:
    """
    Parse CRISPOR offtargets.tsv output file.

    Returns detailed information about each off-target site.
    """
    if not path or not path.exists():
        return []

    offtargets: List[CrisporOffTarget] = []

    with path.open(newline="", encoding="utf-8") as f:
        lines = [l for l in f if not l.startswith("#")]

    if not lines:
        return []

    reader = csv.DictReader(lines, delimiter="\t")

    for r in reader:
        guide_id = _get_str(r, "guideId", "guide", "name")
        chrom = _get_str(r, "chrom", "chr")
        start = _get_int(r, "start", "pos")
        end = _get_int(r, "end")
        strand = _get_str(r, "strand")
        seq = _get_str(r, "seq", "sequence", "otSeq")
        mm = _get_int(r, "mismatchCount", "mm", "mismatches")
        cfd = _get_float(r, "cfdScore", "CFD") or 0.0

        # Mismatch positions (if available)
        mm_pos_str = _get_str(r, "mismatchPos", "mmPos")
        mm_positions = []
        if mm_pos_str:
            try:
                mm_positions = [int(x) for x in mm_pos_str.split(",") if x]
            except ValueError:
                pass

        # Gene information
        gene = _get_str(r, "gene", "geneName") or None
        region = _get_str(r, "region", "geneRegion", "annotation") or None

        offtargets.append(CrisporOffTarget(
            guide_id=guide_id,
            chrom=chrom,
            start=start,
            end=end or start + 23,
            strand=strand,
            sequence=seq.upper(),
            mismatches=mm,
            mismatch_positions=mm_positions,
            cfd_score=cfd,
            gene=gene,
            gene_region=region,
        ))

    return offtargets


def index_guides_by_sequence(guides: List[CrisporGuideRow]) -> Dict[str, CrisporGuideRow]:
    """Create index of guides by sequence for fast lookup."""
    return {g.sequence.upper(): g for g in guides if g.sequence}


def index_guides_by_seq_pam(guides: List[CrisporGuideRow]) -> Dict[tuple, CrisporGuideRow]:
    """Create index of guides by (sequence, PAM) tuple."""
    return {(g.sequence.upper(), g.pam.upper()): g for g in guides if g.sequence}
