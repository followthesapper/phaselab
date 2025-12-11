"""
PhaseLab PAM Scanner: Find PAM sites and extract guide sequences.

Supports multiple Cas systems:
- SpCas9: NGG (default)
- SaCas9: NNGRRT
- Cas12a: TTTV
"""

import re
from typing import List, Dict, Optional
from dataclasses import dataclass


# PAM patterns for different Cas systems
PAM_PATTERNS = {
    "NGG": r"[ACGT]GG",           # SpCas9
    "NAG": r"[ACGT]AG",           # SpCas9 (weaker)
    "NNGRRT": r"[ACGT][ACGT]G[AG][AG]T",  # SaCas9
    "TTTV": r"TTT[ACG]",          # Cas12a (upstream PAM)
    "TTTN": r"TTT[ACGT]",         # Cas12a variant
}


@dataclass
class PAMHit:
    """A PAM site hit with associated guide sequence."""
    position: int          # Position of PAM in sequence (0-based)
    strand: str            # "+" (forward) or "-" (reverse)
    guide: str             # 20bp protospacer sequence
    pam: str               # PAM sequence
    guide_start: int       # Start position of guide
    guide_end: int         # End position of guide


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(b, 'N') for b in reversed(seq))


def find_pam_sites(
    sequence: str,
    pam: str = "NGG",
    guide_length: int = 20,
    both_strands: bool = True,
) -> List[PAMHit]:
    """
    Find all PAM sites in a sequence and extract guide sequences.

    For SpCas9 (NGG): Guide is 20bp upstream of PAM.
    For Cas12a (TTTV): Guide is downstream of PAM.

    Args:
        sequence: DNA sequence to scan (5'->3').
        pam: PAM pattern name (e.g., "NGG", "NNGRRT") or regex.
        guide_length: Length of guide/protospacer (default 20).
        both_strands: Scan both forward and reverse strands.

    Returns:
        List of PAMHit objects.
    """
    sequence = sequence.upper()
    hits = []

    # Get PAM regex pattern
    if pam in PAM_PATTERNS:
        pam_regex = PAM_PATTERNS[pam]
    else:
        pam_regex = pam  # Assume user provided regex

    # Determine if PAM is upstream (Cas12a) or downstream (Cas9)
    pam_upstream = pam.startswith("TTT")  # Cas12a

    # Forward strand
    for match in re.finditer(pam_regex, sequence):
        pam_start = match.start()
        pam_end = match.end()
        pam_seq = match.group()

        if pam_upstream:
            # Cas12a: guide is downstream of PAM
            guide_start = pam_end
            guide_end = guide_start + guide_length
        else:
            # Cas9: guide is upstream of PAM
            guide_start = pam_start - guide_length
            guide_end = pam_start

        # Check bounds
        if guide_start < 0 or guide_end > len(sequence):
            continue

        guide_seq = sequence[guide_start:guide_end]

        hits.append(PAMHit(
            position=pam_start,
            strand="+",
            guide=guide_seq,
            pam=pam_seq,
            guide_start=guide_start,
            guide_end=guide_end,
        ))

    # Reverse strand
    if both_strands:
        rev_seq = reverse_complement(sequence)
        for match in re.finditer(pam_regex, rev_seq):
            pam_start_rev = match.start()
            pam_seq = match.group()

            if pam_upstream:
                guide_start_rev = pam_start_rev + len(pam_seq)
                guide_end_rev = guide_start_rev + guide_length
            else:
                guide_start_rev = pam_start_rev - guide_length
                guide_end_rev = pam_start_rev

            if guide_start_rev < 0 or guide_end_rev > len(rev_seq):
                continue

            guide_seq = rev_seq[guide_start_rev:guide_end_rev]

            # Convert positions back to forward strand coordinates
            # Position on forward strand = len(seq) - 1 - position on reverse
            fwd_pam_pos = len(sequence) - match.end()

            hits.append(PAMHit(
                position=fwd_pam_pos,
                strand="-",
                guide=guide_seq,
                pam=pam_seq,
                guide_start=len(sequence) - guide_end_rev,
                guide_end=len(sequence) - guide_start_rev,
            ))

    return hits


def filter_by_window(
    hits: List[PAMHit],
    tss_position: int,
    window: tuple = (-400, -50),
) -> List[PAMHit]:
    """
    Filter PAM hits to those within a CRISPRa activation window.

    Args:
        hits: List of PAMHit objects.
        tss_position: Position of transcription start site (0-based).
        window: (start, end) relative to TSS (e.g., (-400, -50) for CRISPRa).

    Returns:
        Filtered list of PAMHit objects.
    """
    filtered = []
    for hit in hits:
        # Use guide midpoint for position calculation
        guide_mid = (hit.guide_start + hit.guide_end) // 2
        rel_pos = guide_mid - tss_position

        if window[0] <= rel_pos <= window[1]:
            filtered.append(hit)

    return filtered
