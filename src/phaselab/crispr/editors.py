"""
PhaseLab CRISPR Editors: Base editing and prime editing support.

Extends the CRISPR pipeline to include:
- Adenine base editors (ABE)
- Cytosine base editors (CBE)
- Prime editing (PE)

Version: 0.3.0
"""

import re
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass

from .pam_scan import PAMHit, find_pam_sites, reverse_complement


# Editing windows for base editors (positions from PAM-distal end)
EDITING_WINDOWS = {
    "ABE8e": (4, 8),      # A→G editing window
    "ABE7.10": (4, 7),    # Earlier ABE version
    "BE4": (4, 8),        # C→T editing
    "BE3": (4, 8),        # C→T editing
    "evoAPOBEC1-BE4": (4, 8),
    "Target-AID": (2, 4),  # Narrower window
}


@dataclass
class BaseEditSite:
    """A potential base editing site."""
    guide: str
    pam: str
    position: int
    strand: str

    # Editing info
    edit_type: str          # "A>G" or "C>T"
    edit_positions: List[int]  # Positions of editable bases in guide
    editing_window: Tuple[int, int]

    # Context
    target_codon: Optional[str] = None
    amino_acid_change: Optional[str] = None

    def to_dict(self) -> Dict:
        return {
            "guide": self.guide,
            "pam": self.pam,
            "position": self.position,
            "strand": self.strand,
            "edit_type": self.edit_type,
            "edit_positions": self.edit_positions,
            "editing_window": self.editing_window,
            "target_codon": self.target_codon,
            "amino_acid_change": self.amino_acid_change,
        }


@dataclass
class PrimeEditSite:
    """A potential prime editing site."""
    guide: str              # spacer sequence
    pam: str
    position: int
    strand: str

    # Prime editing template
    pbs_length: int         # Primer binding site length
    rt_template: str        # Reverse transcription template
    edit_distance: int      # Distance from nick to edit

    # Edit details
    edit_type: str          # "substitution", "insertion", "deletion"
    edit_sequence: str      # What's being changed/inserted

    def to_dict(self) -> Dict:
        return {
            "guide": self.guide,
            "pam": self.pam,
            "position": self.position,
            "strand": self.strand,
            "pbs_length": self.pbs_length,
            "rt_template": self.rt_template,
            "edit_distance": self.edit_distance,
            "edit_type": self.edit_type,
            "edit_sequence": self.edit_sequence,
        }


def find_abe_sites(
    sequence: str,
    editor: str = "ABE8e",
    pam: str = "NGG",
) -> List[BaseEditSite]:
    """
    Find adenine base editing sites (A→G conversion).

    Args:
        sequence: DNA sequence to scan
        editor: ABE variant (determines editing window)
        pam: PAM sequence

    Returns:
        List of BaseEditSite objects
    """
    # Get PAM hits
    pam_hits = find_pam_sites(sequence, pam=pam)

    # Get editing window
    window = EDITING_WINDOWS.get(editor, (4, 8))

    sites = []
    for hit in pam_hits:
        guide = hit.guide

        # Find A's in editing window
        # Window positions are from PAM-distal (5') end of guide
        editable_positions = []
        for i in range(window[0] - 1, window[1]):  # 0-indexed
            if i < len(guide) and guide[i] == 'A':
                editable_positions.append(i + 1)  # 1-indexed for output

        if editable_positions:
            sites.append(BaseEditSite(
                guide=guide,
                pam=hit.pam,
                position=hit.position,
                strand=hit.strand,
                edit_type="A>G",
                edit_positions=editable_positions,
                editing_window=window,
            ))

    return sites


def find_cbe_sites(
    sequence: str,
    editor: str = "BE4",
    pam: str = "NGG",
) -> List[BaseEditSite]:
    """
    Find cytosine base editing sites (C→T conversion).

    Args:
        sequence: DNA sequence to scan
        editor: CBE variant (determines editing window)
        pam: PAM sequence

    Returns:
        List of BaseEditSite objects
    """
    pam_hits = find_pam_sites(sequence, pam=pam)
    window = EDITING_WINDOWS.get(editor, (4, 8))

    sites = []
    for hit in pam_hits:
        guide = hit.guide

        editable_positions = []
        for i in range(window[0] - 1, window[1]):
            if i < len(guide) and guide[i] == 'C':
                editable_positions.append(i + 1)

        if editable_positions:
            sites.append(BaseEditSite(
                guide=guide,
                pam=hit.pam,
                position=hit.position,
                strand=hit.strand,
                edit_type="C>T",
                edit_positions=editable_positions,
                editing_window=window,
            ))

    return sites


def design_prime_edit(
    sequence: str,
    edit_position: int,
    edit_type: str,
    edit_sequence: str,
    pbs_length: int = 13,
    rt_extension: int = 10,
    pam: str = "NGG",
) -> List[PrimeEditSite]:
    """
    Design prime editing guides for a specific edit.

    Prime editing uses a pegRNA (prime editing guide RNA) with:
    - Spacer (20bp guide)
    - PBS (primer binding site)
    - RT template (contains the desired edit)

    Args:
        sequence: Target DNA sequence
        edit_position: Position of desired edit (0-based)
        edit_type: "substitution", "insertion", or "deletion"
        edit_sequence: New sequence to introduce
        pbs_length: Primer binding site length (typically 10-17)
        rt_extension: RT template extension past edit
        pam: PAM sequence

    Returns:
        List of PrimeEditSite candidates
    """
    # Find nearby PAM sites (prime editing nicks 3bp upstream of PAM)
    pam_hits = find_pam_sites(sequence, pam=pam)

    candidates = []

    for hit in pam_hits:
        # Nick position is 3bp upstream of PAM
        nick_pos = hit.position - 3

        # Calculate edit distance from nick
        edit_distance = edit_position - nick_pos

        # Prime editing works best for edits within ~30bp of nick
        if abs(edit_distance) > 40:
            continue

        # For + strand, RT template is on the edited strand 3' of nick
        if hit.strand == "+":
            # PBS binds upstream of nick
            pbs_start = nick_pos - pbs_length
            if pbs_start < 0:
                continue

            pbs_seq = reverse_complement(sequence[pbs_start:nick_pos])

            # RT template starts at nick and includes edit
            rt_end = edit_position + len(edit_sequence) + rt_extension
            if rt_end > len(sequence):
                rt_end = len(sequence)

            # Create RT template with edit
            if edit_type == "substitution":
                rt_template = (
                    sequence[nick_pos:edit_position] +
                    edit_sequence +
                    sequence[edit_position + len(edit_sequence):rt_end]
                )
            elif edit_type == "insertion":
                rt_template = (
                    sequence[nick_pos:edit_position] +
                    edit_sequence +
                    sequence[edit_position:rt_end]
                )
            elif edit_type == "deletion":
                # edit_sequence is what to delete
                rt_template = (
                    sequence[nick_pos:edit_position] +
                    sequence[edit_position + len(edit_sequence):rt_end]
                )
            else:
                continue

            candidates.append(PrimeEditSite(
                guide=hit.guide,
                pam=hit.pam,
                position=hit.position,
                strand=hit.strand,
                pbs_length=pbs_length,
                rt_template=rt_template,
                edit_distance=edit_distance,
                edit_type=edit_type,
                edit_sequence=edit_sequence,
            ))

    return candidates


def score_base_edit_site(
    site: BaseEditSite,
    avoid_bystanders: bool = True,
) -> float:
    """
    Score a base editing site for specificity.

    Considers:
    - Number of editable bases (fewer is more specific)
    - Position of editable base in window
    - Bystander edits

    Args:
        site: BaseEditSite to score
        avoid_bystanders: Penalize multiple editable bases

    Returns:
        Score (0-1, higher is better)
    """
    score = 1.0

    # Penalize multiple editable bases (bystander editing)
    n_editable = len(site.edit_positions)
    if avoid_bystanders and n_editable > 1:
        score *= (1.0 / n_editable)

    # Prefer edits in center of window
    window_center = (site.editing_window[0] + site.editing_window[1]) / 2
    for pos in site.edit_positions:
        distance_from_center = abs(pos - window_center)
        window_size = site.editing_window[1] - site.editing_window[0]
        center_score = 1.0 - (distance_from_center / window_size)
        score *= center_score

    return float(score)


def score_prime_edit_site(site: PrimeEditSite) -> float:
    """
    Score a prime editing site.

    Considers:
    - Edit distance from nick (closer is better)
    - RT template length
    - PBS properties

    Args:
        site: PrimeEditSite to score

    Returns:
        Score (0-1, higher is better)
    """
    score = 1.0

    # Penalize long edit distances
    if abs(site.edit_distance) > 15:
        score *= 0.5
    elif abs(site.edit_distance) > 30:
        score *= 0.2

    # Penalize very long RT templates
    rt_len = len(site.rt_template)
    if rt_len > 30:
        score *= (30 / rt_len)

    # PBS GC content (prefer 40-60%)
    pbs_gc = sum(1 for b in site.rt_template[:site.pbs_length] if b in 'GC') / site.pbs_length
    if 0.4 <= pbs_gc <= 0.6:
        score *= 1.0
    else:
        score *= 0.8

    return float(score)


# Convenience functions


def find_correction_sites(
    sequence: str,
    mutation_pos: int,
    mutation_type: str,  # "A>G", "G>A", "C>T", "T>C"
) -> Dict[str, List]:
    """
    Find editing sites to correct a known mutation.

    Args:
        sequence: DNA sequence containing mutation
        mutation_pos: Position of mutation (0-based)
        mutation_type: Type of mutation to correct

    Returns:
        Dictionary with ABE, CBE, and PE candidates
    """
    results = {"ABE": [], "CBE": [], "PE": []}

    # Determine which editors can correct this mutation
    if mutation_type in ["G>A", "g>a"]:
        # Need A→G to restore, use ABE
        abe_sites = find_abe_sites(sequence)
        for site in abe_sites:
            # Check if mutation position is in editing window
            if mutation_pos in range(site.position - 20, site.position):
                results["ABE"].append(site)

    elif mutation_type in ["A>G", "a>g"]:
        # Need G→A equivalent, use CBE on opposite strand
        # (C on sense = G on antisense)
        cbe_sites = find_cbe_sites(sequence)
        for site in cbe_sites:
            if site.strand == "-":
                results["CBE"].append(site)

    elif mutation_type in ["T>C", "t>c"]:
        # Need C→T to restore
        cbe_sites = find_cbe_sites(sequence)
        results["CBE"] = cbe_sites

    elif mutation_type in ["C>T", "c>t"]:
        # Need T→C, use ABE on opposite strand
        abe_sites = find_abe_sites(sequence)
        for site in abe_sites:
            if site.strand == "-":
                results["ABE"].append(site)

    # Prime editing can correct any mutation
    original_base = mutation_type[0]
    pe_sites = design_prime_edit(
        sequence=sequence,
        edit_position=mutation_pos,
        edit_type="substitution",
        edit_sequence=original_base,
    )
    results["PE"] = pe_sites

    return results
