"""
PhaseLab Base Editing: ABE and CBE guide design for single-nucleotide changes.

Implements:
- Cytosine Base Editor (CBE) design - C→T conversions
- Adenine Base Editor (ABE) design - A→G conversions
- Activity window optimization (positions 4-8 for most editors)
- Bystander editing prediction
- Deaminase preference modeling
- IR coherence validation
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass

from .pam_scan import find_pam_sites, PAMHit
from .scoring import (
    gc_content,
    delta_g_santalucia,
    mit_specificity_score,
    cfd_score,
    max_homopolymer_run,
    sequence_complexity,
)
from ..core.coherence import coherence_score, go_no_go
from ..core.hamiltonians import build_grna_hamiltonian


# Base editor activity windows (position from PAM-distal end, 1-indexed)
# Based on Komor et al. 2016 (CBE) and Gaudelli et al. 2017 (ABE)
ACTIVITY_WINDOWS = {
    'BE4': (4, 8),      # Standard CBE, positions 4-8
    'BE4-GAM': (4, 8),  # BE4 with GAM for reduced indels
    'ABE7.10': (4, 7),  # Original ABE
    'ABE8e': (4, 8),    # Enhanced ABE, broader window
    'ABE8.20': (3, 9),  # Wider activity window
    'CBE4max': (4, 8),  # Optimized CBE
    'evoBE4max': (4, 8),
    'APOBEC': (4, 8),   # Generic APOBEC-based
    'TadA': (4, 7),     # TadA-based ABE
}

# Position-specific editing efficiency (relative to peak)
POSITION_EFFICIENCY = {
    1: 0.05, 2: 0.10, 3: 0.30,
    4: 0.85, 5: 1.00, 6: 0.95, 7: 0.80, 8: 0.60,
    9: 0.25, 10: 0.10, 11: 0.05, 12: 0.02,
}

# Sequence context preferences for CBE (APOBEC)
# Format: 5'-N[target C]N-3'
CBE_CONTEXT_PREFERENCE = {
    'TC': 1.2,   # Preferred
    'CC': 1.0,   # Good
    'AC': 0.8,   # Moderate
    'GC': 0.6,   # Less preferred
}

# Sequence context preferences for ABE (TadA)
ABE_CONTEXT_PREFERENCE = {
    'TA': 1.1,   # Slightly preferred
    'CA': 1.0,   # Good
    'AA': 0.9,   # Good
    'GA': 0.7,   # Less preferred
}


@dataclass
class BaseEditConfig:
    """Configuration for base editing guide design."""

    # Editor type
    editor: str = "ABE8e"  # Options: BE4, ABE7.10, ABE8e, etc.

    # PAM settings
    pam: str = "NGG"
    guide_length: int = 20

    # Activity window (will be set based on editor if not specified)
    activity_window: Optional[Tuple[int, int]] = None

    # Target base
    target_base: str = "A"  # "C" for CBE, "A" for ABE

    # Filtering
    min_gc: float = 0.35
    max_gc: float = 0.70
    max_homopolymer: int = 4

    # Bystander editing
    check_bystanders: bool = True
    max_bystanders_in_window: int = 2  # Allow up to 2 other editable bases

    # Scoring weights
    weight_mit: float = 1.0
    weight_cfd: float = 1.0
    weight_position: float = 1.5  # Position in activity window
    weight_context: float = 1.0   # Sequence context
    weight_bystander: float = -0.5  # Negative weight for bystanders
    weight_coherence: float = 1.0

    # Coherence
    compute_coherence: bool = True
    coherence_shots: int = 2000

    # Output
    top_n: int = 10


def get_activity_window(editor: str) -> Tuple[int, int]:
    """Get activity window for a base editor."""
    return ACTIVITY_WINDOWS.get(editor, (4, 8))


def get_target_base(editor: str) -> str:
    """Get target base for an editor (C for CBE, A for ABE)."""
    if 'ABE' in editor or 'TadA' in editor:
        return 'A'
    else:  # CBE, BE4, etc.
        return 'C'


def editing_efficiency_at_position(
    position: int,
    editor: str = "ABE8e",
) -> float:
    """
    Get relative editing efficiency at a position.

    Args:
        position: Position in guide (1-indexed from PAM-distal).
        editor: Base editor type.

    Returns:
        Relative efficiency (0.0 to 1.0).
    """
    window = get_activity_window(editor)

    if window[0] <= position <= window[1]:
        # Inside activity window - use position-specific efficiency
        return POSITION_EFFICIENCY.get(position, 0.5)
    else:
        # Outside window - very low efficiency
        return POSITION_EFFICIENCY.get(position, 0.02)


def sequence_context_score(
    guide_seq: str,
    target_position: int,
    editor: str = "ABE8e",
) -> float:
    """
    Score the sequence context around target base.

    Args:
        guide_seq: Full guide sequence.
        target_position: Position of target base (1-indexed).
        editor: Base editor type.

    Returns:
        Context score (0.5 to 1.2).
    """
    guide_seq = guide_seq.upper()
    idx = target_position - 1  # 0-indexed

    if idx <= 0 or idx >= len(guide_seq):
        return 1.0

    # Get 5' context (base before target)
    context = guide_seq[idx - 1] + guide_seq[idx]

    if get_target_base(editor) == 'C':
        return CBE_CONTEXT_PREFERENCE.get(context, 1.0)
    else:
        return ABE_CONTEXT_PREFERENCE.get(context, 1.0)


def find_bystanders(
    guide_seq: str,
    target_position: int,
    editor: str = "ABE8e",
) -> List[Dict[str, Any]]:
    """
    Find bystander editable bases in the activity window.

    Args:
        guide_seq: Guide sequence.
        target_position: Position of intended edit.
        editor: Base editor type.

    Returns:
        List of bystander positions and predicted efficiencies.
    """
    guide_seq = guide_seq.upper()
    window = get_activity_window(editor)
    target_base = get_target_base(editor)

    bystanders = []

    for pos in range(window[0], window[1] + 1):
        if pos == target_position:
            continue  # Skip intended target

        idx = pos - 1
        if 0 <= idx < len(guide_seq) and guide_seq[idx] == target_base:
            eff = editing_efficiency_at_position(pos, editor)
            context = sequence_context_score(guide_seq, pos, editor)

            bystanders.append({
                'position': pos,
                'base': target_base,
                'efficiency': eff * context,
                'context': guide_seq[max(0, idx-1):idx+2],
            })

    return bystanders


def design_base_edit_guides(
    sequence: str,
    target_position: int,
    target_base: Optional[str] = None,
    config: Optional[BaseEditConfig] = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Design guides for base editing at a specific position.

    Finds guides that place the target base within the editor's
    activity window and ranks by predicted efficiency.

    Args:
        sequence: DNA sequence containing target site.
        target_position: Position of base to edit (0-based in sequence).
        target_base: Target base ("A" for ABE, "C" for CBE). Auto-detected from editor if not specified.
        config: BaseEditConfig with parameters.
        verbose: Print progress messages.

    Returns:
        DataFrame with ranked base editing guide candidates.

    Example:
        >>> from phaselab.crispr import design_base_edit_guides
        >>> guides = design_base_edit_guides(
        ...     sequence=gene_region,
        ...     target_position=100,  # Position of A to edit
        ...     config=BaseEditConfig(editor="ABE8e"),
        ... )
        >>> print(guides[['sequence', 'target_in_window_pos', 'efficiency']].head())
    """
    if config is None:
        config = BaseEditConfig()

    sequence = sequence.upper()

    # Auto-detect target base from editor if not specified
    if target_base is None:
        target_base = config.target_base or get_target_base(config.editor)

    activity_window = config.activity_window or get_activity_window(config.editor)

    if verbose:
        print(f"Designing {config.editor} guides for {target_base} editing at position {target_position}")
        print(f"Activity window: positions {activity_window[0]}-{activity_window[1]}")

    # Verify target base
    if sequence[target_position] != target_base:
        if verbose:
            print(f"Warning: Position {target_position} is {sequence[target_position]}, not {target_base}")

    # Find PAM sites that would position target in activity window
    # For each window position, calculate where PAM should be
    search_start = max(0, target_position - 25)
    search_end = min(len(sequence), target_position + 25)
    search_seq = sequence[search_start:search_end]

    all_hits = find_pam_sites(
        search_seq,
        pam=config.pam,
        guide_length=config.guide_length,
        both_strands=True,
    )

    if verbose:
        print(f"Found {len(all_hits)} PAM sites in search region")

    if not all_hits:
        return _empty_base_edit_df()

    candidates = []

    for hit in all_hits:
        guide_seq = hit.guide

        # Calculate where target base falls in this guide
        # Adjust position back to original coordinates
        abs_guide_start = search_start + hit.guide_start

        if hit.strand == '+':
            target_in_guide = target_position - abs_guide_start + 1  # 1-indexed
        else:
            # For minus strand, positions are reversed
            target_in_guide = (abs_guide_start + config.guide_length) - target_position

        # Check if target is in activity window
        if not (activity_window[0] <= target_in_guide <= activity_window[1]):
            continue

        # Quality filters
        gc = gc_content(guide_seq)
        if gc < config.min_gc or gc > config.max_gc:
            continue

        homo = max_homopolymer_run(guide_seq)
        if homo > config.max_homopolymer:
            continue

        # Editing efficiency
        position_eff = editing_efficiency_at_position(target_in_guide, config.editor)
        context_score = sequence_context_score(guide_seq, target_in_guide, config.editor)

        # Find bystanders
        bystanders = find_bystanders(guide_seq, target_in_guide, config.editor)
        n_bystanders = len(bystanders)

        if config.check_bystanders and n_bystanders > config.max_bystanders_in_window:
            continue

        # Standard scores
        delta_g = delta_g_santalucia(guide_seq)
        mit = mit_specificity_score(guide_seq)
        cfd = cfd_score(guide_seq)

        # IR coherence
        R_bar = None
        go_status = None
        if config.compute_coherence:
            R_bar = _compute_guide_coherence(guide_seq)
            go_status = go_no_go(R_bar)

        # Combined score
        combined = _compute_base_edit_score(
            position_efficiency=position_eff,
            context_score=context_score,
            n_bystanders=n_bystanders,
            mit=mit,
            cfd=cfd,
            R_bar=R_bar,
            config=config,
        )

        candidates.append({
            'sequence': guide_seq,
            'pam': hit.pam,
            'strand': hit.strand,
            'target_in_window_pos': target_in_guide,
            'position_efficiency': round(position_eff, 3),
            'context_score': round(context_score, 3),
            'combined_efficiency': round(position_eff * context_score, 3),
            'n_bystanders': n_bystanders,
            'bystander_positions': [b['position'] for b in bystanders],
            'gc': round(gc, 3),
            'delta_g': round(delta_g, 3),
            'mit_score': round(mit, 1),
            'cfd_score': round(cfd, 1),
            'coherence_R': round(R_bar, 4) if R_bar else None,
            'go_no_go': go_status,
            'combined_score': round(combined, 3),
        })

    if not candidates:
        return _empty_base_edit_df()

    df = pd.DataFrame(candidates)
    df.sort_values(by='combined_score', ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)

    if verbose:
        print(f"Returning {min(len(df), config.top_n)} top candidates")

    return df.head(config.top_n)


def _empty_base_edit_df() -> pd.DataFrame:
    """Return empty DataFrame with base editing columns."""
    return pd.DataFrame(columns=[
        'sequence', 'pam', 'strand', 'target_in_window_pos',
        'position_efficiency', 'context_score', 'combined_efficiency',
        'n_bystanders', 'bystander_positions',
        'gc', 'delta_g', 'mit_score', 'cfd_score',
        'coherence_R', 'go_no_go', 'combined_score',
    ])


def _compute_guide_coherence(guide_seq: str) -> float:
    """Compute IR coherence for guide sequence."""
    try:
        H = build_grna_hamiltonian(guide_seq)
        terms = H.get_terms()
        if not terms:
            return 0.5

        energies = [abs(coeff) for coeff, _ in terms]
        mean_energy = np.mean(energies)
        std_energy = np.std(energies)

        if mean_energy > 0:
            V_phi = std_energy / mean_energy
            R_bar = np.exp(-V_phi / 2)
        else:
            R_bar = 0.5

        return float(np.clip(R_bar, 0, 1))
    except Exception:
        return 0.5


def _compute_base_edit_score(
    position_efficiency: float,
    context_score: float,
    n_bystanders: int,
    mit: float,
    cfd: float,
    R_bar: Optional[float],
    config: BaseEditConfig,
) -> float:
    """Compute combined score for base editing."""
    score = 0.0

    # Editing efficiency (position and context)
    score += config.weight_position * position_efficiency
    score += config.weight_context * context_score

    # Bystander penalty
    score += config.weight_bystander * n_bystanders

    # Specificity
    score += config.weight_mit * (mit / 100.0)
    score += config.weight_cfd * (cfd / 100.0)

    # Coherence
    if R_bar is not None:
        score += config.weight_coherence * R_bar

    return score


def validate_base_edit(
    guide_seq: str,
    target_position: int,
    editor: str = "ABE8e",
) -> Dict[str, Any]:
    """
    Validate a guide for base editing.

    Args:
        guide_seq: Guide sequence.
        target_position: Position of target base in guide (1-indexed).
        editor: Base editor type.

    Returns:
        Validation results dictionary.
    """
    guide_seq = guide_seq.upper()
    window = get_activity_window(editor)
    target_base = get_target_base(editor)

    warnings = []

    # Check position in window
    if not (window[0] <= target_position <= window[1]):
        warnings.append(f"Position {target_position} outside activity window ({window[0]}-{window[1]})")

    # Check target base
    idx = target_position - 1
    if 0 <= idx < len(guide_seq):
        actual_base = guide_seq[idx]
        if actual_base != target_base:
            warnings.append(f"Position {target_position} is {actual_base}, expected {target_base}")

    # Check GC
    gc = gc_content(guide_seq)
    if gc < 0.35 or gc > 0.70:
        warnings.append(f"GC content {gc:.0%} outside optimal range")

    # Check bystanders
    bystanders = find_bystanders(guide_seq, target_position, editor)
    if len(bystanders) > 2:
        warnings.append(f"Found {len(bystanders)} bystander editable bases")

    # Efficiency scores
    position_eff = editing_efficiency_at_position(target_position, editor)
    context = sequence_context_score(guide_seq, target_position, editor)

    R_bar = _compute_guide_coherence(guide_seq)

    return {
        'sequence': guide_seq,
        'editor': editor,
        'target_position': target_position,
        'target_base': target_base,
        'in_activity_window': window[0] <= target_position <= window[1],
        'position_efficiency': position_eff,
        'context_score': context,
        'combined_efficiency': position_eff * context,
        'bystanders': bystanders,
        'n_bystanders': len(bystanders),
        'gc': gc,
        'coherence_R': R_bar,
        'go_no_go': go_no_go(R_bar),
        'warnings': warnings,
        'valid': len(warnings) == 0,
    }


# Convenience functions for specific editors

def design_abe_guides(sequence: str, target_position: int, **kwargs) -> pd.DataFrame:
    """Design ABE (A→G) guides. Shortcut for design_base_edit_guides with ABE8e."""
    config = BaseEditConfig(editor="ABE8e", target_base="A")
    return design_base_edit_guides(sequence, target_position, config=config, **kwargs)


def design_cbe_guides(sequence: str, target_position: int, **kwargs) -> pd.DataFrame:
    """Design CBE (C→T) guides. Shortcut for design_base_edit_guides with BE4."""
    config = BaseEditConfig(editor="BE4", target_base="C")
    return design_base_edit_guides(sequence, target_position, config=config, **kwargs)
