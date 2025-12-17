"""
PhaseLab SURF: Visualization utilities for SURF + coherence analysis.

Provides publication-quality figures for:
- SURF beta profiles with coherence overlay
- Raw vs deconvolved comparison
- Region classification maps
"""

import numpy as np
from typing import Optional, List, Dict, Any, Tuple

from .parser import SURFOutput
from .coherence import SURFCoherenceResult


def plot_surf_coherence(
    result: SURFCoherenceResult,
    figsize: Tuple[float, float] = (14, 10),
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    show_regions: bool = True,
    show_surf_peaks: bool = True,
) -> Any:
    """
    Plot SURF beta profile with coherence overlay.

    Creates a multi-panel figure showing:
    - Top: SURF beta profile with significant regions
    - Middle: Coherence profile
    - Bottom: Local variance with region classification

    Args:
        result: SURFCoherenceResult from compute_surf_coherence.
        figsize: Figure size (width, height).
        title: Optional figure title.
        save_path: Path to save figure (optional).
        show_regions: Whether to show region boundaries.
        show_surf_peaks: Whether to highlight SURF significant peaks.

    Returns:
        Matplotlib figure object.

    Example:
        >>> result = compute_surf_coherence(surf_output)
        >>> fig = plot_surf_coherence(result, save_path='surf_coherence.png')
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
    except ImportError:
        raise ImportError("matplotlib required for visualization")

    fig, axes = plt.subplots(3, 1, figsize=figsize, sharex=True)

    positions = result.surf.positions
    beta = result.surf.beta
    coherence = result.profile.coherence
    variance = result.profile.local_variance
    coh_positions = result.profile.coords

    # Panel 1: SURF beta profile
    ax1 = axes[0]
    ax1.plot(positions, beta, 'b-', linewidth=1, alpha=0.8, label='SURF beta')
    ax1.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax1.set_ylabel('Beta (regulatory signal)', fontsize=11)
    ax1.legend(loc='upper right')

    # Highlight SURF significant regions
    if show_surf_peaks:
        for region in result.surf.regions:
            color = 'green' if region.direction == 'activating' else 'red'
            ax1.axvspan(region.start, region.end, alpha=0.2, color=color)

    # Panel 2: Coherence profile
    ax2 = axes[1]
    ax2.plot(coh_positions, coherence, 'purple', linewidth=1.5, label='Spatial coherence')
    ax2.axhline(y=0.7, color='green', linestyle='--', alpha=0.7, label='Stable threshold')
    ax2.axhline(y=0.4, color='orange', linestyle='--', alpha=0.7, label='Mixed threshold')
    ax2.set_ylabel('Coherence', fontsize=11)
    ax2.set_ylim(0, 1)
    ax2.legend(loc='upper right')

    # Panel 3: Local variance with region classification
    ax3 = axes[2]
    ax3.plot(coh_positions, variance, 'red', linewidth=1, alpha=0.8, label='Local variance')
    ax3.set_ylabel('Local variance', fontsize=11)
    ax3.set_xlabel('Position relative to TSS (bp)', fontsize=11)
    ax3.legend(loc='upper right')

    # Show region classifications
    if show_regions:
        for region in result.regions:
            if region.stability.value == 'stable':
                color = 'green'
                alpha = 0.15
            elif region.stability.value == 'amplifying':
                color = 'red'
                alpha = 0.15
            elif region.stability.value == 'mixed':
                color = 'yellow'
                alpha = 0.1
            else:
                continue

            for ax in axes:
                ax.axvspan(region.start, region.end, alpha=alpha, color=color)

    # Title
    if title:
        fig.suptitle(title, fontsize=14, fontweight='bold')
    else:
        fig.suptitle(
            f"SURF + Coherence Analysis: {result.surf.gene_symbol}\n"
            f"(r = {result.correlation:.3f}, validated = {result.is_validated})",
            fontsize=14,
        )

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')

    return fig


def plot_raw_vs_deconvolved(
    raw_positions: np.ndarray,
    raw_responses: np.ndarray,
    surf: SURFOutput,
    figsize: Tuple[float, float] = (14, 8),
    title: Optional[str] = None,
    save_path: Optional[str] = None,
) -> Any:
    """
    Compare raw screen data with SURF-deconvolved signal.

    Shows how SURF deconvolution cleans up the regulatory signal.

    Args:
        raw_positions: Raw screen positions.
        raw_responses: Raw screen responses.
        surf: SURF output.
        figsize: Figure size.
        title: Optional title.
        save_path: Path to save figure.

    Returns:
        Matplotlib figure object.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise ImportError("matplotlib required for visualization")

    fig, axes = plt.subplots(2, 1, figsize=figsize, sharex=True)

    # Panel 1: Raw data
    ax1 = axes[0]
    ax1.scatter(raw_positions, raw_responses, s=10, alpha=0.5, c='gray', label='Raw sgRNA')
    ax1.axhline(y=0, color='black', linestyle='--', alpha=0.3)
    ax1.set_ylabel('Raw log fold change', fontsize=11)
    ax1.set_title('Raw Screen Data', fontsize=12)
    ax1.legend(loc='upper right')

    # Panel 2: SURF deconvolved
    ax2 = axes[1]
    ax2.plot(surf.positions, surf.beta, 'b-', linewidth=2, label='SURF beta')
    ax2.axhline(y=0, color='black', linestyle='--', alpha=0.3)
    ax2.set_ylabel('Deconvolved beta', fontsize=11)
    ax2.set_xlabel('Position relative to TSS (bp)', fontsize=11)
    ax2.set_title('SURF Deconvolved Signal', fontsize=12)
    ax2.legend(loc='upper right')

    # Highlight SURF significant regions
    for region in surf.regions:
        color = 'green' if region.direction == 'activating' else 'red'
        ax2.axvspan(region.start, region.end, alpha=0.2, color=color)

    if title:
        fig.suptitle(title, fontsize=14, fontweight='bold')
    else:
        fig.suptitle(f"Raw vs SURF Deconvolved: {surf.gene_symbol}", fontsize=14)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')

    return fig


def plot_surf_regions(
    result: SURFCoherenceResult,
    figsize: Tuple[float, float] = (12, 6),
    title: Optional[str] = None,
    save_path: Optional[str] = None,
) -> Any:
    """
    Plot region classification map.

    Shows which regions are stable, mixed, or amplifying,
    with SURF significant regions overlaid.

    Args:
        result: SURFCoherenceResult.
        figsize: Figure size.
        title: Optional title.
        save_path: Path to save figure.

    Returns:
        Matplotlib figure object.
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
        from matplotlib.collections import PatchCollection
    except ImportError:
        raise ImportError("matplotlib required for visualization")

    fig, ax = plt.subplots(figsize=figsize)

    # Get position range
    min_pos = result.surf.positions.min()
    max_pos = result.surf.positions.max()

    # Plot coherence regions as colored bars
    y_base = 0
    region_height = 0.4

    for region in result.regions:
        if region.stability.value == 'stable':
            color = 'green'
            label = 'Stable'
        elif region.stability.value == 'amplifying':
            color = 'red'
            label = 'Amplifying'
        elif region.stability.value == 'mixed':
            color = 'yellow'
            label = 'Mixed'
        else:
            color = 'lightgray'
            label = 'Irrelevant'

        rect = Rectangle(
            (region.start, y_base),
            region.end - region.start,
            region_height,
            facecolor=color,
            edgecolor='black',
            alpha=0.7,
            linewidth=0.5,
        )
        ax.add_patch(rect)

    # Plot SURF significant regions
    y_surf = 0.6
    for region in result.surf.regions:
        color = 'darkgreen' if region.direction == 'activating' else 'darkred'
        rect = Rectangle(
            (region.start, y_surf),
            region.end - region.start,
            region_height,
            facecolor=color,
            edgecolor='black',
            alpha=0.8,
            linewidth=0.5,
        )
        ax.add_patch(rect)

    # Formatting
    ax.set_xlim(min_pos - 50, max_pos + 50)
    ax.set_ylim(-0.1, 1.2)
    ax.set_xlabel('Position relative to TSS (bp)', fontsize=11)
    ax.set_yticks([0.2, 0.8])
    ax.set_yticklabels(['Coherence\nRegions', 'SURF\nPeaks'], fontsize=10)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='green', alpha=0.7, label='Stable'),
        Patch(facecolor='yellow', alpha=0.7, label='Mixed'),
        Patch(facecolor='red', alpha=0.7, label='Amplifying'),
        Patch(facecolor='darkgreen', alpha=0.8, label='SURF Activating'),
        Patch(facecolor='darkred', alpha=0.8, label='SURF Repressing'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=9)

    if title:
        ax.set_title(title, fontsize=14, fontweight='bold')
    else:
        ax.set_title(
            f"Region Classification: {result.surf.gene_symbol}",
            fontsize=14,
        )

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')

    return fig


def create_surf_report(
    result: SURFCoherenceResult,
    output_dir: str,
    gene_symbol: Optional[str] = None,
) -> List[str]:
    """
    Generate a complete visual report for SURF + coherence analysis.

    Creates multiple figures and saves them to output_dir.

    Args:
        result: SURFCoherenceResult.
        output_dir: Directory to save figures.
        gene_symbol: Gene name (uses result.surf.gene_symbol if not provided).

    Returns:
        List of generated file paths.
    """
    from pathlib import Path

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    gene = gene_symbol or result.surf.gene_symbol
    files = []

    # Main coherence plot
    fig = plot_surf_coherence(
        result,
        save_path=output_dir / f"{gene}_surf_coherence.png",
    )
    files.append(str(output_dir / f"{gene}_surf_coherence.png"))

    try:
        import matplotlib.pyplot as plt
        plt.close(fig)
    except:
        pass

    # Region map
    fig = plot_surf_regions(
        result,
        save_path=output_dir / f"{gene}_regions.png",
    )
    files.append(str(output_dir / f"{gene}_regions.png"))

    try:
        import matplotlib.pyplot as plt
        plt.close(fig)
    except:
        pass

    return files
