"""
PhaseLab Landscapes: I/O utilities for loading and exporting landscape data.

Supports common formats from:
- CRISPRa/CRISPRi tiling screens
- TnSeq/RB-TnSeq fitness data
- Deep mutational scanning
- Generic CSV/TSV response matrices
"""

import numpy as np
import json
from pathlib import Path
from typing import Optional, Dict, Any, List, Union
from dataclasses import asdict

from .core import (
    ResponseLandscape,
    CoherenceProfile,
    RegionClassification,
    LandscapeMetrics,
)


def load_tiling_data(
    filepath: Union[str, Path],
    coord_col: str = 'position',
    response_col: str = 'logFC',
    replicate_cols: Optional[List[str]] = None,
    delimiter: str = '\t',
    skip_header: int = 0,
    metadata: Optional[Dict[str, Any]] = None,
) -> ResponseLandscape:
    """
    Load tiling screen data from a tabular file.

    Expects a file with columns for position and response values.
    Optionally can load multiple replicate columns.

    Args:
        filepath: Path to the data file.
        coord_col: Name of the position/coordinate column.
        response_col: Name of the response column (for single replicate).
        replicate_cols: List of column names for replicate responses.
        delimiter: Column delimiter (default tab).
        skip_header: Number of header lines to skip.
        metadata: Additional metadata to attach.

    Returns:
        ResponseLandscape with loaded data.

    Example:
        >>> # Single replicate
        >>> landscape = load_tiling_data(
        ...     'cd69_tiling.tsv',
        ...     coord_col='position',
        ...     response_col='logFC_mean'
        ... )

        >>> # Multiple replicates
        >>> landscape = load_tiling_data(
        ...     'screen_data.tsv',
        ...     coord_col='pos',
        ...     replicate_cols=['rep1_logFC', 'rep2_logFC', 'rep3_logFC']
        ... )
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"Data file not found: {filepath}")

    # Try pandas first (more flexible)
    try:
        import pandas as pd
        df = pd.read_csv(filepath, delimiter=delimiter, skiprows=skip_header)

        coords = df[coord_col].values

        if replicate_cols:
            responses = df[replicate_cols].values
        else:
            responses = df[response_col].values

        coord_labels = [f"{coord_col}_{i}" for i in range(len(coords))]

        file_metadata = {
            'source_file': str(filepath),
            'coord_column': coord_col,
            'response_columns': replicate_cols or [response_col],
            'n_rows': len(df),
        }

        if metadata:
            file_metadata.update(metadata)

        return ResponseLandscape(
            coords=coords,
            responses=responses,
            coord_labels=coord_labels,
            response_labels=replicate_cols or [response_col],
            metadata=file_metadata,
        )

    except ImportError:
        # Fallback to numpy
        data = np.loadtxt(filepath, delimiter=delimiter, skiprows=skip_header + 1)

        # Assume first column is coords, rest are responses
        coords = data[:, 0]
        responses = data[:, 1:] if data.shape[1] > 2 else data[:, 1]

        return ResponseLandscape(
            coords=coords,
            responses=responses,
            metadata=metadata or {'source_file': str(filepath)},
        )


def load_response_matrix(
    filepath: Union[str, Path],
    coords: Optional[np.ndarray] = None,
    transpose: bool = False,
) -> ResponseLandscape:
    """
    Load a response matrix from a file.

    Expects a matrix where rows are positions and columns are replicates
    (or vice versa if transpose=True).

    Args:
        filepath: Path to matrix file (CSV, TSV, or NPY).
        coords: Optional coordinate array. If not provided, uses indices.
        transpose: Whether to transpose the matrix.

    Returns:
        ResponseLandscape with matrix data.
    """
    filepath = Path(filepath)

    if filepath.suffix == '.npy':
        data = np.load(filepath)
    else:
        delimiter = ',' if filepath.suffix == '.csv' else '\t'
        data = np.loadtxt(filepath, delimiter=delimiter)

    if transpose:
        data = data.T

    if coords is None:
        coords = np.arange(data.shape[0])

    responses = data if data.ndim > 1 else data.reshape(-1, 1)

    return ResponseLandscape(
        coords=coords,
        responses=responses,
        metadata={'source_file': str(filepath), 'transposed': transpose},
    )


def load_crispr_surf_output(
    filepath: Union[str, Path],
    signal_col: str = 'beta',
) -> ResponseLandscape:
    """
    Load CRISPR-SURF deconvolution output.

    CRISPR-SURF produces beta profiles representing deconvolved
    regulatory signal.

    Args:
        filepath: Path to CRISPR-SURF output file.
        signal_col: Column name for the deconvolved signal.

    Returns:
        ResponseLandscape with SURF data.
    """
    try:
        import pandas as pd
        df = pd.read_csv(filepath, sep='\t')

        # SURF output typically has chr, start, end, beta columns
        if 'start' in df.columns:
            coords = df['start'].values
        elif 'position' in df.columns:
            coords = df['position'].values
        else:
            coords = np.arange(len(df))

        responses = df[signal_col].values

        return ResponseLandscape(
            coords=coords,
            responses=responses,
            metadata={
                'source': 'CRISPR-SURF',
                'source_file': str(filepath),
                'signal_column': signal_col,
            },
        )

    except ImportError:
        raise ImportError("pandas required for CRISPR-SURF loading")


def export_landscape(
    landscape: ResponseLandscape,
    filepath: Union[str, Path],
    format: str = 'json',
) -> None:
    """
    Export a ResponseLandscape to file.

    Args:
        landscape: Landscape to export.
        filepath: Output path.
        format: Output format ('json', 'tsv', 'csv', 'npy').
    """
    filepath = Path(filepath)

    if format == 'json':
        with open(filepath, 'w') as f:
            json.dump(landscape.to_dict(), f, indent=2)

    elif format in ('tsv', 'csv'):
        delimiter = '\t' if format == 'tsv' else ','

        try:
            import pandas as pd

            data = {'coord': landscape.coords}
            if landscape.responses.ndim == 1:
                data['response'] = landscape.responses
            else:
                for i in range(landscape.n_replicates):
                    col_name = (
                        landscape.response_labels[i]
                        if landscape.response_labels
                        else f'response_{i}'
                    )
                    data[col_name] = landscape.responses[:, i]

            df = pd.DataFrame(data)
            df.to_csv(filepath, sep=delimiter, index=False)

        except ImportError:
            # Fallback to numpy
            if landscape.responses.ndim == 1:
                data = np.column_stack([landscape.coords, landscape.responses])
            else:
                data = np.column_stack([landscape.coords, landscape.responses])
            np.savetxt(filepath, data, delimiter=delimiter)

    elif format == 'npy':
        np.save(filepath, {
            'coords': landscape.coords,
            'responses': landscape.responses,
            'metadata': landscape.metadata,
        })

    else:
        raise ValueError(f"Unknown format: {format}")


def export_regions(
    classification: RegionClassification,
    filepath: Union[str, Path],
    format: str = 'json',
    include_profile: bool = False,
) -> None:
    """
    Export region classification to file.

    Args:
        classification: RegionClassification to export.
        filepath: Output path.
        format: Output format ('json', 'bed', 'tsv').
        include_profile: Whether to include full coherence profile.
    """
    filepath = Path(filepath)

    if format == 'json':
        data = classification.to_dict()
        if include_profile:
            data['profile'] = classification.profile.to_dict()

        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2)

    elif format == 'bed':
        # BED format for genome browsers
        chrom = classification.landscape.metadata.get('chromosome', 'chr1')

        with open(filepath, 'w') as f:
            f.write(f"# Region classification from PhaseLab\n")
            for start, end, cls, score in classification.regions:
                # BED: chrom, start, end, name, score
                name = cls.value
                bed_score = int(score * 1000)  # BED scores are 0-1000
                f.write(f"{chrom}\t{int(start)}\t{int(end)}\t{name}\t{bed_score}\n")

    elif format == 'tsv':
        with open(filepath, 'w') as f:
            f.write("start\tend\tclass\tscore\n")
            for start, end, cls, score in classification.regions:
                f.write(f"{start}\t{end}\t{cls.value}\t{score}\n")

    else:
        raise ValueError(f"Unknown format: {format}")


def export_profile(
    profile: CoherenceProfile,
    filepath: Union[str, Path],
    format: str = 'json',
) -> None:
    """
    Export a CoherenceProfile to file.

    Args:
        profile: CoherenceProfile to export.
        filepath: Output path.
        format: Output format ('json', 'tsv', 'wig').
    """
    filepath = Path(filepath)

    if format == 'json':
        with open(filepath, 'w') as f:
            json.dump(profile.to_dict(), f, indent=2)

    elif format == 'tsv':
        with open(filepath, 'w') as f:
            f.write("coord\tcoherence\tlocal_variance\n")
            for coord, coh, var in zip(profile.coords, profile.coherence, profile.local_variance):
                f.write(f"{coord}\t{coh}\t{var}\n")

    elif format == 'wig':
        # WIG format for genome browsers
        with open(filepath, 'w') as f:
            f.write(f"# Coherence profile from PhaseLab\n")
            f.write(f"variableStep chrom=chr1\n")
            for coord, coh in zip(profile.coords, profile.coherence):
                f.write(f"{int(coord)}\t{coh}\n")

    else:
        raise ValueError(f"Unknown format: {format}")


def load_multiple_landscapes(
    filepaths: List[Union[str, Path]],
    labels: Optional[List[str]] = None,
    **load_kwargs,
) -> Dict[str, ResponseLandscape]:
    """
    Load multiple landscape files into a dictionary.

    Args:
        filepaths: List of file paths.
        labels: Optional labels for each file.
        **load_kwargs: Additional arguments for load_tiling_data.

    Returns:
        Dictionary mapping label → ResponseLandscape.
    """
    if labels is None:
        labels = [Path(f).stem for f in filepaths]

    landscapes = {}
    for filepath, label in zip(filepaths, labels):
        landscapes[label] = load_tiling_data(filepath, **load_kwargs)

    return landscapes
