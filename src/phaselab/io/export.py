"""
PhaseLab Export: File format utilities for CRISPR results.
"""

import json
from datetime import datetime
from pathlib import Path
from typing import Union, Optional
import pandas as pd


def export_guides_json(
    df: pd.DataFrame,
    filepath: Union[str, Path],
    metadata: Optional[dict] = None,
) -> None:
    """
    Export guide results to JSON format.

    Args:
        df: DataFrame from design_guides().
        filepath: Output file path.
        metadata: Optional metadata to include.
    """
    output = {
        'timestamp': datetime.now().isoformat(),
        'n_guides': len(df),
        'metadata': metadata or {},
        'guides': df.to_dict(orient='records'),
    }

    with open(filepath, 'w') as f:
        json.dump(output, f, indent=2, default=str)


def export_guides_fasta(
    df: pd.DataFrame,
    filepath: Union[str, Path],
    include_pam: bool = False,
) -> None:
    """
    Export guide sequences to FASTA format.

    Args:
        df: DataFrame with 'sequence' column.
        filepath: Output file path.
        include_pam: Append PAM to sequence.
    """
    with open(filepath, 'w') as f:
        for idx, row in df.iterrows():
            name = f"guide_{idx+1}_pos{row.get('position', 'NA')}"
            seq = row['sequence']
            if include_pam and 'pam' in row:
                seq += row['pam']
            f.write(f">{name}\n{seq}\n")


def export_crispor_batch(
    sequence: str,
    filepath: Union[str, Path],
    name: str = "promoter",
) -> None:
    """
    Export promoter sequence for CRISPOR batch submission.

    Args:
        sequence: Promoter DNA sequence.
        filepath: Output file path.
        name: Sequence name for FASTA header.
    """
    with open(filepath, 'w') as f:
        f.write(f">{name}\n{sequence}\n")
