"""
CRISPOR Client: Interface to CRISPOR CLI.

Handles running CRISPOR on sequences and parsing results.
"""

from __future__ import annotations

import subprocess
import tempfile
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List
import os


@dataclass
class CrisporConfig:
    """Configuration for CRISPOR client."""

    crispor_path: Path  # Path to crispor directory (containing crispor.py)
    genome: str = "hg38"
    genome_dir: Optional[Path] = None  # Path to genome index files
    pam: str = "NGG"
    max_mismatches: int = 4
    python_bin: str = "python3"
    timeout: int = 300  # seconds

    def __post_init__(self):
        self.crispor_path = Path(self.crispor_path)
        if self.genome_dir:
            self.genome_dir = Path(self.genome_dir)


@dataclass
class CrisporOutput:
    """Output from a CRISPOR run."""

    guides_tsv: Path
    offtargets_tsv: Optional[Path]
    stdout: str
    stderr: str
    success: bool
    output_dir: Path


class CrisporClient:
    """
    Client for running CRISPOR analyses.

    CRISPOR provides:
    - On-target activity scoring (Doench 2016, Moreno-Mateos)
    - Off-target enumeration and scoring (MIT, CFD)
    - Genome-wide specificity analysis

    Example:
        config = CrisporConfig(crispor_path="/opt/crispor")
        client = CrisporClient(config)

        result = client.score_sequence(
            sequence="ATCGATCG...",
            name="my_target"
        )
    """

    def __init__(self, config: CrisporConfig):
        self.config = config
        self._validate_config()

    def _validate_config(self):
        """Validate CRISPOR installation."""
        crispor_py = self.config.crispor_path / "crispor.py"
        if not crispor_py.exists():
            # Try alternate locations
            alt_paths = [
                self.config.crispor_path / "bin" / "crispor.py",
                self.config.crispor_path / "crispor" / "crispor.py",
            ]
            for alt in alt_paths:
                if alt.exists():
                    self.config.crispor_path = alt.parent
                    return
            raise FileNotFoundError(
                f"crispor.py not found at {crispor_py} or alternate locations. "
                f"Please verify CRISPOR installation path."
            )

    @property
    def crispor_script(self) -> Path:
        return self.config.crispor_path / "crispor.py"

    def score_sequence(
        self,
        sequence: str,
        name: str = "target",
        include_offtargets: bool = True,
        output_dir: Optional[Path] = None,
    ) -> CrisporOutput:
        """
        Run CRISPOR on a DNA sequence.

        Args:
            sequence: DNA sequence to analyze (will find all PAM sites)
            name: Name for the sequence in output
            include_offtargets: Whether to compute off-targets (slower but important)
            output_dir: Where to save results (temp dir if None)

        Returns:
            CrisporOutput with paths to result files
        """
        # Clean sequence
        sequence = sequence.upper().replace("\n", "").replace(" ", "")

        # Create output directory
        if output_dir is None:
            output_dir = Path(tempfile.mkdtemp(prefix="phaselab_crispor_"))
        else:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

        # Write input FASTA
        fasta_path = output_dir / "input.fa"
        fasta_path.write_text(f">{name}\n{sequence}\n")

        # Output files
        guides_tsv = output_dir / "guides.tsv"
        offtargets_tsv = output_dir / "offtargets.tsv" if include_offtargets else None

        # Build command
        cmd = [
            self.config.python_bin,
            str(self.crispor_script),
            self.config.genome,
            str(fasta_path),
            str(guides_tsv),
            "--pam", self.config.pam,
        ]

        if include_offtargets and offtargets_tsv:
            cmd.extend(["--offtargets", str(offtargets_tsv)])
            cmd.extend(["--mm", str(self.config.max_mismatches)])

        if self.config.genome_dir:
            cmd.extend(["--genomeDir", str(self.config.genome_dir)])

        # Run CRISPOR
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=self.config.timeout,
                cwd=str(self.config.crispor_path),
            )
            success = result.returncode == 0
        except subprocess.TimeoutExpired:
            return CrisporOutput(
                guides_tsv=guides_tsv,
                offtargets_tsv=offtargets_tsv,
                stdout="",
                stderr="CRISPOR timed out",
                success=False,
                output_dir=output_dir,
            )
        except Exception as e:
            return CrisporOutput(
                guides_tsv=guides_tsv,
                offtargets_tsv=offtargets_tsv,
                stdout="",
                stderr=str(e),
                success=False,
                output_dir=output_dir,
            )

        return CrisporOutput(
            guides_tsv=guides_tsv,
            offtargets_tsv=offtargets_tsv if offtargets_tsv and offtargets_tsv.exists() else None,
            stdout=result.stdout,
            stderr=result.stderr,
            success=success,
            output_dir=output_dir,
        )

    def score_guides(
        self,
        guides: List[str],
        include_offtargets: bool = True,
        output_dir: Optional[Path] = None,
    ) -> CrisporOutput:
        """
        Score specific guide sequences.

        Args:
            guides: List of 20bp guide sequences (without PAM)
            include_offtargets: Whether to compute off-targets
            output_dir: Where to save results

        Returns:
            CrisporOutput with paths to result files
        """
        # Build FASTA with each guide + NGG
        fasta_lines = []
        for i, guide in enumerate(guides):
            guide = guide.upper().replace(" ", "")
            # Add NGG PAM for CRISPOR to recognize
            seq_with_pam = guide + "NGG"
            fasta_lines.append(f">guide_{i+1}")
            fasta_lines.append(seq_with_pam)

        fasta_text = "\n".join(fasta_lines) + "\n"

        # Create output directory
        if output_dir is None:
            output_dir = Path(tempfile.mkdtemp(prefix="phaselab_crispor_"))
        else:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

        # Write input FASTA
        fasta_path = output_dir / "input.fa"
        fasta_path.write_text(fasta_text)

        # Run with noGenome mode for quick scoring if no off-targets needed
        guides_tsv = output_dir / "guides.tsv"
        offtargets_tsv = output_dir / "offtargets.tsv" if include_offtargets else None

        cmd = [
            self.config.python_bin,
            str(self.crispor_script),
        ]

        if include_offtargets:
            cmd.append(self.config.genome)
        else:
            cmd.append("noGenome")

        cmd.extend([
            str(fasta_path),
            str(guides_tsv),
            "--pam", self.config.pam,
        ])

        if include_offtargets and offtargets_tsv:
            cmd.extend(["--offtargets", str(offtargets_tsv)])
            cmd.extend(["--mm", str(self.config.max_mismatches)])

        if self.config.genome_dir:
            cmd.extend(["--genomeDir", str(self.config.genome_dir)])

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=self.config.timeout,
                cwd=str(self.config.crispor_path),
            )
            success = result.returncode == 0
        except Exception as e:
            return CrisporOutput(
                guides_tsv=guides_tsv,
                offtargets_tsv=offtargets_tsv,
                stdout="",
                stderr=str(e),
                success=False,
                output_dir=output_dir,
            )

        return CrisporOutput(
            guides_tsv=guides_tsv,
            offtargets_tsv=offtargets_tsv if offtargets_tsv and offtargets_tsv.exists() else None,
            stdout=result.stdout,
            stderr=result.stderr,
            success=success,
            output_dir=output_dir,
        )

    @staticmethod
    def is_available(crispor_path: Path) -> bool:
        """Check if CRISPOR is available at the given path."""
        crispor_py = Path(crispor_path) / "crispor.py"
        return crispor_py.exists()
