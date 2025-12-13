"""
PhaseLab CRISPOR Integration Module.

Provides automated off-target validation using a local CRISPOR installation.
CRISPOR is a separate tool that must be installed independently due to:
- Large genome index files (~6GB per species)
- External dependencies (BWA, scikit-learn with specific model versions)

This module wraps CRISPOR to provide:
- Automated guide validation
- Off-target counting by mismatch distance
- MIT/CFD specificity score validation
- Integration with PhaseLab's composite scoring

Installation:
    # Clone CRISPOR
    git clone https://github.com/maximilianh/crisporWebsite /path/to/crispor

    # Download genome (e.g., hg38)
    cd /path/to/crispor
    python tools/crisporDownloadGenome hg38

    # Configure PhaseLab
    from phaselab.crispr.crispor_integration import CrisporValidator
    validator = CrisporValidator(crispor_path="/path/to/crispor")
"""

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


@dataclass
class CrisporConfig:
    """Configuration for CRISPOR integration."""

    # Path to CRISPOR installation directory
    crispor_path: str

    # Python executable to use (should have CRISPOR dependencies)
    python_executable: Optional[str] = None

    # Genome to use (e.g., "hg38", "mm10", "danRer11")
    genome: str = "hg38"

    # Timeout for CRISPOR commands (seconds)
    timeout: int = 300

    # Whether to compute efficiency scores (slower, requires trained models)
    compute_efficiency: bool = False

    def __post_init__(self):
        # Validate paths
        crispor_py = Path(self.crispor_path) / "crispor.py"
        if not crispor_py.exists():
            raise FileNotFoundError(
                f"CRISPOR not found at {self.crispor_path}. "
                f"Expected to find crispor.py at {crispor_py}"
            )

        genome_dir = Path(self.crispor_path) / "genomes" / self.genome
        if not genome_dir.exists():
            raise FileNotFoundError(
                f"Genome {self.genome} not found at {genome_dir}. "
                f"Run: python {self.crispor_path}/tools/crisporDownloadGenome {self.genome}"
            )

        # Use system python if not specified
        if self.python_executable is None:
            self.python_executable = "python3"


class CrisporValidator:
    """
    Validates guide RNAs using local CRISPOR installation.

    Provides genome-wide off-target search and specificity scoring.

    Example:
        >>> validator = CrisporValidator("/path/to/crispor", genome="hg38")
        >>> results = validator.validate_guides(guides, sequence)
        >>> for guide in results:
        ...     print(f"{guide['sequence']}: MIT={guide['mit_specificity']}")
    """

    def __init__(
        self,
        crispor_path: str,
        genome: str = "hg38",
        python_executable: Optional[str] = None,
        compute_efficiency: bool = False,
    ):
        """
        Initialize CRISPOR validator.

        Args:
            crispor_path: Path to CRISPOR installation directory.
            genome: Genome identifier (e.g., "hg38", "mm10").
            python_executable: Python to use for running CRISPOR.
            compute_efficiency: Whether to compute efficiency scores.
        """
        self.config = CrisporConfig(
            crispor_path=crispor_path,
            genome=genome,
            python_executable=python_executable,
            compute_efficiency=compute_efficiency,
        )
        self._validate_installation()

    def _validate_installation(self):
        """Check that CRISPOR is properly installed."""
        # Check for required files
        required = [
            "crispor.py",
            f"genomes/{self.config.genome}/{self.config.genome}.2bit",
            f"genomes/{self.config.genome}/{self.config.genome}.fa.bwt",
        ]

        missing = []
        for f in required:
            path = Path(self.config.crispor_path) / f
            if not path.exists():
                missing.append(str(path))

        if missing:
            raise RuntimeError(
                f"CRISPOR installation incomplete. Missing: {missing}"
            )

        logger.info(f"CRISPOR validator initialized for genome {self.config.genome}")

    def validate_sequence(
        self,
        sequence: str,
        sequence_name: str = "input",
    ) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        """
        Run CRISPOR on a sequence and return all guides with off-target data.

        Args:
            sequence: Input DNA sequence (promoter region).
            sequence_name: Name for the sequence (used in output).

        Returns:
            (guides, off_targets)
            - guides: List of guide dicts with MIT/CFD scores
            - off_targets: List of off-target hit dicts
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            # Write input FASTA
            input_fa = Path(tmpdir) / "input.fa"
            with open(input_fa, 'w') as f:
                f.write(f">{sequence_name}\n{sequence}\n")

            # Output files
            output_tsv = Path(tmpdir) / "guides.tsv"
            ot_tsv = Path(tmpdir) / "offtargets.tsv"

            # Build command
            cmd = [
                self.config.python_executable,
                str(Path(self.config.crispor_path) / "crispor.py"),
                self.config.genome,
                str(input_fa),
                str(output_tsv),
                "-o", str(ot_tsv),
            ]

            if not self.config.compute_efficiency:
                cmd.append("--noEffScores")

            # Run CRISPOR
            logger.info(f"Running CRISPOR: {' '.join(cmd)}")
            try:
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=self.config.timeout,
                    cwd=self.config.crispor_path,
                )

                if result.returncode != 0:
                    logger.error(f"CRISPOR failed: {result.stderr}")
                    raise RuntimeError(f"CRISPOR error: {result.stderr}")

            except subprocess.TimeoutExpired:
                raise RuntimeError(
                    f"CRISPOR timed out after {self.config.timeout}s"
                )

            # Parse results
            guides = self._parse_guide_tsv(output_tsv)
            off_targets = self._parse_offtarget_tsv(ot_tsv) if ot_tsv.exists() else []

            return guides, off_targets

    def _parse_guide_tsv(self, filepath: Path) -> List[Dict[str, Any]]:
        """Parse CRISPOR guide output TSV.

        IMPORTANT: Guides with MIT=0, CFD=0, OTs=0 and "NotEnoughFlankSeq" are
        treated as INVALID/UNSCORABLE, not as "perfect guides with no off-targets".
        """
        guides = []

        if not filepath.exists():
            return guides

        with open(filepath, 'r') as f:
            header = None
            for line in f:
                if line.startswith('#'):
                    header = line[1:].strip().split('\t')
                    continue
                if not header:
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue

                mit = self._safe_float(parts[3])
                cfd = self._safe_float(parts[4])
                ot_count = int(parts[5]) if parts[5] else 0

                # Check for "NotEnoughFlankSeq" in any efficiency columns
                has_not_enough_flank = any(
                    'NotEnoughFlankSeq' in p for p in parts[6:] if isinstance(p, str)
                )

                # Detect invalid/unscorable guides (MIT=0, CFD=0, OTs=0)
                # These are NOT perfect guides - CRISPOR couldn't score them
                is_unscorable = (
                    (mit == 0 or mit is None) and
                    (cfd == 0 or cfd is None) and
                    ot_count == 0
                )

                guide = {
                    'guide_id': parts[1],
                    'sequence': parts[2][:20],  # 20bp without PAM
                    'sequence_with_pam': parts[2],  # Full 23bp
                    'mit_specificity': mit,
                    'cfd_specificity': cfd,
                    'offtarget_count': ot_count,
                    'not_enough_flank': has_not_enough_flank,
                    'is_unscorable': is_unscorable,
                }
                guides.append(guide)

        return guides

    def _parse_offtarget_tsv(self, filepath: Path) -> List[Dict[str, Any]]:
        """Parse CRISPOR off-target TSV and aggregate by guide."""
        from collections import defaultdict

        ot_by_guide = defaultdict(lambda: {0: 0, 1: 0, 2: 0, 3: 0, 4: 0})

        if not filepath.exists():
            return []

        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('seqId'):
                    continue  # Skip header

                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue

                guide_seq = parts[2][:20]
                mm_count = int(parts[5])

                if mm_count <= 4:
                    ot_by_guide[guide_seq][mm_count] += 1

        # Convert to list format
        return [
            {'sequence': seq, **{f'ot_{k}mm': v for k, v in counts.items()}}
            for seq, counts in ot_by_guide.items()
        ]

    @staticmethod
    def _safe_float(value: str) -> Optional[float]:
        """Parse float, returning None for invalid values."""
        if value in ('-1', '', 'NotEnoughFlankSeq'):
            return None
        try:
            return float(value)
        except ValueError:
            return None

    def validate_guides(
        self,
        guides: List[Dict[str, Any]],
        sequence: str,
        sequence_name: str = "input",
    ) -> List[Dict[str, Any]]:
        """
        Validate a list of guide dictionaries against CRISPOR.

        IMPORTANT: Only re-ranks within the PhaseLab candidate set.
        CRISPOR validation constrains but doesn't override PhaseLab's
        biological window logic.

        Args:
            guides: List of guide dicts (from design_guides()).
            sequence: The full sequence used for guide design.
            sequence_name: Name for the sequence.

        Returns:
            Guides with CRISPOR validation data merged in.
            Fields added:
            - mit_specificity: CRISPOR MIT score (None if unscorable)
            - cfd_specificity: CRISPOR CFD score (None if unscorable)
            - ot_0mm, ot_1mm, etc.: Off-target counts by mismatch
            - crispor_validated: True if CRISPOR could score this guide
            - is_unscorable: True if MIT=0, CFD=0, OTs=0 (invalid data)
        """
        # Run CRISPOR on the sequence
        crispor_guides, off_targets = self.validate_sequence(sequence, sequence_name)

        # Index by sequence
        crispor_index = {g['sequence'].upper(): g for g in crispor_guides}
        ot_index = {g['sequence'].upper(): g for g in off_targets}

        # Merge data into input guides
        validated = []
        for guide in guides:
            seq = guide.get('sequence', '').upper()
            result = {**guide}

            if seq in crispor_index:
                cg = crispor_index[seq]

                # Check if this guide is unscorable (MIT=0, CFD=0, OTs=0)
                is_unscorable = cg.get('is_unscorable', False)

                if is_unscorable:
                    # INVALID: Treat as unvalidated, not as "perfect"
                    result['mit_specificity'] = None
                    result['cfd_specificity'] = None
                    result['crispor_validated'] = False
                    result['is_unscorable'] = True
                    result['unscorable_reason'] = "CRISPOR returned MIT=0, CFD=0, OTs=0 (NotEnoughFlankSeq or unmappable)"
                else:
                    # Valid CRISPOR data
                    result['mit_specificity'] = cg.get('mit_specificity')
                    result['cfd_specificity'] = cg.get('cfd_specificity')
                    result['crispor_validated'] = True
                    result['is_unscorable'] = False

                result['crispor_ot_count'] = cg.get('offtarget_count', 0)

                # Add off-target breakdown if available
                if seq in ot_index:
                    ot = ot_index[seq]
                    result['ot_0mm'] = ot.get('ot_0mm', 0)
                    result['ot_1mm'] = ot.get('ot_1mm', 0)
                    result['ot_2mm'] = ot.get('ot_2mm', 0)
                    result['ot_3mm'] = ot.get('ot_3mm', 0)
                    result['ot_4mm'] = ot.get('ot_4mm', 0)
                else:
                    result['ot_0mm'] = 0
                    result['ot_1mm'] = 0
                    result['ot_2mm'] = 0
                    result['ot_3mm'] = 0
                    result['ot_4mm'] = 0

            else:
                # Not found in CRISPOR output at all
                result['crispor_validated'] = False
                result['is_unscorable'] = True
                result['unscorable_reason'] = "Guide not found in CRISPOR output"
                result['mit_specificity'] = None
                result['cfd_specificity'] = None
                result['ot_0mm'] = 0
                result['ot_1mm'] = 0
                result['ot_2mm'] = 0
                result['ot_3mm'] = 0
                result['ot_4mm'] = 0

            validated.append(result)

        return validated


def setup_crispor(
    install_path: str = "~/.phaselab/crispor",
    genome: str = "hg38",
) -> str:
    """
    Download and set up CRISPOR for local use.

    This is a convenience function that:
    1. Clones CRISPOR from GitHub
    2. Downloads the specified genome index

    Note: Genome downloads are large (~6GB for hg38).

    Args:
        install_path: Where to install CRISPOR.
        genome: Genome to download (default: hg38).

    Returns:
        Path to CRISPOR installation.

    Example:
        >>> crispor_path = setup_crispor()
        >>> validator = CrisporValidator(crispor_path)
    """
    install_path = Path(install_path).expanduser()

    # Clone CRISPOR if not present
    if not (install_path / "crispor.py").exists():
        logger.info(f"Cloning CRISPOR to {install_path}")
        install_path.parent.mkdir(parents=True, exist_ok=True)

        subprocess.run([
            "git", "clone",
            "https://github.com/maximilianh/crisporWebsite",
            str(install_path)
        ], check=True)

    # Download genome if not present
    genome_dir = install_path / "genomes" / genome
    if not (genome_dir / f"{genome}.2bit").exists():
        logger.info(f"Downloading genome {genome} (this may take a while)...")

        download_script = install_path / "tools" / "crisporDownloadGenome"
        if download_script.exists():
            subprocess.run([
                "python3", str(download_script), genome
            ], cwd=str(install_path), check=True)
        else:
            raise RuntimeError(
                f"Genome download script not found at {download_script}"
            )

    return str(install_path)


# Convenience function for one-shot validation
def validate_with_crispor(
    guides: List[Dict[str, Any]],
    sequence: str,
    crispor_path: str,
    genome: str = "hg38",
) -> List[Dict[str, Any]]:
    """
    One-shot guide validation with CRISPOR.

    Args:
        guides: List of guide dicts from design_guides().
        sequence: The sequence used for guide design.
        crispor_path: Path to CRISPOR installation.
        genome: Genome identifier.

    Returns:
        Guides with CRISPOR validation data.
    """
    validator = CrisporValidator(crispor_path, genome=genome)
    return validator.validate_guides(guides, sequence)
