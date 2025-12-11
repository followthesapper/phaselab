"""
PhaseLab Tissue Models: Pre-configured chromatin models for common tissues.

Provides tissue-specific chromatin accessibility models optimized for
different therapeutic targets.
"""

import numpy as np
from typing import Dict, List, Optional
from .accessibility import ChromatinModel, ChromatinScore, ChromatinState


# Available tissue models
AVAILABLE_TISSUES = [
    "brain",
    "brain_cortex",
    "brain_hippocampus",
    "liver",
    "blood",
    "blood_pbmc",
    "blood_t_cell",
    "muscle",
    "heart",
    "lung",
    "kidney",
    "generic",
]


def get_tissue_model(tissue: str) -> ChromatinModel:
    """
    Get pre-configured chromatin model for a tissue.

    Args:
        tissue: Tissue name from AVAILABLE_TISSUES

    Returns:
        ChromatinModel configured for the tissue
    """
    tissue = tissue.lower()

    if tissue.startswith("brain"):
        return _get_brain_model(tissue)
    elif tissue.startswith("blood"):
        return _get_blood_model(tissue)
    elif tissue == "liver":
        return _get_liver_model()
    elif tissue == "muscle":
        return _get_muscle_model()
    elif tissue == "heart":
        return _get_heart_model()
    elif tissue == "lung":
        return _get_lung_model()
    elif tissue == "kidney":
        return _get_kidney_model()
    else:
        return ChromatinModel(
            name="generic",
            tissue="generic"
        )


def _get_brain_model(tissue: str) -> ChromatinModel:
    """Get brain-specific chromatin model."""

    # Brain region specific settings
    if "cortex" in tissue:
        cell_type = "cortex"
        # Cortex has high neuronal accessibility
        dnase_weight = 0.45
    elif "hippocampus" in tissue:
        cell_type = "hippocampus"
        dnase_weight = 0.40
    else:
        cell_type = "mixed"
        dnase_weight = 0.40

    return ChromatinModel(
        name=f"brain_{cell_type}",
        tissue="brain",
        cell_type=cell_type,
        dnase_weight=dnase_weight,
        atac_weight=0.30,
        histone_weight=0.30,
        # Brain has more restrictive accessibility
        open_threshold=0.65,
        closed_threshold=0.35
    )


def _get_blood_model(tissue: str) -> ChromatinModel:
    """Get blood/immune cell chromatin model."""

    if "pbmc" in tissue:
        cell_type = "PBMC"
    elif "t_cell" in tissue:
        cell_type = "T_cell"
    else:
        cell_type = "mixed_blood"

    return ChromatinModel(
        name=f"blood_{cell_type}",
        tissue="blood",
        cell_type=cell_type,
        dnase_weight=0.35,
        atac_weight=0.40,  # ATAC works well in blood
        histone_weight=0.25,
        # Blood cells have more accessible chromatin
        open_threshold=0.60,
        closed_threshold=0.25
    )


def _get_liver_model() -> ChromatinModel:
    """Get liver chromatin model."""
    return ChromatinModel(
        name="liver_hepatocyte",
        tissue="liver",
        cell_type="hepatocyte",
        dnase_weight=0.40,
        atac_weight=0.30,
        histone_weight=0.30,
        open_threshold=0.70,
        closed_threshold=0.30
    )


def _get_muscle_model() -> ChromatinModel:
    """Get skeletal muscle chromatin model."""
    return ChromatinModel(
        name="muscle_skeletal",
        tissue="muscle",
        cell_type="skeletal",
        dnase_weight=0.40,
        atac_weight=0.35,
        histone_weight=0.25,
        # Muscle has more heterogeneous accessibility
        open_threshold=0.75,
        closed_threshold=0.35
    )


def _get_heart_model() -> ChromatinModel:
    """Get cardiac muscle chromatin model."""
    return ChromatinModel(
        name="heart_cardiomyocyte",
        tissue="heart",
        cell_type="cardiomyocyte",
        dnase_weight=0.45,
        atac_weight=0.30,
        histone_weight=0.25,
        open_threshold=0.70,
        closed_threshold=0.35
    )


def _get_lung_model() -> ChromatinModel:
    """Get lung epithelial chromatin model."""
    return ChromatinModel(
        name="lung_epithelial",
        tissue="lung",
        cell_type="epithelial",
        dnase_weight=0.40,
        atac_weight=0.35,
        histone_weight=0.25
    )


def _get_kidney_model() -> ChromatinModel:
    """Get kidney chromatin model."""
    return ChromatinModel(
        name="kidney_tubular",
        tissue="kidney",
        cell_type="tubular",
        dnase_weight=0.40,
        atac_weight=0.30,
        histone_weight=0.30
    )


# Convenience functions for common tissues


def brain_accessibility(
    chrom: str,
    position: int,
    region: str = "cortex"
) -> ChromatinScore:
    """
    Compute brain chromatin accessibility.

    Optimized for neurological disease targets like:
    - SCN2A (autism)
    - RAI1 (Smith-Magenis Syndrome)
    - Other CNS genes

    Args:
        chrom: Chromosome
        position: Genomic position
        region: Brain region (cortex, hippocampus)

    Returns:
        ChromatinScore with brain-specific accessibility
    """
    model = get_tissue_model(f"brain_{region}")
    return model.score_position(chrom, position)


def liver_accessibility(
    chrom: str,
    position: int
) -> ChromatinScore:
    """
    Compute liver chromatin accessibility.

    Useful for metabolic disease targets.

    Args:
        chrom: Chromosome
        position: Genomic position

    Returns:
        ChromatinScore for liver
    """
    model = get_tissue_model("liver")
    return model.score_position(chrom, position)


def blood_accessibility(
    chrom: str,
    position: int,
    cell_type: str = "pbmc"
) -> ChromatinScore:
    """
    Compute blood cell chromatin accessibility.

    Useful for:
    - Hematological disease targets
    - Ex vivo gene therapy (CAR-T, etc.)

    Args:
        chrom: Chromosome
        position: Genomic position
        cell_type: Blood cell type (pbmc, t_cell)

    Returns:
        ChromatinScore for blood cells
    """
    model = get_tissue_model(f"blood_{cell_type}")
    return model.score_position(chrom, position)


# Disease-specific accessibility functions


def scn2a_accessibility(position: int) -> ChromatinScore:
    """
    Compute chromatin accessibility for SCN2A locus.

    SCN2A is expressed in brain neurons, so uses cortex model.

    Args:
        position: Position within SCN2A promoter

    Returns:
        ChromatinScore for SCN2A
    """
    # SCN2A is on chromosome 2
    return brain_accessibility("chr2", 165239414 + position, region="cortex")


def rai1_accessibility(position: int) -> ChromatinScore:
    """
    Compute chromatin accessibility for RAI1 locus.

    RAI1 is expressed in brain, uses cortex model.

    Args:
        position: Position within RAI1 promoter

    Returns:
        ChromatinScore for RAI1
    """
    # RAI1 is on chromosome 17
    return brain_accessibility("chr17", 17681458 + position, region="cortex")
