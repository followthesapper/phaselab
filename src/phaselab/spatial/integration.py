"""
PhaseLab Spatial: Integration with external tools and data sources.

Integrates spatial coherence analysis with:
- CRISPOR (off-target scoring)
- Chromatin accessibility (ATAC-seq, DNase-seq)
- BioGRID (gene interaction context)
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
from pathlib import Path

from .regulatory import (
    RegulatoryLandscape,
    RegulatoryRegion,
    classify_regulatory_regions,
    compute_regulatory_coherence,
)
from .targeting import rank_guides_in_regions
from ..landscapes.core import StabilityClass


@dataclass
class IntegratedAnalysis:
    """
    Complete integrated analysis result.

    Combines:
    - Spatial coherence (region classification)
    - CRISPOR (guide-level scoring)
    - Chromatin (accessibility context)
    - BioGRID (gene network context)
    """
    gene_symbol: str
    regions: List[RegulatoryRegion]
    guides: List[Dict[str, Any]]
    chromatin_context: Optional[Dict[str, Any]] = None
    network_context: Optional[Dict[str, Any]] = None
    validation: Optional[Dict[str, Any]] = None
    recommendations: List[str] = field(default_factory=list)

    @property
    def n_stable_regions(self) -> int:
        """Number of stable regions."""
        return sum(1 for r in self.regions if r.is_safe)

    @property
    def n_safe_guides(self) -> int:
        """Number of guides in stable regions."""
        return sum(1 for g in self.guides if g.get('in_stable_region', False))

    @property
    def top_guides(self) -> List[Dict[str, Any]]:
        """Top-ranked guides in stable regions."""
        safe_guides = [g for g in self.guides if g.get('in_stable_region', False)]
        return sorted(safe_guides, key=lambda x: x.get('rank', 999))[:10]

    def summary(self) -> str:
        """Generate human-readable summary."""
        lines = [
            "=" * 70,
            f"INTEGRATED ANALYSIS: {self.gene_symbol}",
            "=" * 70,
            "",
            "SPATIAL COHERENCE:",
            f"  Total regions: {len(self.regions)}",
            f"  Stable regions: {self.n_stable_regions}",
            "",
            "GUIDE ANALYSIS:",
            f"  Total guides: {len(self.guides)}",
            f"  Guides in stable regions: {self.n_safe_guides}",
        ]

        if self.chromatin_context:
            lines.extend([
                "",
                "CHROMATIN CONTEXT:",
                f"  Accessibility: {self.chromatin_context.get('accessibility_status', 'unknown')}",
            ])

        if self.network_context:
            lines.extend([
                "",
                "NETWORK CONTEXT:",
                f"  BioGRID interactions: {self.network_context.get('n_interactions', 'unknown')}",
            ])

        if self.recommendations:
            lines.extend(["", "RECOMMENDATIONS:"])
            for rec in self.recommendations[:5]:
                lines.append(f"  - {rec}")

        lines.append("=" * 70)
        return "\n".join(lines)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'gene_symbol': self.gene_symbol,
            'regions': [r.to_dict() for r in self.regions],
            'n_guides': len(self.guides),
            'n_safe_guides': self.n_safe_guides,
            'top_guides': self.top_guides,
            'chromatin_context': self.chromatin_context,
            'network_context': self.network_context,
            'validation': self.validation,
            'recommendations': self.recommendations,
        }


def integrate_with_crispor(
    landscape: RegulatoryLandscape,
    crispor_guides: List[Dict[str, Any]],
    regions: Optional[List[RegulatoryRegion]] = None,
    position_key: str = 'position',
) -> List[Dict[str, Any]]:
    """
    Integrate CRISPOR guide data with spatial coherence regions.

    CRISPOR provides:
    - Off-target counts (0mm, 1mm, 2mm, 3mm, 4mm)
    - MIT specificity score
    - CFD score
    - Doench on-target score

    Spatial coherence provides:
    - Region stability classification
    - Expected variance reduction
    - Amplification warnings

    Args:
        landscape: RegulatoryLandscape for context.
        crispor_guides: List of guides from CRISPOR.
        regions: Pre-computed regions (computed if not provided).
        position_key: Key for TSS-relative position in guide data.

    Returns:
        List of guides with spatial coherence annotations.

    Example:
        >>> # After running CRISPOR
        >>> from phaselab.integrations.crispor import parse_guides_tsv
        >>> crispor_guides = parse_guides_tsv('crispor_output.tsv')
        >>>
        >>> # Integrate with spatial coherence
        >>> integrated = integrate_with_crispor(landscape, crispor_guides)
        >>> safe_guides = [g for g in integrated if g['in_stable_region']]
    """
    if regions is None:
        regions = classify_regulatory_regions(landscape)

    # Rank guides within regions
    ranked = rank_guides_in_regions(
        crispor_guides,
        regions,
        position_key=position_key,
        score_key='mit_specificity',
        higher_is_better=True,
    )

    # Add integration metadata
    for guide in ranked:
        guide['integration_source'] = 'crispor'
        guide['spatial_coherence_validated'] = True

        # Add warnings for guides in amplifying regions
        if guide.get('region_stability') == 'amplifying':
            guide['warnings'] = guide.get('warnings', [])
            guide['warnings'].append(
                "CAUTION: Guide is in an amplifying region - results may be unpredictable"
            )

    return ranked


def integrate_with_chromatin(
    landscape: RegulatoryLandscape,
    regions: List[RegulatoryRegion],
    atac_data: Optional[np.ndarray] = None,
    dnase_data: Optional[np.ndarray] = None,
    histone_marks: Optional[Dict[str, np.ndarray]] = None,
) -> Dict[str, Any]:
    """
    Integrate chromatin accessibility data with region classification.

    Chromatin state affects CRISPRa/CRISPRi efficacy:
    - Open chromatin → better dCas9 access
    - Active marks (H3K27ac) → transcriptionally active
    - Closed chromatin → may need chromatin remodelers

    Args:
        landscape: RegulatoryLandscape for coordinate context.
        regions: Classified regulatory regions.
        atac_data: ATAC-seq signal (same coordinates as landscape).
        dnase_data: DNase-seq signal.
        histone_marks: Dict of histone ChIP-seq signals (e.g., 'H3K27ac').

    Returns:
        Dictionary with chromatin integration results.
    """
    results = {
        'gene': landscape.gene_symbol,
        'regions_analyzed': len(regions),
        'chromatin_data_available': {
            'atac': atac_data is not None,
            'dnase': dnase_data is not None,
            'histones': list(histone_marks.keys()) if histone_marks else [],
        },
        'region_chromatin': [],
    }

    for region in regions:
        region_result = {
            'start': region.start,
            'end': region.end,
            'stability': region.stability.value,
        }

        # Get region mask
        mask = (landscape.positions >= region.start) & (landscape.positions <= region.end)

        if atac_data is not None and len(atac_data) == len(landscape.positions):
            region_atac = atac_data[mask]
            region_result['mean_atac'] = float(np.nanmean(region_atac))
            region_result['atac_status'] = (
                'OPEN' if np.nanmean(region_atac) > 0.5 else 'CLOSED'
            )

        if dnase_data is not None and len(dnase_data) == len(landscape.positions):
            region_dnase = dnase_data[mask]
            region_result['mean_dnase'] = float(np.nanmean(region_dnase))

        if histone_marks:
            for mark, data in histone_marks.items():
                if len(data) == len(landscape.positions):
                    region_data = data[mask]
                    region_result[f'mean_{mark}'] = float(np.nanmean(region_data))

        results['region_chromatin'].append(region_result)

    # Overall accessibility status
    if atac_data is not None:
        stable_regions = [r for r in results['region_chromatin']
                        if r['stability'] == 'stable']
        if stable_regions:
            stable_atac = np.mean([r.get('mean_atac', 0) for r in stable_regions])
            results['accessibility_status'] = 'FAVORABLE' if stable_atac > 0.5 else 'CHECK_ACCESS'
            results['stable_region_mean_atac'] = stable_atac

    return results


def integrate_with_biogrid(
    gene_symbol: str,
    regions: List[RegulatoryRegion],
    biogrid_data: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Integrate BioGRID gene interaction data.

    BioGRID provides context about:
    - Gene interaction partners
    - CRISPR screen phenotypes (ORCS)
    - Genetic interactions

    Args:
        gene_symbol: Gene to look up.
        regions: Classified regulatory regions.
        biogrid_data: Pre-loaded BioGRID data (or fetched if None).

    Returns:
        Dictionary with network context.
    """
    # Placeholder for BioGRID integration
    # In production, this would fetch from BioGRID API or local database

    results = {
        'gene': gene_symbol,
        'n_interactions': 0,
        'interaction_types': {},
        'crispr_phenotypes': [],
        'recommendations': [],
    }

    if biogrid_data:
        # Process BioGRID data
        results['n_interactions'] = biogrid_data.get('interaction_count', 0)
        results['interaction_types'] = biogrid_data.get('types', {})

        # Add ORCS CRISPR screen data if available
        if 'orcs' in biogrid_data:
            results['crispr_phenotypes'] = biogrid_data['orcs']

    # Generate recommendations based on network context
    if results['n_interactions'] > 100:
        results['recommendations'].append(
            f"{gene_symbol} is a highly connected gene ({results['n_interactions']} interactions). "
            f"Consider potential off-target network effects."
        )

    return results


def run_integrated_analysis(
    landscape: RegulatoryLandscape,
    crispor_guides: Optional[List[Dict[str, Any]]] = None,
    atac_data: Optional[np.ndarray] = None,
    biogrid_data: Optional[Dict[str, Any]] = None,
    window: int = 50,
) -> IntegratedAnalysis:
    """
    Run complete integrated analysis combining all data sources.

    This is the top-level integration function that combines:
    1. Spatial coherence (region classification)
    2. CRISPOR (guide scoring)
    3. Chromatin (accessibility)
    4. BioGRID (network context)

    Args:
        landscape: RegulatoryLandscape to analyze.
        crispor_guides: CRISPOR guide data (optional).
        atac_data: ATAC-seq accessibility data (optional).
        biogrid_data: BioGRID interaction data (optional).
        window: Window size for coherence computation.

    Returns:
        IntegratedAnalysis with complete results.

    Example:
        >>> analysis = run_integrated_analysis(
        ...     landscape,
        ...     crispor_guides=guides,
        ...     atac_data=atac,
        ... )
        >>> print(analysis.summary())
        >>> safe_guides = analysis.top_guides
    """
    # Step 1: Classify regions
    regions = classify_regulatory_regions(landscape, window=window)
    profile = compute_regulatory_coherence(landscape, window=window)

    # Step 2: Integrate CRISPOR guides
    if crispor_guides:
        guides = integrate_with_crispor(landscape, crispor_guides, regions)
    else:
        guides = []

    # Step 3: Integrate chromatin
    chromatin_context = None
    if atac_data is not None:
        chromatin_context = integrate_with_chromatin(landscape, regions, atac_data=atac_data)

    # Step 4: Integrate BioGRID
    network_context = None
    if biogrid_data:
        network_context = integrate_with_biogrid(landscape.gene_symbol, regions, biogrid_data)

    # Step 5: Generate recommendations
    recommendations = _generate_recommendations(
        landscape, regions, guides, chromatin_context, network_context, profile
    )

    # Step 6: Validation summary
    validation = {
        'coherence_variance_correlation': profile.correlation,
        'p_value': profile.p_value,
        'is_validated': profile.is_validated,
        'variance_reduction_estimate': profile.variance_reduction_estimate,
    }

    return IntegratedAnalysis(
        gene_symbol=landscape.gene_symbol,
        regions=regions,
        guides=guides,
        chromatin_context=chromatin_context,
        network_context=network_context,
        validation=validation,
        recommendations=recommendations,
    )


def _generate_recommendations(
    landscape: RegulatoryLandscape,
    regions: List[RegulatoryRegion],
    guides: List[Dict[str, Any]],
    chromatin: Optional[Dict[str, Any]],
    network: Optional[Dict[str, Any]],
    profile,
) -> List[str]:
    """Generate recommendations based on integrated analysis."""
    recommendations = []

    # Region-based recommendations
    stable = [r for r in regions if r.is_safe]
    amplifying = [r for r in regions if r.is_dangerous]

    if stable:
        best = max(stable, key=lambda r: r.coherence_score)
        recommendations.append(
            f"TARGET REGION: [{best.start}, {best.end}] relative to TSS - "
            f"{best.variance_reduction:.0%} expected variance reduction"
        )

    if amplifying:
        recommendations.append(
            f"AVOID: {len(amplifying)} amplifying region(s) detected - "
            f"these show MYC-like behavior"
        )

    # Validation status
    if profile.is_validated:
        recommendations.append(
            f"VALIDATED: Spatial coherence model confirmed (r={profile.correlation:.3f})"
        )
    else:
        recommendations.append(
            f"CAUTION: Spatial coherence not strongly validated for this gene"
        )

    # Guide recommendations
    if guides:
        safe_guides = [g for g in guides if g.get('in_stable_region', False)]
        if safe_guides:
            recommendations.append(
                f"GUIDES: {len(safe_guides)} guides available in stable regions"
            )

            # Top guide
            if safe_guides[0].get('sequence'):
                recommendations.append(
                    f"TOP GUIDE: {safe_guides[0]['sequence']} "
                    f"(MIT={safe_guides[0].get('mit_specificity', 'N/A')})"
                )

    # Chromatin recommendations
    if chromatin and chromatin.get('accessibility_status') == 'CHECK_ACCESS':
        recommendations.append(
            "CHROMATIN: Stable regions show low accessibility - consider dCas9-p300 or similar"
        )

    return recommendations
