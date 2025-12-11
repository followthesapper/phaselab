"""
PhaseLab Report Generator: Generate comprehensive validation reports.

Implements:
- Guide ranking summary reports
- IR coherence validation reports
- Therapeutic assessment reports
- Multi-format export (JSON, HTML, Markdown)
- CRISPOR compatibility links

Version: 0.5.0
"""

import json
from typing import Dict, List, Optional, Any, Union
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


@dataclass
class ReportMetadata:
    """Report metadata."""
    title: str
    generated_at: str
    phaselab_version: str = "0.5.0"
    target_gene: Optional[str] = None
    target_disease: Optional[str] = None
    author: Optional[str] = None


@dataclass
class GuideReport:
    """Individual guide report."""
    rank: int
    sequence: str
    position: int
    strand: str

    # Core scores
    combined_score: float
    coherence_R: float
    go_no_go: str

    # Detailed scores
    gc_content: float
    delta_g: float
    mit_score: float
    cfd_score: float

    # Optional application-specific scores
    cut_efficiency: Optional[float] = None
    repression_efficiency: Optional[float] = None
    activation_potential: Optional[float] = None

    # Warnings
    warnings: List[str] = None


def generate_guide_report(
    guides: List[Dict[str, Any]],
    gene: str,
    application: str = "CRISPRa",
    top_n: int = 10,
) -> Dict[str, Any]:
    """
    Generate comprehensive guide ranking report.

    Args:
        guides: List of guide dictionaries from design pipeline
        gene: Target gene symbol
        application: CRISPR application (CRISPRa, CRISPRi, knockout, etc.)
        top_n: Number of top guides to include

    Returns:
        Report dictionary

    Example:
        >>> report = generate_guide_report(guides, "RAI1", "CRISPRa")
        >>> export_report(report, "/path/to/report.json")
    """
    report = {
        'metadata': {
            'title': f'{application} Guide Design Report - {gene}',
            'generated_at': datetime.now().isoformat(),
            'phaselab_version': '0.5.0',
            'target_gene': gene,
            'application': application,
        },
        'summary': {
            'total_candidates_evaluated': len(guides),
            'top_n_reported': min(top_n, len(guides)),
            'go_count': sum(1 for g in guides if g.get('go_no_go') == 'GO'),
            'no_go_count': sum(1 for g in guides if g.get('go_no_go') == 'NO-GO'),
        },
        'guides': [],
        'quality_metrics': {},
    }

    # Process top guides
    for i, guide in enumerate(guides[:top_n]):
        guide_report = {
            'rank': i + 1,
            'sequence': guide.get('sequence', ''),
            'position': guide.get('position', guide.get('cds_position', 0)),
            'strand': guide.get('strand', '+'),
            'combined_score': round(guide.get('combined_score', 0), 3),
            'coherence_R': round(guide.get('coherence_R', 0), 4),
            'go_no_go': guide.get('go_no_go', 'UNKNOWN'),
            'gc_content': round(guide.get('gc', 0), 3),
            'delta_g': round(guide.get('delta_g', 0), 2),
            'mit_score': round(guide.get('mit_score', 0), 1),
            'cfd_score': round(guide.get('cfd_score', 0), 1),
        }

        # Application-specific scores
        if 'cut_efficiency' in guide:
            guide_report['cut_efficiency'] = round(guide['cut_efficiency'], 3)
        if 'repression_efficiency' in guide:
            guide_report['repression_efficiency'] = round(guide['repression_efficiency'], 3)

        report['guides'].append(guide_report)

    # Quality metrics summary
    go_guides = [g for g in guides if g.get('go_no_go') == 'GO']
    if go_guides:
        report['quality_metrics'] = {
            'mean_coherence': round(
                sum(g.get('coherence_R', 0) for g in go_guides) / len(go_guides), 4
            ),
            'mean_combined_score': round(
                sum(g.get('combined_score', 0) for g in go_guides) / len(go_guides), 3
            ),
            'gc_range': (
                round(min(g.get('gc', 0) for g in go_guides), 2),
                round(max(g.get('gc', 0) for g in go_guides), 2),
            ),
        }

    return report


def generate_therapeutic_report(
    guides: List[Dict[str, Any]],
    gene: str,
    disease: str,
    therapeutic_window: Dict[str, float],
    dosage_results: Optional[Dict] = None,
) -> Dict[str, Any]:
    """
    Generate therapeutic assessment report.

    Args:
        guides: Guide candidates
        gene: Target gene
        disease: Disease name
        therapeutic_window: Therapeutic expression window
        dosage_results: Optional dosage optimization results

    Returns:
        Therapeutic report dictionary
    """
    report = {
        'metadata': {
            'title': f'Therapeutic Assessment Report - {disease}',
            'generated_at': datetime.now().isoformat(),
            'phaselab_version': '0.5.0',
            'target_gene': gene,
            'disease': disease,
        },
        'therapeutic_window': {
            'baseline_expression': therapeutic_window.get('baseline', 0.5),
            'therapeutic_min': therapeutic_window.get('therapeutic_min', 0.7),
            'therapeutic_max': therapeutic_window.get('therapeutic_max', 1.2),
            'optimal': therapeutic_window.get('optimal', 0.85),
        },
        'recommended_guides': [],
        'therapeutic_assessment': {},
    }

    # Filter to GO guides
    go_guides = [g for g in guides if g.get('go_no_go') == 'GO']

    for guide in go_guides[:5]:
        guide_entry = {
            'sequence': guide.get('sequence', ''),
            'position': guide.get('position', 0),
            'coherence_R': round(guide.get('coherence_R', 0), 4),
            'predicted_fold_change': round(guide.get('estimated_fold_change', 1.5), 2),
        }

        # Check if achieves therapeutic level
        baseline = therapeutic_window.get('baseline', 0.5)
        fold_change = guide.get('estimated_fold_change', 1.5)
        achieved = baseline * fold_change

        guide_entry['predicted_expression'] = round(achieved, 2)
        guide_entry['in_therapeutic_window'] = (
            therapeutic_window.get('therapeutic_min', 0.7) <=
            achieved <=
            therapeutic_window.get('therapeutic_max', 1.2)
        )

        report['recommended_guides'].append(guide_entry)

    # Dosage results if provided
    if dosage_results:
        report['dosage_assessment'] = {
            'estimated_fold_change': round(dosage_results.get('estimated_fold_change', 0), 2),
            'achieved_expression': round(dosage_results.get('achieved_expression', 0), 2),
            'in_therapeutic_window': dosage_results.get('in_therapeutic_window', False),
            'recommendation': dosage_results.get('recommendation', ''),
        }

    return report


def generate_coherence_report(
    coherence_values: List[float],
    source: str = "guide_design",
    hardware_validation: Optional[Dict] = None,
) -> Dict[str, Any]:
    """
    Generate IR coherence validation report.

    Args:
        coherence_values: List of R̄ values
        source: Source of coherence measurements
        hardware_validation: Optional IBM Quantum hardware data

    Returns:
        Coherence report dictionary
    """
    import numpy as np
    from ..core.constants import E_MINUS_2

    coherence_array = np.array(coherence_values)

    report = {
        'metadata': {
            'title': 'IR Coherence Validation Report',
            'generated_at': datetime.now().isoformat(),
            'phaselab_version': '0.5.0',
            'source': source,
        },
        'coherence_statistics': {
            'count': len(coherence_values),
            'mean': round(float(np.mean(coherence_array)), 4),
            'std': round(float(np.std(coherence_array)), 4),
            'min': round(float(np.min(coherence_array)), 4),
            'max': round(float(np.max(coherence_array)), 4),
            'median': round(float(np.median(coherence_array)), 4),
        },
        'go_no_go_analysis': {
            'threshold': round(E_MINUS_2, 4),
            'go_count': int(np.sum(coherence_array > E_MINUS_2)),
            'no_go_count': int(np.sum(coherence_array <= E_MINUS_2)),
            'go_fraction': round(
                float(np.mean(coherence_array > E_MINUS_2)), 3
            ),
        },
        'classification': {
            'excellent': int(np.sum(coherence_array >= 0.9)),
            'good': int(np.sum((coherence_array >= 0.7) & (coherence_array < 0.9))),
            'marginal': int(np.sum((coherence_array >= E_MINUS_2) & (coherence_array < 0.7))),
            'unreliable': int(np.sum(coherence_array < E_MINUS_2)),
        },
    }

    # Hardware validation if provided
    if hardware_validation:
        report['hardware_validation'] = {
            'backend': hardware_validation.get('backend', 'unknown'),
            'hardware_R': round(hardware_validation.get('hardware_R', 0), 4),
            'simulator_R': round(hardware_validation.get('simulator_R', 0), 4),
            'agreement': hardware_validation.get('agreement', 'unknown'),
        }

    return report


def export_report(
    report: Dict[str, Any],
    output_path: str,
    format: str = "json",
) -> str:
    """
    Export report to file.

    Args:
        report: Report dictionary
        output_path: Output file path
        format: Output format (json, markdown, html)

    Returns:
        Path to written file
    """
    output_path = Path(output_path)

    if format == "json":
        with open(output_path, 'w') as f:
            json.dump(report, f, indent=2, default=str)

    elif format == "markdown":
        md_content = _report_to_markdown(report)
        with open(output_path, 'w') as f:
            f.write(md_content)

    elif format == "html":
        html_content = _report_to_html(report)
        with open(output_path, 'w') as f:
            f.write(html_content)

    else:
        raise ValueError(f"Unknown format: {format}")

    logger.info(f"Report exported to {output_path}")
    return str(output_path)


def _report_to_markdown(report: Dict[str, Any]) -> str:
    """Convert report to Markdown format."""
    lines = []

    # Title
    metadata = report.get('metadata', {})
    lines.append(f"# {metadata.get('title', 'PhaseLab Report')}")
    lines.append("")
    lines.append(f"**Generated:** {metadata.get('generated_at', 'N/A')}")
    lines.append(f"**PhaseLab Version:** {metadata.get('phaselab_version', 'N/A')}")

    if metadata.get('target_gene'):
        lines.append(f"**Target Gene:** {metadata['target_gene']}")
    if metadata.get('disease'):
        lines.append(f"**Disease:** {metadata['disease']}")

    lines.append("")

    # Summary section
    if 'summary' in report:
        lines.append("## Summary")
        lines.append("")
        summary = report['summary']
        lines.append(f"- Total candidates: {summary.get('total_candidates_evaluated', 'N/A')}")
        lines.append(f"- GO guides: {summary.get('go_count', 'N/A')}")
        lines.append(f"- NO-GO guides: {summary.get('no_go_count', 'N/A')}")
        lines.append("")

    # Guides table
    if 'guides' in report and report['guides']:
        lines.append("## Top Guides")
        lines.append("")
        lines.append("| Rank | Sequence | Score | R̄ | GO/NO-GO |")
        lines.append("|------|----------|-------|-----|----------|")

        for guide in report['guides']:
            lines.append(
                f"| {guide['rank']} | `{guide['sequence']}` | "
                f"{guide['combined_score']:.3f} | {guide['coherence_R']:.3f} | "
                f"{guide['go_no_go']} |"
            )
        lines.append("")

    # Therapeutic window
    if 'therapeutic_window' in report:
        lines.append("## Therapeutic Window")
        lines.append("")
        tw = report['therapeutic_window']
        lines.append(f"- Baseline: {tw.get('baseline_expression', 'N/A')}")
        lines.append(f"- Therapeutic minimum: {tw.get('therapeutic_min', 'N/A')}")
        lines.append(f"- Therapeutic maximum: {tw.get('therapeutic_max', 'N/A')}")
        lines.append(f"- Optimal: {tw.get('optimal', 'N/A')}")
        lines.append("")

    # Coherence statistics
    if 'coherence_statistics' in report:
        lines.append("## Coherence Statistics")
        lines.append("")
        stats = report['coherence_statistics']
        lines.append(f"- Mean R̄: {stats.get('mean', 'N/A')}")
        lines.append(f"- Std R̄: {stats.get('std', 'N/A')}")
        lines.append(f"- Range: [{stats.get('min', 'N/A')}, {stats.get('max', 'N/A')}]")
        lines.append("")

    # Footer
    lines.append("---")
    lines.append("*Generated by PhaseLab - Informational Relativity Framework*")

    return "\n".join(lines)


def _report_to_html(report: Dict[str, Any]) -> str:
    """Convert report to HTML format."""
    metadata = report.get('metadata', {})

    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>{metadata.get('title', 'PhaseLab Report')}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        h1 {{ color: #2c3e50; }}
        h2 {{ color: #34495e; border-bottom: 1px solid #eee; padding-bottom: 5px; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #3498db; color: white; }}
        tr:nth-child(even) {{ background-color: #f2f2f2; }}
        .go {{ color: #27ae60; font-weight: bold; }}
        .no-go {{ color: #e74c3c; font-weight: bold; }}
        .sequence {{ font-family: monospace; background: #f8f9fa; padding: 2px 4px; }}
        .metadata {{ color: #7f8c8d; font-size: 0.9em; }}
    </style>
</head>
<body>
    <h1>{metadata.get('title', 'PhaseLab Report')}</h1>
    <p class="metadata">
        Generated: {metadata.get('generated_at', 'N/A')}<br>
        PhaseLab Version: {metadata.get('phaselab_version', 'N/A')}
    </p>
"""

    # Summary
    if 'summary' in report:
        summary = report['summary']
        html += f"""
    <h2>Summary</h2>
    <ul>
        <li>Total candidates evaluated: {summary.get('total_candidates_evaluated', 'N/A')}</li>
        <li>GO guides: <span class="go">{summary.get('go_count', 'N/A')}</span></li>
        <li>NO-GO guides: <span class="no-go">{summary.get('no_go_count', 'N/A')}</span></li>
    </ul>
"""

    # Guides table
    if 'guides' in report and report['guides']:
        html += """
    <h2>Top Guides</h2>
    <table>
        <tr>
            <th>Rank</th>
            <th>Sequence</th>
            <th>Combined Score</th>
            <th>R̄</th>
            <th>GO/NO-GO</th>
            <th>GC%</th>
            <th>ΔG</th>
        </tr>
"""
        for guide in report['guides']:
            go_class = 'go' if guide['go_no_go'] == 'GO' else 'no-go'
            html += f"""
        <tr>
            <td>{guide['rank']}</td>
            <td class="sequence">{guide['sequence']}</td>
            <td>{guide['combined_score']:.3f}</td>
            <td>{guide['coherence_R']:.4f}</td>
            <td class="{go_class}">{guide['go_no_go']}</td>
            <td>{guide['gc_content']:.1%}</td>
            <td>{guide['delta_g']:.1f}</td>
        </tr>
"""
        html += "    </table>\n"

    # Coherence stats
    if 'coherence_statistics' in report:
        stats = report['coherence_statistics']
        html += f"""
    <h2>Coherence Statistics</h2>
    <ul>
        <li>Mean R̄: {stats.get('mean', 'N/A')}</li>
        <li>Standard Deviation: {stats.get('std', 'N/A')}</li>
        <li>Range: [{stats.get('min', 'N/A')}, {stats.get('max', 'N/A')}]</li>
    </ul>
"""

    # Footer
    html += """
    <hr>
    <p class="metadata"><em>Generated by PhaseLab - Informational Relativity Framework</em></p>
</body>
</html>
"""

    return html


def generate_crispor_link(
    sequence: str,
    genome: str = "hg38",
) -> str:
    """
    Generate CRISPOR compatibility link for guide validation.

    Args:
        sequence: Target sequence (with context)
        genome: Genome assembly

    Returns:
        CRISPOR URL for the sequence
    """
    import urllib.parse

    base_url = "http://crispor.tefor.net/crispor.py"
    params = {
        'batchId': '',
        'pos': '',
        'org': genome,
        'seq': sequence[:500],  # CRISPOR has length limit
        'pam': 'NGG',
    }

    query = urllib.parse.urlencode(params)
    return f"{base_url}?{query}"
