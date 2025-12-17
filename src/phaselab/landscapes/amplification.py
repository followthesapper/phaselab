"""
PhaseLab Landscapes: Amplification detection.

Detects MYC-like amplifying behavior where perturbations amplify variance
rather than producing predictable effects. These regions should be AVOIDED
for therapeutic targeting.

The key finding from E213:
    MYC showed POSITIVE coherence-variance correlation (+0.45)
    while CD69, IL2RA, GATA1 showed NEGATIVE correlation.
    This is the hallmark of super-enhancer/amplifying regulatory architecture.
"""

import numpy as np
from scipy import stats
from typing import Optional, List, Tuple, Dict, Any
from dataclasses import dataclass
from enum import Enum

from .core import ResponseLandscape, CoherenceProfile
from .coherence import compute_spatial_coherence


class AmplificationType(Enum):
    """Types of amplifying behavior detected."""
    NONE = "none"
    GLOBAL = "global"  # Entire landscape amplifies
    LOCAL = "local"  # Specific regions amplify
    THRESHOLD = "threshold"  # Amplification above certain response levels


@dataclass
class AmplificationResult:
    """
    Result of amplification detection analysis.

    Attributes:
        is_amplifying: Whether significant amplification was detected.
        amplification_type: Type of amplification (NONE, GLOBAL, LOCAL, THRESHOLD).
        global_correlation: Coherence-variance correlation for entire landscape.
        amplifying_regions: List of (start, end, score) for local amplifying regions.
        risk_level: Overall risk assessment (LOW, MEDIUM, HIGH, CRITICAL).
        recommendation: Human-readable recommendation.
    """
    is_amplifying: bool
    amplification_type: AmplificationType
    global_correlation: float
    global_p_value: float
    amplifying_regions: List[Tuple[float, float, float]]
    risk_level: str
    recommendation: str

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            'is_amplifying': self.is_amplifying,
            'amplification_type': self.amplification_type.value,
            'global_correlation': self.global_correlation,
            'global_p_value': self.global_p_value,
            'amplifying_regions': [
                {'start': s, 'end': e, 'score': sc}
                for s, e, sc in self.amplifying_regions
            ],
            'risk_level': self.risk_level,
            'recommendation': self.recommendation,
        }


class AmplificationDetector:
    """
    Detector for MYC-like amplifying regulatory behavior.

    This class provides comprehensive amplification detection
    including global, local, and threshold-based analysis.

    Example:
        >>> detector = AmplificationDetector(window=50)
        >>> result = detector.analyze(landscape)
        >>> if result.is_amplifying:
        ...     print(f"WARNING: {result.recommendation}")
        ...     for start, end, score in result.amplifying_regions:
        ...         print(f"  Amplifying region: {start}-{end}")
    """

    def __init__(
        self,
        window: int = 50,
        global_threshold: float = 0.2,
        local_threshold: float = 0.3,
        significance_level: float = 0.05,
    ):
        """
        Initialize the amplification detector.

        Args:
            window: Window size for coherence computation.
            global_threshold: Positive correlation threshold for global amplification.
            local_threshold: Positive correlation threshold for local amplification.
            significance_level: P-value threshold for significance.
        """
        self.window = window
        self.global_threshold = global_threshold
        self.local_threshold = local_threshold
        self.significance_level = significance_level

    def analyze(
        self,
        landscape: ResponseLandscape,
        profile: Optional[CoherenceProfile] = None,
    ) -> AmplificationResult:
        """
        Perform comprehensive amplification analysis.

        Args:
            landscape: Input ResponseLandscape.
            profile: Pre-computed CoherenceProfile (computed if not provided).

        Returns:
            AmplificationResult with full analysis.
        """
        # Compute profile if needed
        if profile is None:
            profile = compute_spatial_coherence(landscape, window=self.window)

        # Global analysis
        global_corr = profile.correlation
        global_p = profile.p_value
        is_global_amplifying = (
            global_corr > self.global_threshold and
            global_p < self.significance_level
        )

        # Local analysis
        local_regions = self._detect_local_amplification(profile)

        # Threshold analysis
        threshold_regions = self._detect_threshold_amplification(landscape, profile)

        # Combine all amplifying regions
        all_regions = local_regions + threshold_regions

        # Determine amplification type
        if is_global_amplifying:
            amp_type = AmplificationType.GLOBAL
        elif local_regions:
            amp_type = AmplificationType.LOCAL
        elif threshold_regions:
            amp_type = AmplificationType.THRESHOLD
        else:
            amp_type = AmplificationType.NONE

        # Risk assessment
        risk_level = self._assess_risk(global_corr, len(all_regions), is_global_amplifying)

        # Generate recommendation
        recommendation = self._generate_recommendation(
            amp_type, risk_level, global_corr, all_regions
        )

        return AmplificationResult(
            is_amplifying=amp_type != AmplificationType.NONE,
            amplification_type=amp_type,
            global_correlation=global_corr,
            global_p_value=global_p,
            amplifying_regions=all_regions,
            risk_level=risk_level,
            recommendation=recommendation,
        )

    def _detect_local_amplification(
        self,
        profile: CoherenceProfile,
    ) -> List[Tuple[float, float, float]]:
        """Detect locally amplifying regions."""
        regions = []

        # Sliding window for local correlation
        local_window = max(20, len(profile.coords) // 10)

        for i in range(0, len(profile.coords) - local_window, local_window // 2):
            local_coh = profile.coherence[i:i + local_window]
            local_var = profile.local_variance[i:i + local_window]

            if len(local_coh) < 10:
                continue

            # Remove NaN
            mask = ~(np.isnan(local_coh) | np.isnan(local_var))
            if np.sum(mask) < 10:
                continue

            local_corr, local_p = stats.pearsonr(local_coh[mask], local_var[mask])

            if local_corr > self.local_threshold and local_p < 0.1:
                start = profile.coords[i]
                end_idx = min(i + local_window - 1, len(profile.coords) - 1)
                end = profile.coords[end_idx]
                regions.append((float(start), float(end), float(local_corr)))

        return regions

    def _detect_threshold_amplification(
        self,
        landscape: ResponseLandscape,
        profile: CoherenceProfile,
    ) -> List[Tuple[float, float, float]]:
        """
        Detect amplification above response thresholds.

        Some systems only amplify when responses exceed certain levels.
        """
        regions = []

        # Split by response magnitude
        responses = landscape.mean_response
        high_response_mask = np.abs(responses) > np.percentile(np.abs(responses), 75)

        # Map back to profile coordinates
        if len(profile.coords) < 20:
            return []

        # Find coordinates with high responses
        landscape_sorted = landscape.sort_by_coords()
        coords = landscape_sorted.coords if landscape_sorted.coords.ndim == 1 else landscape_sorted.coords[:, 0]

        # Check correlation in high-response regions only
        high_indices = []
        for i, coord in enumerate(profile.coords):
            # Find nearest landscape coordinate
            nearest_idx = np.argmin(np.abs(coords - coord))
            if high_response_mask[nearest_idx]:
                high_indices.append(i)

        if len(high_indices) > 15:
            high_coh = profile.coherence[high_indices]
            high_var = profile.local_variance[high_indices]

            mask = ~(np.isnan(high_coh) | np.isnan(high_var))
            if np.sum(mask) > 10:
                corr, p = stats.pearsonr(high_coh[mask], high_var[mask])

                if corr > self.local_threshold and p < 0.1:
                    # Mark high-response regions as potentially amplifying
                    start = profile.coords[min(high_indices)]
                    end = profile.coords[max(high_indices)]
                    regions.append((float(start), float(end), float(corr)))

        return regions

    def _assess_risk(
        self,
        global_corr: float,
        n_regions: int,
        is_global: bool,
    ) -> str:
        """Assess overall risk level."""
        if is_global and global_corr > 0.4:
            return "CRITICAL"
        elif is_global:
            return "HIGH"
        elif n_regions > 2:
            return "MEDIUM"
        elif n_regions > 0:
            return "LOW"
        else:
            return "MINIMAL"

    def _generate_recommendation(
        self,
        amp_type: AmplificationType,
        risk_level: str,
        global_corr: float,
        regions: List[Tuple[float, float, float]],
    ) -> str:
        """Generate human-readable recommendation."""
        if amp_type == AmplificationType.NONE:
            return (
                "No significant amplification detected. "
                "Standard spatial coherence analysis is valid for this landscape."
            )

        if amp_type == AmplificationType.GLOBAL:
            return (
                f"CRITICAL: Global amplifying behavior detected (r=+{global_corr:.2f}). "
                f"This landscape resembles MYC-like super-enhancer architecture. "
                f"Perturbations will likely amplify variance rather than produce "
                f"predictable effects. STRONGLY recommend alternative targeting strategy."
            )

        if amp_type == AmplificationType.LOCAL:
            region_str = ", ".join(f"[{s:.0f}-{e:.0f}]" for s, e, _ in regions[:3])
            return (
                f"Local amplifying regions detected: {region_str}. "
                f"Avoid these specific regions. Stable regions elsewhere "
                f"may still be valid targets."
            )

        if amp_type == AmplificationType.THRESHOLD:
            return (
                f"Threshold-dependent amplification: high-response regions show "
                f"amplifying behavior. Consider titrating perturbation strength "
                f"to stay below amplification threshold."
            )

        return "Unknown amplification pattern. Exercise caution."


def amplification_score(
    landscape: ResponseLandscape,
    window: int = 50,
) -> float:
    """
    Compute a single amplification score for the landscape.

    Higher scores indicate more amplifying behavior.
    Score > 0.2 suggests MYC-like architecture.

    Args:
        landscape: Input ResponseLandscape.
        window: Window size for coherence computation.

    Returns:
        Amplification score (coherence-variance correlation).
    """
    profile = compute_spatial_coherence(landscape, window=window)
    return max(0, profile.correlation)  # Only positive correlations count


def is_amplifying(
    landscape: ResponseLandscape,
    window: int = 50,
    threshold: float = 0.2,
) -> bool:
    """
    Quick check for amplifying behavior.

    Args:
        landscape: Input ResponseLandscape.
        window: Window size for coherence computation.
        threshold: Positive correlation threshold.

    Returns:
        True if landscape shows amplifying behavior.
    """
    profile = compute_spatial_coherence(landscape, window=window)
    return profile.correlation > threshold and profile.p_value < 0.05
