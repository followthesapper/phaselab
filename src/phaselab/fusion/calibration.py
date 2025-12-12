"""
Score calibration for evidence fusion.

Provides calibration curves to convert raw predictor scores
to well-calibrated probabilities.
"""

from dataclasses import dataclass, field
from typing import Optional, Dict, List, Tuple, Callable
import logging

import numpy as np
from scipy import interpolate

logger = logging.getLogger(__name__)


@dataclass
class CalibrationCurve:
    """
    Calibration curve for a single predictor.

    Attributes
    ----------
    name : str
        Predictor name
    raw_scores : np.ndarray
        Raw score bins
    calibrated_scores : np.ndarray
        Calibrated probability for each bin
    n_samples : int
        Number of samples used for calibration
    brier_score : float
        Brier score (lower is better)
    """
    name: str
    raw_scores: np.ndarray
    calibrated_scores: np.ndarray
    n_samples: int = 0
    brier_score: float = 0.0

    def __post_init__(self):
        # Create interpolation function
        self._interp = interpolate.interp1d(
            self.raw_scores,
            self.calibrated_scores,
            kind="linear",
            bounds_error=False,
            fill_value=(self.calibrated_scores[0], self.calibrated_scores[-1]),
        )

    def calibrate(self, raw_score: float) -> float:
        """
        Convert raw score to calibrated probability.

        Parameters
        ----------
        raw_score : float
            Raw predictor score

        Returns
        -------
        float
            Calibrated probability
        """
        return float(np.clip(self._interp(raw_score), 0.0, 1.0))

    def calibrate_batch(self, raw_scores: np.ndarray) -> np.ndarray:
        """Calibrate multiple scores."""
        return np.clip(self._interp(raw_scores), 0.0, 1.0)

    @classmethod
    def identity(cls, name: str) -> "CalibrationCurve":
        """Create identity (no-op) calibration curve."""
        scores = np.linspace(0, 1, 11)
        return cls(
            name=name,
            raw_scores=scores,
            calibrated_scores=scores,
        )

    @classmethod
    def from_bins(
        cls,
        name: str,
        bin_edges: List[float],
        observed_rates: List[float],
        n_samples: int = 0,
    ) -> "CalibrationCurve":
        """
        Create calibration curve from binned observations.

        Parameters
        ----------
        name : str
            Predictor name
        bin_edges : List[float]
            Edges of score bins (n+1 values)
        observed_rates : List[float]
            Observed positive rate in each bin (n values)
        n_samples : int
            Total samples used

        Returns
        -------
        CalibrationCurve
        """
        # Use bin centers as raw scores
        raw_scores = [
            (bin_edges[i] + bin_edges[i + 1]) / 2
            for i in range(len(bin_edges) - 1)
        ]

        # Add endpoints
        raw_scores = [0.0] + raw_scores + [1.0]
        calibrated = [observed_rates[0]] + list(observed_rates) + [observed_rates[-1]]

        return cls(
            name=name,
            raw_scores=np.array(raw_scores),
            calibrated_scores=np.array(calibrated),
            n_samples=n_samples,
        )


class Calibrator:
    """
    Calibrator for multiple predictors.

    Manages calibration curves for different predictors and
    provides calibrated score aggregation.

    Parameters
    ----------
    method : str
        Calibration method: "platt" (logistic), "isotonic", "beta"

    Examples
    --------
    >>> calibrator = Calibrator()
    >>>
    >>> # Fit calibration from validation data
    >>> calibrator.fit(
    ...     predictor_name="DeepCRISPR",
    ...     predicted_scores=predictions,
    ...     actual_outcomes=outcomes,
    ... )
    >>>
    >>> # Calibrate new predictions
    >>> calibrated = calibrator.calibrate("DeepCRISPR", raw_score=0.75)
    """

    # Default calibration curves based on published data
    DEFAULT_CURVES: Dict[str, Tuple[List[float], List[float]]] = {
        "DeepCRISPR": (
            [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            [0.05, 0.15, 0.35, 0.55, 0.75, 0.90],
        ),
        "DeepSpCas9": (
            [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            [0.08, 0.18, 0.38, 0.58, 0.78, 0.88],
        ),
        "RuleBased": (
            [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            [0.10, 0.25, 0.40, 0.55, 0.70, 0.80],
        ),
    }

    def __init__(self, method: str = "isotonic"):
        self.method = method
        self._curves: Dict[str, CalibrationCurve] = {}
        self._load_defaults()

    def _load_defaults(self) -> None:
        """Load default calibration curves."""
        for name, (raw, cal) in self.DEFAULT_CURVES.items():
            self._curves[name] = CalibrationCurve(
                name=name,
                raw_scores=np.array(raw),
                calibrated_scores=np.array(cal),
            )

    def fit(
        self,
        predictor_name: str,
        predicted_scores: np.ndarray,
        actual_outcomes: np.ndarray,
        n_bins: int = 10,
    ) -> CalibrationCurve:
        """
        Fit calibration curve from validation data.

        Parameters
        ----------
        predictor_name : str
            Name of predictor
        predicted_scores : np.ndarray
            Predicted scores
        actual_outcomes : np.ndarray
            Binary actual outcomes (0 or 1)
        n_bins : int
            Number of calibration bins

        Returns
        -------
        CalibrationCurve
            Fitted calibration curve
        """
        predicted_scores = np.asarray(predicted_scores)
        actual_outcomes = np.asarray(actual_outcomes)

        if len(predicted_scores) != len(actual_outcomes):
            raise ValueError("Predictions and outcomes must have same length")

        if self.method == "isotonic":
            curve = self._fit_isotonic(
                predictor_name, predicted_scores, actual_outcomes
            )
        elif self.method == "platt":
            curve = self._fit_platt(
                predictor_name, predicted_scores, actual_outcomes
            )
        else:
            # Default to binned calibration
            curve = self._fit_binned(
                predictor_name, predicted_scores, actual_outcomes, n_bins
            )

        self._curves[predictor_name] = curve
        return curve

    def _fit_binned(
        self,
        name: str,
        scores: np.ndarray,
        outcomes: np.ndarray,
        n_bins: int,
    ) -> CalibrationCurve:
        """Fit binned calibration curve."""
        bin_edges = np.linspace(0, 1, n_bins + 1)
        observed_rates = []

        for i in range(n_bins):
            mask = (scores >= bin_edges[i]) & (scores < bin_edges[i + 1])
            if mask.sum() > 0:
                rate = outcomes[mask].mean()
            else:
                # Use midpoint if no samples
                rate = (bin_edges[i] + bin_edges[i + 1]) / 2
            observed_rates.append(rate)

        # Calculate Brier score
        calibrated = self._apply_binned(scores, bin_edges, observed_rates)
        brier = float(np.mean((calibrated - outcomes) ** 2))

        return CalibrationCurve.from_bins(
            name=name,
            bin_edges=list(bin_edges),
            observed_rates=observed_rates,
            n_samples=len(scores),
        )

    def _apply_binned(
        self,
        scores: np.ndarray,
        bin_edges: np.ndarray,
        rates: List[float],
    ) -> np.ndarray:
        """Apply binned calibration."""
        result = np.zeros_like(scores)
        for i, rate in enumerate(rates):
            mask = (scores >= bin_edges[i]) & (scores < bin_edges[i + 1])
            result[mask] = rate
        return result

    def _fit_isotonic(
        self,
        name: str,
        scores: np.ndarray,
        outcomes: np.ndarray,
    ) -> CalibrationCurve:
        """Fit isotonic regression calibration."""
        try:
            from sklearn.isotonic import IsotonicRegression
            ir = IsotonicRegression(out_of_bounds="clip")
            ir.fit(scores, outcomes)

            # Create curve from isotonic transform
            x = np.linspace(0, 1, 101)
            y = ir.predict(x)

            return CalibrationCurve(
                name=name,
                raw_scores=x,
                calibrated_scores=y,
                n_samples=len(scores),
            )

        except ImportError:
            logger.warning("sklearn not available, using binned calibration")
            return self._fit_binned(name, scores, outcomes, 10)

    def _fit_platt(
        self,
        name: str,
        scores: np.ndarray,
        outcomes: np.ndarray,
    ) -> CalibrationCurve:
        """Fit Platt scaling (logistic calibration)."""
        try:
            from sklearn.linear_model import LogisticRegression
            lr = LogisticRegression()
            lr.fit(scores.reshape(-1, 1), outcomes)

            x = np.linspace(0, 1, 101)
            y = lr.predict_proba(x.reshape(-1, 1))[:, 1]

            return CalibrationCurve(
                name=name,
                raw_scores=x,
                calibrated_scores=y,
                n_samples=len(scores),
            )

        except ImportError:
            logger.warning("sklearn not available, using binned calibration")
            return self._fit_binned(name, scores, outcomes, 10)

    def calibrate(
        self,
        predictor_name: str,
        raw_score: float,
    ) -> float:
        """
        Calibrate a single score.

        Parameters
        ----------
        predictor_name : str
            Name of predictor
        raw_score : float
            Raw prediction score

        Returns
        -------
        float
            Calibrated probability
        """
        if predictor_name not in self._curves:
            # Use identity if no curve
            return raw_score

        return self._curves[predictor_name].calibrate(raw_score)

    def calibrate_batch(
        self,
        predictor_name: str,
        raw_scores: np.ndarray,
    ) -> np.ndarray:
        """Calibrate multiple scores."""
        if predictor_name not in self._curves:
            return raw_scores

        return self._curves[predictor_name].calibrate_batch(raw_scores)

    def get_curve(self, predictor_name: str) -> Optional[CalibrationCurve]:
        """Get calibration curve for predictor."""
        return self._curves.get(predictor_name)

    @property
    def available_predictors(self) -> List[str]:
        """List of predictors with calibration curves."""
        return list(self._curves.keys())

    def reliability_diagram(
        self,
        predictor_name: str,
        predicted_scores: np.ndarray,
        actual_outcomes: np.ndarray,
        n_bins: int = 10,
    ) -> Dict[str, np.ndarray]:
        """
        Generate data for reliability diagram.

        Parameters
        ----------
        predictor_name : str
            Predictor name
        predicted_scores : np.ndarray
            Predicted scores
        actual_outcomes : np.ndarray
            Actual binary outcomes
        n_bins : int
            Number of bins

        Returns
        -------
        Dict[str, np.ndarray]
            Data for plotting: bin_centers, observed_rates, predicted_rates, counts
        """
        bin_edges = np.linspace(0, 1, n_bins + 1)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        observed_rates = []
        predicted_rates = []
        counts = []

        for i in range(n_bins):
            mask = (predicted_scores >= bin_edges[i]) & (predicted_scores < bin_edges[i + 1])
            count = mask.sum()
            counts.append(count)

            if count > 0:
                observed_rates.append(actual_outcomes[mask].mean())
                predicted_rates.append(predicted_scores[mask].mean())
            else:
                observed_rates.append(np.nan)
                predicted_rates.append(np.nan)

        return {
            "bin_centers": bin_centers,
            "observed_rates": np.array(observed_rates),
            "predicted_rates": np.array(predicted_rates),
            "counts": np.array(counts),
        }
