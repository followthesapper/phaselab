"""
PhaseLab Constants: Universal constants from Informational Relativity framework.
"""

import numpy as np

# The e^-2 threshold: fundamental GO/NO-GO boundary
# Below this, quantum simulations are unreliable
E_MINUS_2 = np.exp(-2)  # ≈ 0.1353

# 4π² universal constant (appears in IR spectral analysis)
FOUR_PI_SQUARED = 4 * np.pi ** 2  # ≈ 39.478

# Phase thresholds for classification
COHERENCE_THRESHOLDS = {
    "excellent": 0.9,
    "good": 0.7,
    "marginal": E_MINUS_2,
    "unreliable": 0.0,
}

# Circadian model constants
CIRCADIAN_PERIOD_HOURS = 24.0
OMEGA_CIRCADIAN = 2 * np.pi / CIRCADIAN_PERIOD_HOURS  # rad/hour
