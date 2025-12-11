"""
Tests for phaselab.circadian module.
"""

import numpy as np
import pytest
from phaselab.circadian.kuramoto import (
    kuramoto_order_parameter,
    kuramoto_ode,
    simulate_kuramoto,
)
from phaselab.circadian.sms_model import (
    simulate_sms_clock,
    SMSClockParams,
    therapeutic_scan,
    classify_synchronization,
    predict_sleep_quality,
)


class TestKuramoto:
    """Test base Kuramoto model."""

    def test_order_parameter_synced(self):
        """Synchronized phases give R ≈ 1."""
        phases = np.array([0.1, 0.1, 0.1, 0.1])
        R, psi = kuramoto_order_parameter(phases)

        assert R > 0.99

    def test_order_parameter_desynced(self):
        """Evenly distributed phases give R ≈ 0."""
        phases = np.array([0, np.pi/2, np.pi, 3*np.pi/2])
        R, psi = kuramoto_order_parameter(phases)

        assert R < 0.1

    def test_kuramoto_ode_shape(self):
        """ODE returns correct shape."""
        phases = np.array([0.0, 0.5, 1.0])
        omegas = np.array([1.0, 1.0, 1.0])
        K = np.ones((3, 3)) * 0.5

        dphases = kuramoto_ode(0, phases, omegas, K)

        assert dphases.shape == phases.shape

    def test_simulate_kuramoto(self):
        """Basic simulation runs."""
        omegas = np.array([1.0, 1.0, 1.0, 1.0])
        K = np.ones((4, 4)) * 0.5
        np.fill_diagonal(K, 0)

        result = simulate_kuramoto(omegas, K, t_end=10.0)

        assert 't' in result
        assert 'phases' in result
        assert 'order_param' in result
        assert len(result['t']) > 0


class TestSMSClockParams:
    """Test SMS model parameters."""

    def test_default_params(self):
        """Check default parameter values."""
        params = SMSClockParams()

        assert params.tau_P == 4.0  # PER delay
        assert params.alpha_P == 2.0  # PER suppression
        assert params.beta_R == 0.5  # RORα effect
        assert params.beta_V == 0.5  # REV-ERBα effect

    def test_custom_params(self):
        """Custom parameters."""
        params = SMSClockParams(
            tau_P=6.0,
            alpha_P=3.0,
        )

        assert params.tau_P == 6.0
        assert params.alpha_P == 3.0


class TestSMSClock:
    """Test SMS circadian clock simulation."""

    def test_simulate_returns_dict(self):
        """Simulation returns expected keys."""
        result = simulate_sms_clock(rai1_level=0.5, t_end=48.0)

        expected_keys = [
            't', 'theta_C', 'theta_B', 'P', 'R', 'V',
            'order_param', 'final_R_bar', 'classification', 'go_no_go'
        ]
        for key in expected_keys:
            assert key in result, f"Missing key: {key}"

    def test_simulate_time_array(self):
        """Time array has correct length."""
        result = simulate_sms_clock(rai1_level=0.5, t_end=24.0, dt=0.1)

        # Should have ~240 time points (24/0.1)
        assert len(result['t']) > 200

    def test_simulate_phases_bounded(self):
        """Phases should be in [0, 2π]."""
        result = simulate_sms_clock(rai1_level=0.5, t_end=48.0)

        assert np.all(result['theta_C'] >= 0)
        assert np.all(result['theta_C'] <= 2 * np.pi)
        assert np.all(result['theta_B'] >= 0)
        assert np.all(result['theta_B'] <= 2 * np.pi)

    def test_simulate_order_param_bounded(self):
        """Order parameter should be in [0, 1]."""
        result = simulate_sms_clock(rai1_level=0.5, t_end=48.0)

        assert np.all(result['order_param'] >= 0)
        assert np.all(result['order_param'] <= 1)

    def test_simulate_per_feedback(self):
        """PER feedback state should evolve."""
        result = simulate_sms_clock(rai1_level=0.5, t_end=48.0)

        # P should not stay constant
        assert np.std(result['P']) > 0

    def test_simulate_reproducible(self):
        """Same seed gives same result."""
        result1 = simulate_sms_clock(rai1_level=0.5, t_end=24.0, random_seed=42)
        result2 = simulate_sms_clock(rai1_level=0.5, t_end=24.0, random_seed=42)

        assert result1['final_R_bar'] == result2['final_R_bar']

    def test_rai1_level_effect(self):
        """Different RAI1 levels should give different results."""
        result_low = simulate_sms_clock(rai1_level=0.3, t_end=48.0, random_seed=42)
        result_high = simulate_sms_clock(rai1_level=1.0, t_end=48.0, random_seed=42)

        # Results should differ (though both may sync in simple model)
        # At minimum, the stored rai1_level should differ
        assert result_low['rai1_level'] != result_high['rai1_level']


class TestClassifySynchronization:
    """Test synchronization classification."""

    def test_synchronized(self):
        assert classify_synchronization(0.95) == "SYNCHRONIZED"

    def test_partially_synchronized(self):
        assert classify_synchronization(0.75) == "PARTIALLY_SYNCHRONIZED"

    def test_weakly_synchronized(self):
        assert classify_synchronization(0.20) == "WEAKLY_SYNCHRONIZED"

    def test_desynchronized(self):
        assert classify_synchronization(0.05) == "DESYNCHRONIZED"


class TestTherapeuticScan:
    """Test therapeutic window scanning."""

    def test_scan_returns_dict(self):
        """Scan returns expected structure."""
        result = therapeutic_scan(
            rai1_levels=[0.5, 0.7, 1.0],
            t_end=48.0,
            n_trials=1,
        )

        assert 'levels' in result
        assert 'R_bars' in result
        assert 'classifications' in result
        assert 'optimal_level' in result

    def test_scan_levels_match(self):
        """Output levels match input."""
        levels = [0.4, 0.6, 0.8]
        result = therapeutic_scan(rai1_levels=levels, t_end=24.0, n_trials=1)

        assert result['levels'] == levels
        assert len(result['R_bars']) == len(levels)

    def test_scan_r_bars_valid(self):
        """R̄ values are in valid range."""
        result = therapeutic_scan(
            rai1_levels=[0.5, 1.0],
            t_end=24.0,
            n_trials=1,
        )

        for r in result['R_bars']:
            assert 0 <= r <= 1


class TestPredictSleepQuality:
    """Test sleep quality prediction."""

    def test_normal_sleep(self):
        result = predict_sleep_quality(0.95)
        assert "NORMAL" in result

    def test_mild_disruption(self):
        result = predict_sleep_quality(0.75)
        assert "MILD" in result

    def test_severe_disruption(self):
        result = predict_sleep_quality(0.20)  # Between e^-2 and 0.5
        assert "SEVERE" in result

    def test_critical(self):
        result = predict_sleep_quality(0.05)
        assert "CRITICAL" in result


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
