"""
Tests for phaselab.core module.
"""

import numpy as np
import pytest
from phaselab.core.coherence import (
    coherence_score,
    go_no_go,
    phase_variance,
    classify_coherence,
    compare_sim_hardware,
)
from phaselab.core.constants import E_MINUS_2, FOUR_PI_SQUARED
from phaselab.core.hamiltonians import (
    PauliHamiltonian,
    build_pauli_hamiltonian,
    build_grna_hamiltonian,
)


class TestConstants:
    """Test IR constants."""

    def test_e_minus_2(self):
        """e^-2 should be approximately 0.1353."""
        assert abs(E_MINUS_2 - 0.1353352832366127) < 1e-10

    def test_four_pi_squared(self):
        """4π² should be approximately 39.478."""
        assert abs(FOUR_PI_SQUARED - 39.4784176043574) < 1e-10


class TestCoherence:
    """Test coherence score calculations."""

    def test_perfect_sync_phases(self):
        """Identical phases should give R̄ = 1."""
        phases = np.array([0.5, 0.5, 0.5, 0.5])
        R_bar = coherence_score(phases, mode='phases')
        assert abs(R_bar - 1.0) < 1e-10

    def test_opposite_phases(self):
        """Opposite phases should give R̄ = 0."""
        phases = np.array([0.0, np.pi, 0.0, np.pi])
        R_bar = coherence_score(phases, mode='phases')
        assert R_bar < 0.01

    def test_random_phases(self):
        """Random phases should give intermediate R̄."""
        np.random.seed(42)
        phases = np.random.uniform(0, 2 * np.pi, 100)
        R_bar = coherence_score(phases, mode='phases')
        assert 0 < R_bar < 1

    def test_variance_mode(self):
        """R̄ = exp(-V_φ/2) relationship."""
        V_phi = 1.0
        R_bar = coherence_score(V_phi, mode='variance')
        expected = np.exp(-0.5)
        assert abs(R_bar - expected) < 1e-10

    def test_zero_variance(self):
        """Zero variance should give R̄ = 1."""
        R_bar = coherence_score(0.0, mode='variance')
        assert abs(R_bar - 1.0) < 1e-10


class TestGoNoGo:
    """Test GO/NO-GO classification."""

    def test_go_above_threshold(self):
        """R̄ > e^-2 should be GO."""
        assert go_no_go(0.5) == "GO"
        assert go_no_go(0.9) == "GO"
        assert go_no_go(0.14) == "GO"

    def test_no_go_below_threshold(self):
        """R̄ < e^-2 should be NO-GO."""
        assert go_no_go(0.1) == "NO-GO"
        assert go_no_go(0.05) == "NO-GO"
        assert go_no_go(0.0) == "NO-GO"

    def test_none_is_no_go(self):
        """None should be NO-GO."""
        assert go_no_go(None) == "NO-GO"

    def test_boundary(self):
        """Test boundary value."""
        # Just above threshold
        assert go_no_go(E_MINUS_2 + 0.001) == "GO"
        # Just below threshold
        assert go_no_go(E_MINUS_2 - 0.001) == "NO-GO"


class TestClassifyCoherence:
    """Test coherence classification."""

    def test_excellent(self):
        assert classify_coherence(0.95) == "EXCELLENT"

    def test_good(self):
        assert classify_coherence(0.75) == "GOOD"

    def test_marginal(self):
        assert classify_coherence(0.2) == "MARGINAL"

    def test_unreliable(self):
        assert classify_coherence(0.1) == "UNRELIABLE"


class TestPhaseVariance:
    """Test phase variance calculation."""

    def test_zero_variance_synced(self):
        """Synchronized phases should have low variance."""
        phases = np.array([1.0, 1.0, 1.0])
        V_phi = phase_variance(phases)
        assert V_phi < 0.01

    def test_high_variance_desynced(self):
        """Desynchronized phases should have high variance."""
        phases = np.array([0, np.pi/2, np.pi, 3*np.pi/2])
        V_phi = phase_variance(phases)
        assert V_phi > 1.0


class TestCompareSimHardware:
    """Test simulator-hardware comparison."""

    def test_excellent_agreement(self):
        diff, agreement = compare_sim_hardware(0.85, 0.86)
        assert agreement == "EXCELLENT"
        assert diff < 0.02

    def test_good_agreement(self):
        diff, agreement = compare_sim_hardware(0.85, 0.90)
        assert agreement == "GOOD"

    def test_poor_agreement(self):
        diff, agreement = compare_sim_hardware(0.85, 0.50)
        assert agreement == "POOR"


class TestHamiltonians:
    """Test Hamiltonian builders."""

    def test_pauli_hamiltonian_creation(self):
        """Test basic Hamiltonian creation."""
        H = PauliHamiltonian(4)
        H.add_term(1.0, "ZZII")
        H.add_term(-0.5, "XXII")

        terms = H.get_terms()
        assert len(terms) == 2
        assert terms[0] == (1.0, "ZZII")

    def test_pauli_hamiltonian_wrong_length(self):
        """Wrong Pauli string length should raise error."""
        H = PauliHamiltonian(4)
        with pytest.raises(ValueError):
            H.add_term(1.0, "ZZI")  # Only 3 chars, need 4

    def test_build_pauli_hamiltonian(self):
        """Test generic Hamiltonian builder."""
        H = build_pauli_hamiltonian(
            n_qubits=4,
            zz_terms=[(0, 1, 0.5)],
            z_terms=[(2, -0.3)],
        )
        terms = H.get_terms()
        assert len(terms) == 2

    def test_build_grna_hamiltonian(self):
        """Test gRNA Hamiltonian builder."""
        guide = "ATCGATCGATCGATCGATCG"  # 20bp
        H = build_grna_hamiltonian(guide)

        assert H.n_qubits == 20
        assert len(H.terms) > 0

    def test_grna_hamiltonian_gc_effect(self):
        """Higher GC should affect Hamiltonian."""
        low_gc = "AAAAAAAAAAAAAAAAAAAA"
        high_gc = "GCGCGCGCGCGCGCGCGCGC"

        H_low = build_grna_hamiltonian(low_gc)
        H_high = build_grna_hamiltonian(high_gc)

        # They should be different
        assert len(H_low.terms) != len(H_high.terms) or \
               H_low.terms[0].coefficient != H_high.terms[0].coefficient


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
