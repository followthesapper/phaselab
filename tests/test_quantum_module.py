"""
Tests for PhaseLab v0.6.0 quantum module (ATLAS-Q integration).
"""

import numpy as np
import pytest

from phaselab.quantum import is_atlas_q_available, get_atlas_q_version


class TestQuantumModuleAvailability:
    """Test basic module availability."""

    def test_is_atlas_q_available_returns_bool(self):
        """is_atlas_q_available should return a boolean."""
        result = is_atlas_q_available()
        assert isinstance(result, bool)

    def test_get_atlas_q_version_returns_string(self):
        """get_atlas_q_version should return a string."""
        result = get_atlas_q_version()
        assert isinstance(result, str)

    def test_atlas_q_version_format(self):
        """If ATLAS-Q is available, version should be semver-like."""
        if is_atlas_q_available():
            version = get_atlas_q_version()
            # Should contain at least one dot (e.g., "0.6.3")
            assert "." in version


class TestCoherenceModule:
    """Test quantum coherence calculations."""

    def test_compute_coherence_from_expectations_basic(self):
        """Test coherence from expectation values."""
        from phaselab.quantum.coherence import compute_coherence_from_expectations

        # High coherence case (expectations near 1)
        expectations = np.array([0.9, 0.85, 0.88, 0.92, 0.87])
        result = compute_coherence_from_expectations(expectations)

        assert 0 <= result.R_bar <= 1
        assert result.V_phi >= 0
        assert result.n_measurements == 5
        assert bool(result.is_go) is True  # Should be above e^-2

    def test_compute_coherence_from_expectations_low_coherence(self):
        """Test lower coherence case with varied expectations."""
        from phaselab.quantum.coherence import compute_coherence_from_expectations

        # Expectations spanning a wider range = lower coherence than clustered values
        # Note: uniform [-1,1] maps via arccos to [0, pi], giving moderate coherence
        # For truly low coherence, we need expectations that map to diverse phases
        expectations = np.array([0.9, -0.9, 0.5, -0.5, 0.0, 0.3, -0.3])
        result = compute_coherence_from_expectations(expectations)

        assert 0 <= result.R_bar <= 1
        # Varied expectations should have lower coherence than clustered ones
        assert result.R_bar < 0.9  # Not perfectly coherent

    def test_compute_coherence_empty_input(self):
        """Test handling of empty input."""
        from phaselab.quantum.coherence import compute_coherence_from_expectations

        result = compute_coherence_from_expectations(np.array([]))
        assert result.R_bar == 0.0
        assert result.is_go is False

    def test_coherence_native_vs_atlas(self):
        """Test that native and ATLAS-Q give same results."""
        from phaselab.quantum.coherence import compute_coherence_from_expectations

        expectations = np.array([0.8, 0.75, 0.82, 0.78, 0.85])

        native = compute_coherence_from_expectations(expectations, use_atlas_q=False)

        if is_atlas_q_available():
            atlas = compute_coherence_from_expectations(expectations, use_atlas_q=True)
            # Results should be very close
            assert abs(native.R_bar - atlas.R_bar) < 0.01

    def test_compute_coherence_from_phases(self):
        """Test coherence from phase angles."""
        from phaselab.quantum.coherence import compute_coherence_from_phases

        # All phases near 0 = high coherence
        phases = np.array([0.1, 0.05, 0.08, 0.12, 0.07])
        result = compute_coherence_from_phases(phases)

        assert result.R_bar > 0.95

    def test_compare_coherence(self):
        """Test coherence comparison function."""
        from phaselab.quantum.coherence import (
            compute_coherence_from_expectations,
            compare_coherence,
        )

        exp1 = np.array([0.9, 0.85, 0.88])
        exp2 = np.array([0.88, 0.84, 0.87])

        result1 = compute_coherence_from_expectations(exp1)
        result2 = compute_coherence_from_expectations(exp2)

        diff, agreement = compare_coherence(result1, result2)
        assert diff >= 0
        assert agreement in ["EXCELLENT", "GOOD", "MODERATE", "POOR"]


class TestGroupingModule:
    """Test Hamiltonian grouping."""

    def test_pauli_commutes(self):
        """Test Pauli commutativity check."""
        from phaselab.quantum.grouping import pauli_commutes

        # Same Paulis commute
        assert pauli_commutes("XX", "XX") is True
        assert pauli_commutes("ZZ", "ZZ") is True

        # Identity commutes with everything
        assert pauli_commutes("XI", "IX") is True
        assert pauli_commutes("II", "ZZ") is True

        # XY and YX anti-commute at both positions (2 positions = even, so they commute!)
        # For non-commuting, need odd number of anti-commuting positions
        assert pauli_commutes("XY", "YX") is True  # Anti-commute at both → even → commute
        assert pauli_commutes("XZ", "ZI") is False  # Anti-commute at position 0 only → odd → don't commute

        # Check multi-qubit
        assert pauli_commutes("ZZII", "IZZI") is True

    def test_ir_hamiltonian_grouping_basic(self):
        """Test basic grouping functionality."""
        from phaselab.quantum.grouping import ir_hamiltonian_grouping

        coefficients = np.array([1.0, -0.5, 0.3, -0.2])
        pauli_strings = ["ZZII", "IZZI", "XXII", "IXXI"]

        result = ir_hamiltonian_grouping(
            coefficients,
            pauli_strings,
            total_shots=1000,
            max_group_size=3,
        )

        # Should have groups
        assert len(result.groups) > 0

        # All terms should be in exactly one group
        all_terms = set()
        for group in result.groups:
            for idx in group:
                assert idx not in all_terms
                all_terms.add(idx)
        assert all_terms == set(range(4))

        # Shots should sum to total
        assert sum(result.shots_per_group) == 1000

    def test_ir_hamiltonian_grouping_without_paulis(self):
        """Test grouping without Pauli strings."""
        from phaselab.quantum.grouping import ir_hamiltonian_grouping

        coefficients = np.array([1.0, -0.5, 0.3])

        result = ir_hamiltonian_grouping(
            coefficients,
            pauli_strings=None,
            total_shots=500,
        )

        assert len(result.groups) > 0
        assert result.method in ["simple_commuting", "atlas_q:vra_coherence"]


class TestVQEModule:
    """Test VQE optimization."""

    def test_vqe_config_defaults(self):
        """Test VQE configuration defaults."""
        from phaselab.quantum.vqe import VQEConfig

        config = VQEConfig()
        assert config.ansatz == "hardware_efficient"
        assert config.n_layers == 2
        assert config.enable_coherence is True

    def test_run_vqe_basic(self):
        """Test basic VQE run."""
        from phaselab.quantum.vqe import run_vqe, VQEConfig

        hamiltonian_terms = [
            (1.0, "ZZ"),
            (-0.5, "XX"),
            (0.3, "YY"),
        ]

        config = VQEConfig(max_iterations=10, verbose=False)
        result = run_vqe(hamiltonian_terms, n_qubits=2, config=config)

        assert isinstance(result.energy, float)
        assert len(result.optimal_parameters) > 0
        assert result.method in ["simple", "atlas_q"]

    def test_run_vqe_empty_hamiltonian(self):
        """Test VQE with empty Hamiltonian."""
        from phaselab.quantum.vqe import run_vqe

        result = run_vqe([], n_qubits=2)
        assert result.energy == 0.0
        assert result.converged is True


class TestBackendModule:
    """Test backend selection."""

    def test_get_available_backends(self):
        """Test listing available backends."""
        from phaselab.quantum.backend import get_available_backends

        backends = get_available_backends()
        assert len(backends) > 0

        # At least Qiskit Aer should be available
        names = [b.name for b in backends]
        assert "qiskit_aer" in names

    def test_get_optimal_backend_small(self):
        """Test backend selection for small circuits."""
        from phaselab.quantum.backend import get_optimal_backend, BackendType

        backend = get_optimal_backend(n_qubits=4, is_clifford=False)
        assert isinstance(backend, BackendType)

    def test_get_optimal_backend_large(self):
        """Test backend selection for large circuits."""
        from phaselab.quantum.backend import get_optimal_backend, BackendType

        backend = get_optimal_backend(n_qubits=30, is_clifford=False)
        # Should select MPS for large circuits
        assert backend in [BackendType.ATLAS_MPS, BackendType.QISKIT_AER]

    def test_get_optimal_backend_clifford(self):
        """Test backend selection for Clifford circuits."""
        from phaselab.quantum.backend import get_optimal_backend, BackendType

        if is_atlas_q_available():
            backend = get_optimal_backend(n_qubits=10, is_clifford=True)
            assert backend == BackendType.ATLAS_STABILIZER


class TestGPUModule:
    """Test GPU acceleration."""

    def test_check_gpu(self):
        """Test GPU info check."""
        from phaselab.quantum.gpu import check_gpu

        info = check_gpu()
        assert isinstance(info.available, bool)
        assert isinstance(info.device_name, str)

    def test_batch_coherence_cpu(self):
        """Test batch coherence on CPU."""
        from phaselab.quantum.gpu import batch_coherence_gpu

        batches = [
            np.array([0.9, 0.85, 0.88]),
            np.array([0.7, 0.65, 0.72]),
        ]

        R_bars = batch_coherence_gpu(batches)
        assert len(R_bars) == 2
        assert all(0 <= r <= 1 for r in R_bars)

    def test_batch_coherence_empty(self):
        """Test batch coherence with empty input."""
        from phaselab.quantum.gpu import batch_coherence_gpu

        R_bars = batch_coherence_gpu([])
        assert R_bars == []
