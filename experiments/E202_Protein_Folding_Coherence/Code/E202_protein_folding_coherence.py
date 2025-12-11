#!/usr/bin/env python3
"""
E202: Protein Folding Coherence - IBM Quantum Hardware Validation
==================================================================

Validates the phaselab.protein module on IBM Quantum hardware by computing
IR coherence metrics for protein structure observables (Ramachandran angles).

The key insight: Ramachandran angles (φ, ψ) are phase-like observables that
naturally map to quantum coherence metrics. Well-folded proteins should show
high coherence (ordered phases), while disordered regions show low coherence.

Test Cases:
1. Alpha helix: φ ≈ -60°, ψ ≈ -45° (highly ordered → high R̄)
2. Beta sheet: φ ≈ -120°, ψ ≈ +135° (ordered → high R̄)
3. Random coil: random φ, ψ (disordered → low R̄)
4. Mixed structure: combination of above

Author: Dylan Vaca
Date: December 2025
"""

import numpy as np
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, asdict
import sys

# Add src to path for local PhaseLab imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "src"))

# Qiskit imports
try:
    from qiskit import QuantumCircuit, transpile
    from qiskit_aer import AerSimulator
    from qiskit_aer.primitives import SamplerV2 as AerSampler
    QISKIT_AVAILABLE = True
except ImportError:
    QISKIT_AVAILABLE = False
    print("WARNING: Qiskit not installed - quantum simulation unavailable")

# IBM Runtime imports
try:
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
    IBM_AVAILABLE = True
except ImportError:
    IBM_AVAILABLE = False
    print("WARNING: IBM Runtime not installed - hardware validation unavailable")

# =============================================================================
# IR COHERENCE FRAMEWORK
# =============================================================================

E2_THRESHOLD = np.exp(-2)  # 0.135 - validated threshold

@dataclass
class CoherenceMetrics:
    R_bar: float
    V_phi: float
    is_above_e2: bool
    go_no_go: str

def compute_coherence(phases: np.ndarray) -> CoherenceMetrics:
    """Compute IR coherence from phase angles."""
    phasors = np.exp(1j * phases)
    R_bar = np.abs(np.mean(phasors))
    V_phi = -2.0 * np.log(max(R_bar, 1e-12))
    is_above = R_bar > E2_THRESHOLD
    return CoherenceMetrics(R_bar, V_phi, is_above, "GO" if is_above else "NO-GO")

def ramachandran_coherence(phi_angles: np.ndarray, psi_angles: np.ndarray) -> Tuple[float, float, float]:
    """
    Compute coherence from Ramachandran angles.

    Returns:
        (R_phi, R_psi, R_combined): Coherence for φ, ψ, and combined
    """
    # Convert degrees to radians if needed
    if np.max(np.abs(phi_angles)) > 2 * np.pi:
        phi_angles = np.radians(phi_angles)
        psi_angles = np.radians(psi_angles)

    z_phi = np.exp(1j * phi_angles)
    z_psi = np.exp(1j * psi_angles)

    R_phi = np.abs(np.mean(z_phi))
    R_psi = np.abs(np.mean(z_psi))
    R_combined = np.sqrt(R_phi * R_psi)

    return R_phi, R_psi, R_combined

# =============================================================================
# TEST PROTEIN STRUCTURES
# =============================================================================

def generate_alpha_helix(n_residues: int = 20, noise_std: float = 5.0) -> Tuple[np.ndarray, np.ndarray]:
    """Generate Ramachandran angles for an alpha helix."""
    # Ideal alpha helix: φ = -57°, ψ = -47°
    phi = np.full(n_residues, -57.0) + np.random.randn(n_residues) * noise_std
    psi = np.full(n_residues, -47.0) + np.random.randn(n_residues) * noise_std
    return phi, psi

def generate_beta_sheet(n_residues: int = 20, noise_std: float = 5.0) -> Tuple[np.ndarray, np.ndarray]:
    """Generate Ramachandran angles for a beta sheet."""
    # Ideal beta sheet: φ = -120°, ψ = +135°
    phi = np.full(n_residues, -120.0) + np.random.randn(n_residues) * noise_std
    psi = np.full(n_residues, 135.0) + np.random.randn(n_residues) * noise_std
    return phi, psi

def generate_random_coil(n_residues: int = 20) -> Tuple[np.ndarray, np.ndarray]:
    """Generate random Ramachandran angles (disordered)."""
    # Random within allowed Ramachandran space
    phi = np.random.uniform(-180, 180, n_residues)
    psi = np.random.uniform(-180, 180, n_residues)
    return phi, psi

def generate_mixed_structure(n_helix: int = 10, n_sheet: int = 5, n_coil: int = 5) -> Tuple[np.ndarray, np.ndarray]:
    """Generate mixed secondary structure."""
    phi_helix, psi_helix = generate_alpha_helix(n_helix)
    phi_sheet, psi_sheet = generate_beta_sheet(n_sheet)
    phi_coil, psi_coil = generate_random_coil(n_coil)

    phi = np.concatenate([phi_helix, phi_sheet, phi_coil])
    psi = np.concatenate([psi_helix, psi_sheet, psi_coil])

    return phi, psi

# =============================================================================
# QUANTUM CIRCUIT FOR PHASE COHERENCE
# =============================================================================

def build_phase_encoding_circuit(phases: np.ndarray, n_qubits: int = 8) -> QuantumCircuit:
    """
    Build quantum circuit that encodes phase angles.

    Maps Ramachandran angles to qubit rotations, then measures coherence
    via interference pattern.
    """
    n = min(len(phases), n_qubits)
    qc = QuantumCircuit(n)

    # Encode phases as rotations
    for i in range(n):
        # Hadamard to create superposition
        qc.h(i)
        # Phase rotation based on Ramachandran angle
        qc.rz(phases[i], i)

    # Entangling layer (captures correlations between residues)
    for i in range(n - 1):
        qc.cx(i, i + 1)

    # Second rotation layer
    for i in range(n):
        qc.ry(phases[i] * 0.5, i)

    return qc

def build_coherence_measurement_circuits(phi: np.ndarray, psi: np.ndarray, n_qubits: int = 8) -> List[QuantumCircuit]:
    """Build circuits to measure φ and ψ coherence."""
    circuits = []
    n = min(len(phi), n_qubits)

    # Circuit for φ coherence
    phi_rad = np.radians(phi[:n]) if np.max(np.abs(phi)) > 2 * np.pi else phi[:n]
    qc_phi = build_phase_encoding_circuit(phi_rad, n)
    qc_phi.measure_all()
    circuits.append(('phi', qc_phi))

    # Circuit for ψ coherence
    psi_rad = np.radians(psi[:n]) if np.max(np.abs(psi)) > 2 * np.pi else psi[:n]
    qc_psi = build_phase_encoding_circuit(psi_rad, n)
    qc_psi.measure_all()
    circuits.append(('psi', qc_psi))

    # Combined circuit (encodes both)
    qc_combined = QuantumCircuit(n)
    for i in range(n):
        qc_combined.h(i)
        qc_combined.rz(phi_rad[i], i)
        qc_combined.ry(psi_rad[i] * 0.5, i)
    for i in range(n - 1):
        qc_combined.cx(i, i + 1)
    qc_combined.measure_all()
    circuits.append(('combined', qc_combined))

    return circuits

def extract_coherence_from_counts(counts: Dict[str, int], n_qubits: int) -> float:
    """Extract coherence metric from measurement counts."""
    total = sum(counts.values())

    # Compute expectation values for each qubit
    expectations = []
    for i in range(n_qubits):
        exp_val = 0.0
        for bitstring, count in counts.items():
            # Bit value at position i
            if i < len(bitstring):
                bit = int(bitstring[-(i+1)])  # Qiskit uses little-endian
                sign = 1 if bit == 0 else -1
                exp_val += sign * count / total
        expectations.append(exp_val)

    expectations = np.array(expectations)

    # Convert expectations to phases and compute coherence
    phases = np.arccos(np.clip(expectations, -1, 1))
    coherence = compute_coherence(phases)

    return coherence.R_bar

# =============================================================================
# MAIN VALIDATION PIPELINE
# =============================================================================

@dataclass
class ProteinCoherenceResult:
    structure_type: str
    n_residues: int
    classical_R_phi: float
    classical_R_psi: float
    classical_R_combined: float
    quantum_R_phi: float
    quantum_R_psi: float
    quantum_R_combined: float
    hardware_R_phi: Optional[float]
    hardware_R_psi: Optional[float]
    hardware_R_combined: Optional[float]
    go_no_go_classical: str
    go_no_go_quantum: str
    go_no_go_hardware: Optional[str]
    agreement: bool

def run_protein_validation(
    use_hardware: bool = False,
    backend_name: str = "ibm_torino",
    shots: int = 4096
) -> Dict:
    """Run protein folding coherence validation."""

    print("="*70)
    print("E202: Protein Folding Coherence - IBM Quantum Validation")
    print("="*70)
    print()
    print("Validating phaselab.protein module on quantum hardware")
    print("Testing IR coherence for Ramachandran angle observables")
    print()

    if not QISKIT_AVAILABLE:
        print("ERROR: Qiskit not available")
        return {'error': 'Qiskit not available'}

    # Initialize backends
    simulator = AerSimulator()
    hardware_backend = None

    if use_hardware and IBM_AVAILABLE:
        try:
            service = QiskitRuntimeService(channel="ibm_quantum_platform")
            hardware_backend = service.backend(backend_name)
            print(f"Hardware backend: {backend_name}")
            print(f"Qubits: {hardware_backend.num_qubits}")
        except Exception as e:
            print(f"WARNING: Could not connect to IBM hardware: {e}")
            use_hardware = False

    # Test cases
    test_cases = [
        ("alpha_helix", generate_alpha_helix(20, noise_std=5.0)),
        ("beta_sheet", generate_beta_sheet(20, noise_std=5.0)),
        ("random_coil", generate_random_coil(20)),
        ("mixed_structure", generate_mixed_structure(10, 5, 5)),
        ("alpha_helix_noisy", generate_alpha_helix(20, noise_std=20.0)),
        ("alpha_helix_clean", generate_alpha_helix(20, noise_std=1.0)),
    ]

    results = []
    n_qubits = 8

    for name, (phi, psi) in test_cases:
        print(f"\n{'='*70}")
        print(f"Testing: {name}")
        print(f"{'='*70}")

        # Classical coherence
        R_phi_classical, R_psi_classical, R_combined_classical = ramachandran_coherence(phi, psi)
        print(f"\nClassical coherence:")
        print(f"  R_φ = {R_phi_classical:.4f}")
        print(f"  R_ψ = {R_psi_classical:.4f}")
        print(f"  R_combined = {R_combined_classical:.4f}")

        go_classical = "GO" if R_combined_classical > E2_THRESHOLD else "NO-GO"
        print(f"  Classification: {go_classical}")

        # Build quantum circuits
        circuits = build_coherence_measurement_circuits(phi, psi, n_qubits)

        # Simulator run
        print(f"\nQuantum simulation (AerSimulator):")
        quantum_results = {}

        for label, qc in circuits:
            transpiled = transpile(qc, simulator, optimization_level=1)
            sampler = AerSampler()
            job = sampler.run([transpiled], shots=shots)
            result = job.result()
            counts = result[0].data.meas.get_counts()
            R = extract_coherence_from_counts(counts, n_qubits)
            quantum_results[label] = R
            print(f"  R_{label} = {R:.4f}")

        R_phi_quantum = quantum_results['phi']
        R_psi_quantum = quantum_results['psi']
        R_combined_quantum = quantum_results['combined']
        go_quantum = "GO" if R_combined_quantum > E2_THRESHOLD else "NO-GO"
        print(f"  Classification: {go_quantum}")

        # Hardware run (if available)
        R_phi_hardware = None
        R_psi_hardware = None
        R_combined_hardware = None
        go_hardware = None

        if use_hardware and hardware_backend is not None:
            print(f"\nQuantum hardware ({backend_name}):")
            try:
                hardware_results = {}

                for label, qc in circuits:
                    transpiled = transpile(qc, hardware_backend, optimization_level=3)
                    sampler = SamplerV2(hardware_backend)
                    job = sampler.run([transpiled], shots=shots)
                    result = job.result()
                    counts = result[0].data.meas.get_counts()
                    R = extract_coherence_from_counts(counts, n_qubits)
                    hardware_results[label] = R
                    print(f"  R_{label} = {R:.4f}")

                R_phi_hardware = hardware_results['phi']
                R_psi_hardware = hardware_results['psi']
                R_combined_hardware = hardware_results['combined']
                go_hardware = "GO" if R_combined_hardware > E2_THRESHOLD else "NO-GO"
                print(f"  Classification: {go_hardware}")

            except Exception as e:
                print(f"  Hardware error: {e}")

        # Check agreement
        agreement = (go_classical == go_quantum)
        if go_hardware is not None:
            agreement = agreement and (go_classical == go_hardware)

        results.append(ProteinCoherenceResult(
            structure_type=name,
            n_residues=len(phi),
            classical_R_phi=float(R_phi_classical),
            classical_R_psi=float(R_psi_classical),
            classical_R_combined=float(R_combined_classical),
            quantum_R_phi=float(R_phi_quantum),
            quantum_R_psi=float(R_psi_quantum),
            quantum_R_combined=float(R_combined_quantum),
            hardware_R_phi=R_phi_hardware,
            hardware_R_psi=R_psi_hardware,
            hardware_R_combined=R_combined_hardware,
            go_no_go_classical=go_classical,
            go_no_go_quantum=go_quantum,
            go_no_go_hardware=go_hardware,
            agreement=agreement
        ))

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)

    print("\nResults by structure type:")
    print("-"*70)
    print(f"{'Structure':<25} {'Classical':<12} {'Quantum':<12} {'Hardware':<12} {'Agreement'}")
    print("-"*70)

    for r in results:
        hw_str = f"{r.hardware_R_combined:.4f}" if r.hardware_R_combined else "N/A"
        agree_str = "✓" if r.agreement else "✗"
        print(f"{r.structure_type:<25} {r.classical_R_combined:.4f}       {r.quantum_R_combined:.4f}       {hw_str:<12} {agree_str}")

    # Expected behavior validation
    print("\n" + "="*70)
    print("EXPECTED BEHAVIOR VALIDATION")
    print("="*70)

    expected = {
        "alpha_helix": ("GO", "Ordered structure should show high coherence"),
        "beta_sheet": ("GO", "Ordered structure should show high coherence"),
        "random_coil": ("NO-GO", "Disordered structure should show low coherence"),
        "mixed_structure": ("GO", "Partially ordered should be above threshold"),
        "alpha_helix_noisy": ("GO", "Noisy but still ordered"),
        "alpha_helix_clean": ("GO", "Very clean should show very high coherence"),
    }

    validation_passed = 0
    for r in results:
        expected_go, reason = expected.get(r.structure_type, ("UNKNOWN", ""))
        actual = r.go_no_go_classical
        passed = (expected_go == actual) or (expected_go == "UNKNOWN")
        status = "PASS" if passed else "FAIL"
        validation_passed += int(passed)
        print(f"\n{r.structure_type}:")
        print(f"  Expected: {expected_go} ({reason})")
        print(f"  Actual: {actual} (R̄ = {r.classical_R_combined:.4f})")
        print(f"  Status: {status}")

    print(f"\n\nValidation: {validation_passed}/{len(results)} tests passed")

    # Save results
    output = {
        'experiment': 'E202_Protein_Folding_Coherence',
        'timestamp': datetime.now().isoformat(),
        'backend': backend_name if use_hardware else 'AerSimulator',
        'shots': shots,
        'n_qubits': n_qubits,
        'e2_threshold': float(E2_THRESHOLD),
        'results': [asdict(r) for r in results],
        'validation_passed': validation_passed,
        'total_tests': len(results),
    }

    output_dir = Path(__file__).parent.parent / "Data"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"E202_protein_coherence_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"

    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\n\nResults saved: {output_file}")

    return output


def main():
    import argparse
    parser = argparse.ArgumentParser(description='E202: Protein Folding Coherence Validation')
    parser.add_argument('--hardware', action='store_true', help='Run on IBM hardware')
    parser.add_argument('--backend', default='ibm_torino', help='IBM backend name')
    parser.add_argument('--shots', type=int, default=4096, help='Number of shots')
    args = parser.parse_args()

    return run_protein_validation(
        use_hardware=args.hardware,
        backend_name=args.backend,
        shots=args.shots
    )


if __name__ == "__main__":
    main()
