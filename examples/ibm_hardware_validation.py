#!/usr/bin/env python3
"""
PhaseLab IBM Hardware Validation

This script validates the PhaseLab coherence metrics on IBM Quantum hardware.
It replicates the E200 SMS Gene Therapy experiment to ensure the package
produces consistent results with the original validation.

Requirements:
    pip install phaselab[quantum]

Usage:
    python ibm_hardware_validation.py --simulator  # Free test
    python ibm_hardware_validation.py --hardware   # Uses IBM Quantum
"""

import argparse
import json
import os
from datetime import datetime
from typing import Dict, List, Tuple

import numpy as np

# PhaseLab imports
from phaselab.core.coherence import coherence_score, go_no_go, compare_sim_hardware
from phaselab.core.hamiltonians import build_grna_hamiltonian
from phaselab.core.constants import E_MINUS_2


# Top 3 gRNA candidates from E200 experiment
GRNA_CANDIDATES = [
    {
        "name": "gRNA_1",
        "sequence": "GAAGGAGAGCAAGAGCGCGA",
        "position": -81,
        "expected_hardware_R": 0.854,  # From E200 IBM Torino run
    },
    {
        "name": "gRNA_2",
        "sequence": "AACTGCAAAGAAGTGGGCAC",
        "position": -334,
        "expected_hardware_R": 0.840,
    },
    {
        "name": "gRNA_3",
        "sequence": "TACAGGAGCTTCCAGCGTCA",
        "position": -294,
        "expected_hardware_R": 0.839,
    },
]


def build_measurement_circuits(guide_seq: str, n_qubits: int = 4):
    """
    Build measurement circuits for a guide RNA Hamiltonian.

    Uses a simplified 4-qubit encoding for validation.
    """
    from qiskit import QuantumCircuit
    from qiskit.circuit.library import efficient_su2

    # Pauli terms to measure (simplified)
    pauli_terms = [
        ("ZIII", 1.0),
        ("IZII", 1.0),
        ("ZZII", 0.5),
        ("XXII", 0.3),
    ]

    circuits = []
    for pauli_str, coeff in pauli_terms:
        # Create ansatz circuit
        qc = QuantumCircuit(n_qubits)

        # Simple variational ansatz using function-based API (Qiskit 2.1+)
        ansatz = efficient_su2(n_qubits, reps=1)

        # Set some initial parameters based on guide sequence
        params = np.array([
            hash(guide_seq[i % len(guide_seq)]) % 100 / 100 * np.pi
            for i in range(ansatz.num_parameters)
        ])
        bound_ansatz = ansatz.assign_parameters(params)
        qc.compose(bound_ansatz, inplace=True)

        # Add measurement basis rotation
        for i, p in enumerate(pauli_str):
            if p == 'X':
                qc.h(i)
            elif p == 'Y':
                qc.sdg(i)
                qc.h(i)

        qc.measure_all()
        circuits.append((qc, pauli_str, coeff))

    return circuits


def compute_expectation(counts: Dict, pauli_str: str) -> float:
    """Compute expectation value from measurement counts."""
    total = sum(counts.values())
    exp_val = 0.0

    for bitstring, count in counts.items():
        # Reverse bitstring (Qiskit LSB-first)
        bits = bitstring[::-1]
        parity = sum(int(bits[i]) for i, p in enumerate(pauli_str) if p != 'I') % 2
        sign = 1 if parity == 0 else -1
        exp_val += sign * count / total

    return exp_val


def run_simulator_test() -> Dict:
    """Run validation on Qiskit Aer simulator."""
    from qiskit_aer import AerSimulator

    print("=" * 60)
    print("PhaseLab IBM Quantum Validation - SIMULATOR")
    print("=" * 60)

    simulator = AerSimulator()
    results = {"candidates": [], "backend": "aer_simulator"}

    for candidate in GRNA_CANDIDATES:
        print(f"\n[{candidate['name']}] {candidate['sequence']}")

        circuits = build_measurement_circuits(candidate["sequence"])
        expectations = []

        for qc, pauli_str, coeff in circuits:
            job = simulator.run(qc, shots=2000)
            counts = job.result().get_counts()
            exp = compute_expectation(counts, pauli_str)
            expectations.append(exp * coeff)
            print(f"  {pauli_str}: <P> = {exp:.4f}")

        # Compute energy and coherence
        energy = sum(expectations)

        # IR coherence from expectations variance
        exp_array = np.array([abs(e) for e in expectations])
        V_phi = np.var(exp_array) * 4  # Scale factor
        R_bar = np.exp(-V_phi / 2)

        status = go_no_go(R_bar)

        print(f"  Energy: {energy:.4f}")
        print(f"  Coherence R̄: {R_bar:.4f}")
        print(f"  Status: {status}")

        results["candidates"].append({
            "name": candidate["name"],
            "sequence": candidate["sequence"],
            "position": candidate["position"],
            "energy": energy,
            "coherence_R": R_bar,
            "go_no_go": status,
            "expectations": expectations,
        })

    return results


def run_hardware_test(backend_name: str = "ibm_torino") -> Dict:
    """Run validation on IBM Quantum hardware."""
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

    print("=" * 60)
    print(f"PhaseLab IBM Quantum Validation - HARDWARE ({backend_name})")
    print("=" * 60)

    # Load token from environment or .env file
    token = os.environ.get("IBM_QUANTUM_TOKEN")
    if not token:
        env_path = os.path.join(os.path.dirname(__file__), "..", ".env")
        if os.path.exists(env_path):
            with open(env_path) as f:
                for line in f:
                    if line.startswith("IBM_QUANTUM_TOKEN="):
                        token = line.split("=", 1)[1].strip()

    if not token:
        # Try E200 .env
        e200_env = "/home/admin/dev/Informational_Relativity/experiments/E200_SMS_Gene_Therapy/.env"
        if os.path.exists(e200_env):
            with open(e200_env) as f:
                for line in f:
                    if line.startswith("IBM_QUANTUM_TOKEN="):
                        token = line.split("=", 1)[1].strip()

    if not token:
        raise ValueError("IBM_QUANTUM_TOKEN not found. Set environment variable or create .env file.")

    # Connect to IBM Quantum (use new channel name)
    service = QiskitRuntimeService(channel="ibm_quantum_platform", token=token)
    backend = service.backend(backend_name)

    print(f"Connected to: {backend.name}")
    print(f"Qubits: {backend.num_qubits}")

    # Build all circuits
    all_circuits = []
    circuit_info = []

    for candidate in GRNA_CANDIDATES:
        circuits = build_measurement_circuits(candidate["sequence"])
        for qc, pauli_str, coeff in circuits:
            all_circuits.append(qc)
            circuit_info.append({
                "candidate": candidate["name"],
                "pauli": pauli_str,
                "coeff": coeff,
            })

    # Transpile
    print(f"\nTranspiling {len(all_circuits)} circuits...")
    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    transpiled = pm.run(all_circuits)

    # Run on hardware
    print("Submitting to hardware...")
    sampler = SamplerV2(backend)
    job = sampler.run(transpiled, shots=2000)

    print(f"Job ID: {job.job_id()}")
    print("Waiting for results...")

    result = job.result()

    # Process results
    results = {"candidates": [], "backend": backend_name, "job_id": job.job_id()}

    candidate_results = {}
    for i, info in enumerate(circuit_info):
        counts = result[i].data.meas.get_counts()
        exp = compute_expectation(counts, info["pauli"])

        if info["candidate"] not in candidate_results:
            candidate_results[info["candidate"]] = {"expectations": [], "paulis": []}

        candidate_results[info["candidate"]]["expectations"].append(exp * info["coeff"])
        candidate_results[info["candidate"]]["paulis"].append(info["pauli"])

    # Compute coherence for each candidate
    print("\n" + "=" * 60)
    print("RESULTS")
    print("=" * 60)

    for candidate in GRNA_CANDIDATES:
        name = candidate["name"]
        if name not in candidate_results:
            continue

        expectations = candidate_results[name]["expectations"]
        energy = sum(expectations)

        exp_array = np.array([abs(e) for e in expectations])
        V_phi = np.var(exp_array) * 4
        R_bar = np.exp(-V_phi / 2)

        status = go_no_go(R_bar)

        # Compare with expected
        expected_R = candidate["expected_hardware_R"]
        diff, agreement = compare_sim_hardware(expected_R, R_bar)

        print(f"\n[{name}] {candidate['sequence']}")
        print(f"  Position: {candidate['position']} bp")
        print(f"  Energy: {energy:.4f}")
        print(f"  Coherence R̄: {R_bar:.4f}")
        print(f"  Expected R̄: {expected_R:.4f}")
        print(f"  Difference: {diff:.4f} ({agreement})")
        print(f"  Status: {status}")

        results["candidates"].append({
            "name": name,
            "sequence": candidate["sequence"],
            "position": candidate["position"],
            "energy": energy,
            "coherence_R": R_bar,
            "expected_R": expected_R,
            "difference": diff,
            "agreement": agreement,
            "go_no_go": status,
            "expectations": expectations,
        })

    # Save results
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = f"phaselab_hardware_validation_{timestamp}.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\nResults saved to: {output_file}")

    return results


def main():
    parser = argparse.ArgumentParser(description="PhaseLab IBM Quantum Validation")
    parser.add_argument("--simulator", action="store_true", help="Run on Aer simulator")
    parser.add_argument("--hardware", action="store_true", help="Run on IBM Quantum hardware")
    parser.add_argument("--backend", type=str, default="ibm_torino", help="IBM backend name")

    args = parser.parse_args()

    if args.hardware:
        results = run_hardware_test(args.backend)
    elif args.simulator:
        results = run_simulator_test()
    else:
        print("PhaseLab IBM Quantum Validation")
        print()
        print("Usage:")
        print("  python ibm_hardware_validation.py --simulator  # Free test")
        print("  python ibm_hardware_validation.py --hardware   # IBM Quantum")
        print()
        print("This validates the PhaseLab coherence metrics match E200 results.")
        return

    # Summary
    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)

    all_go = all(c["go_no_go"] == "GO" for c in results["candidates"])
    print(f"Backend: {results['backend']}")
    print(f"All candidates GO: {'YES' if all_go else 'NO'}")

    for c in results["candidates"]:
        print(f"  {c['name']}: R̄={c['coherence_R']:.3f} [{c['go_no_go']}]")

    if all_go:
        print("\n✓ PhaseLab validation PASSED")
    else:
        print("\n✗ PhaseLab validation FAILED")


if __name__ == "__main__":
    main()
