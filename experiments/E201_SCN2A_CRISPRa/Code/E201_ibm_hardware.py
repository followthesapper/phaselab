#!/usr/bin/env python3
"""
E201: SCN2A IBM Quantum Hardware Validation
============================================

Run the top SCN2A CRISPRa guide candidates on real IBM Quantum hardware
to validate the coherence metric predictions from simulator.

Same protocol as E200 (RAI1/SMS) validation on IBM Torino.

Author: Dylan Vaca
Date: December 2025
"""

import os
import json
import numpy as np
from pathlib import Path
from datetime import datetime
from dotenv import load_dotenv

# Load environment variables
env_path = Path(__file__).parent.parent / ".env"
if env_path.exists():
    load_dotenv(env_path)

# IBM Quantum imports
try:
    from qiskit import QuantumCircuit, transpile
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
    IBM_AVAILABLE = True
except ImportError:
    IBM_AVAILABLE = False
    print("WARNING: qiskit-ibm-runtime not installed")

# =============================================================================
# TOP SCN2A GUIDE CANDIDATES (from simulator results)
# =============================================================================

TOP_GUIDES = [
    {
        "sequence": "TTCCACTTTTGACCAGGAGA",
        "position": -64,
        "gc": 0.45,
        "chromatin": "OPEN (DNase HS)",
        "simulator_coherence": 0.944,
    },
    {
        "sequence": "AGATGGTTCCACTTTTGACC",
        "position": -70,
        "gc": 0.45,
        "chromatin": "OPEN (DNase HS)",
        "simulator_coherence": 0.945,
    },
    {
        "sequence": "GCTGACTGCTACATAGCCAA",
        "position": -104,
        "gc": 0.50,
        "chromatin": "OPEN (DNase HS)",
        "simulator_coherence": 0.945,
    },
]

# SCN2A promoter (subset for target region)
SCN2A_PROMOTER = """
CAAATCTGACAAAGACAGTATTTCAGACCCGTTTAAGACCTAGAACTCACCATGAGTTCTAAAATTGGTT
CTCAGCACCATGGACAGCGTTACTGCAATAGGAAATTAAAGATCGATTTGGCCCCAAATTAAAATGGTGT
TGTAAAAAAGTGGGAGAAAAAAAATGCCTATCCTTTTACTTCAAATTTTAAAAAAATGATCCTGGCCTTC
ACAACTGTTCATAAGAAGAATAATTAATTAAACAAACATATATTGAGAACATCATATGCTCAGTAAAACT
TTGATTCTATAAATGGTGTCATTTACCAAATGGATTCTTTTGACAATTTAATTTTCTCTTATCTCTCTAA
GAAGATGTAACTACACACTATAGTATACTACTACAATTATCAAATTTCATGTTGCATGTAACTTGTCGTC
TGTATTTTTGTAGTTAGATTAGATTAACTAAAGATTTTTCAAGTTTGCCTTTAAGTCATTTAATTTTCCT
GCCTTATCTTTAACCTTTCAACATTCCTCCAAACAATAGCAACACAAGTGTTATGTGTTAACTTCTCTAG
TGACAAAAACTTATACTTCTCCACAAAGAGATGTGATGTTCATTATCAATAAGCTTGACATCTAAAATTG
TTTTATAGGAGATACATATTACTTTTTCAGATGGTATATAAAGTTAAATAAATCTTAAGTTTTCAATGAT
GGGAAAAGCTTCCATTTAGTTTAAACATAATGTAAAGAAATTTGAATCCCCAAAATAGAATTATAATTCT
AAAAATTCATACTATAATTCTTCTTAAATGTTTAAATTACAGTTAATTAAAGTAGTTGATTTCAAATAGA
GTGGAATTATGGGCTGTACATCATTTAATTTTATGTGCTGACTGCTACATAGCCAAAGGAACGTGAATTA
AGATGGTTCCACTTTTGACCAGGAGATGGAGCTGTCATGTAAGATGCTGCCTTTATTTATTTATTTTTCT
AATTTAGCATGCTGTTTTCTAACAGACATTGGGTACCATCGAATGACTGTCAGAACAGAAAGCTAAGGCA
""".replace('\n', '').replace(' ', '')

TSS_INDEX = 1000

# =============================================================================
# IR COHERENCE
# =============================================================================

E2_THRESHOLD = np.exp(-2)

def compute_coherence(measurements):
    """Compute IR coherence from measurement outcomes."""
    phases = np.arccos(np.clip(measurements, -1, 1))
    phasors = np.exp(1j * phases)
    R_bar = np.abs(np.mean(phasors))
    return float(R_bar)

# =============================================================================
# HAMILTONIAN BUILDING
# =============================================================================

NN_PARAMS = {
    'AA': -1.00, 'AT': -0.88, 'AG': -1.28, 'AC': -1.44,
    'TA': -0.58, 'TT': -1.00, 'TG': -1.45, 'TC': -1.30,
    'GA': -1.30, 'GT': -1.44, 'GG': -1.84, 'GC': -2.24,
    'CA': -1.45, 'CT': -1.28, 'CG': -2.17, 'CC': -1.84,
}

def complement(base):
    return {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}.get(base, 'N')

def build_hamiltonian(grna, target, n_qubits=4):
    """Build simplified Hamiltonian for hardware validation."""
    coeffs = []
    paulis = []
    n = n_qubits

    # Use first n bases
    grna_short = grna[:n]
    target_short = target[:n]

    # ZZ interactions
    for i in range(n - 1):
        dinuc = grna_short[i:i+2]
        base_energy = NN_PARAMS.get(dinuc, -1.0)

        match_i = grna_short[i] == complement(target_short[i])
        match_j = grna_short[i+1] == complement(target_short[i+1])

        if match_i and match_j:
            energy = base_energy
        elif match_i or match_j:
            energy = base_energy * 0.3
        else:
            energy = abs(base_energy) * 0.5

        pauli = ['I'] * n
        pauli[i] = 'Z'
        pauli[i+1] = 'Z'
        paulis.append(''.join(pauli))
        coeffs.append(energy)

    # Single Z terms
    for i in range(n):
        pauli = ['I'] * n
        pauli[i] = 'Z'
        paulis.append(''.join(pauli))
        if grna_short[i] == complement(target_short[i]):
            coeffs.append(-0.5)
        else:
            coeffs.append(0.3)

    return np.array(coeffs), paulis

def build_ansatz(n_qubits=4):
    """Build EfficientSU2-style ansatz."""
    qc = QuantumCircuit(n_qubits)
    np.random.seed(42)

    for i in range(n_qubits):
        qc.ry(np.random.randn() * 0.5, i)

    for i in range(n_qubits - 1):
        qc.cx(i, i + 1)

    for i in range(n_qubits):
        qc.ry(np.random.randn() * 0.5, i)

    return qc

# =============================================================================
# HARDWARE VALIDATION
# =============================================================================

def run_hardware_validation(backend_name="ibm_torino", shots=2000):
    """Run top guides on IBM Quantum hardware."""

    print("="*70)
    print("E201: SCN2A IBM Quantum Hardware Validation")
    print("="*70)
    print()

    if not IBM_AVAILABLE:
        print("ERROR: qiskit-ibm-runtime not available")
        return None

    # Initialize service
    token = os.environ.get("IBM_QUANTUM_TOKEN")
    if not token:
        print("ERROR: IBM_QUANTUM_TOKEN not set")
        print("Set it in experiments/E201_SCN2A_CRISPRa/.env")
        return None

    print(f"Connecting to IBM Quantum...")
    service = QiskitRuntimeService(
        channel="ibm_quantum_platform",
        token=token
    )

    print(f"Using backend: {backend_name}")
    backend = service.backend(backend_name)
    print(f"Backend status: {backend.status().status_msg}")
    print()

    results = {
        "experiment": "E201_SCN2A_Hardware",
        "timestamp": datetime.now().isoformat(),
        "backend": backend_name,
        "shots": shots,
        "guides": []
    }

    # Process each guide
    for idx, guide_info in enumerate(TOP_GUIDES):
        guide_seq = guide_info["sequence"]
        print(f"\n[{idx+1}/{len(TOP_GUIDES)}] Processing: {guide_seq}")
        print(f"    Position: {guide_info['position']} bp from TSS")
        print(f"    GC: {guide_info['gc']:.0%}")
        print(f"    Chromatin: {guide_info['chromatin']}")

        # Get target region (using simplified approach)
        target_region = "A" * 20  # Placeholder - in production, extract from promoter

        # Build Hamiltonian
        coeffs, paulis = build_hamiltonian(guide_seq, target_region)
        ansatz = build_ansatz()

        # Build measurement circuits
        circuits = []
        for pauli_str in paulis:
            qc = ansatz.copy()
            for i, p in enumerate(pauli_str):
                if p == 'X':
                    qc.h(i)
                elif p == 'Y':
                    qc.sdg(i)
                    qc.h(i)
            qc.measure_all()
            circuits.append(qc)

        print(f"    Transpiling {len(circuits)} circuits...")
        transpiled = transpile(circuits, backend, optimization_level=1)

        print(f"    Submitting to {backend_name}...")
        sampler = SamplerV2(backend)
        job = sampler.run(transpiled, shots=shots)

        print(f"    Job ID: {job.job_id()}")
        print(f"    Waiting for results...")

        result = job.result()

        # Compute expectations
        expectations = []
        for i, pauli_str in enumerate(paulis):
            pub_result = result[i]
            counts = pub_result.data.meas.get_counts()

            exp_val = 0.0
            total = sum(counts.values())
            for bitstring, count in counts.items():
                parity = sum(int(bitstring[j]) for j, p in enumerate(pauli_str) if p != 'I') % 2
                sign = 1 if parity == 0 else -1
                exp_val += sign * count / total
            expectations.append(exp_val)

        expectations = np.array(expectations)

        # Compute coherence
        hardware_coherence = compute_coherence(expectations)
        quantum_energy = float(np.sum(coeffs * expectations))

        # Compare with simulator
        sim_coherence = guide_info["simulator_coherence"]
        difference = abs(hardware_coherence - sim_coherence) / sim_coherence * 100

        go_status = "GO" if hardware_coherence > E2_THRESHOLD else "NO-GO"

        print(f"\n    RESULTS:")
        print(f"    Simulator R̄:  {sim_coherence:.3f}")
        print(f"    Hardware R̄:   {hardware_coherence:.3f}")
        print(f"    Difference:    {difference:.1f}%")
        print(f"    Quantum E:     {quantum_energy:.2f}")
        print(f"    Status:        {go_status}")

        results["guides"].append({
            "sequence": guide_seq,
            "position": guide_info["position"],
            "gc_content": guide_info["gc"],
            "chromatin": guide_info["chromatin"],
            "simulator_coherence": sim_coherence,
            "hardware_coherence": hardware_coherence,
            "difference_percent": difference,
            "quantum_energy": quantum_energy,
            "go_no_go": go_status,
            "job_id": job.job_id()
        })

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)

    all_go = all(g["go_no_go"] == "GO" for g in results["guides"])
    avg_diff = np.mean([g["difference_percent"] for g in results["guides"]])

    print(f"\nAll guides GO: {all_go}")
    print(f"Average sim-hardware difference: {avg_diff:.1f}%")

    for g in results["guides"]:
        print(f"\n{g['sequence']}")
        print(f"  Sim: {g['simulator_coherence']:.3f} → HW: {g['hardware_coherence']:.3f} [{g['go_no_go']}]")

    # Save results
    output_dir = Path(__file__).parent.parent / "Data"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"E201_hardware_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"

    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved: {output_file}")

    return results


def main():
    return run_hardware_validation()


if __name__ == "__main__":
    main()
