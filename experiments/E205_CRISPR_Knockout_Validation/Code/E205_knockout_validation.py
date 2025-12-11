"""
E205: CRISPR Knockout Module - IBM Quantum Hardware Validation

Validates the phaselab.crispr.knockout module on IBM Quantum hardware.
Tests cutting efficiency prediction reliability using IR coherence metrics.

Author: Dylan Vaca
Date: December 2025
"""

import sys
import os
import json
from datetime import datetime
import numpy as np

# Add parent paths
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'src'))

# Try imports
try:
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import SparsePauliOp
    from qiskit_aer import AerSimulator
    QISKIT_AVAILABLE = True
except ImportError:
    QISKIT_AVAILABLE = False
    print("Warning: Qiskit not available, using classical simulation only")

try:
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
    IBM_AVAILABLE = True
except ImportError:
    IBM_AVAILABLE = False
    print("Warning: IBM Runtime not available")


# ============================================================================
# KNOCKOUT SCORING FUNCTIONS
# ============================================================================

def gc_content(seq):
    """Calculate GC content."""
    seq = seq.upper()
    return sum(1 for b in seq if b in 'GC') / len(seq) if seq else 0


def cut_efficiency_score(guide_seq):
    """Predict cutting efficiency using position-specific model."""
    guide_seq = guide_seq.upper()
    if len(guide_seq) != 20:
        return 0.5

    score = 0.5

    # Position-specific weights (Rule Set 2 approximation)
    position_weights = {
        1: {'T': -0.1}, 2: {'T': -0.1, 'C': -0.05}, 4: {'A': 0.05},
        13: {'G': 0.1}, 14: {'G': 0.1, 'C': 0.05}, 15: {'G': 0.1},
        16: {'G': 0.15, 'C': 0.1}, 17: {'C': 0.1},
        18: {'G': 0.1, 'A': 0.05}, 19: {'G': 0.1, 'C': 0.1},
        20: {'G': 0.2, 'A': -0.1},
    }

    for i, base in enumerate(guide_seq):
        pos = i + 1
        if pos in position_weights:
            score += position_weights[pos].get(base, 0)

    gc = gc_content(guide_seq)
    if 0.45 <= gc <= 0.65:
        score += 0.1
    elif gc < 0.30 or gc > 0.80:
        score -= 0.15

    if 'TTTT' in guide_seq:
        score -= 0.2

    if guide_seq[-1] == 'G':
        score += 0.1

    return float(np.clip(score, 0, 1))


def frameshift_probability(cds_position, exon_length):
    """Estimate frameshift probability from cut position."""
    relative_pos = cds_position / exon_length if exon_length > 0 else 0.5

    if relative_pos < 0.25:
        base_prob = 0.90
    elif relative_pos < 0.50:
        base_prob = 0.75
    elif relative_pos < 0.75:
        base_prob = 0.60
    else:
        base_prob = 0.45

    return base_prob * 0.66


# ============================================================================
# QUANTUM CIRCUITS
# ============================================================================

def create_knockout_coherence_circuit(guide_seq, n_qubits=5):
    """
    Create quantum circuit to validate knockout guide coherence.

    Encodes cutting efficiency factors as phase rotations.
    """
    qc = QuantumCircuit(n_qubits)

    # Hadamard superposition
    for i in range(n_qubits):
        qc.h(i)

    # Encode guide properties as rotations
    gc = gc_content(guide_seq)
    cut_eff = cut_efficiency_score(guide_seq)

    # GC content encoding
    qc.rz(gc * np.pi, 0)

    # Cutting efficiency encoding
    qc.rz(cut_eff * np.pi, 1)

    # Position-specific encoding (seed region importance)
    seed_region = guide_seq[-12:] if len(guide_seq) >= 12 else guide_seq
    seed_gc = gc_content(seed_region)
    qc.rz(seed_gc * np.pi, 2)

    # PAM-proximal base encoding
    if len(guide_seq) >= 20:
        pam_proximal = guide_seq[-3:]
        g_count = pam_proximal.count('G')
        qc.rz(g_count * np.pi / 3, 3)

    # Entangling layer for correlations
    for i in range(n_qubits - 1):
        qc.cx(i, i + 1)

    # Final rotations
    qc.ry(cut_eff * np.pi / 2, 4)

    qc.measure_all()
    return qc


def compute_coherence_from_counts(counts, n_qubits):
    """Compute IR coherence R̄ from measurement counts."""
    total = sum(counts.values())
    phases = []

    for bitstring, count in counts.items():
        # Convert bitstring to phase
        bitstring = bitstring.replace(' ', '')
        value = int(bitstring, 2)
        phase = 2 * np.pi * value / (2 ** n_qubits)
        phases.extend([phase] * count)

    if not phases:
        return 0.5

    z = np.mean(np.exp(1j * np.array(phases)))
    R_bar = np.abs(z)
    return float(R_bar)


# ============================================================================
# TEST SCENARIOS
# ============================================================================

TEST_KNOCKOUT_GUIDES = {
    'high_efficiency': {
        'sequence': 'GCGACTGCTACATAGCCAGG',  # High cut efficiency expected
        'description': 'Optimized knockout guide with G at position 20, good GC',
        'expected_efficiency': 'HIGH',
    },
    'moderate_efficiency': {
        'sequence': 'ATCGATCGATCGATCGATCG',
        'description': 'Moderate guide with even GC distribution',
        'expected_efficiency': 'MODERATE',
    },
    'low_efficiency': {
        'sequence': 'TTTTAAAACCCCGGGGAAAA',
        'description': 'Poor guide with poly-T, extreme GC regions',
        'expected_efficiency': 'LOW',
    },
    'seed_optimized': {
        'sequence': 'AACTGACGGCGCTAGCCGTG',
        'description': 'Guide optimized for seed region (G-rich in positions 13-20)',
        'expected_efficiency': 'HIGH',
    },
    'pam_proximal_optimal': {
        'sequence': 'ACTGATCGATCGATCGATGG',
        'description': 'Guide with GG at PAM-proximal positions',
        'expected_efficiency': 'MODERATE_HIGH',
    },
}


# ============================================================================
# VALIDATION
# ============================================================================

def run_knockout_validation(
    use_hardware: bool = True,
    backend_name: str = "ibm_torino",
    shots: int = 4096,
):
    """
    Run knockout module validation on IBM Quantum hardware.
    """
    print("=" * 60)
    print("E205: CRISPR Knockout Module Validation")
    print("=" * 60)
    print(f"Timestamp: {datetime.now().isoformat()}")
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

    results = {
        'experiment': 'E205_CRISPR_Knockout_Validation',
        'timestamp': datetime.now().isoformat(),
        'backend': backend_name if use_hardware else 'simulator',
        'shots': shots,
        'guides': {},
    }

    print("\n" + "-" * 40)
    print("Testing Knockout Guides")
    print("-" * 40)

    for name, guide_info in TEST_KNOCKOUT_GUIDES.items():
        seq = guide_info['sequence']
        print(f"\n[{name}]")
        print(f"  Sequence: {seq}")
        print(f"  Description: {guide_info['description']}")

        # Classical metrics
        gc = gc_content(seq)
        cut_eff = cut_efficiency_score(seq)
        fs_prob = frameshift_probability(50, 300)  # Example position

        print(f"  GC content: {gc:.1%}")
        print(f"  Cut efficiency: {cut_eff:.3f}")
        print(f"  Frameshift prob: {fs_prob:.3f}")

        # Create and run circuit
        n_qubits = 5
        qc = create_knockout_coherence_circuit(seq, n_qubits)

        # Simulator run
        from qiskit import transpile
        qc_sim = transpile(qc, simulator)
        sim_result = simulator.run(qc_sim, shots=shots).result()
        sim_counts = sim_result.get_counts()
        sim_R = compute_coherence_from_counts(sim_counts, n_qubits)

        print(f"  Simulator R̄: {sim_R:.4f}")

        hw_R = None
        if use_hardware and hardware_backend:
            try:
                qc_hw = transpile(qc, hardware_backend, optimization_level=3)
                sampler = SamplerV2(mode=hardware_backend)
                job = sampler.run([qc_hw], shots=shots)
                hw_result = job.result()

                # Extract counts from SamplerV2 result
                pub_result = hw_result[0]
                hw_counts = {}
                for bits, count in pub_result.data.meas.get_counts().items():
                    hw_counts[bits] = count

                hw_R = compute_coherence_from_counts(hw_counts, n_qubits)
                print(f"  Hardware R̄: {hw_R:.4f}")
            except Exception as e:
                print(f"  Hardware error: {e}")

        # GO/NO-GO classification
        E_MINUS_2 = np.exp(-2)
        R_final = hw_R if hw_R is not None else sim_R
        go_status = "GO" if R_final > E_MINUS_2 else "NO-GO"

        print(f"  Status: {go_status}")

        results['guides'][name] = {
            'sequence': seq,
            'gc': gc,
            'cut_efficiency': cut_eff,
            'frameshift_prob': fs_prob,
            'simulator_R': sim_R,
            'hardware_R': hw_R,
            'go_no_go': go_status,
            'expected': guide_info['expected_efficiency'],
        }

    # Summary
    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)

    all_go = all(g['go_no_go'] == 'GO' for g in results['guides'].values())
    results['all_validated'] = all_go

    for name, data in results['guides'].items():
        R_val = data['hardware_R'] if data['hardware_R'] else data['simulator_R']
        print(f"  {name}: R̄={R_val:.4f} [{data['go_no_go']}]")

    print(f"\nOverall: {'✓ VALIDATED' if all_go else '✗ ISSUES DETECTED'}")

    # Save results
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'Data')
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_file = os.path.join(output_dir, f'E205_knockout_validation_{timestamp}.json')

    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to: {output_file}")

    return results


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-hardware', action='store_true', help='Skip hardware validation')
    parser.add_argument('--backend', default='ibm_torino', help='IBM backend name')
    parser.add_argument('--shots', type=int, default=4096, help='Number of shots')
    args = parser.parse_args()

    run_knockout_validation(
        use_hardware=not args.no_hardware,
        backend_name=args.backend,
        shots=args.shots,
    )
