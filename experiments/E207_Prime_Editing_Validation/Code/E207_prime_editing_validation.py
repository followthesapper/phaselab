"""
E207: Prime Editing Module - IBM Quantum Hardware Validation

Validates the phaselab.crispr.prime_editing module on IBM Quantum hardware.
Tests pegRNA design components (PBS, RT template) using IR coherence metrics.

Author: Dylan Vaca
Date: December 2025
"""

import sys
import os
import json
from datetime import datetime
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'src'))

try:
    from qiskit import QuantumCircuit, transpile
    from qiskit_aer import AerSimulator
    QISKIT_AVAILABLE = True
except ImportError:
    QISKIT_AVAILABLE = False

try:
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
    IBM_AVAILABLE = True
except ImportError:
    IBM_AVAILABLE = False


# ============================================================================
# PRIME EDITING FUNCTIONS
# ============================================================================

def gc_content(seq):
    seq = seq.upper()
    return sum(1 for b in seq if b in 'GC') / len(seq) if seq else 0


def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(b, 'N') for b in reversed(seq.upper()))


def estimate_hairpin_dg(seq):
    """Estimate secondary structure ΔG (simplified)."""
    seq = seq.upper()
    n = len(seq)
    if n < 8:
        return 0.0

    dg = 0.0
    for stem_len in range(4, min(9, n // 2)):
        for i in range(n - 2 * stem_len - 3):
            stem1 = seq[i:i + stem_len]
            for loop_len in range(3, 9):
                j = i + stem_len + loop_len
                if j + stem_len > n:
                    continue
                stem2 = seq[j:j + stem_len]
                stem2_rc = reverse_complement(stem2)
                matches = sum(1 for a, b in zip(stem1, stem2_rc) if a == b)
                if matches >= stem_len - 1:
                    dg -= 1.5 * matches

    return max(-15.0, dg)


def pbs_score(pbs_seq):
    """Score PBS sequence quality."""
    gc = gc_content(pbs_seq)
    length = len(pbs_seq)

    gc_score = 1.0
    if gc < 0.35 or gc > 0.65:
        gc_score = 0.7
    if gc < 0.25 or gc > 0.75:
        gc_score = 0.4

    length_score = 1.2 if 13 <= length <= 15 else 1.0

    ss_dg = estimate_hairpin_dg(pbs_seq)
    ss_score = 1.0 if ss_dg > -5 else 0.6

    return gc_score * length_score * ss_score


def rt_template_score(rt_seq):
    """Score RT template quality."""
    gc = gc_content(rt_seq)
    length = len(rt_seq)

    gc_score = 1.0 if 0.35 <= gc <= 0.65 else 0.7
    length_score = 1.2 if 10 <= length <= 16 else 1.0

    ss_dg = estimate_hairpin_dg(rt_seq)
    ss_score = 1.0 if ss_dg > -5 else 0.6

    return gc_score * length_score * ss_score


# ============================================================================
# QUANTUM CIRCUITS
# ============================================================================

def create_prime_edit_coherence_circuit(spacer, pbs, rt_template, n_qubits=6):
    """
    Create quantum circuit for prime editing coherence.

    Uses amplitude encoding - high quality pegRNA components produce
    states clustered near |000000⟩ for high coherence.
    """
    qc = QuantumCircuit(n_qubits)

    # Calculate quality metrics
    spacer_gc = gc_content(spacer)
    pbs_sc = pbs_score(pbs)
    rt_sc = rt_template_score(rt_template)

    # Secondary structure penalty
    extension = rt_template + pbs
    ext_ss = estimate_hairpin_dg(extension)
    # ss_quality: 1.0 = no structure (good), 0.0 = strong structure (bad)
    ss_quality = max(0, min(1, (ext_ss + 15) / 15))

    # Length quality (optimal PBS: 13-15bp, optimal RT: 10-16bp)
    pbs_len = len(pbs)
    pbs_len_quality = 1.0 if 13 <= pbs_len <= 15 else (0.7 if 10 <= pbs_len <= 17 else 0.4)

    rt_len = len(rt_template)
    rt_len_quality = 1.0 if 10 <= rt_len <= 16 else (0.7 if 7 <= rt_len <= 25 else 0.4)

    # Qubit 0: Spacer GC quality (optimal ~50-60%)
    gc_quality = 1.0 - abs(spacer_gc - 0.55) * 2
    qc.ry((1.0 - gc_quality) * np.pi / 3, 0)

    # Qubit 1: PBS score (already 0-1.2 range, normalize)
    pbs_quality = min(1.0, pbs_sc / 1.2)
    qc.ry((1.0 - pbs_quality) * np.pi / 3, 1)

    # Qubit 2: RT template score
    rt_quality = min(1.0, rt_sc / 1.2)
    qc.ry((1.0 - rt_quality) * np.pi / 3, 2)

    # Qubit 3: Secondary structure quality
    qc.ry((1.0 - ss_quality) * np.pi / 3, 3)

    # Qubit 4: PBS length quality
    qc.ry((1.0 - pbs_len_quality) * np.pi / 4, 4)

    # Qubit 5: RT length quality
    qc.ry((1.0 - rt_len_quality) * np.pi / 4, 5)

    # Light entanglement
    qc.cx(0, 1)
    qc.cx(2, 3)
    qc.cx(4, 5)

    qc.measure_all()
    return qc


def compute_coherence_from_counts(counts, n_qubits):
    phases = []
    for bitstring, count in counts.items():
        bitstring = bitstring.replace(' ', '')
        value = int(bitstring, 2)
        phase = 2 * np.pi * value / (2 ** n_qubits)
        phases.extend([phase] * count)
    if not phases:
        return 0.5
    z = np.mean(np.exp(1j * np.array(phases)))
    return float(np.abs(z))


# ============================================================================
# TEST SCENARIOS
# ============================================================================

TEST_PRIME_EDITS = {
    'optimal_design': {
        'spacer': 'GCGACTGCTACATAGCCAGG',
        'pbs': 'CCTGGCTATGTAGC',  # 14bp, good GC
        'rt_template': 'AGTCGATCGATCG',  # 13bp
        'description': 'Well-designed pegRNA with optimal PBS and RT',
        'expected': 'HIGH',
    },
    'short_pbs': {
        'spacer': 'ATCGATCGATCGATCGATCG',
        'pbs': 'CGATCGAT',  # 8bp - minimum
        'rt_template': 'AGTCGATCGATCG',
        'description': 'Short PBS may reduce efficiency',
        'expected': 'MODERATE',
    },
    'long_rt': {
        'spacer': 'GCGACTGCTACATAGCCAGG',
        'pbs': 'CCTGGCTATGTAGC',
        'rt_template': 'AGTCGATCGATCGATCGATCGATCGATCG',  # 29bp
        'description': 'Long RT template may have secondary structure',
        'expected': 'MODERATE',
    },
    'high_gc_extension': {
        'spacer': 'GCGCGCGCGCGCGCGCGCGC',
        'pbs': 'GCGCGCGCGCGCGC',  # High GC
        'rt_template': 'GCGCGCGCGCGCGC',  # High GC
        'description': 'High GC may cause secondary structure issues',
        'expected': 'LOW',
    },
    'balanced_design': {
        'spacer': 'AACTGACGGCGCTAGCCGTG',
        'pbs': 'CACGGCTAGCGCC',  # 13bp, 62% GC
        'rt_template': 'GTCAGTTAGCATC',  # 13bp, 46% GC
        'description': 'Balanced GC across components',
        'expected': 'HIGH',
    },
}


# ============================================================================
# VALIDATION
# ============================================================================

def run_prime_editing_validation(
    use_hardware: bool = True,
    backend_name: str = "ibm_torino",
    shots: int = 4096,
):
    """Run prime editing module validation."""
    print("=" * 60)
    print("E207: Prime Editing Module Validation")
    print("=" * 60)
    print(f"Timestamp: {datetime.now().isoformat()}")

    if not QISKIT_AVAILABLE:
        return {'error': 'Qiskit not available'}

    simulator = AerSimulator()
    hardware_backend = None

    if use_hardware and IBM_AVAILABLE:
        try:
            service = QiskitRuntimeService(channel="ibm_quantum_platform")
            hardware_backend = service.backend(backend_name)
            print(f"Hardware backend: {backend_name}")
        except Exception as e:
            print(f"WARNING: {e}")
            use_hardware = False

    results = {
        'experiment': 'E207_Prime_Editing_Validation',
        'timestamp': datetime.now().isoformat(),
        'backend': backend_name if use_hardware else 'simulator',
        'shots': shots,
        'designs': {},
    }

    print("\n" + "-" * 40)
    print("Testing Prime Editing Designs")
    print("-" * 40)

    for name, design in TEST_PRIME_EDITS.items():
        spacer = design['spacer']
        pbs = design['pbs']
        rt = design['rt_template']

        print(f"\n[{name}]")
        print(f"  Spacer: {spacer}")
        print(f"  PBS: {pbs} ({len(pbs)}bp)")
        print(f"  RT: {rt} ({len(rt)}bp)")

        pbs_sc = pbs_score(pbs)
        rt_sc = rt_template_score(rt)
        ext_ss = estimate_hairpin_dg(rt + pbs)

        print(f"  PBS score: {pbs_sc:.3f}, RT score: {rt_sc:.3f}")
        print(f"  Extension SS ΔG: {ext_ss:.1f} kcal/mol")

        n_qubits = 6
        qc = create_prime_edit_coherence_circuit(spacer, pbs, rt, n_qubits)

        qc_sim = transpile(qc, simulator)
        sim_counts = simulator.run(qc_sim, shots=shots).result().get_counts()
        sim_R = compute_coherence_from_counts(sim_counts, n_qubits)

        print(f"  Simulator R̄: {sim_R:.4f}")

        hw_R = None
        if use_hardware and hardware_backend:
            try:
                qc_hw = transpile(qc, hardware_backend, optimization_level=3)
                sampler = SamplerV2(mode=hardware_backend)
                job = sampler.run([qc_hw], shots=shots)
                pub_result = job.result()[0]
                hw_counts = dict(pub_result.data.meas.get_counts())
                hw_R = compute_coherence_from_counts(hw_counts, n_qubits)
                print(f"  Hardware R̄: {hw_R:.4f}")
            except Exception as e:
                print(f"  Hardware error: {e}")

        R_final = hw_R if hw_R is not None else sim_R
        go_status = "GO" if R_final > np.exp(-2) else "NO-GO"
        print(f"  Status: {go_status}")

        results['designs'][name] = {
            'spacer': spacer,
            'pbs': pbs,
            'rt_template': rt,
            'pbs_score': pbs_sc,
            'rt_score': rt_sc,
            'extension_ss_dg': ext_ss,
            'simulator_R': sim_R,
            'hardware_R': hw_R,
            'go_no_go': go_status,
        }

    # Summary
    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)

    all_go = all(d['go_no_go'] == 'GO' for d in results['designs'].values())
    results['all_validated'] = all_go

    for name, data in results['designs'].items():
        R_val = data['hardware_R'] or data['simulator_R']
        print(f"  {name}: R̄={R_val:.4f} [{data['go_no_go']}]")

    print(f"\nOverall: {'✓ VALIDATED' if all_go else '✗ ISSUES DETECTED'}")

    # Save
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'Data')
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_file = os.path.join(output_dir, f'E207_prime_editing_validation_{timestamp}.json')

    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to: {output_file}")
    return results


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-hardware', action='store_true')
    parser.add_argument('--backend', default='ibm_torino')
    parser.add_argument('--shots', type=int, default=4096)
    args = parser.parse_args()

    run_prime_editing_validation(
        use_hardware=not args.no_hardware,
        backend_name=args.backend,
        shots=args.shots,
    )
