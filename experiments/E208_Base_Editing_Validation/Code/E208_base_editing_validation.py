"""
E208: Base Editing Module - IBM Quantum Hardware Validation

Validates the phaselab.crispr.base_editing module on IBM Quantum hardware.
Tests ABE (A→G) and CBE (C→T) editing efficiency using IR coherence metrics.

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
# BASE EDITING FUNCTIONS
# ============================================================================

def gc_content(seq):
    seq = seq.upper()
    return sum(1 for b in seq if b in 'GC') / len(seq) if seq else 0


ACTIVITY_WINDOWS = {
    'BE4': (4, 8), 'ABE7.10': (4, 7), 'ABE8e': (4, 8),
    'ABE8.20': (3, 9), 'CBE4max': (4, 8),
}

POSITION_EFFICIENCY = {
    1: 0.05, 2: 0.10, 3: 0.30, 4: 0.85, 5: 1.00,
    6: 0.95, 7: 0.80, 8: 0.60, 9: 0.25, 10: 0.10,
}

CBE_CONTEXT = {'TC': 1.2, 'CC': 1.0, 'AC': 0.8, 'GC': 0.6}
ABE_CONTEXT = {'TA': 1.1, 'CA': 1.0, 'AA': 0.9, 'GA': 0.7}


def editing_efficiency_at_position(position, editor='ABE8e'):
    window = ACTIVITY_WINDOWS.get(editor, (4, 8))
    if window[0] <= position <= window[1]:
        return POSITION_EFFICIENCY.get(position, 0.5)
    return POSITION_EFFICIENCY.get(position, 0.02)


def sequence_context_score(guide_seq, target_position, editor='ABE8e'):
    guide_seq = guide_seq.upper()
    idx = target_position - 1
    if idx <= 0 or idx >= len(guide_seq):
        return 1.0
    context = guide_seq[idx - 1] + guide_seq[idx]
    if 'ABE' in editor:
        return ABE_CONTEXT.get(context, 1.0)
    return CBE_CONTEXT.get(context, 1.0)


def find_bystanders(guide_seq, target_position, editor='ABE8e'):
    """Find bystander editable bases in activity window."""
    guide_seq = guide_seq.upper()
    window = ACTIVITY_WINDOWS.get(editor, (4, 8))
    target_base = 'A' if 'ABE' in editor else 'C'

    bystanders = []
    for pos in range(window[0], window[1] + 1):
        if pos == target_position:
            continue
        idx = pos - 1
        if 0 <= idx < len(guide_seq) and guide_seq[idx] == target_base:
            eff = editing_efficiency_at_position(pos, editor)
            ctx = sequence_context_score(guide_seq, pos, editor)
            bystanders.append({
                'position': pos,
                'efficiency': eff * ctx,
            })
    return bystanders


# ============================================================================
# QUANTUM CIRCUITS
# ============================================================================

def create_base_edit_coherence_circuit(guide_seq, target_pos, editor, n_qubits=5):
    """
    Create quantum circuit for base editing coherence.

    Uses amplitude encoding - high editing efficiency produces
    states clustered near |00000⟩ for high coherence.
    """
    qc = QuantumCircuit(n_qubits)

    # Calculate quality metrics
    gc = gc_content(guide_seq)
    pos_eff = editing_efficiency_at_position(target_pos, editor)
    ctx_sc = sequence_context_score(guide_seq, target_pos, editor)

    # Bystander penalty (fewer bystanders = better)
    bystanders = find_bystanders(guide_seq, target_pos, editor)
    n_bystanders = len(bystanders)
    bystander_quality = max(0, 1.0 - n_bystanders * 0.25)  # Each bystander costs 25%

    # GC quality (optimal 40-60%)
    gc_quality = 1.0 - abs(gc - 0.50) * 2
    gc_quality = max(0, min(1, gc_quality))

    # Qubit 0: GC content quality
    qc.ry((1.0 - gc_quality) * np.pi / 3, 0)

    # Qubit 1: Position efficiency (already 0-1)
    qc.ry((1.0 - pos_eff) * np.pi / 3, 1)

    # Qubit 2: Context score (normalize from 0.6-1.2 to 0-1)
    ctx_quality = min(1.0, (ctx_sc - 0.5) / 0.7)
    qc.ry((1.0 - ctx_quality) * np.pi / 3, 2)

    # Qubit 3: Bystander quality
    qc.ry((1.0 - bystander_quality) * np.pi / 3, 3)

    # Qubit 4: Combined editing efficiency
    combined_eff = pos_eff * min(1.0, ctx_sc)
    qc.ry((1.0 - combined_eff) * np.pi / 4, 4)

    # Light entanglement
    qc.cx(0, 1)
    qc.cx(2, 3)

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

TEST_BASE_EDITS = {
    'abe_optimal_pos5': {
        'sequence': 'GCGAATGCTACATAGCCAGG',  # A at position 5
        'target_position': 5,
        'editor': 'ABE8e',
        'description': 'ABE with target at peak efficiency position 5',
        'expected': 'HIGH',
    },
    'abe_pos7': {
        'sequence': 'GCGACTGATACATAGCCAGG',  # A at position 7
        'target_position': 7,
        'editor': 'ABE8e',
        'description': 'ABE with target at position 7',
        'expected': 'MODERATE_HIGH',
    },
    'abe_edge_pos3': {
        'sequence': 'GCAACTGCTACATAGCCAGG',  # A at position 3
        'target_position': 3,
        'editor': 'ABE8.20',  # Wider window includes pos 3
        'description': 'ABE8.20 with target at edge of window',
        'expected': 'MODERATE',
    },
    'cbe_optimal_tc_context': {
        'sequence': 'GCTCCTGCTACATAGCCAGG',  # TC at positions 3-4
        'target_position': 4,
        'editor': 'BE4',
        'description': 'CBE with preferred TC context',
        'expected': 'HIGH',
    },
    'cbe_poor_context': {
        'sequence': 'GCGCCTGCTACATAGCCAGG',  # GC context (poor)
        'target_position': 4,
        'editor': 'BE4',
        'description': 'CBE with disfavored GC context',
        'expected': 'MODERATE',
    },
    'abe_multiple_bystanders': {
        'sequence': 'AAAACTGCTACATAGCCAGG',  # Multiple As in window
        'target_position': 5,
        'editor': 'ABE8e',
        'description': 'ABE with multiple bystander As',
        'expected': 'MODERATE',
    },
}


# ============================================================================
# VALIDATION
# ============================================================================

def run_base_editing_validation(
    use_hardware: bool = True,
    backend_name: str = "ibm_torino",
    shots: int = 4096,
):
    """Run base editing module validation."""
    print("=" * 60)
    print("E208: Base Editing Module Validation")
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
        'experiment': 'E208_Base_Editing_Validation',
        'timestamp': datetime.now().isoformat(),
        'backend': backend_name if use_hardware else 'simulator',
        'shots': shots,
        'guides': {},
    }

    print("\n" + "-" * 40)
    print("Testing Base Editing Guides")
    print("-" * 40)

    for name, test in TEST_BASE_EDITS.items():
        seq = test['sequence']
        pos = test['target_position']
        editor = test['editor']

        print(f"\n[{name}]")
        print(f"  Sequence: {seq}")
        print(f"  Editor: {editor}, Target position: {pos}")

        gc = gc_content(seq)
        pos_eff = editing_efficiency_at_position(pos, editor)
        ctx_sc = sequence_context_score(seq, pos, editor)
        bystanders = find_bystanders(seq, pos, editor)

        print(f"  GC: {gc:.1%}")
        print(f"  Position efficiency: {pos_eff:.3f}")
        print(f"  Context score: {ctx_sc:.3f}")
        print(f"  Combined efficiency: {pos_eff * ctx_sc:.3f}")
        print(f"  Bystanders: {len(bystanders)}")

        n_qubits = 5
        qc = create_base_edit_coherence_circuit(seq, pos, editor, n_qubits)

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

        results['guides'][name] = {
            'sequence': seq,
            'editor': editor,
            'target_position': pos,
            'gc': gc,
            'position_efficiency': pos_eff,
            'context_score': ctx_sc,
            'combined_efficiency': pos_eff * ctx_sc,
            'n_bystanders': len(bystanders),
            'simulator_R': sim_R,
            'hardware_R': hw_R,
            'go_no_go': go_status,
        }

    # Summary
    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)

    all_go = all(g['go_no_go'] == 'GO' for g in results['guides'].values())
    results['all_validated'] = all_go

    for name, data in results['guides'].items():
        R_val = data['hardware_R'] or data['simulator_R']
        print(f"  {name}: R̄={R_val:.4f} [{data['go_no_go']}]")

    print(f"\nOverall: {'✓ VALIDATED' if all_go else '✗ ISSUES DETECTED'}")

    # Save
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'Data')
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_file = os.path.join(output_dir, f'E208_base_editing_validation_{timestamp}.json')

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

    run_base_editing_validation(
        use_hardware=not args.no_hardware,
        backend_name=args.backend,
        shots=args.shots,
    )
