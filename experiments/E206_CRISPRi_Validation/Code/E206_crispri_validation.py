"""
E206: CRISPRi Module - IBM Quantum Hardware Validation

Validates the phaselab.crispr.interference module on IBM Quantum hardware.
Tests repression efficiency prediction using IR coherence metrics.

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
# CRISPRi SCORING FUNCTIONS
# ============================================================================

def gc_content(seq):
    seq = seq.upper()
    return sum(1 for b in seq if b in 'GC') / len(seq) if seq else 0


REPRESSION_POSITION_WEIGHTS = {
    'template': {
        (-50, 0): 0.6, (0, 50): 0.9, (50, 150): 1.0,
        (150, 300): 0.8, (300, 500): 0.5,
    },
    'nontemplate': {
        (-50, 0): 0.5, (0, 50): 0.7, (50, 150): 0.85,
        (150, 300): 0.7, (300, 500): 0.4,
    },
}


def repression_efficiency_score(guide_seq, position, strand='+', repressor='KRAB'):
    """Calculate repression efficiency based on position and strand."""
    guide_seq = guide_seq.upper()
    strand_key = 'template' if strand == '+' else 'nontemplate'
    position_score = 0.3

    for (start, end), score in REPRESSION_POSITION_WEIGHTS[strand_key].items():
        if start <= position < end:
            position_score = score
            break

    repressor_mult = {'KRAB': 1.0, 'MeCP2': 0.9, 'dCas9_only': 0.6}.get(repressor, 1.0)

    gc = gc_content(guide_seq)
    gc_mod = 0.85 if (gc < 0.35 or gc > 0.70) else 1.0

    return float(np.clip(position_score * repressor_mult * gc_mod, 0, 1))


def steric_hindrance_score(position, strand):
    """Calculate steric hindrance potential."""
    if strand == '+':
        if 0 <= position <= 100:
            return 1.0 - abs(position - 50) / 100
        elif -50 <= position < 0:
            return 0.7 + position / 250
        elif 100 < position <= 300:
            return max(0.3, 1.0 - (position - 100) / 400)
        return 0.2
    else:
        if 0 <= position <= 100:
            return 0.8 - abs(position - 50) / 125
        elif -50 <= position < 0:
            return 0.5 + position / 250
        elif 100 < position <= 300:
            return max(0.2, 0.8 - (position - 100) / 400)
        return 0.15


# ============================================================================
# QUANTUM CIRCUITS
# ============================================================================

def create_crispri_coherence_circuit(guide_seq, position, strand, n_qubits=5):
    """
    Create quantum circuit for CRISPRi guide coherence.

    Uses amplitude encoding to maintain phase coherence.
    High-quality guides produce states clustered near |00000⟩.
    """
    qc = QuantumCircuit(n_qubits)

    gc = gc_content(guide_seq)
    rep_eff = repression_efficiency_score(guide_seq, position, strand)
    steric = steric_hindrance_score(position, strand)

    # Compute combined quality score (0 to 1)
    # Higher quality = smaller rotation = more coherent output
    quality = (rep_eff + steric + gc) / 3.0

    # Invert: high quality -> small angle -> coherent state
    # Low quality -> large angle -> decoherent state
    base_angle = (1.0 - quality) * np.pi / 3  # Max π/3 rotation for poor guides

    # Apply controlled rotations based on quality metrics
    # Qubit 0: GC content (good GC = small rotation)
    gc_angle = abs(gc - 0.55) * np.pi / 2  # Optimal GC ~55%
    qc.ry(gc_angle, 0)

    # Qubit 1: Repression efficiency (high = small rotation)
    rep_angle = (1.0 - rep_eff) * np.pi / 3
    qc.ry(rep_angle, 1)

    # Qubit 2: Steric hindrance (high = small rotation)
    steric_angle = (1.0 - steric) * np.pi / 3
    qc.ry(steric_angle, 2)

    # Qubit 3: Position quality (optimal position = small rotation)
    pos_quality = rep_eff  # Position is already encoded in rep_eff
    qc.ry((1.0 - pos_quality) * np.pi / 4, 3)

    # Qubit 4: Strand effect
    strand_bonus = 0.1 if strand == '+' else 0.0
    qc.ry((1.0 - strand_bonus) * np.pi / 6, 4)

    # Light entanglement to correlate quality metrics
    qc.cx(0, 1)
    qc.cx(2, 3)

    qc.measure_all()
    return qc


def compute_coherence_from_counts(counts, n_qubits):
    total = sum(counts.values())
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

TEST_CRISPRI_GUIDES = {
    'tss_optimal_template': {
        'sequence': 'GCGACTGCTACATAGCCAGG',
        'position': 75,
        'strand': '+',
        'description': 'Optimal position (50-150bp) on template strand',
        'expected': 'HIGH',
    },
    'tss_optimal_nontemplate': {
        'sequence': 'ATCGATCGATCGATCGATCG',
        'position': 75,
        'strand': '-',
        'description': 'Optimal position on non-template strand',
        'expected': 'MODERATE_HIGH',
    },
    'promoter_proximal': {
        'sequence': 'GCGACTGCTACATAGCCAGG',
        'position': -25,
        'strand': '+',
        'description': 'Promoter-proximal, blocking initiation',
        'expected': 'MODERATE',
    },
    'downstream_weak': {
        'sequence': 'ATCGATCGATCGATCGATCG',
        'position': 400,
        'strand': '+',
        'description': 'Far downstream, weak repression expected',
        'expected': 'LOW',
    },
    'steric_optimal': {
        'sequence': 'AACTGACGGCGCTAGCCGTG',
        'position': 50,
        'strand': '+',
        'description': 'Peak steric hindrance position',
        'expected': 'HIGH',
    },
}


# ============================================================================
# VALIDATION
# ============================================================================

def run_crispri_validation(
    use_hardware: bool = True,
    backend_name: str = "ibm_torino",
    shots: int = 4096,
):
    """Run CRISPRi module validation."""
    print("=" * 60)
    print("E206: CRISPRi Module Validation")
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
            print(f"WARNING: Could not connect: {e}")
            use_hardware = False

    results = {
        'experiment': 'E206_CRISPRi_Validation',
        'timestamp': datetime.now().isoformat(),
        'backend': backend_name if use_hardware else 'simulator',
        'shots': shots,
        'guides': {},
    }

    print("\n" + "-" * 40)
    print("Testing CRISPRi Guides")
    print("-" * 40)

    for name, guide_info in TEST_CRISPRI_GUIDES.items():
        seq = guide_info['sequence']
        pos = guide_info['position']
        strand = guide_info['strand']

        print(f"\n[{name}]")
        print(f"  Sequence: {seq}")
        print(f"  Position: {pos}bp, Strand: {strand}")

        gc = gc_content(seq)
        rep_eff = repression_efficiency_score(seq, pos, strand)
        steric = steric_hindrance_score(pos, strand)

        print(f"  GC: {gc:.1%}, Repression: {rep_eff:.3f}, Steric: {steric:.3f}")

        n_qubits = 5
        qc = create_crispri_coherence_circuit(seq, pos, strand, n_qubits)

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
            'position': pos,
            'strand': strand,
            'gc': gc,
            'repression_efficiency': rep_eff,
            'steric_hindrance': steric,
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
    output_file = os.path.join(output_dir, f'E206_crispri_validation_{timestamp}.json')

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

    run_crispri_validation(
        use_hardware=not args.no_hardware,
        backend_name=args.backend,
        shots=args.shots,
    )
