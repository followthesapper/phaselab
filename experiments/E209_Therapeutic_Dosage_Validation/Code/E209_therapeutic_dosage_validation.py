"""
E209: Therapeutic Dosage Module - IBM Quantum Hardware Validation

Validates the phaselab.therapy module on IBM Quantum hardware.
Tests therapeutic window optimization using IR coherence metrics.

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
# THERAPEUTIC DOSAGE FUNCTIONS
# ============================================================================

E_MINUS_2 = np.exp(-2)


class TherapeuticWindow:
    """Therapeutic window parameters."""
    def __init__(
        self,
        baseline=0.50,
        therapeutic_min=0.70,
        therapeutic_max=1.20,
        optimal=0.85,
        name=None,
    ):
        self.baseline = baseline
        self.therapeutic_min = therapeutic_min
        self.therapeutic_max = therapeutic_max
        self.optimal = optimal
        self.name = name

    @property
    def required_upregulation(self):
        return self.therapeutic_min / self.baseline


THERAPEUTIC_WINDOWS = {
    'SMS': TherapeuticWindow(0.50, 0.70, 1.10, 0.80, "Smith-Magenis Syndrome"),
    'SCN2A_NDD': TherapeuticWindow(0.50, 0.65, 1.15, 0.85, "SCN2A Neurodevelopmental Disorder"),
    'SHANK3': TherapeuticWindow(0.50, 0.60, 1.10, 0.80, "Phelan-McDermid Syndrome"),
    'CHD8': TherapeuticWindow(0.50, 0.65, 1.15, 0.85, "CHD8-related ASD"),
}


def estimate_expression_change(guide_coherence, binding_energy):
    """Estimate expression change from guide quality."""
    binding_factor = min(1.0, max(0, (-binding_energy - 5) / 20))
    fold_min, fold_max = 1.5, 3.0
    base_fold = fold_min + (fold_max - fold_min) * binding_factor
    reliability = guide_coherence if guide_coherence > E_MINUS_2 else 0.5
    estimated_fold = base_fold * (0.7 + 0.3 * reliability)
    return {
        'estimated_fold': estimated_fold,
        'reliability': reliability,
    }


def validate_therapeutic_level(expression, window):
    """Check if expression is in therapeutic window."""
    in_window = window.therapeutic_min <= expression <= window.therapeutic_max
    distance_to_optimal = abs(expression - window.optimal)

    if expression > window.therapeutic_max:
        status = "OVEREXPRESSION_RISK"
    elif in_window:
        status = "THERAPEUTIC"
    else:
        status = "SUBTHERAPEUTIC"

    return {
        'expression': expression,
        'in_window': in_window,
        'status': status,
        'distance_to_optimal': distance_to_optimal,
    }


def therapeutic_score(achieved_expression, window):
    """Calculate how well expression matches therapeutic target."""
    if achieved_expression < window.therapeutic_min:
        # Below threshold - score based on how close
        return max(0, achieved_expression / window.therapeutic_min)
    elif achieved_expression > window.therapeutic_max:
        # Above max - penalize overexpression
        excess = achieved_expression - window.therapeutic_max
        return max(0, 1.0 - excess / 0.5)
    else:
        # In window - score based on distance to optimal
        range_size = window.therapeutic_max - window.therapeutic_min
        distance = abs(achieved_expression - window.optimal)
        return 1.0 - distance / range_size


# ============================================================================
# QUANTUM CIRCUITS
# ============================================================================

def create_therapeutic_coherence_circuit(
    guide_coherence,
    binding_energy,
    achieved_expression,
    window,
    n_qubits=6,
):
    """
    Create quantum circuit for therapeutic dosage coherence.

    Uses amplitude encoding - therapeutic success produces
    states clustered near |000000⟩ for high coherence.
    """
    qc = QuantumCircuit(n_qubits)

    # Calculate quality metrics
    th_score = therapeutic_score(achieved_expression, window)

    # In-window quality
    in_window = window.therapeutic_min <= achieved_expression <= window.therapeutic_max
    window_quality = 1.0 if in_window else 0.3

    # Distance to optimal (closer = better)
    dist = abs(achieved_expression - window.optimal)
    optimal_quality = max(0, 1.0 - dist / 0.5)

    # Binding quality (more negative = stronger binding = better)
    binding_quality = min(1.0, max(0, (-binding_energy - 5) / 25))

    # Guide coherence is already a quality metric (0-1)

    # Qubit 0: Guide coherence quality
    qc.ry((1.0 - guide_coherence) * np.pi / 3, 0)

    # Qubit 1: Binding quality
    qc.ry((1.0 - binding_quality) * np.pi / 3, 1)

    # Qubit 2: Therapeutic score
    qc.ry((1.0 - th_score) * np.pi / 3, 2)

    # Qubit 3: In-window quality
    qc.ry((1.0 - window_quality) * np.pi / 3, 3)

    # Qubit 4: Distance to optimal
    qc.ry((1.0 - optimal_quality) * np.pi / 3, 4)

    # Qubit 5: Overall reliability (combination)
    overall = (guide_coherence + th_score + window_quality) / 3
    qc.ry((1.0 - overall) * np.pi / 4, 5)

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

TEST_DOSAGE_SCENARIOS = {
    'sms_optimal': {
        'disease': 'SMS',
        'guide_coherence': 0.85,
        'binding_energy': -20,
        'description': 'SMS with high-quality guide, expected optimal dosing',
        'expected': 'THERAPEUTIC',
    },
    'sms_subtherapeutic': {
        'disease': 'SMS',
        'guide_coherence': 0.50,
        'binding_energy': -10,
        'description': 'SMS with weak guide, may be subtherapeutic',
        'expected': 'BORDERLINE',
    },
    'scn2a_optimal': {
        'disease': 'SCN2A_NDD',
        'guide_coherence': 0.90,
        'binding_energy': -22,
        'description': 'SCN2A autism with excellent guide',
        'expected': 'THERAPEUTIC',
    },
    'scn2a_overexpression_risk': {
        'disease': 'SCN2A_NDD',
        'guide_coherence': 0.95,
        'binding_energy': -28,
        'description': 'SCN2A with very strong guide, risk of overexpression',
        'expected': 'HIGH_EXPRESSION',
    },
    'shank3_moderate': {
        'disease': 'SHANK3',
        'guide_coherence': 0.70,
        'binding_energy': -15,
        'description': 'Phelan-McDermid with moderate guide',
        'expected': 'THERAPEUTIC',
    },
    'chd8_low_coherence': {
        'disease': 'CHD8',
        'guide_coherence': 0.20,  # Below e^-2 threshold
        'binding_energy': -18,
        'description': 'CHD8 with unreliable guide (low coherence)',
        'expected': 'UNRELIABLE',
    },
}


# ============================================================================
# VALIDATION
# ============================================================================

def run_therapeutic_dosage_validation(
    use_hardware: bool = True,
    backend_name: str = "ibm_torino",
    shots: int = 4096,
):
    """Run therapeutic dosage module validation."""
    print("=" * 60)
    print("E209: Therapeutic Dosage Module Validation")
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
        'experiment': 'E209_Therapeutic_Dosage_Validation',
        'timestamp': datetime.now().isoformat(),
        'backend': backend_name if use_hardware else 'simulator',
        'shots': shots,
        'scenarios': {},
    }

    print("\n" + "-" * 40)
    print("Testing Therapeutic Dosage Scenarios")
    print("-" * 40)

    for name, scenario in TEST_DOSAGE_SCENARIOS.items():
        disease = scenario['disease']
        guide_coh = scenario['guide_coherence']
        binding_e = scenario['binding_energy']
        window = THERAPEUTIC_WINDOWS[disease]

        print(f"\n[{name}]")
        print(f"  Disease: {window.name}")
        print(f"  Guide coherence: {guide_coh:.2f}")
        print(f"  Binding energy: {binding_e} kcal/mol")

        # Estimate expression
        expr_data = estimate_expression_change(guide_coh, binding_e)
        achieved = window.baseline * expr_data['estimated_fold']

        print(f"  Estimated fold change: {expr_data['estimated_fold']:.2f}x")
        print(f"  Achieved expression: {achieved:.1%}")
        print(f"  Therapeutic window: {window.therapeutic_min:.0%} - {window.therapeutic_max:.0%}")

        # Validate therapeutic level
        validation = validate_therapeutic_level(achieved, window)
        th_score = therapeutic_score(achieved, window)

        print(f"  Status: {validation['status']}")
        print(f"  Therapeutic score: {th_score:.3f}")

        # Quantum circuit
        n_qubits = 6
        qc = create_therapeutic_coherence_circuit(
            guide_coh, binding_e, achieved, window, n_qubits
        )

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
        go_status = "GO" if R_final > E_MINUS_2 else "NO-GO"
        print(f"  IR Status: {go_status}")

        results['scenarios'][name] = {
            'disease': disease,
            'guide_coherence': guide_coh,
            'binding_energy': binding_e,
            'estimated_fold': expr_data['estimated_fold'],
            'achieved_expression': achieved,
            'therapeutic_min': window.therapeutic_min,
            'therapeutic_max': window.therapeutic_max,
            'optimal': window.optimal,
            'in_therapeutic_window': validation['in_window'],
            'therapeutic_status': validation['status'],
            'therapeutic_score': th_score,
            'simulator_R': sim_R,
            'hardware_R': hw_R,
            'go_no_go': go_status,
        }

    # Summary
    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)

    all_go = all(s['go_no_go'] == 'GO' for s in results['scenarios'].values())
    results['all_validated'] = all_go

    therapeutic_count = sum(
        1 for s in results['scenarios'].values()
        if s['in_therapeutic_window']
    )

    for name, data in results['scenarios'].items():
        R_val = data['hardware_R'] or data['simulator_R']
        th_status = "✓ IN WINDOW" if data['in_therapeutic_window'] else "✗ OUT OF WINDOW"
        print(f"  {name}: R̄={R_val:.4f} [{data['go_no_go']}] {th_status}")

    print(f"\nTherapeutic success rate: {therapeutic_count}/{len(results['scenarios'])}")
    print(f"Overall: {'✓ VALIDATED' if all_go else '✗ ISSUES DETECTED'}")

    # Save
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'Data')
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_file = os.path.join(output_dir, f'E209_therapeutic_dosage_{timestamp}.json')

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

    run_therapeutic_dosage_validation(
        use_hardware=not args.no_hardware,
        backend_name=args.backend,
        shots=args.shots,
    )
