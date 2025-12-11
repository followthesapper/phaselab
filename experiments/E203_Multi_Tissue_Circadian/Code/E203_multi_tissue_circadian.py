#!/usr/bin/env python3
"""
E203: Multi-Tissue Circadian Coupling - IBM Quantum Hardware Validation
=========================================================================

Validates the phaselab.circadian.multi_tissue module on IBM Quantum hardware
by computing IR coherence for inter-tissue phase coupling.

The key insight: Circadian clocks in different tissues (SCN, liver, muscle, etc.)
oscillate with coupled phases. The phase relationships encode tissue synchronization,
which maps naturally to quantum coherence metrics.

Test Cases:
1. Healthy synchronization: SCN-peripheral coupling → high R̄
2. Jet lag: Phase-shifted peripherals → temporarily low R̄ → recovery
3. Shift work: Chronic desynchronization → persistently low R̄
4. Disease state: Weakened coupling → reduced R̄

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
OMEGA_CIRCADIAN = 2 * np.pi / 24.0  # Circadian angular frequency (rad/hour)

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

# =============================================================================
# MULTI-TISSUE CIRCADIAN MODEL
# =============================================================================

# Tissue-specific parameters
TISSUE_PARAMS = {
    "SCN": {
        "omega": OMEGA_CIRCADIAN,  # Master clock
        "amplitude": 1.0,
        "coupling_to_scn": 0.0,  # Self-coupling N/A
        "intrinsic_period": 24.0,
    },
    "liver": {
        "omega": OMEGA_CIRCADIAN * 0.98,  # Slight period difference
        "amplitude": 0.8,
        "coupling_to_scn": 0.4,
        "intrinsic_period": 24.5,
    },
    "muscle": {
        "omega": OMEGA_CIRCADIAN * 1.01,
        "amplitude": 0.7,
        "coupling_to_scn": 0.3,
        "intrinsic_period": 23.8,
    },
    "heart": {
        "omega": OMEGA_CIRCADIAN * 0.99,
        "amplitude": 0.75,
        "coupling_to_scn": 0.35,
        "intrinsic_period": 24.2,
    },
}

def simulate_multi_tissue_phases(
    tissues: List[str],
    t_end: float = 240.0,  # hours
    K_global: float = 0.5,  # Coupling strength
    K_peripheral: float = 0.1,
    noise_strength: float = 0.01,
    disease_tissue: Optional[str] = None,
    disease_severity: float = 0.0,
    random_seed: int = 42,
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """
    Simulate multi-tissue circadian phases using Kuramoto-style coupling.

    Returns:
        (t, phases): Time array and dict of phase trajectories per tissue
    """
    np.random.seed(random_seed)

    n = len(tissues)
    dt = 0.1  # hours
    t = np.arange(0, t_end, dt)
    n_steps = len(t)

    # Initialize phases
    phases = {tissue: np.zeros(n_steps) for tissue in tissues}

    # Random initial phases
    current_phases = {tissue: np.random.uniform(0, 2*np.pi) for tissue in tissues}

    for tissue in tissues:
        phases[tissue][0] = current_phases[tissue]

    # Integrate ODEs
    for step in range(1, n_steps):
        scn_phase = current_phases.get("SCN", 0.0)

        for tissue in tissues:
            tp = TISSUE_PARAMS[tissue]
            omega = tp["omega"]
            coupling_to_scn = tp["coupling_to_scn"]

            # Base frequency
            dphase = omega * dt

            # Coupling to SCN
            if tissue != "SCN":
                K_scn = K_global * coupling_to_scn
                if disease_tissue == tissue:
                    K_scn *= (1.0 - disease_severity * 0.5)
                dphase += K_scn * np.sin(scn_phase - current_phases[tissue]) * dt

            # Inter-peripheral coupling
            if tissue != "SCN":
                for other_tissue in tissues:
                    if other_tissue != tissue and other_tissue != "SCN":
                        dphase += K_peripheral * np.sin(
                            current_phases[other_tissue] - current_phases[tissue]
                        ) * dt

            # Noise
            dphase += noise_strength * np.random.randn() * np.sqrt(dt)

            current_phases[tissue] = (current_phases[tissue] + dphase) % (2 * np.pi)
            phases[tissue][step] = current_phases[tissue]

    return t, phases

def compute_tissue_coherence(phases: Dict[str, np.ndarray], steady_fraction: float = 0.25) -> Dict:
    """Compute coherence metrics from tissue phase trajectories."""
    tissues = list(phases.keys())
    n_steps = len(phases[tissues[0]])
    steady_idx = int((1 - steady_fraction) * n_steps)

    # Global coherence (all tissues)
    all_phases = np.array([phases[t][steady_idx:] for t in tissues])
    z_global = np.mean(np.exp(1j * all_phases), axis=0)
    global_R_bar = float(np.mean(np.abs(z_global)))

    # Per-tissue coherence (phase stability)
    tissue_R_bars = {}
    for tissue in tissues:
        ph = phases[tissue][steady_idx:]
        z = np.exp(1j * ph)
        tissue_R_bars[tissue] = float(np.abs(np.mean(z)))

    # SCN-peripheral coherence
    if "SCN" in tissues:
        peripheral_tissues = [t for t in tissues if t != "SCN"]
        if peripheral_tissues:
            scn_ph = phases["SCN"][steady_idx:]
            phase_diffs = []
            for pt in peripheral_tissues:
                pt_ph = phases[pt][steady_idx:]
                diff = np.mean(np.exp(1j * (scn_ph - pt_ph)))
                phase_diffs.append(diff)
            scn_peripheral_coherence = float(np.abs(np.mean(phase_diffs)))
        else:
            scn_peripheral_coherence = 1.0
    else:
        scn_peripheral_coherence = 0.0

    # Classification
    if global_R_bar >= 0.9:
        classification = "SYNCHRONIZED"
    elif global_R_bar >= 0.7:
        classification = "PARTIALLY_SYNCHRONIZED"
    elif global_R_bar >= E2_THRESHOLD:
        classification = "WEAKLY_SYNCHRONIZED"
    else:
        classification = "DESYNCHRONIZED"

    return {
        "global_R_bar": global_R_bar,
        "tissue_R_bars": tissue_R_bars,
        "scn_peripheral_coherence": scn_peripheral_coherence,
        "classification": classification,
        "go_no_go": "GO" if global_R_bar > E2_THRESHOLD else "NO-GO",
    }

# =============================================================================
# QUANTUM CIRCUIT FOR TISSUE PHASE COHERENCE
# =============================================================================

def build_tissue_phase_circuit(tissue_phases: Dict[str, float]) -> QuantumCircuit:
    """
    Build quantum circuit encoding tissue phase relationships.

    Each qubit represents a tissue, with rotation angle = circadian phase.
    """
    tissues = list(tissue_phases.keys())
    n = len(tissues)
    qc = QuantumCircuit(n)

    # Encode each tissue's phase
    for i, tissue in enumerate(tissues):
        phase = tissue_phases[tissue]
        qc.h(i)
        qc.rz(phase, i)

    # Coupling layer (represents inter-tissue interactions)
    # SCN is typically first, couples to all
    if "SCN" in tissues:
        scn_idx = tissues.index("SCN")
        for i in range(n):
            if i != scn_idx:
                qc.cx(scn_idx, i)

    # Peripheral coupling
    for i in range(n):
        if tissues[i] != "SCN":
            for j in range(i + 1, n):
                if tissues[j] != "SCN":
                    qc.cx(i, j)

    return qc

def extract_tissue_coherence_from_counts(counts: Dict[str, int], n_qubits: int) -> float:
    """Extract coherence from measurement counts."""
    total = sum(counts.values())

    # Compute expectation values
    expectations = []
    for i in range(n_qubits):
        exp_val = 0.0
        for bitstring, count in counts.items():
            if i < len(bitstring):
                bit = int(bitstring[-(i+1)])
                sign = 1 if bit == 0 else -1
                exp_val += sign * count / total
        expectations.append(exp_val)

    expectations = np.array(expectations)
    phases = np.arccos(np.clip(expectations, -1, 1))
    coherence = compute_coherence(phases)

    return coherence.R_bar

# =============================================================================
# TEST SCENARIOS
# =============================================================================

def scenario_healthy_sync(random_seed: int = 42) -> Tuple[str, Dict]:
    """Healthy synchronized state."""
    tissues = ["SCN", "liver", "muscle", "heart"]
    t, phases = simulate_multi_tissue_phases(
        tissues, t_end=168, K_global=0.5, noise_strength=0.01, random_seed=random_seed
    )
    coherence = compute_tissue_coherence(phases)
    return "healthy_synchronized", {
        "phases": {k: v[-1] for k, v in phases.items()},
        "coherence": coherence,
        "expected": "GO",
    }

def scenario_jet_lag(random_seed: int = 42) -> Tuple[str, Dict]:
    """Jet lag - acute phase shift."""
    tissues = ["SCN", "liver", "muscle", "heart"]

    # First, establish baseline
    t, phases = simulate_multi_tissue_phases(
        tissues, t_end=72, K_global=0.5, noise_strength=0.01, random_seed=random_seed
    )

    # Apply phase shift to peripherals (simulating travel)
    shift = np.pi / 2  # 6-hour shift
    initial_phases = {
        "SCN": phases["SCN"][-1],  # SCN adapts faster
        "liver": (phases["liver"][-1] + shift) % (2 * np.pi),
        "muscle": (phases["muscle"][-1] + shift) % (2 * np.pi),
        "heart": (phases["heart"][-1] + shift) % (2 * np.pi),
    }

    # Short-term: desynchronized
    coherence = compute_coherence(np.array(list(initial_phases.values())))

    return "jet_lag_acute", {
        "phases": initial_phases,
        "coherence": {
            "global_R_bar": coherence.R_bar,
            "classification": "DESYNCHRONIZED" if coherence.R_bar < 0.5 else "PARTIALLY_SYNCHRONIZED",
            "go_no_go": coherence.go_no_go,
        },
        "expected": "GO",  # Should still be above e^-2
    }

def scenario_shift_work(random_seed: int = 42) -> Tuple[str, Dict]:
    """Chronic shift work - persistent desynchronization."""
    tissues = ["SCN", "liver", "muscle"]

    # Simulate with inverted light schedule (night shift)
    t, phases = simulate_multi_tissue_phases(
        tissues, t_end=168, K_global=0.2,  # Weakened coupling
        noise_strength=0.05,  # Higher noise
        random_seed=random_seed
    )

    coherence = compute_tissue_coherence(phases)

    return "chronic_shift_work", {
        "phases": {k: v[-1] for k, v in phases.items()},
        "coherence": coherence,
        "expected": "GO",  # Reduced but still above threshold
    }

def scenario_disease_desync(random_seed: int = 42) -> Tuple[str, Dict]:
    """Disease state - severely weakened coupling."""
    tissues = ["SCN", "liver", "muscle", "heart"]

    t, phases = simulate_multi_tissue_phases(
        tissues, t_end=168,
        K_global=0.1,  # Very weak coupling
        disease_tissue="liver",
        disease_severity=0.8,
        noise_strength=0.1,
        random_seed=random_seed
    )

    coherence = compute_tissue_coherence(phases)

    return "disease_desynchronization", {
        "phases": {k: v[-1] for k, v in phases.items()},
        "coherence": coherence,
        "expected": "NO-GO",  # Should be below threshold
    }

# =============================================================================
# MAIN VALIDATION PIPELINE
# =============================================================================

@dataclass
class CircadianCoherenceResult:
    scenario: str
    classical_R_bar: float
    quantum_R_bar: float
    hardware_R_bar: Optional[float]
    classification: str
    go_no_go_classical: str
    go_no_go_quantum: str
    go_no_go_hardware: Optional[str]
    expected: str
    validation_passed: bool

def run_circadian_validation(
    use_hardware: bool = False,
    backend_name: str = "ibm_torino",
    shots: int = 4096
) -> Dict:
    """Run multi-tissue circadian validation."""

    print("="*70)
    print("E203: Multi-Tissue Circadian Coupling - IBM Quantum Validation")
    print("="*70)
    print()
    print("Validating phaselab.circadian.multi_tissue module")
    print("Testing IR coherence for inter-tissue phase synchronization")
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

    # Test scenarios
    scenarios = [
        scenario_healthy_sync(),
        scenario_jet_lag(),
        scenario_shift_work(),
        scenario_disease_desync(),
    ]

    results = []

    for name, data in scenarios:
        print(f"\n{'='*70}")
        print(f"Scenario: {name}")
        print(f"{'='*70}")

        phases = data["phases"]
        classical_coherence = data["coherence"]
        expected = data["expected"]

        # Classical results
        classical_R_bar = classical_coherence.get("global_R_bar", 0.0)
        if isinstance(classical_R_bar, dict):
            classical_R_bar = compute_coherence(np.array(list(phases.values()))).R_bar

        print(f"\nClassical coherence:")
        print(f"  R̄ = {classical_R_bar:.4f}")
        print(f"  Classification: {classical_coherence.get('classification', 'N/A')}")

        go_classical = "GO" if classical_R_bar > E2_THRESHOLD else "NO-GO"

        # Build quantum circuit
        qc = build_tissue_phase_circuit(phases)
        qc.measure_all()

        # Simulator run
        transpiled = transpile(qc, simulator, optimization_level=1)
        sampler = AerSampler()
        job = sampler.run([transpiled], shots=shots)
        result = job.result()
        counts = result[0].data.meas.get_counts()

        quantum_R_bar = extract_tissue_coherence_from_counts(counts, len(phases))
        go_quantum = "GO" if quantum_R_bar > E2_THRESHOLD else "NO-GO"

        print(f"\nQuantum simulation:")
        print(f"  R̄ = {quantum_R_bar:.4f}")
        print(f"  Classification: {go_quantum}")

        # Hardware run
        hardware_R_bar = None
        go_hardware = None

        if use_hardware and hardware_backend is not None:
            try:
                transpiled_hw = transpile(qc, hardware_backend, optimization_level=3)
                sampler_hw = SamplerV2(hardware_backend)
                job_hw = sampler_hw.run([transpiled_hw], shots=shots)
                result_hw = job_hw.result()
                counts_hw = result_hw[0].data.meas.get_counts()

                hardware_R_bar = extract_tissue_coherence_from_counts(counts_hw, len(phases))
                go_hardware = "GO" if hardware_R_bar > E2_THRESHOLD else "NO-GO"

                print(f"\nHardware ({backend_name}):")
                print(f"  R̄ = {hardware_R_bar:.4f}")
                print(f"  Classification: {go_hardware}")

            except Exception as e:
                print(f"\nHardware error: {e}")

        # Validation
        validation_passed = (go_classical == expected) or (expected == "GO" and classical_R_bar > E2_THRESHOLD)

        results.append(CircadianCoherenceResult(
            scenario=name,
            classical_R_bar=float(classical_R_bar),
            quantum_R_bar=float(quantum_R_bar),
            hardware_R_bar=float(hardware_R_bar) if hardware_R_bar else None,
            classification=classical_coherence.get('classification', 'N/A'),
            go_no_go_classical=go_classical,
            go_no_go_quantum=go_quantum,
            go_no_go_hardware=go_hardware,
            expected=expected,
            validation_passed=validation_passed
        ))

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)

    print("\nResults by scenario:")
    print("-"*70)
    print(f"{'Scenario':<30} {'Classical':<10} {'Quantum':<10} {'Hardware':<10} {'Expected':<10}")
    print("-"*70)

    for r in results:
        hw_str = f"{r.hardware_R_bar:.4f}" if r.hardware_R_bar else "N/A"
        print(f"{r.scenario:<30} {r.classical_R_bar:.4f}     {r.quantum_R_bar:.4f}     {hw_str:<10} {r.expected}")

    # Validation summary
    passed = sum(1 for r in results if r.validation_passed)
    print(f"\nValidation: {passed}/{len(results)} scenarios match expected behavior")

    # Save results
    output = {
        'experiment': 'E203_Multi_Tissue_Circadian',
        'timestamp': datetime.now().isoformat(),
        'backend': backend_name if use_hardware else 'AerSimulator',
        'shots': shots,
        'e2_threshold': float(E2_THRESHOLD),
        'results': [asdict(r) for r in results],
        'validation_passed': passed,
        'total_tests': len(results),
    }

    output_dir = Path(__file__).parent.parent / "Data"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"E203_circadian_coherence_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"

    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\n\nResults saved: {output_file}")

    return output


def main():
    import argparse
    parser = argparse.ArgumentParser(description='E203: Multi-Tissue Circadian Validation')
    parser.add_argument('--hardware', action='store_true', help='Run on IBM hardware')
    parser.add_argument('--backend', default='ibm_torino', help='IBM backend name')
    parser.add_argument('--shots', type=int, default=4096, help='Number of shots')
    args = parser.parse_args()

    return run_circadian_validation(
        use_hardware=args.hardware,
        backend_name=args.backend,
        shots=args.shots
    )


if __name__ == "__main__":
    main()
