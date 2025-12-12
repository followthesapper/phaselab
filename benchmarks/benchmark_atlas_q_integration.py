#!/usr/bin/env python3
"""
PhaseLab v0.6.0 Benchmark: ATLAS-Q Integration Performance

This benchmark compares:
1. Coherence calculation: Heuristic vs ATLAS-Q circular stats
2. Measurement grouping: Simple vs IR-enhanced grouping
3. VQE performance: Simple fallback vs ATLAS-Q backend

Run with: python benchmarks/benchmark_atlas_q_integration.py
"""

import json
import time
from dataclasses import asdict, dataclass
from datetime import datetime
from typing import Dict, List, Tuple

import numpy as np

# Check what's available
print("=" * 60)
print("PhaseLab v0.6.0 ATLAS-Q Integration Benchmark")
print("=" * 60)

try:
    import phaselab
    print(f"✓ PhaseLab version: {phaselab.__version__}")
except ImportError:
    print("✗ PhaseLab not installed")
    exit(1)

try:
    import atlas_q
    ATLAS_Q_AVAILABLE = True
    print(f"✓ ATLAS-Q version: {atlas_q.__version__}")
except ImportError:
    ATLAS_Q_AVAILABLE = False
    print("✗ ATLAS-Q not installed - will test fallback paths only")

try:
    import torch
    GPU_AVAILABLE = torch.cuda.is_available()
    if GPU_AVAILABLE:
        print(f"✓ GPU available: {torch.cuda.get_device_name(0)}")
    else:
        print("○ No GPU available - CPU only")
except ImportError:
    GPU_AVAILABLE = False
    print("○ PyTorch not available")

print("=" * 60)


@dataclass
class BenchmarkResult:
    """Result from a single benchmark."""
    name: str
    method: str
    time_ms: float
    iterations: int
    value: float
    extra: Dict = None


def benchmark_coherence_calculation():
    """Benchmark coherence calculation methods."""
    print("\n--- Benchmark 1: Coherence Calculation ---")

    from phaselab.core.coherence import coherence_score
    from phaselab.quantum.coherence import (
        compute_coherence_from_expectations,
        compute_coherence_from_hamiltonian,
    )

    results = []

    # Test data: simulated Pauli expectations
    np.random.seed(42)
    n_terms = 100
    expectations = 0.8 + 0.15 * np.random.randn(n_terms)
    expectations = np.clip(expectations, -1, 1)

    # Benchmark 1a: Original heuristic (from Hamiltonian structure)
    coefficients = np.random.randn(n_terms)

    n_iter = 1000
    start = time.perf_counter()
    for _ in range(n_iter):
        result = compute_coherence_from_hamiltonian(coefficients, use_atlas_q=False)
    elapsed = (time.perf_counter() - start) * 1000

    results.append(BenchmarkResult(
        name="coherence_heuristic",
        method="native_heuristic",
        time_ms=elapsed / n_iter,
        iterations=n_iter,
        value=result.R_bar,
        extra={"V_phi": result.V_phi, "is_go": result.is_go}
    ))
    print(f"  Heuristic (from coefficients): {elapsed/n_iter:.4f} ms, R̄={result.R_bar:.4f}")

    # Benchmark 1b: Real circular stats (native)
    start = time.perf_counter()
    for _ in range(n_iter):
        result = compute_coherence_from_expectations(expectations, use_atlas_q=False)
    elapsed = (time.perf_counter() - start) * 1000

    results.append(BenchmarkResult(
        name="coherence_circular_native",
        method="native_circular",
        time_ms=elapsed / n_iter,
        iterations=n_iter,
        value=result.R_bar,
        extra={"V_phi": result.V_phi, "is_go": result.is_go}
    ))
    print(f"  Circular stats (native):        {elapsed/n_iter:.4f} ms, R̄={result.R_bar:.4f}")

    # Benchmark 1c: ATLAS-Q circular stats (if available)
    if ATLAS_Q_AVAILABLE:
        start = time.perf_counter()
        for _ in range(n_iter):
            result = compute_coherence_from_expectations(expectations, use_atlas_q=True)
        elapsed = (time.perf_counter() - start) * 1000

        results.append(BenchmarkResult(
            name="coherence_circular_atlas",
            method="atlas_q",
            time_ms=elapsed / n_iter,
            iterations=n_iter,
            value=result.R_bar,
            extra={"V_phi": result.V_phi, "is_go": result.is_go}
        ))
        print(f"  Circular stats (ATLAS-Q):       {elapsed/n_iter:.4f} ms, R̄={result.R_bar:.4f}")

    return results


def benchmark_measurement_grouping():
    """Benchmark Hamiltonian grouping methods."""
    print("\n--- Benchmark 2: Measurement Grouping ---")

    from phaselab.quantum.grouping import ir_hamiltonian_grouping

    results = []

    # Test data: mock Hamiltonian
    np.random.seed(42)
    n_terms = 50
    n_qubits = 8

    coefficients = np.random.randn(n_terms)
    pauli_ops = ['I', 'X', 'Y', 'Z']
    pauli_strings = [
        ''.join(np.random.choice(pauli_ops, n_qubits))
        for _ in range(n_terms)
    ]

    n_iter = 100
    total_shots = 10000

    # Benchmark 2a: Simple grouping (fallback)
    start = time.perf_counter()
    for _ in range(n_iter):
        result = ir_hamiltonian_grouping(
            coefficients, pauli_strings,
            total_shots=total_shots,
            use_atlas_q=False
        )
    elapsed = (time.perf_counter() - start) * 1000

    results.append(BenchmarkResult(
        name="grouping_simple",
        method="simple_commuting",
        time_ms=elapsed / n_iter,
        iterations=n_iter,
        value=result.variance_reduction,
        extra={
            "n_groups": len(result.groups),
            "method": result.method,
        }
    ))
    print(f"  Simple grouping:    {elapsed/n_iter:.4f} ms, {len(result.groups)} groups, {result.variance_reduction:.2f}× reduction")

    # Benchmark 2b: IR-enhanced grouping (ATLAS-Q)
    if ATLAS_Q_AVAILABLE:
        start = time.perf_counter()
        for _ in range(n_iter):
            result = ir_hamiltonian_grouping(
                coefficients, pauli_strings,
                total_shots=total_shots,
                use_atlas_q=True
            )
        elapsed = (time.perf_counter() - start) * 1000

        results.append(BenchmarkResult(
            name="grouping_ir_enhanced",
            method="atlas_q_vra",
            time_ms=elapsed / n_iter,
            iterations=n_iter,
            value=result.variance_reduction,
            extra={
                "n_groups": len(result.groups),
                "method": result.method,
            }
        ))
        print(f"  IR-enhanced (ATLAS-Q): {elapsed/n_iter:.4f} ms, {len(result.groups)} groups, {result.variance_reduction:.2f}× reduction")

    return results


def benchmark_vqe():
    """Benchmark VQE optimization methods."""
    print("\n--- Benchmark 3: VQE Optimization ---")

    from phaselab.quantum.vqe import run_vqe, VQEConfig

    results = []

    # Test Hamiltonian
    np.random.seed(42)
    n_qubits = 4
    n_terms = 10

    coefficients = np.random.randn(n_terms)
    pauli_ops = ['I', 'X', 'Y', 'Z']
    hamiltonian_terms = [
        (coefficients[i], ''.join(np.random.choice(pauli_ops, n_qubits)))
        for i in range(n_terms)
    ]

    config = VQEConfig(
        n_layers=2,
        max_iterations=20,
        verbose=False,
    )

    n_iter = 5  # VQE is slow

    # Benchmark 3a: Simple VQE (fallback)
    start = time.perf_counter()
    for _ in range(n_iter):
        result = run_vqe(hamiltonian_terms, n_qubits, config, use_atlas_q=False)
    elapsed = (time.perf_counter() - start) * 1000

    results.append(BenchmarkResult(
        name="vqe_simple",
        method="simple",
        time_ms=elapsed / n_iter,
        iterations=n_iter,
        value=result.energy,
        extra={
            "converged": result.converged,
            "n_iterations": result.n_iterations,
            "coherence_R": result.coherence.R_bar if result.coherence else None,
        }
    ))
    print(f"  Simple VQE:     {elapsed/n_iter:.1f} ms, E={result.energy:.4f}, R̄={result.coherence.R_bar if result.coherence else 'N/A':.4f}")

    # Benchmark 3b: ATLAS-Q VQE (if available)
    if ATLAS_Q_AVAILABLE:
        start = time.perf_counter()
        for _ in range(n_iter):
            result = run_vqe(hamiltonian_terms, n_qubits, config, use_atlas_q=True)
        elapsed = (time.perf_counter() - start) * 1000

        results.append(BenchmarkResult(
            name="vqe_atlas",
            method="atlas_q",
            time_ms=elapsed / n_iter,
            iterations=n_iter,
            value=result.energy,
            extra={
                "converged": result.converged,
                "n_iterations": result.n_iterations,
                "coherence_R": result.coherence.R_bar if result.coherence else None,
            }
        ))
        print(f"  ATLAS-Q VQE:    {elapsed/n_iter:.1f} ms, E={result.energy:.4f}, R̄={result.coherence.R_bar if result.coherence else 'N/A':.4f}")

    return results


def benchmark_backend_selection():
    """Benchmark backend auto-selection."""
    print("\n--- Benchmark 4: Backend Selection ---")

    from phaselab.quantum.backend import get_available_backends, get_optimal_backend, BackendType

    results = []

    backends = get_available_backends()
    print(f"  Available backends: {len(backends)}")
    for b in backends:
        print(f"    - {b.name}: {b.description}")

    # Test backend selection
    test_cases = [
        (4, False, "small_circuit"),
        (15, False, "medium_circuit"),
        (30, False, "large_circuit"),
        (10, True, "clifford_circuit"),
    ]

    for n_qubits, is_clifford, name in test_cases:
        backend = get_optimal_backend(n_qubits, is_clifford)
        results.append(BenchmarkResult(
            name=f"backend_selection_{name}",
            method="auto",
            time_ms=0,
            iterations=1,
            value=0,
            extra={
                "n_qubits": n_qubits,
                "is_clifford": is_clifford,
                "selected_backend": backend.value,
            }
        ))
        print(f"  {name} ({n_qubits}q, clifford={is_clifford}): {backend.value}")

    return results


def benchmark_gpu():
    """Benchmark GPU operations (if available)."""
    print("\n--- Benchmark 5: GPU Acceleration ---")

    from phaselab.quantum.gpu import check_gpu, batch_coherence_gpu

    results = []

    gpu_info = check_gpu()
    print(f"  GPU available: {gpu_info.available}")
    if gpu_info.available:
        print(f"  Device: {gpu_info.device_name}")
        print(f"  Memory: {gpu_info.memory_total / 1e9:.1f} GB total, {gpu_info.memory_free / 1e9:.1f} GB free")

    # Benchmark batch coherence
    np.random.seed(42)
    n_batches = 100
    n_terms = 50

    batches = [
        0.8 + 0.15 * np.random.randn(n_terms)
        for _ in range(n_batches)
    ]
    batches = [np.clip(b, -1, 1) for b in batches]

    n_iter = 10

    start = time.perf_counter()
    for _ in range(n_iter):
        R_bars = batch_coherence_gpu(batches)
    elapsed = (time.perf_counter() - start) * 1000

    results.append(BenchmarkResult(
        name="batch_coherence",
        method="gpu" if gpu_info.available else "cpu",
        time_ms=elapsed / n_iter,
        iterations=n_iter,
        value=np.mean(R_bars),
        extra={
            "n_batches": n_batches,
            "gpu_available": gpu_info.available,
        }
    ))
    print(f"  Batch coherence ({n_batches} batches): {elapsed/n_iter:.2f} ms, mean R̄={np.mean(R_bars):.4f}")

    return results


def main():
    """Run all benchmarks and save results."""
    all_results = []

    # Run benchmarks
    all_results.extend(benchmark_coherence_calculation())
    all_results.extend(benchmark_measurement_grouping())
    all_results.extend(benchmark_vqe())
    all_results.extend(benchmark_backend_selection())
    all_results.extend(benchmark_gpu())

    # Summary
    print("\n" + "=" * 60)
    print("BENCHMARK SUMMARY")
    print("=" * 60)

    def make_json_serializable(obj):
        """Convert numpy types to Python native types."""
        if isinstance(obj, dict):
            return {k: make_json_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, (list, tuple)):
            return [make_json_serializable(v) for v in obj]
        elif isinstance(obj, np.bool_):
            return bool(obj)
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    summary = {
        "timestamp": datetime.now().isoformat(),
        "phaselab_version": phaselab.__version__,
        "atlas_q_available": bool(ATLAS_Q_AVAILABLE),
        "atlas_q_version": atlas_q.__version__ if ATLAS_Q_AVAILABLE else None,
        "gpu_available": bool(GPU_AVAILABLE),
        "results": [
            make_json_serializable({
                "name": r.name,
                "method": r.method,
                "time_ms": r.time_ms,
                "value": r.value,
                "extra": r.extra,
            })
            for r in all_results
        ]
    }

    # Key comparisons
    print("\nKey Performance Metrics:")

    # Coherence comparison
    heuristic = next((r for r in all_results if r.name == "coherence_heuristic"), None)
    circular = next((r for r in all_results if r.name == "coherence_circular_native"), None)
    atlas = next((r for r in all_results if r.name == "coherence_circular_atlas"), None)

    if heuristic and circular:
        print(f"\n  Coherence Accuracy:")
        print(f"    Heuristic R̄:     {heuristic.value:.4f}")
        print(f"    Circular R̄:      {circular.value:.4f}")
        if atlas:
            print(f"    ATLAS-Q R̄:       {atlas.value:.4f}")

    # Grouping comparison
    simple = next((r for r in all_results if r.name == "grouping_simple"), None)
    enhanced = next((r for r in all_results if r.name == "grouping_ir_enhanced"), None)

    if simple and enhanced:
        improvement = enhanced.value / simple.value if simple.value > 0 else 0
        print(f"\n  Measurement Grouping:")
        print(f"    Simple variance reduction:  {simple.value:.2f}×")
        print(f"    IR-enhanced reduction:      {enhanced.value:.2f}×")
        print(f"    Improvement factor:         {improvement:.2f}×")

    # Save results
    output_file = "benchmarks/benchmark_results.json"
    with open(output_file, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\n✓ Results saved to {output_file}")
    print("=" * 60)

    return summary


if __name__ == "__main__":
    main()
