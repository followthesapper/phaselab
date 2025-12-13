#!/usr/bin/env python3
"""
Atlas-Q Performance Benchmark Suite.

This script benchmarks the actual Atlas-Q value-add paths to determine
where performance gains are achievable.

The key insight: Atlas-Q only wins when:
1. Processing large expectation value arrays
2. Running many coherence computations in bulk
3. Using optimized backends (GPU/MPS) via ATLAS-Q

Current guide coherence uses heuristic path (Hamiltonian coefficient variance)
which is already fast and doesn't benefit from Atlas-Q acceleration.

Usage:
    python tests/bench_atlas_q.py
"""

import time
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import numpy as np


def benchmark_coherence_from_expectations():
    """
    Benchmark compute_coherence_from_expectations with varying array sizes.

    This is the path that CAN benefit from Atlas-Q optimization.
    """
    from phaselab.quantum.coherence import compute_coherence_from_expectations
    from phaselab.quantum import is_atlas_q_available

    print("\n" + "=" * 70)
    print("BENCHMARK: compute_coherence_from_expectations")
    print("This is the path that can benefit from Atlas-Q acceleration")
    print("=" * 70)

    atlas_available = is_atlas_q_available()
    print(f"\nAtlas-Q Available: {atlas_available}")

    # Test sizes representing realistic workloads
    sizes = [32, 128, 512, 2048, 8192, 32768]
    n_iterations = 100

    results = []

    print(f"\n{'Size':<10} {'Atlas-Q (ms)':<15} {'Native (ms)':<15} {'Speedup':<10} {'Method':<15}")
    print("-" * 65)

    for size in sizes:
        # Generate random expectation values in [-1, 1]
        expectations = np.random.uniform(-1, 1, size)

        # Benchmark Atlas-Q path
        start = time.perf_counter()
        for _ in range(n_iterations):
            result_aq = compute_coherence_from_expectations(expectations, use_atlas_q=True)
        atlas_time = (time.perf_counter() - start) * 1000 / n_iterations

        # Benchmark Native path
        start = time.perf_counter()
        for _ in range(n_iterations):
            result_native = compute_coherence_from_expectations(expectations, use_atlas_q=False)
        native_time = (time.perf_counter() - start) * 1000 / n_iterations

        # Calculate speedup
        speedup = native_time / atlas_time if atlas_time > 0 else 0

        # Verify methods match
        method = result_aq.method

        results.append({
            'size': size,
            'atlas_time_ms': atlas_time,
            'native_time_ms': native_time,
            'speedup': speedup,
            'method': method,
        })

        print(f"{size:<10} {atlas_time:<15.3f} {native_time:<15.3f} {speedup:<10.2f}x {method:<15}")

    # Summary
    print("\n" + "-" * 65)
    if atlas_available:
        avg_speedup = np.mean([r['speedup'] for r in results])
        print(f"Average speedup with Atlas-Q: {avg_speedup:.2f}x")
    else:
        print("Atlas-Q not available - using native fallback")

    return results


def benchmark_guide_coherence():
    """
    Benchmark guide coherence computation (heuristic path).

    This shows that guide coherence doesn't benefit from Atlas-Q
    because it uses Hamiltonian coefficient variance as a proxy.
    """
    from phaselab.crispr.coherence_utils import (
        compute_guide_coherence,
        compute_coherence_batch,
        get_coherence_eligibility_info,
    )

    print("\n" + "=" * 70)
    print("BENCHMARK: compute_guide_coherence (CRISPR guides)")
    print("This path uses HEURISTIC - does NOT benefit from Atlas-Q")
    print("=" * 70)

    # Show eligibility info
    info = get_coherence_eligibility_info()
    print(f"\nMethod: {info['method']}")
    print(f"Atlas-Q Available: {info['atlas_q_available']}")
    print(f"Acceleration Active: {info['acceleration_active']}")
    print(f"Reason: {info['reason']}")

    # Generate random guide sequences
    def random_guide(length=20):
        return ''.join(np.random.choice(['A', 'T', 'C', 'G']) for _ in range(length))

    batch_sizes = [10, 50, 100, 500]

    print(f"\n{'Batch Size':<12} {'Total Time (ms)':<18} {'Per Guide (ms)':<15}")
    print("-" * 45)

    for batch_size in batch_sizes:
        guides = [random_guide() for _ in range(batch_size)]

        start = time.perf_counter()
        results = compute_coherence_batch(guides)
        total_time = (time.perf_counter() - start) * 1000
        per_guide = total_time / batch_size

        print(f"{batch_size:<12} {total_time:<18.2f} {per_guide:<15.3f}")

    print("\nNote: These times are already fast because heuristic coherence")
    print("is computationally cheap. Atlas-Q won't help here.")


def benchmark_where_atlas_q_matters():
    """
    Demonstrate scenarios where Atlas-Q WOULD provide speedup.

    This requires calling the expectation-based path directly.
    """
    from phaselab.quantum.coherence import compute_coherence_from_expectations
    from phaselab.quantum import is_atlas_q_available

    print("\n" + "=" * 70)
    print("WHEN ATLAS-Q MATTERS")
    print("=" * 70)

    print("""
Atlas-Q provides speedup in these scenarios:

1. VQE Optimization
   - When running VQE with many measurement shots
   - Each shot produces expectation values that need coherence

2. Large-Scale Hamiltonian Simulation
   - When simulating Hamiltonians with many terms
   - Coherence from 1000+ expectation values

3. Batch Quantum Circuit Analysis
   - Analyzing multiple quantum circuits
   - Each circuit produces expectation vectors

4. Hardware Validation
   - Processing real IBM Quantum data
   - Large arrays of measurement outcomes

ATLAS-Q does NOT help with:
- Guide coherence (uses heuristic proxy)
- Small expectation arrays (<100 elements)
- Classical CRISPOR integration
- Off-target scoring
""")

    if is_atlas_q_available():
        print("\nAtlas-Q IS available in this environment.")
        print("For maximum benefit, use compute_coherence_from_expectations()")
        print("with large expectation value arrays (2000+ elements).")
    else:
        print("\nAtlas-Q is NOT available in this environment.")
        print("Install with: pip install atlas-q")


def verify_method_selection():
    """
    CI-friendly test that verifies method selection works correctly.
    """
    from phaselab.quantum.coherence import compute_coherence_from_expectations
    from phaselab.quantum import is_atlas_q_available

    print("\n" + "=" * 70)
    print("METHOD SELECTION VERIFICATION")
    print("=" * 70)

    expectations = np.random.uniform(-1, 1, 100)

    # Test with use_atlas_q=True
    result_true = compute_coherence_from_expectations(expectations, use_atlas_q=True)

    # Test with use_atlas_q=False
    result_false = compute_coherence_from_expectations(expectations, use_atlas_q=False)

    print(f"\nuse_atlas_q=True  -> method: {result_true.method}")
    print(f"use_atlas_q=False -> method: {result_false.method}")

    # Verify results are consistent
    r_bar_diff = abs(result_true.R_bar - result_false.R_bar)
    print(f"\nR̄ difference: {r_bar_diff:.6f}")

    # Assertions for CI
    passed = True

    # Method should be 'native' when atlas_q=False
    if result_false.method != "native":
        print(f"FAIL: Expected 'native' when use_atlas_q=False, got {result_false.method}")
        passed = False
    else:
        print("PASS: use_atlas_q=False returns 'native' method")

    # When atlas_q=True, method depends on availability
    if is_atlas_q_available():
        if result_true.method != "atlas_q":
            print(f"FAIL: Expected 'atlas_q' when available, got {result_true.method}")
            passed = False
        else:
            print("PASS: use_atlas_q=True returns 'atlas_q' method when available")
    else:
        if result_true.method != "native":
            print(f"FAIL: Expected 'native' fallback when Atlas-Q unavailable, got {result_true.method}")
            passed = False
        else:
            print("PASS: use_atlas_q=True falls back to 'native' when Atlas-Q unavailable")

    # Results should be consistent (within numerical precision)
    if r_bar_diff > 0.01:
        print(f"FAIL: R̄ values differ too much: {r_bar_diff}")
        passed = False
    else:
        print("PASS: R̄ values are consistent between methods")

    return passed


def main():
    """Run all benchmarks."""
    print("=" * 70)
    print("ATLAS-Q PERFORMANCE BENCHMARK SUITE")
    print("PhaseLab v0.6.1")
    print("=" * 70)

    # Run benchmarks
    benchmark_coherence_from_expectations()
    benchmark_guide_coherence()
    benchmark_where_atlas_q_matters()

    # Run verification
    passed = verify_method_selection()

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    if passed:
        print("\n✓ All method selection tests PASSED")
    else:
        print("\n✗ Some tests FAILED")
        return 1

    print("""
KEY TAKEAWAYS:

1. Guide coherence uses HEURISTIC path
   - Fast by default (~0.05-0.2ms per guide)
   - Does NOT benefit from Atlas-Q

2. Atlas-Q helps with LARGE expectation arrays
   - compute_coherence_from_expectations() with 2000+ elements
   - VQE optimization loops
   - Hardware data processing

3. For CRISPR guide design:
   - Coherence is already fast
   - Bottlenecks are likely in CRISPOR integration or off-target enumeration
   - Consider parallelization for batch guide scoring
""")

    return 0


if __name__ == "__main__":
    sys.exit(main())
