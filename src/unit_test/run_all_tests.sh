#!/bin/bash
# Test runner script for clustering program unit tests
# Compiles and runs all unit tests
# Usage: bash run_all_tests.sh or from src: make test
# @version 1.0

set -e  # Exit on first error

echo "========================================="
echo "Building Unit Tests"
echo "========================================="
echo ""

TESTS=(
    "test_maths"
    "test_array"
    "test_threshold"
    "test_pdb"
    "test_clustering"
    "test_assoc"
    "test_read_input"
)

SRC_DIR="../"

# Build all Fortran test binaries before executing them.
if ! make -C "$SRC_DIR" build_fortran_tests >/dev/null; then
    echo "[FAIL] Could not build Fortran test binaries."
    exit 1
fi

echo ""
echo "========================================="
echo "Running All Tests"
echo "========================================="
echo ""

total_passed=0
total_failed=0

# Run all tests
for test_name in "${TESTS[@]}"; do
    echo ""
    echo "Running $test_name..."
    echo "----------"
    
    if ./"$test_name"; then
        total_passed=$((total_passed + 1))
        echo ""
        echo "[OK] $test_name PASSED"
    else
        total_failed=$((total_failed + 1))
        echo ""
        echo "[FAIL] $test_name FAILED"
    fi
done

echo ""
echo "========================================="
echo "OVERALL TEST SUMMARY"
echo "========================================="
echo "Tests passed: $total_passed/${#TESTS[@]}"
echo "Tests failed: $total_failed/${#TESTS[@]}"
echo "========================================="

echo ""
echo "Cleaning generated test files..."
make -C "$SRC_DIR" cleanup_test_artifacts

if [ $total_failed -eq 0 ]; then
    echo "[OK] All tests passed successfully!"
    exit 0
else
    echo "[FAIL] Some tests failed. Review output above."
    exit 1
fi
