#!/bin/bash
# Test runner script for clustering program unit tests
# Compiles and runs all unit tests
# Usage: bash run_all_tests.sh or from src: make test

set -e  # Exit on first error

echo "========================================="
echo "Building Unit Tests"
echo "========================================="
echo ""

TESTS=(
    "test_maths"
    "test_matrix" 
    "test_threshold"
    "test_pdb"
    "test_clustering"
    "test_assoc"
    "test_read_input"
)

SRC_DIR="../src"

COMPILE_COMMANDS=(
    "gfortran -O0 -Wall -llapack -lblas -fopenmp -march=native -funroll-loops -fno-fast-math -o test_maths test_maths.f90 $SRC_DIR/maths.o $SRC_DIR/mod_pdb.o"
    "gfortran -O0 -Wall -llapack -lblas -fopenmp -march=native -funroll-loops -fno-fast-math -o test_matrix test_matrix.f90 $SRC_DIR/mod_matrix.o $SRC_DIR/maths.o $SRC_DIR/mod_pdb.o"
    "gfortran -O0 -Wall -llapack -lblas -fopenmp -march=native -funroll-loops -fno-fast-math -o test_threshold test_threshold.f90 $SRC_DIR/mod_threshold.o $SRC_DIR/maths.o"
    "gfortran -O0 -Wall -llapack -lblas -fopenmp -march=native -funroll-loops -fno-fast-math -o test_pdb test_pdb.f90 $SRC_DIR/mod_pdb.o"
    "gfortran -O0 -Wall -llapack -lblas -fopenmp -march=native -funroll-loops -fno-fast-math -o test_clustering test_clustering.f90 $SRC_DIR/mod_clust_algorithm.o $SRC_DIR/read_input.o $SRC_DIR/mod_pdb.o $SRC_DIR/mod_assoc.o $SRC_DIR/maths.o"
    "gfortran -O0 -Wall -llapack -lblas -fopenmp -march=native -funroll-loops -fno-fast-math -o test_assoc test_assoc.f90 $SRC_DIR/mod_assoc.o"
    "gfortran -O0 -Wall -llapack -lblas -fopenmp -march=native -funroll-loops -fno-fast-math -o test_read_input test_read_input.f90 $SRC_DIR/read_input.o $SRC_DIR/mod_pdb.o $SRC_DIR/mod_assoc.o"
)

# Compile all tests
for i in "${!TESTS[@]}"; do
    test_name="${TESTS[$i]}"
    compile_cmd="${COMPILE_COMMANDS[$i]}"
    
    echo "Compiling $test_name..."
    if eval "$compile_cmd"; then
        echo "  [OK] Compilation successful"
    else
        echo "  [FAIL] Compilation failed"
        exit 1
    fi
done

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

if [ $total_failed -eq 0 ]; then
    echo "[OK] All tests passed successfully!"
    exit 0
else
    echo "[FAIL] Some tests failed. Review output above."
    exit 1
fi
