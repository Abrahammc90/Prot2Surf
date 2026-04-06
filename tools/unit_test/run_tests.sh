#!/bin/bash

################################################################################
#
# Unit Test Runner for Clustering Analysis Tools
#
# This script runs the unit test suite for all Python tools in the clustering
# analysis pipeline using pytest.
#
# Usage:
#     ./run_tests.sh [OPTIONS]
#
# Options:
#     -h, --help              Show this help message
#     -a, --all               Run all tests with coverage report
#     -v, --verbose           Run tests in verbose mode
#     -c, --coverage          Generate HTML coverage report
#     -f, --fast              Run fast tests only (skip slow tests)
#     -k PATTERN              Run tests matching PATTERN
#     --cluster               Run only cluster_kon tests
#     --combine               Run only combine tests
#     --plot-angle            Run only angle plotting tests
#     --plot-dist             Run only distance plotting tests
#     --utilities             Run only utility tests
#     --failfast              Stop on first failure
#
# Examples:
#     ./run_tests.sh                    # Run all tests
#     ./run_tests.sh -v                 # Run tests verbosely
#     ./run_tests.sh -c                 # Generate coverage report
#     ./run_tests.sh --cluster -v       # Run cluster tests verbosely
#     ./run_tests.sh -k "NAM_algorithm" # Run tests matching pattern
#
# Author: Abraham Muñiz-Chicharro
# Version: 1.0
#
################################################################################

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Script configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARENT_DIR="$(dirname "$SCRIPT_DIR")"
VENV_DIR="venv"

# Default options
VERBOSE=false
COVERAGE=false
FAST=false
TEST_PATTERN=""
FAILFAST=false
SPECIFIC_TESTS=""

################################################################################
# Function: Print usage information
################################################################################
print_usage() {
    sed -n '/^# Usage:/,/^# Author:/p' "$0" | sed 's/^# //'
}

################################################################################
# Function: Print error message
################################################################################
print_error() {
    echo -e "${RED}ERROR: $1${NC}" >&2
}

################################################################################
# Function: Print info message
################################################################################
print_info() {
    echo -e "${BLUE}INFO: $1${NC}"
}

################################################################################
# Function: Print success message
################################################################################
print_success() {
    echo -e "${GREEN}SUCCESS: $1${NC}"
}

################################################################################
# Function: Print warning message
################################################################################
print_warning() {
    echo -e "${YELLOW}WARNING: $1${NC}"
}

################################################################################
# Function: Check if pytest is installed
################################################################################
check_pytest() {
    if ! command -v pytest &> /dev/null; then
        print_error "pytest is not installed"
        echo "Install pytest with: pip install pytest pytest-cov"
        return 1
    fi
    return 0
}

################################################################################
# Function: Check if numpy is installed
################################################################################
check_numpy() {
    if ! python3 -c "import numpy" 2>/dev/null; then
        print_warning "numpy is not installed - some tests may fail"
        echo "Install numpy with: pip install numpy"
    fi
}

################################################################################
# Function: Print test environment information
################################################################################
print_environment_info() {
    echo -e "\n${BLUE}========== Test Environment ==========${NC}"
    echo "Python version: $(python3 --version 2>&1)"
    echo "Pytest version: $(pytest --version 2>&1 | head -1)"
    echo "Test directory: $SCRIPT_DIR"
    echo "Working directory: $(pwd)"
    echo ""
}

################################################################################
# Function: Build pytest command
################################################################################
build_pytest_command() {
    local cmd="pytest"
    
    # Add test path
    cmd="$cmd ."
    
    # Add verbose flag
    if [[ "$VERBOSE" == true ]]; then
        cmd="$cmd -v"
    fi
    
    # Add coverage flag
    if [[ "$COVERAGE" == true ]]; then
        cmd="$cmd --cov=.. --cov-report=html --cov-report=term-missing"
    fi
    
    # Add fast mode (skip slow tests)
    if [[ "$FAST" == true ]]; then
        cmd="$cmd -m 'not slow'"
    fi
    
    # Add test pattern matching
    if [[ -n "$TEST_PATTERN" ]]; then
        cmd="$cmd -k '$TEST_PATTERN'"
    fi
    
    # Add specific test file
    if [[ -n "$SPECIFIC_TESTS" ]]; then
        cmd="$cmd $SPECIFIC_TESTS"
    fi
    
    # Add failfast flag
    if [[ "$FAILFAST" == true ]]; then
        cmd="$cmd -x"
    fi
    
    # Add output options
    cmd="$cmd --tb=short"
    
    echo "$cmd"
}

################################################################################
# Function: Run the tests
################################################################################
run_tests() {
    local pytest_cmd=$(build_pytest_command)
    
    print_info "Running tests with command:"
    echo "  $pytest_cmd"
    echo ""
    
    cd "$SCRIPT_DIR"
    
    # Execute pytest command
    if eval "$pytest_cmd"; then
        return 0
    else
        return 1
    fi
}

################################################################################
# Function: Display test summary
################################################################################
display_test_summary() {
    local exit_code="$1"
    if [[ "$exit_code" -eq 0 ]]; then
        print_success "All tests passed!"
        echo ""
        
        if [[ "$COVERAGE" == true ]]; then
            print_info "Coverage report generated in: htmlcov/index.html"
        fi
    else
        print_error "Some tests failed"
        return 1
    fi
}

################################################################################
# Function: Parse command-line arguments
################################################################################
parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                print_usage
                exit 0
                ;;
            -a|--all)
                VERBOSE=true
                COVERAGE=true
                shift
                ;;
            -v|--verbose)
                VERBOSE=true
                shift
                ;;
            -c|--coverage)
                COVERAGE=true
                shift
                ;;
            -f|--fast)
                FAST=true
                shift
                ;;
            -k)
                TEST_PATTERN="$2"
                shift 2
                ;;
            -x|--failfast)
                FAILFAST=true
                shift
                ;;
            --cluster)
                SPECIFIC_TESTS="test_cluster_kon.py"
                shift
                ;;
            --combine)
                SPECIFIC_TESTS="test_combine.py"
                shift
                ;;
            --plot-angle)
                SPECIFIC_TESTS="test_plot_angle.py"
                shift
                ;;
            --plot-dist)
                SPECIFIC_TESTS="test_plot_dist.py"
                shift
                ;;
            --utilities)
                SPECIFIC_TESTS="test_utilities.py"
                shift
                ;;
            *)
                print_error "Unknown option: $1"
                print_usage
                exit 1
                ;;
        esac
    done
}

################################################################################
# Function: Main execution
################################################################################
main() {
    # Print header
    echo ""
    echo -e "${BLUE}╔════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║  Unit Test Runner - Clustering Analysis Tools            ║${NC}"
    echo -e "${BLUE}╚════════════════════════════════════════════════════════════╝${NC}"
    echo ""
    
    # Parse arguments
    parse_arguments "$@"
    
    # Check dependencies
    print_info "Checking dependencies..."
    if ! check_pytest; then
        exit 1
    fi
    check_numpy
    echo ""
    
    # Print environment info
    print_environment_info
    
    # Run tests
    echo -e "${BLUE}========== Running Tests ==========${NC}"
    if run_tests; then
        exit_code=0
    else
        exit_code=1
    fi
    
    echo ""
    echo -e "${BLUE}========== Test Results ==========${NC}"
    display_test_summary "$exit_code"
    
    exit $exit_code
}

# Run main function with all arguments
main "$@"
