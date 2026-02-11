#!/bin/bash

################################################################################
#
# Quick Test Runner - Shortcut for running all tests
#
# This script is a quick wrapper around run_tests.sh with common use cases.
#
# Usage:
#     ./test.sh [quick|verbose|coverage|pattern]
#
# Commands:
#     (no args)   - Run all tests quickly
#     quick       - Run tests without verbose output
#     verbose     - Run tests with verbose output
#     coverage    - Run tests with coverage report
#     fast        - Run only fast tests (skip marked slow tests)
#     pattern     - Interactive pattern selection
#
# Examples:
#     ./test.sh                    # Quick test run
#     ./test.sh verbose            # Verbose test run
#     ./test.sh coverage           # Generate coverage report
#
# Author: Abraham Muñiz-Chicharro
# Version: 1.0
#
################################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_TESTS="$SCRIPT_DIR/run_tests.sh"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

# Check if run_tests.sh exists
if [[ ! -f "$RUN_TESTS" ]]; then
    echo -e "${RED}ERROR: run_tests.sh not found${NC}" >&2
    exit 1
fi

# Determine which test command to run
case "${1:-}" in
    quick)
        "$RUN_TESTS"
        ;;
    verbose|-v)
        "$RUN_TESTS" -v
        ;;
    coverage|-c)
        "$RUN_TESTS" -c
        ;;
    all|-a)
        "$RUN_TESTS" -a
        ;;
    fast|-f)
        "$RUN_TESTS" -f
        ;;
    cluster)
        "$RUN_TESTS" --cluster -v
        ;;
    combine)
        "$RUN_TESTS" --combine -v
        ;;
    plot-angle)
        "$RUN_TESTS" --plot-angle -v
        ;;
    plot-dist)
        "$RUN_TESTS" --plot-dist -v
        ;;
    utilities)
        "$RUN_TESTS" --utilities -v
        ;;
    help|-h|--help)
        sed -n '/^# Usage:/,/^# Author:/p' "$0" | sed 's/^# //'
        exit 0
        ;;
    *)
        # Default: run tests normally
        "$RUN_TESTS"
        ;;
esac
