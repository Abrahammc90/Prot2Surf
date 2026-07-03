"""
Unit tests for regioselectivity.py - Main function and argument parsing.

@version 1.0
"""
import pytest
import sys
from regioselectivity import main

def test_main_runs(monkeypatch):
    # Simulate command-line arguments for main()
    monkeypatch.setattr(sys, 'argv', ['regioselectivity.py', '-complexes', 'file1', '-references', 'ref1'])
    try:
        main()
    except SystemExit:
        pass  # main() may call sys.exit()
    # This test checks that main() runs without crashing for minimal args
