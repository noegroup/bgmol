"""
Unit and regression test for the ommsystems package.
"""

# Import package, test suite, and other packages as needed
import ommsystems
import pytest
import sys

def test_ommsystems_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "ommsystems" in sys.modules
