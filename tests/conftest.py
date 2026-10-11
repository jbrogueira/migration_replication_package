# tests/conftest.py
# The package is run from its root (config.py at top level, relative output paths):
# put the package root on sys.path and chdir there before any model import.
import os
import sys

PKG_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PKG_ROOT not in sys.path:
    sys.path.insert(0, PKG_ROOT)
os.chdir(PKG_ROOT)
