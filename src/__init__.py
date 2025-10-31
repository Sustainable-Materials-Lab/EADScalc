"""
EADScalc - Elemental Analysis and XPS Data Calculator for Cellulose

This package provides tools for calculating the degree of substitution (DS) 
for modified cellulose samples based on elemental analysis (CHNS) data and 
XPS (X-ray Photoelectron Spectroscopy) data.

Modules:
    EADScalc: Calculate DS from elemental analysis (CHNS) data
    XPSDScalc: Calculate DS from XPS atomic concentration data
"""

__version__ = "1.1.0"
__author__ = "Sam Eyley"

# Import main CLI functions for convenience
from .EADScalc import cli as eads_cli
from .XPSDScalc import *

__all__ = ['eads_cli']
