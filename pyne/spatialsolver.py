"""This module accesses the AHOTN, DGFEM and SCT-STEP flavored nuetron transport solvers.
"""
import re
import sys
from warnings import warn
from pyne.utils import VnVWarning

import numpy as np

if sys.version_info[0] > 2:
  basestring = str

warn(__name__ + " is not yet V&V compliant.", VnVWarning)

#from spatialsolver.solver import solver as internalsolver
import spatial.solver

def test():
	print("TEST")

def solve(inputdict_unchecked): 
				spatial.solver.solve(inputdict_unchecked)
        #spatialsolver.solver.solve(inputdict_unchecked)
