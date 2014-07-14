"""This module accesses the AHOTN, DGFEM and SCT-STEP flavored nuetron transport solvers.
"""
from __future__ import division

import sys
import os

#Solver imports
#sys.path.append("../fortran/spatial_solvers_3d/source")
#sys.path.append("../../fortran/spatial_solvers_3d/source")

#from spatialsolvers.solver import solve as internalsolver
import spatialsolvers.solver

def solve(inputdict_unchecked):
        solve(inputdict_unchecked)
        print("Solver called!")
