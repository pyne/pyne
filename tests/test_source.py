"""Source tests"""
from pyne.source import PointSource
import os
import warnings

from unittest import TestCase
import nose

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in

from pyne.utils import QAWarning
warnings.simplefilter("ignore", QAWarning)


###############################################################################
# tests the constructor
def test_point_source():
    pt_src = PointSource(1, 2, 3, 4, 5, 6, 12, "p", 0.5)
    assert_equal(pt_src.x, 1)
    assert_equal(pt_src.y, 2)
    assert_equal(pt_src.z, 3)
    assert_equal(pt_src.i, 4)
    assert_equal(pt_src.j, 5)
    assert_equal(pt_src.k, 6)
    assert_equal(pt_src.particle, "p")
    assert_equal(pt_src.weight, 0.5)


################################################################################
# test write particle for mcnp
def test_point_source_mcnp():
    pt_src = PointSource(1, 2, 3, 4, 5, 6, 12, "Proton", 0.5)
    assert_equal(
        "SDEF POS=1 2 3\n     VEC=4 5 6 DIR=1\n     ERG=12\n     WGT=0.5\n     PAR=h", pt_src.mcnp(6))
    
    pt_src = PointSource(1, 2, 3, 4, 5, 6, 12, "Proton")
    assert_equal(
        "SDEF POS=1 2 3\n     VEC=4 5 6 DIR=1\n     ERG=12\n     PAR=h", pt_src.mcnp(6))

    pt_src = PointSource(1, 2, 3, 4, 5, 6, 12)
    assert_equal(
        "SDEF POS=1 2 3\n     VEC=4 5 6 DIR=1\n     ERG=12", pt_src.mcnp(6))
    
    pt_src = PointSource(1, 2, 3, 4, 5, 6)
    assert_equal(
        "SDEF POS=1 2 3\n     VEC=4 5 6 DIR=1", pt_src.mcnp(6))
    
    pt_src = PointSource(1, 2, 3)
    assert_equal(
        "SDEF POS=1 2 3", pt_src.mcnp(6))

# Run as script
#
if __name__ == "__main__":
    nose.runmodule()
