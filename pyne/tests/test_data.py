"""PyNE nuclear data tests"""
import os
import math

import nose 
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, assert_in

import pyne
from pyne import data
import numpy as np

# These tests require nuc_data
if not os.path.isfile(pyne.nuc_data):
    raise RuntimeError("Tests require nuc_data.h5.  Please run nuc_data_make.")

def test_atomic_mass():
    o16 = [15.99491461957, 16.0]
    u235 = [235.043930131, 235.0]
    am242m = [242.059549364, 242.0]

    # zzaam form
    assert_in(data.atomic_mass(80160),  o16)
    assert_in(data.atomic_mass(922350), u235)
    assert_in(data.atomic_mass(952421), am242m)


def test_b_coherent():
    assert_equal(data.b_coherent('H1'), -3.7406E-13 + 0j)
    assert_equal(data.b_coherent(491150), 4.01E-13 - 5.62E-15j)


def test_b_incoherent():
    assert_equal(data.b_incoherent('PD105'), -2.6E-13 + 0j)
    assert_equal(data.b_incoherent(621490), 3.14E-12 - 1.03E-12j)


def test_b():
    bc = data.b_coherent(621490)
    bi = data.b_incoherent('SM149')
    assert_equal(data.b('SM149'), math.sqrt(abs(bc)**2 + abs(bi)**2))


def test_half_life():
    assert_equal(data.half_life('H1'), np.inf)
    assert_equal(data.half_life(922351), 1560.0)    


def test_decay_const():
    assert_equal(data.decay_const('H1'), 0.0)
    assert_equal(data.decay_const(922351), np.log(2.0)/1560.0)    


def test_branch_ratio():
    assert_equal(data.branch_ratio('H1', 'H1'), 1.0)
    assert_equal(data.branch_ratio(922351, 922350), 1.0)
    assert_equal(data.branch_ratio(922351, 922360), 0.0)
    assert_equal(data.branch_ratio(611460, 621460), 0.34)


def test_state_energy():
    assert_equal(data.state_energy('H1'), 0.0)
    assert_equal(data.state_energy(922351), 7.65e-5)


def test_decay_children():
    assert_equal(data.decay_children('H1'), set())
    assert_equal(data.decay_children(922351), set([922350]))
    assert_equal(data.decay_children(611460), set([601460, 621460]))

if __name__ == "__main__":
    nose.runmodule()

