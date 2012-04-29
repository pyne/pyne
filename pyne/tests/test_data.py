"""PyNE nuclear data tests"""
import math

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, assert_in

from pyne import data
import numpy as np

def test_atomic_mass():
    o16 = [15.99491461957, 16.0]
    u235 = [235.043931368, 235.0]
    am242m = [242.059550625, 242.0]

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


if __name__ == "__main__":
    nose.main()

