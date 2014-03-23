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


def test_q_val():
    assert_equal(data.q_val(110240001), 0.473)
    assert_equal(data.q_val('H1'), 0.0)
    assert_equal(data.q_val(92235), 4.674)

def test_gamma_frac():
    assert_equal(data.gamma_frac('H1'), 0.0)
    assert_equal(data.gamma_frac(92235), 0.036)
    assert_equal(data.gamma_frac(110240001), 0.998)

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

def test_wims_fpyield():
    assert_equal(data.fpyield('Th-232', 'Eu-154'), 2.2000E-13)
    assert_equal(data.fpyield(962440000, 611480001), 1.3800E-06)

def test_nds_fpyield():
    assert_equal(data.fpyield('Th-232', 'Eu-154', 3), 9.6000E-8)
    assert_equal(data.fpyield('Th-232', 'Eu-154', 3, True), 3.8000E-8)

def test_half_life():
    assert_equal(data.half_life('H1'), np.inf)
    assert_equal(data.half_life(922350001), 1560.0)    


def test_decay_const():
    assert_equal(data.decay_const('H1'), 0.0)
    assert_equal(data.decay_const(922350001), np.log(2.0)/1560.0)    

def test_branch_ratio():
    assert_equal(data.branch_ratio('H1', 'H1'), 1.0)
    assert_equal(data.branch_ratio(922350001, 922350000), 1.0)
    assert_equal(data.branch_ratio(922350001, 922360000), 0.0)
    assert_equal(data.branch_ratio(611460000, 621460000), 0.34)

def test_state_energy():
    assert_equal(data.state_energy('H1'), 0.0)
    assert_equal(data.state_energy(922350001), 7.65e-5)


def test_decay_children():
    assert_equal(data.decay_children('H1'), set())
    assert_equal(data.decay_children(922350001), set([922350000]))
    assert_equal(data.decay_children(611460000), set([601460000, 621460000]))
    assert_equal(data.decay_children('O16'), set())
    assert_equal(data.decay_children('80161'), set([60120000, 80160000]))


def test_abundance_by_z_for_soundness():
    for vs in data.abundance_by_z.values():
        if vs:
            assert(abs(1-sum([v[1] for v in vs]))<1e-12)

def test_constants():
    cases = [
        (3.14159265359, data.pi),
        (6.0221415e+23, data.N_A), 
        (1e24, data.barns_per_cm2),
        (1e-24, data.cm2_per_barn),
        (24.0 * 3600.0, data.sec_per_day),
        ]
    for exp, obs in cases:
        yield assert_equal, exp, obs

if __name__ == "__main__":
    nose.runmodule()

