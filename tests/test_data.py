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
    assert_equal(data.decay_children('80166'), set([60120000, 80160000]))


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


def test_metastable_id():
    assert_equal(data.metastable_id(430990000, 1), 430990002)


def test_decay_half_life():
    assert_equal(data.decay_half_life(551370000, 561370000), [(949252608.0, 2840184.0)])


def test_decay_half_life_byparent():
    assert_equal(data.decay_half_life_byparent(551370000), [(949252608.0, 2840184.0)])


def test_decay_branch_ratio_byparent():
    assert_equal(data.decay_branch_ratio_byparent(551370000), [1.0])


def test_decay_photon_branch_ratio_byparent():
    assert_equal(data.decay_photon_branch_ratio_byparent(551370000), [(1.0, 0.0)])


def test_decay_beta_branch_ratio_byparent():
    assert_equal(data.decay_beta_branch_ratio_byparent(551370000), [(1.0, 0.0)])


def test_gamma_energy():
    assert_equal(data.gamma_energy(551370000), [(283.5, 0.1), (661.657, 0.003)])


def test_gamma_photon_intensity():
    assert_equal(data.gamma_photon_intensity(551370000), [(0.00058, 8e-05), (85.1, 0.2)])


def test_gamma_conversion_intensity():
    assert_equal(data.gamma_conversion_intensity(551370000), [(np.nan, np.nan), (0.1124, np.nan)])

def test_gamma_total_intensity():
    assert_equal(data.gamma_total_intensity(561370000), [(100.0, np.nan)]) 

def test_gamma_from_to_byparent():
    assert_equal(data.gamma_from_to_byparent(551370000), [(561370001, 561370000), (561370002, 561370000)])

def test_gamma_from_to_byen():
    assert_equal(data.gamma_from_to_byen(661.65, 0.1),
                 [(390990016, 390990005),
                 (822040062, 822040024),
                 (902290055, 902290000),
                 (400880011, 400880004),
                 (551310023, 551310009),
                 (0, 0),
                 (431070028, 431070020),
                 (972490039, 972490003),
                 (0, 0),
                 (380930068, 380930050),
                 (561370002, 561370000),
                 (561370002, 561370000),
                 (621520071, 621520020),
                 (621540026, 621540006),
                 (781810026, 781810000),
                 (791930069, 791930033),
                 (541390033, 541390028)])


def test_gamma_parent():
    assert_equal(data.gamma_parent(661.65, 0.1),
                 [611510000, 651500000, 380990000, 832040000, 892290000,
                  410880000, 561310000, 771830000, 962480000, 992530000,
                  661550000, 370930000, 551370000, 561370000, 611520000,
                  611540000, 791810000, 801930000, 982520000])


if __name__ == "__main__":
    nose.runmodule()

