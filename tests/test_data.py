"""PyNE nuclear data tests"""
import os
import math
import warnings

import nose
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, assert_in
import numpy as np
import numpy.testing as npt

from pyne.utils import VnVWarning

warnings.simplefilter("ignore", VnVWarning)

import pyne
from pyne import data


# These tests require nuc_data
if not os.path.isfile(pyne.nuc_data):
    raise RuntimeError("Tests require nuc_data.h5.  Please run nuc_data_make.")


def test_atomic_mass():
    o16 = [15.99491461957, 16.0]
    u235 = [235.043930131, 235.0]
    am242m = [242.059549364, 242.0]

    # zzaam form
    assert_in(data.atomic_mass(80160), o16)
    assert_in(data.atomic_mass(922350), u235)
    assert_in(data.atomic_mass(952421), am242m)


def test_q_val():
    assert_equal(data.q_val(110240001), 0.473)
    assert_equal(data.q_val('H1'), 0.0)
    assert_equal(data.q_val(92235), 4.674)


def test_simple_xs():
    assert_equal(data.simple_xs(922350000, 'tot', 'thermal'), 698.2)
    assert_equal(data.simple_xs('u235', 'elastic', 'thermal'), 15.04)
    assert_equal(data.simple_xs(922350000, b'gamma', 'thermal'), 98.81)
    assert_equal(data.simple_xs(922350000, 'fission', 'thermal'), 584.4)
    assert_equal(data.simple_xs(922350000, 'tot', 'thermal'), 698.2)

    assert_equal(data.simple_xs(922350000, 'tot', 'thermal_maxwell_ave'), 608.4)
    assert_equal(data.simple_xs(922350000, 'absorption', 'resonance_integral'), 411.1)
    assert_equal(data.simple_xs(922350000, 'tot', 'fourteen_MeV'), 5.865)
    assert_equal(data.simple_xs(922350000, 'tot', 'fission_spectrum_ave'), 7.705)


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
    assert_equal(data.b('SM149'), math.sqrt(abs(bc) ** 2 + abs(bi) ** 2))


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
    assert_equal(data.decay_const(922350001), np.log(2.0) / 1560.0)


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
    assert_equal(data.decay_children('80166', False), set([60120000, 80160000]))


def test_abundance_by_z_for_soundness():
    for vs in data.abundance_by_z.values():
        if vs:
            assert (abs(1 - sum([v[1] for v in vs])) < 1e-12)


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
    assert_equal(data.decay_half_life(551370000, 561370000), (949252608.0,
                                                              2840184.0))


def test_decay_half_life_byparent():
    assert_equal(data.decay_half_life_byparent(551370000), [(949252608.0,
                                                             2840184.0)])


def test_decay_branch_ratio():
    assert_equal(data.decay_branch_ratio(551370000, 561370000), 1.0)


def test_decay_photon_branch_ratio():
    npt.assert_array_almost_equal(
        data.decay_photon_branch_ratio(551370000, 561370000), (1.0, np.nan))


def test_decay_beta_branch_ratio():
    npt.assert_array_almost_equal(
        data.decay_beta_branch_ratio(551370000, 561370000), (1.0, np.nan))


def test_decay_branch_ratio_byparent():
    assert_equal(data.decay_branch_ratio_byparent(551370000), [1.0])


def test_decay_photon_branch_ratio_byparent():
    npt.assert_array_almost_equal(
        data.decay_photon_branch_ratio_byparent(551370000), [(1.0, np.nan)])


def test_decay_beta_branch_ratio_byparent():
    npt.assert_array_almost_equal(
        data.decay_beta_branch_ratio_byparent(551370000), [(1.0, np.nan)])


def test_gamma_energy():
    assert_equal(data.gamma_energy(551370000), [(283.5, 0.1), (661.657, 0.003)])


def test_gamma_photon_intensity():
    assert_equal(data.gamma_photon_intensity(551370000), [(0.00058, 8e-05),
                                                          (85.1, 0.2)])


def test_gamma_conversion_intensity():
    npt.assert_array_almost_equal(data.gamma_conversion_intensity(551370000),
                                  [(np.nan, np.nan), (0.1124, np.nan)])


def test_gamma_total_intensity():
    npt.assert_array_almost_equal(data.gamma_total_intensity(561370000),
                                  [(100.0, np.nan)])


def test_gamma_from_to_byparent():
    assert_equal(data.gamma_from_to_byparent(551370000),
                 [(561370001, 561370000), (561370002, 561370000)])


def test_gamma_from_to_byen():
    assert_equal(data.gamma_from_to_byen(661.65, 0.1),
                 [(621510087, 621510015),
                  (641500021, 641500006),
                  (390990016, 390990005),
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


def test_alpha_energy():
    assert_equal(data.alpha_energy(952410000),
                 [4758.0, 4800.0, 4834.0, 4889.0, 4956.0, 4962.0, 4964.0,
                  5004.0, 5055.0, 5068.0, 5089.0, 5096.0, 5106.0, 5114.0,
                  5137.0, 5155.0, 5178.0, 5182.0, 5192.0, 5217.0, 5223.0,
                  5232.0, 5244.0, 5279.0, 5322.0, 5388.0, 5416.5, 5442.8,
                  5469.0, 5485.56, 5511.5, 5544.5])


def test_alpha_intensity():
    npt.assert_array_almost_equal(data.alpha_intensity(952410000),
                                  [1e-05, 8.6e-05, 0.0007, np.nan, np.nan,
                                   np.nan, np.nan,
                                   0.0001, np.nan, 0.00014, 0.0004, 0.0004,
                                   np.nan, 0.0004,
                                   0.00032, 0.0007, 0.0003, 0.0009, 0.0006,
                                   1e-05, 0.0013,
                                   np.nan, 0.0024, 0.0005, 0.015, 1.66, 0.01,
                                   13.1, 0.04, 84.8,
                                   0.225, 0.37])


def test_alpha_parent():
    assert_equal(data.alpha_parent(5322.0, 0.1), [952410000, 972490000])


def test_alpha_child_byen():
    assert_equal(data.alpha_child_byen(5322.0, 0.1), [932370008, 952450000])


def test_alpha_child_byparent():
    assert_equal(data.alpha_child_byparent(952410000),
                 [932370047, 932370043, 932370041, 932370038, 932370034,
                  932370033, 932370032, 932370031, 932370028, 932370027,
                  932370026, 932370024, 932370023, 932370022, 932370021,
                  932370020, 932370019, 932370018, 932370017, 932370015,
                  932370014, 932370013, 932370012, 932370009, 932370008,
                  932370006, 932370005, 932370004, 932370003, 932370002,
                  932370001, 932370000])


def test_beta_endpoint_energy():
    npt.assert_array_almost_equal(data.beta_endpoint_energy(551370000),
                                  [514.03, np.nan, 1176.])


def test_beta_average_energy():
    assert_equal(data.beta_average_energy(551370000), [174.32, 334.65, 416.26])


def test_beta_intensity():
    assert_equal(data.beta_intensity(551370000), [94.7, 0.00058, 5.3])


def test_beta_parent():
    assert_equal(data.beta_parent(1000, 0.1), [310760000, 410990001, 441070000])


def test_beta_child_byen():
    assert_equal(data.beta_child_byen(1000, 0.1),
                 [320760084, 420990050, 451070007])


def test_beta_child_byparent():
    assert_equal(data.beta_child_byparent(551370000),
                 [561370002, 561370001, 561370000])


def test_ecbp_endpoint_energy():
    assert_equal(data.ecbp_endpoint_energy(170320000)[-1], 11600.0)


def test_ecbp_average_energy():
    assert_equal(data.ecbp_average_energy(110220000), [215.54, 835.0])


def test_ec_intensity():
    npt.assert_array_almost_equal(data.ec_intensity(110220000), [9.618, np.nan])


def test_beta_plus_intensity():
    assert_equal(data.beta_plus_intensity(110220000), [90.326, 0.056])


def test_ecbp_parent():
    assert_equal(data.ecbp_parent(215.54, 0.5),
                 [110220000, 541230000, 571330000])


def test_ecbp_child_byen():
    assert_equal(data.ecbp_child_byen(215.54, 0.5),
                 [100220001, 531230020, 561330006])


def test_ecbp_child_byparent():
    assert_equal(data.ecbp_child_byparent(110220000), [100220001, 100220000])


def test_id_from_level():
    assert_equal(data.id_from_level(811920000, 445, 'X'), 811920010)
    assert_equal(data.id_from_level(561370000, 662), 561370002)


def test_xray_data():
    npt.assert_almost_equal(data.calculate_xray_data(551370000, 0.1, 0.1),
                            [(30.9728, 0.04693557573523976),
                             (30.6251, 0.025406227145485287),
                             (35.0, 0.017058197119274962),
                             (4.29, 0.019708)])


def test_gamma_xray():
    npt.assert_almost_equal(data.gamma_xrays(551370000),
                            [[(32.1936, 0.0), (31.8171, 0.0),
                              (36.4, 0.0), (4.47, 0.0)],
                             [(32.1936, 23759153.763765693),
                              (31.8171, 12896468.662972018),
                              (36.4, 8749697.07326229),
                              (4.47, 5927514.212400001)]])


def test_ecbp_xray():
    npt.assert_almost_equal(data.ecbp_xrays(110220000),
                            [[(1.041, 0.11273075771609505),
                              (1.041, 0.05669229805542418),
                              (1.07, 0.0014231536684807533),
                              (np.nan, 0.0)],
                             [(1.041, 1.273768717251104e-05),
                              (1.041, 6.4057828790558e-06),
                              (1.07, 1.6080514843315978e-07),
                              (np.nan, 0.0)]])


def test_gamma_photon_intensity_byen():
    npt.assert_almost_equal(data.gamma_photon_intensity_byen(661.657, 0.05),
                            [(0.08, 0.017),
                             (16.0, 3.0),
                             (85.1, 0.2),
                             (89.9, 0.14),
                             (1.5, 0.1),
                             (0.14, np.nan),
                             (160.0, 24.0),
                             (0.32, 0.1),
                             (5.0, np.nan)])


if __name__ == "__main__":
    nose.runmodule()

