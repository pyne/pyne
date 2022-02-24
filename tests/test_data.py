"""PyNE nuclear data tests"""
import os
import math
import warnings

import nose
from nose.tools import assert_equal, assert_in, assert_true
import numpy as np
import numpy.testing as npt

from pyne.utils import QAWarning

warnings.simplefilter("ignore", QAWarning)

import pyne
from pyne import data, nucname

from pyne import utils

if utils.use_warnings():
    utils.toggle_warnings()

# These tests require nuc_data
if not os.path.isfile(pyne.nuc_data):
    raise RuntimeError("Tests require nuc_data.h5.  Please run nuc_data_make.")


def test_atomic_mass():
    o16 = [15.9949146196, 16.0]
    u235 = [235.04392819, 235.0]
    am242m = [242.059547428, 242.0]

    # zzaam form
    assert_in(data.atomic_mass(80160), o16)
    assert_in(data.atomic_mass(922350), u235)
    assert_in(data.atomic_mass(952421), am242m)


def test_natural_abund_excited_state():
    # initialize natural_abund_map
    gnd = 902320000
    excited = gnd + 1
    data.natural_abund(gnd)
    # excited state should not be in the map yet
    assert_equal(data.natural_abund_map.get(excited), None)
    nabund = data.natural_abund(excited)
    assert_equal(nabund, data.natural_abund_map.get(excited))


def test_q_val():
    assert_equal(data.q_val(110240001), 0.473)
    assert_equal(data.q_val("H1"), 0.0)
    assert_equal(data.q_val(92235), 4.674)


def test_simple_xs():
    assert_equal(data.simple_xs(922350000, "tot", "thermal"), 698.2)
    assert_equal(data.simple_xs("u235", "elastic", "thermal"), 15.04)
    assert_equal(data.simple_xs(922350000, b"gamma", "thermal"), 98.81)
    assert_equal(data.simple_xs(922350000, "fission", "thermal"), 584.4)
    assert_equal(data.simple_xs(922350000, "tot", "thermal"), 698.2)

    assert_equal(data.simple_xs(922350000, "tot", "thermal_maxwell_ave"), 608.4)
    assert_equal(data.simple_xs(922350000, "absorption", "resonance_integral"), 411.1)
    assert_equal(data.simple_xs(922350000, "tot", "fourteen_MeV"), 5.865)
    assert_equal(data.simple_xs(922350000, "tot", "fission_spectrum_ave"), 7.705)


def test_gamma_frac():
    assert_equal(data.gamma_frac("H1"), 0.0)
    assert_equal(data.gamma_frac(92235), 0.036)
    assert_equal(data.gamma_frac(110240001), 0.998)


def test_ext_air_dose():
    assert_equal(data.ext_air_dose(40100000), 1.49e-10)
    assert_equal(data.ext_air_dose("H3", 0), 4.41e-12)
    assert_true(math.isnan(data.ext_air_dose(25054, 1)))


def test_ext_soil_dose():
    assert_equal(data.ext_soil_dose(40100000, 0), 0.537)
    assert_equal(data.ext_soil_dose("H3", 2), 3.49e-8)
    assert_equal(data.ext_soil_dose(25054, 1), 9590.0)


def test_ingest_dose():
    assert_equal(data.ingest_dose(40100000), 4.66e-6)
    assert_equal(data.ingest_dose("H3", 2), 6.12e-8)
    assert_equal(data.ingest_dose(25054, 1), 2.7e-6)


def test_inhale_dose():
    assert_equal(data.inhale_dose(40100000), 0.000354)
    assert_equal(data.inhale_dose("H3", 2), 9.02e-8)
    assert_equal(data.inhale_dose(25054, 1), 6.4e-6)


def test_b_coherent():
    assert_equal(data.b_coherent("H1"), -3.7406e-13 + 0j)
    assert_equal(data.b_coherent(491150), 4.01e-13 - 5.62e-15j)


def test_b_incoherent():
    assert_equal(data.b_incoherent("PD105"), -2.6e-13 + 0j)
    assert_equal(data.b_incoherent(621490), 3.14e-12 - 1.03e-12j)


def test_b():
    bc = data.b_coherent(621490)
    bi = data.b_incoherent("SM149")
    assert_equal(data.b("SM149"), math.sqrt(abs(bc) ** 2 + abs(bi) ** 2))


def test_wims_fpyield():
    assert_equal(data.fpyield("Th-232", "Eu-154"), 2.2000e-13)
    assert_equal(data.fpyield(962440000, 611480001), 1.3800e-06)


def test_nds_fpyield():
    assert_equal(data.fpyield("Th-232", "Eu-154", 3), 2.79e-07)
    assert_equal(data.fpyield("Th-232", "Eu-154", 3, True), 9.3e-08)


def test_half_life():
    assert_equal(data.half_life("H1"), np.inf)
    assert_equal(data.half_life(922350001), 1560.0)
    assert_equal(data.half_life("Eu151"), np.inf)


def test_decay_const():
    assert_equal(data.decay_const("H1"), 0.0)
    assert_equal(data.decay_const(922350001), np.log(2.0) / 1560.0)


def test_branch_ratio():
    assert_equal(data.branch_ratio("H1", "H1"), 1.0)
    assert_equal(data.branch_ratio(922350001, 922350000), 1.0)
    assert_equal(data.branch_ratio(922350001, 922360000), 0.0)
    assert_equal(data.branch_ratio(611460000, 621460000), 0.34299999999999997)

    children = data.decay_children("U235")
    for child in children:
        obs = data.branch_ratio("U235", child)
        assert_true(obs >= 0.0 and obs <= 1.0)

    # There was a bug with metastable ids being dropped prematurely,
    # which would then lead to the branch ratio for the ground state being reported
    # Obviously, this was bad, so here is a test.
    assert_equal(data.branch_ratio("Se86", "Br86M"), 0.0)

    # Not all isomeric transitions have a 100% branch ratio
    assert_equal(data.branch_ratio(932400001, 932400000), 0.0012)


def test_state_energy():
    assert_equal(data.state_energy("H1"), 0.0)
    assert_equal(data.state_energy(922350001), 7.6e-5)


def test_decay_children():
    assert_equal(data.decay_children("H1"), set())
    assert_equal(data.decay_children(922350001), set([922350000]))
    assert_equal(data.decay_children(611460000), set([601460000, 621460000]))
    assert_equal(data.decay_children("O16"), set())
    assert_equal(data.decay_children("80166", False), set([60120000, 80160000]))
    # Spontaneous fission case
    assert_equal(
        data.decay_children("U-235"),
        set(
            [
                360830000,
                420950000,
                430990000,
                441010000,
                441030000,
                441060000,
                451030000,
                451050000,
                461050000,
                461070000,
                461080000,
                471090000,
                481130000,
                491150000,
                511250000,
                521270000,
                531270000,
                531350000,
                541310000,
                541340000,
                541350000,
                541360000,
                551330000,
                551340000,
                551350000,
                551370000,
                601430000,
                601450000,
                611470000,
                611480000,
                611480001,
                611490000,
                621470000,
                621480000,
                621490000,
                621500000,
                621510000,
                621520000,
                631510000,
                631520000,
                631530000,
                631540000,
                631550000,
                641540000,
                641550000,
                641560000,
                641570000,
                641580000,
                661600000,
                661610000,
                661620000,
                661630000,
                661640000,
                671650000,
                681660000,
                681670000,
                902310000,
            ]
        ),
    )


def test_abundance_by_z_for_soundness():
    for vs in data.abundance_by_z.values():
        if vs:
            assert abs(1 - sum([v[1] for v in vs])) < 1e-12


def test_constants():
    cases = [
        (3.14159265359, data.pi),
        (6.0221415e23, data.N_A),
        (1e24, data.barns_per_cm2),
        (1e-24, data.cm2_per_barn),
        (24.0 * 3600.0, data.sec_per_day),
    ]
    for exp, obs in cases:
        yield assert_equal, exp, obs


def test_metastable_id():
    assert_equal(data.metastable_id(430990000, 1), 430990002)
    assert_equal(data.metastable_id(310720000, 1), 310720002)
    assert_equal(data.metastable_id(451080000, 1), 451080004)
    assert_equal(data.metastable_id(611360000, 0), 611360000)
    assert_equal(data.metastable_id(611360000, 1), 611360001)


def test_decay_half_life():
    assert_equal(data.decay_half_life(551370000, 561370000), (949252608.0, 2840184.0))


def test_decay_half_life_byparent():
    assert_equal(data.decay_half_life_byparent(551370000), [(949252608.0, 2840184.0)])


def test_decay_branch_ratio():
    npt.assert_array_almost_equal(
        data.decay_branch_ratio(551370000, 561370000), (1.0, np.nan)
    )


def test_decay_photon_branch_ratio():
    npt.assert_array_almost_equal(
        data.decay_photon_branch_ratio(551370000, 561370000), (1.0, np.nan)
    )


def test_decay_beta_branch_ratio():
    npt.assert_array_almost_equal(
        data.decay_beta_branch_ratio(551370000, 561370000), (1.0, np.nan)
    )


def test_decay_branch_ratio_byparent():
    assert_equal(data.decay_branch_ratio_byparent(551370000), [1.0])


def test_decay_photon_branch_ratio_byparent():
    npt.assert_array_almost_equal(
        data.decay_photon_branch_ratio_byparent(551370000), [(1.0, np.nan)]
    )


def test_decay_beta_branch_ratio_byparent():
    npt.assert_array_almost_equal(
        data.decay_beta_branch_ratio_byparent(551370000), [(1.0, np.nan)]
    )


def test_gamma_energy():
    assert_equal(data.gamma_energy(551370000), [(283.5, 0.1), (661.657, 0.003)])


def test_gamma_energy_byen():
    npt.assert_equal(
        data.gamma_energy_byen(103.5, 0.05),
        [
            (103.5, 0.4),
            (103.5, np.nan),
            (103.5, 0.1),
            (103.5, 0.3),
            (103.5, 0.2),
            (103.5, 0.04),
            (103.5, np.nan),
            (103.519, 0.004),
            (103.519, 0.004),
            (103.54, 0.08),
        ],
    )


def test_gamma_parent_child():
    assert_equal(
        data.gamma_parent_child(103.5, 0.05),
        [
            (-811910000, 801910000),
            (380780000, 370780000),
            (691520018, 691520000),
            (781860000, 771860000),
            (902250000, 882210000),
            (942420000, 922380000),
            (982520000, 521320000),
            (511320000, 521320000),
            (511320001, 521320000),
            (671700000, 681700000),
        ],
    )


def test_gamma_child_byen():
    assert_equal(
        data.gamma_child_byen(103.5, 0.1),
        [
            391020000,
            521320000,
            731720000,
            791870000,
            791870000,
            832120000,
            832120000,
            922380000,
            982500000,
            801910000,
            370780000,
            691520000,
            771860000,
            882210000,
            922380000,
            521320000,
            521320000,
            521320000,
            681700000,
            741800000,
            741800000,
            591330000,
            591330000,
            741800000,
            792000000,
            872210000,
            611560000,
            1002560000,
        ],
    )


def test_gamma_child_byparent():
    assert_equal(data.gamma_child_byparent(551370000), [561370000, 561370000])


def test_gamma_photon_intensity():
    assert_equal(
        data.gamma_photon_intensity(551370000), [(0.00058, 8e-05), (85.1, 0.2)]
    )


def test_gamma_conversion_intensity():
    npt.assert_array_almost_equal(
        data.gamma_conversion_intensity(551370000), [(np.nan, np.nan), (0.1124, np.nan)]
    )


def test_gamma_total_intensity():
    npt.assert_array_almost_equal(
        data.gamma_total_intensity(561370002), [(100.0, np.nan)]
    )


def test_gamma_from_to_byparent():
    assert_equal(
        data.gamma_from_to_byparent(551370000),
        [(561370001, 561370000), (561370002, 561370000)],
    )


def test_gamma_from_to_byen():
    assert_equal(
        data.gamma_from_to_byen(661.65, 0.1),
        [
            (621510087, 621510015),
            (641500021, 641500006),
            (390990016, 390990005),
            (822040062, 822040024),
            (902290055, 902290001),
            (400880011, 400880004),
            (400880011, 400880004),
            (551310023, 551310009),
            (0, 0),
            (431070028, 431070020),
            (972490039, 972490003),
            (0, 0),
            (380930068, 380930050),
            (561370002, 561370000),
            (561370002, 561370000),
            (621520096, 621520019),
            (621540026, 621540006),
            (621540026, 621540006),
            (781810026, 781810000),
            (791930069, 791930033),
            (431060030, 431060027),
        ],
    )


def test_gamma_parent():
    assert_equal(
        data.gamma_parent(661.65, 0.1),
        [
            611510000,
            651500000,
            380990000,
            832040000,
            892290000,
            410880000,
            410880001,
            561310000,
            771830000,
            962480000,
            992530000,
            661550000,
            370930000,
            551370000,
            561370002,
            611520000,
            611540000,
            611540001,
            791810000,
            801930003,
            982520000,
        ],
    )


def test_alpha_energy():
    assert_equal(
        data.alpha_energy(952410000),
        [
            4758.0,
            4800.0,
            4834.0,
            4889.0,
            4956.0,
            4962.0,
            4964.0,
            5004.0,
            5055.0,
            5068.0,
            5089.0,
            5096.0,
            5106.0,
            5114.0,
            5137.0,
            5155.0,
            5178.0,
            5182.0,
            5192.0,
            5217.0,
            5223.0,
            5232.0,
            5244.0,
            5279.0,
            5322.0,
            5388.0,
            5416.5,
            5442.8,
            5469.0,
            5485.56,
            5511.5,
            5544.5,
        ],
    )


def test_alpha_intensity():
    npt.assert_array_almost_equal(
        data.alpha_intensity(952410000),
        [
            1e-05,
            8.6e-05,
            0.0007,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.0001,
            np.nan,
            0.00014,
            0.0004,
            0.0004,
            np.nan,
            0.0004,
            0.00032,
            0.0007,
            0.0003,
            0.0009,
            0.0006,
            1e-05,
            0.0013,
            np.nan,
            0.0024,
            0.0005,
            0.015,
            1.66,
            0.01,
            13.1,
            0.04,
            84.8,
            0.225,
            0.37,
        ],
    )


def test_alpha_parent():
    assert_equal(data.alpha_parent(5322.0, 0.1), [952410000, 972490000])


def test_alpha_child_byen():
    assert_equal(data.alpha_child_byen(5322.0, 0.1), [932370008, -952450000])


def test_alpha_child_byparent():
    assert_equal(
        data.alpha_child_byparent(952410000),
        [
            932370047,
            932370043,
            932370041,
            932370038,
            932370034,
            932370033,
            932370032,
            932370031,
            932370028,
            932370027,
            932370026,
            932370025,
            932370023,
            932370022,
            932370021,
            932370020,
            932370019,
            932370018,
            932370017,
            932370015,
            932370014,
            932370013,
            932370012,
            932370009,
            932370008,
            932370006,
            932370005,
            932370004,
            932370003,
            932370002,
            932370001,
            932370000,
        ],
    )
    assert_equal(
        data.alpha_child_byparent(922350000),
        [
            902310038,
            -902310000,
            902310028,
            902310023,
            902310020,
            902310018,
            -902310000,
            902310016,
            902310014,
            902310013,
            902310012,
            902310007,
            902310005,
            902310004,
            902310003,
            902310002,
            902310001,
            902310000,
        ],
    )


def test_beta_endpoint_energy():
    npt.assert_array_almost_equal(
        data.beta_endpoint_energy(551370000), [514.03, np.nan, 1176.0]
    )


def test_beta_average_energy():
    assert_equal(data.beta_average_energy(551370000), [174.32, 334.65, 416.26])


def test_beta_intensity():
    assert_equal(data.beta_intensity(551370000), [94.7, 0.00058, 5.3])


def test_beta_parent():
    assert_equal(data.beta_parent(1000, 0.1), [310760000, 410990001, 441070000])


def test_beta_child_byen():
    assert_equal(data.beta_child_byen(1000, 0.1), [320760084, 420990050, 451070007])


def test_beta_child_byparent():
    assert_equal(data.beta_child_byparent(551370000), [561370002, 561370001, 561370000])


def test_ecbp_endpoint_energy():
    assert_equal(data.ecbp_endpoint_energy(170320000)[-1], 11600.0)


def test_ecbp_average_energy():
    assert_equal(data.ecbp_average_energy(110220000), [215.54, 835.0])


def test_ec_intensity():
    npt.assert_array_almost_equal(data.ec_intensity(110220000), [9.618, np.nan])


def test_beta_plus_intensity():
    assert_equal(data.beta_plus_intensity(110220000), [90.326, 0.056])


def test_ecbp_parent():
    assert_equal(
        data.ecbp_parent(215.54, 0.5),
        [110220000, 340690000, 340700000, 541230000, 571330000, 601390000],
    )


def test_ecbp_child_byen():
    assert_equal(
        data.ecbp_child_byen(215.54, 0.5),
        [100220001, -330690000, 330700031, 531230020, 561330006, 591390015],
    )


def test_ecbp_child_byparent():
    assert_equal(data.ecbp_child_byparent(110220000), [100220001, 100220000])


def test_id_from_level():
    assert_equal(data.id_from_level(811920000, 445, "X"), 811920010)
    assert_equal(data.id_from_level(561370000, 662), 561370002)


def test_xray_data():
    npt.assert_almost_equal(
        data.calculate_xray_data(551370000, 0.1, 0.1),
        [
            (30.9728, 0.04693557573523976),
            (30.6251, 0.025406227145485287),
            (35.0, 0.017058197119274962),
            (4.29, 0.019708),
        ],
    )


def test_gamma_xray():
    npt.assert_almost_equal(
        data.gamma_xrays(551370000),
        [
            (32.1936, 3.667054764126558),
            (31.8171, 1.9904773259678956),
            (36.4, 1.3504529099055453),
            (4.47, 0.9148692519999999),
        ],
    )


def test_ecbp_xray():
    npt.assert_almost_equal(
        data.ecbp_xrays(120230000),
        [(1.041, 0.0015155635216184057), (1.07, 1.2730733581594526e-05)],
    )


def test_gamma_photon_intensity_byen():
    npt.assert_almost_equal(
        data.gamma_photon_intensity_byen(661.657, 0.05),
        [
            (0.08, 0.017),
            (16.0, 3.0),
            (85.1, 0.2),
            (89.9, 0.14),
            (1.5, 0.1),
            (0.27, np.nan),
            (0.14, np.nan),
            (160.0, 24.0),
            (0.32, 0.1),
            (np.nan, np.nan),
        ],
    )


# Tests associated with "special cases" from decaygen.py


def test_special_branches():
    special_branches = {
        (320770001, 320770000): 0.19,
        (360850001, 360850000): 0.212,
        (451040000, 441040000): 0.0045,
        (451040000, 461040000): 0.9955,
        (461110001, 461110000): 0.73,
        (471100001, 471100000): 0.0133,
        (491190001, 491190000): 0.044,
        (511260001, 511260000): 0.14,
        (521270000, 531270000): 1.0,
        (521290001, 531290000): 0.36,
        (711770001, 711770000): 0.214,
        (842110001, 842110000): 0.00016,
    }
    for item in special_branches:
        assert_equal(
            data.decay_branch_ratio(
                nucname.id_to_state_id(item[0]), nucname.id_to_state_id(item[1])
            )[0],
            special_branches[item],
        )


def test_special_children():
    special_children = {
        451040000: set([441040000, 461040000]),
        521270000: set([531270000]),
    }
    for item in special_children:
        assert_equal(
            set(data.decay_data_children(nucname.id_to_state_id(item))),
            special_children[item],
        )


if __name__ == "__main__":
    nose.runmodule()
