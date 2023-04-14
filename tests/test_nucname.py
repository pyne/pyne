"""nucname tests"""
from __future__ import unicode_literals, division
from unittest import TestCase
import pytest
import warnings

from pyne.utils import QAWarning

warnings.simplefilter("ignore", QAWarning)
from pyne import nucname


def test_name_zz():
    assert nucname.name_zz["He"] == 2
    assert nucname.name_zz["U"] == 92


def test_zz_name():
    assert nucname.zz_name[1] == "H"
    assert nucname.zz_name[94] == "Pu"


def test_LAN():
    assert (
        nucname.LAN ==
        set(
            [
                "Ce",
                "Dy",
                "Er",
                "Eu",
                "Gd",
                "Ho",
                "La",
                "Lu",
                "Nd",
                "Pm",
                "Pr",
                "Sm",
                "Tb",
                "Tm",
                "Yb",
            ]
        ))


def test_ACT():
    assert (
        nucname.ACT ==
        set(
            [
                "Ac",
                "Th",
                "Pa",
                "U",
                "Np",
                "Pu",
                "Am",
                "Cm",
                "Bk",
                "Cf",
                "Es",
                "Fm",
                "Md",
                "No",
                "Lr",
            ]
        ))


def test_TRU():
    assert (
        nucname.TRU ==
        set(
            [
                "Am",
                "Bh",
                "Bk",
                "Cf",
                "Cm",
                "Db",
                "Ds",
                "Es",
                "Fm",
                "Hs",
                "Lr",
                "Md",
                "Mt",
                "No",
                "Np",
                "Pu",
                "Rf",
                "Rg",
                "Sg",
                "Fl",
                "Lv",
                "Cn",
            ]
        ))


def test_MA():
    assert (
        nucname.MA == set(["Np", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"]))


def test_FP():
    assert (
        nucname.FP ==
        set(
            [
                "Ag",
                "Al",
                "Ar",
                "As",
                "At",
                "Au",
                "B",
                "Ba",
                "Be",
                "Bi",
                "Br",
                "C",
                "Ca",
                "Cd",
                "Ce",
                "Cl",
                "Co",
                "Cr",
                "Cs",
                "Cu",
                "Dy",
                "Er",
                "Eu",
                "F",
                "Fe",
                "Fr",
                "Ga",
                "Gd",
                "Ge",
                "H",
                "He",
                "Hf",
                "Hg",
                "Ho",
                "I",
                "In",
                "Ir",
                "K",
                "Kr",
                "La",
                "Li",
                "Lu",
                "Mg",
                "Mn",
                "Mo",
                "N",
                "Na",
                "Nb",
                "Nd",
                "Ne",
                "Ni",
                "O",
                "Os",
                "P",
                "Pb",
                "Pd",
                "Pm",
                "Po",
                "Pr",
                "Pt",
                "Ra",
                "Rb",
                "Re",
                "Rh",
                "Rn",
                "Ru",
                "S",
                "Sb",
                "Sc",
                "Se",
                "Si",
                "Sm",
                "Sn",
                "Sr",
                "Ta",
                "Tb",
                "Tc",
                "Te",
                "Ti",
                "Tl",
                "Tm",
                "V",
                "W",
                "Xe",
                "Y",
                "Yb",
                "Zn",
                "Zr",
            ]
        ))


def test_lan():
    assert (
        nucname.lan == set([57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71]))


def test_act():
    assert (
        nucname.act ==
        set([89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103]))


def test_tru():
    assert (
        nucname.tru ==
        set(
            [
                93,
                94,
                95,
                96,
                97,
                98,
                99,
                100,
                101,
                102,
                103,
                104,
                105,
                106,
                107,
                108,
                109,
                110,
                111,
                112,
                114,
                116,
            ]
        ))


def test_ma():
    assert nucname.ma == set([93, 95, 96, 97, 98, 99, 100, 101, 102, 103])


def test_fp():
    assert (
        nucname.fp ==
        set(
            [
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10,
                11,
                12,
                13,
                14,
                15,
                16,
                17,
                18,
                19,
                20,
                21,
                22,
                23,
                24,
                25,
                26,
                27,
                28,
                29,
                30,
                31,
                32,
                33,
                34,
                35,
                36,
                37,
                38,
                39,
                40,
                41,
                42,
                43,
                44,
                45,
                46,
                47,
                48,
                49,
                50,
                51,
                52,
                53,
                54,
                55,
                56,
                57,
                58,
                59,
                60,
                61,
                62,
                63,
                64,
                65,
                66,
                67,
                68,
                69,
                70,
                71,
                72,
                73,
                74,
                75,
                76,
                77,
                78,
                79,
                80,
                81,
                82,
                83,
                84,
                85,
                86,
                87,
                88,
            ]
        ))


def check_cases(f, x, exp):
    obs = f(x)
    assert exp == obs


cases = [
    20040,
    "he4",
    "Cm-244",
    "PU239",
    "AM242M",
    2004,
    95642,
    95242,
    92636,
    95942,
    "Am-242m",
    "he",
    "U",
    "Np",
    "4he",
    "244CM",
    "239Pu",
    "242AM",
    40020,
    2440961,
    2390940,
    2420950,
    92,
]

caseids = [
    20040000,
    20040000,
    962440000,
    942390000,
    952420001,
    20040000,
    952420000,
    952420001,
    922360001,
    952420004,
    952420001,
    20000000,
    920000000,
    930000000,
    20040000,
    962440000,
    942390000,
    952420000,
    20040000,
    962440001,
    942390000,
    952420000,
    920000000,
]


def test_id():
    assert nucname.id(20040) == 20040000

    assert nucname.id("he4") == 20040000
    assert nucname.id("Cm-244") == 962440000
    assert nucname.id("PU239") == 942390000
    assert nucname.id("AM242M") == 952420001

    assert nucname.id(2004) == 20040000
    assert nucname.id(95642) == 952420000
    assert nucname.id(95242) == 952420001
    assert nucname.id(92636) == 922360001
    assert nucname.id(95942) == 952420004

    assert nucname.id("Am-242m") == 952420001

    assert nucname.id("he") == 20000000
    assert nucname.id("U") == 920000000
    assert nucname.id("Np") == 930000000
    assert nucname.id("Cl") == 170000000

    assert nucname.id("4he") == 20040000
    assert nucname.id("244CM") == 962440000
    assert nucname.id("239Pu") == 942390000
    assert nucname.id("242AM") == 952420000

    assert nucname.id(40020) == 20040000
    assert nucname.id(2440961) == 962440001
    assert nucname.id(2390940) == 942390000
    assert nucname.id(2420950) == 952420000
    assert nucname.id(92) == 920000000

    assert nucname.id("94-Pu-239") == nucname.id("Pu-239")
    assert nucname.id("95-Am-242m") == nucname.id("Am-242m")
    assert nucname.id("94-Pu-239") == nucname.id("Pu-239")
    assert nucname.id("95-Am-242") == nucname.id("Am-242")

    pytest.raises(RuntimeError, nucname.id, "0-H-1")


def test_name():
    assert nucname.name(942390) == "Pu239"
    assert nucname.name(952421) == "Am242M"

    assert nucname.name("PU239") == "Pu239"

    assert nucname.name(94239) == "Pu239"
    assert nucname.name(95242) == "Am242M"
    assert nucname.name(95642) == "Am242"
    assert nucname.name(92636) == "U236M"

    assert nucname.name("Am-242m") == "Am242M"

    assert nucname.name(40020) == "He4"
    assert nucname.name(2440961) == "Cm244M"
    assert nucname.name(2390940) == "Pu239"
    assert nucname.name(2420950) == "Am242"


@pytest.mark.parametrize("case, exp", zip(cases,[
        2,
        2,
        96,
        94,
        95,
        2,
        95,
        95,
        92,
        95,
        95,
        2,
        92,
        93,
        2,
        96,
        94,
        95,
        2,
        96,
        94,
        95,
        92,
    ]))
def test_znum(case,exp):
    check_cases(nucname.znum, case, exp)


@pytest.mark.parametrize("case, exp",zip(cases,[
        4,
        4,
        244,
        239,
        242,
        4,
        242,
        242,
        236,
        242,
        242,
        0,
        0,
        0,
        4,
        244,
        239,
        242,
        4,
        244,
        239,
        242,
        0,
    ]))
def test_anum(case, exp):
    check_cases(nucname.anum, case, exp)


@pytest.mark.parametrize("case, exp",zip(cases,
    [0, 0, 0, 0, 1, 0, 0, 1, 1, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]))
def test_snum(case, exp):
    check_cases(nucname.snum, case, exp)


def test_zzaaam():
    assert nucname.zzaaam(20040) == 20040

    assert nucname.zzaaam("he4") == 20040
    assert nucname.zzaaam("Cm-244") == 962440
    assert nucname.zzaaam("PU239") == 942390
    assert nucname.zzaaam("AM242M") == 952421

    assert nucname.zzaaam(2004) == 20040
    assert nucname.zzaaam(95642) == 952420
    assert nucname.zzaaam(95242) == 952421
    assert nucname.zzaaam(92636) == 922361
    assert nucname.zzaaam(95942) == 952424

    assert nucname.zzaaam("Am-242m") == 952421

    assert nucname.zzaaam("he") == 20000
    assert nucname.zzaaam("U") == 920000
    assert nucname.zzaaam("Np") == 930000

    assert nucname.zzaaam("4he") == 20040
    assert nucname.zzaaam("244CM") == 962440
    assert nucname.zzaaam("239Pu") == 942390
    assert nucname.zzaaam("242AM") == 952420

    assert nucname.zzaaam(40020) == 20040
    assert nucname.zzaaam(2440961) == 962441
    assert nucname.zzaaam(2390940) == 942390
    assert nucname.zzaaam(2420950) == 952420


@pytest.mark.parametrize("val, id",set(zip([
        20040,
        20040,
        962440,
        942390,
        952421,
        20040,
        952420,
        952421,
        922361,
        952424,
        952421,
        20000,
        920000,
        930000,
        20040,
        962440,
        942390,
        952420,
        20040,
        962441,
        942390,
        952420,
        920000,
    ],caseids)))
def test_zzaaam_to_id(val, id):
    if val is not None:
        check_cases(nucname.zzaaam_to_id, val, id)


def test_zzzaaa():

    assert nucname.zzzaaa(20040) == 2004
    assert nucname.zzzaaa("he4") == 2004
    assert nucname.zzzaaa("Cm-244") == 96244
    assert nucname.zzzaaa("PU239") == 94239
    assert nucname.zzzaaa("AM242M") == 95242

    assert nucname.zzzaaa(2004) == 2004
    assert nucname.zzzaaa(95642) == 95242
    assert nucname.zzzaaa(95242) == 95242
    assert nucname.zzzaaa(92636) == 92236
    assert nucname.zzzaaa(95942) == 95242

    assert nucname.zzzaaa("Am-242m") == 95242

    assert nucname.zzzaaa("he") == 2000
    assert nucname.zzzaaa("U") == 92000
    assert nucname.zzzaaa("Np") == 93000

    assert nucname.zzzaaa("4he") == 2004
    assert nucname.zzzaaa("244CM") == 96244
    assert nucname.zzzaaa("239Pu") == 94239
    assert nucname.zzzaaa("242AM") == 95242

    assert nucname.zzzaaa(40020) == 2004
    assert nucname.zzzaaa(2440961) == 96244
    assert nucname.zzzaaa(2390940) == 94239
    assert nucname.zzzaaa(2420950) == 95242


def test_zzzaaa_to_id():

    assert nucname.zzzaaa_to_id(2004) == nucname.id(20040)
    assert nucname.zzzaaa_to_id(96244) == nucname.id("Cm-244")
    assert nucname.zzzaaa_to_id(94239) == nucname.id("PU239")

    assert nucname.zzzaaa_to_id(95242) == nucname.id(95642)
    # Added /10*10 to remove the state information from id (ZZZAAA carries no state
    # information, defaults to zero when converting ZZZAAA back to ID)
    assert nucname.zzzaaa_to_id(95242) == (nucname.id(95942) // 10) * 10
    assert nucname.zzzaaa_to_id(95242) == (nucname.id("AM-242m") // 10) * 10

    assert nucname.zzzaaa_to_id(2000) == nucname.id("he")
    assert nucname.zzzaaa_to_id(92000) == nucname.id("U")
    assert nucname.zzzaaa_to_id(93000) == nucname.id("Np")

    assert nucname.zzzaaa_to_id(2004) == nucname.id(40020)


def test_zzllaaam():
    assert nucname.zzllaaam(942390) == "94-Pu-239"
    assert nucname.zzllaaam(952421) == "95-Am-242m"

    assert nucname.zzllaaam("Pu-239") == "94-Pu-239"

    assert nucname.zzllaaam(94239) == "94-Pu-239"
    assert nucname.zzllaaam(95642) == "95-Am-242"
    assert nucname.zzllaaam(95242) == "95-Am-242m"
    assert nucname.zzllaaam(92636) == "92-U-236m"

    assert nucname.zzllaaam(2390940) == "94-Pu-239"
    assert nucname.zzllaaam(2420951) == "95-Am-242m"


def test_zzllaaam_to_id():
    assert nucname.zzllaaam_to_id("94-Pu-239") == nucname.id("Pu-239")
    assert nucname.zzllaaam_to_id("95-Am-242m") == nucname.id("Am-242m")

    assert nucname.zzllaaam_to_id("94-Pu-239") == nucname.id("Pu-239")
    assert nucname.zzllaaam_to_id("95-Am-242") == nucname.id("Am-242")
    assert nucname.zzllaaam_to_id("95-Am-242m") == nucname.id("Am-242m")


def test_mcnp():
    assert nucname.mcnp(10010) == 1001
    assert nucname.mcnp(952421) == 95242
    assert nucname.mcnp(952420) == 95642
    assert nucname.mcnp(922361) == 92636

    assert nucname.mcnp("H1") == 1001
    assert nucname.mcnp("AM242M") == 95242
    assert nucname.mcnp("AM242") == 95642
    assert nucname.mcnp("U236M") == 92636

    assert nucname.mcnp(1001) == 1001
    assert nucname.mcnp(95642) == 95642

    assert nucname.mcnp("Am-242") == 95642
    assert nucname.mcnp("Am-242m") == 95242
    assert nucname.mcnp("U-236m") == 92636

    assert nucname.mcnp(10010) == 1001
    assert nucname.mcnp(2420951) == 95242
    assert nucname.mcnp(2420950) == 95642
    assert nucname.mcnp(2360921) == 92636


@pytest.mark.parametrize("val, id",set(zip([
        2004,
        2004,
        96244,
        94239,
        95242,
        2004,
        95642,
        95242,
        92636,
        95942,
        95242,
        2000,
        92000,
        93000,
        2004,
        96244,
        94239,
        95642,
        2004,
        96644,
        94239,
        95642,
        92000,
    ], caseids)))
def test_mcnp_to_id(val, id):
    if val is not None:
        check_cases(nucname.mcnp_to_id, val, id)

    # tests for invalid inputs
    pytest.raises(RuntimeError, nucname.mcnp_to_id, 92)


def test_openmc():
    assert nucname.openmc(10010) == "H1"
    assert nucname.openmc(952421) == "Am242_m1"
    assert nucname.openmc(952420) == "Am242"
    assert nucname.openmc(922361) == "U236_m1"

    assert nucname.openmc("H1") == "H1"
    assert nucname.openmc("AM242M") == "Am242_m1"
    assert nucname.openmc("AM242") == "Am242"
    assert nucname.openmc("U236M") == "U236_m1"

    assert nucname.openmc(1001) == "H1"
    assert nucname.openmc(95642) == "Am242"

    assert nucname.openmc("Am-242") == "Am242"
    assert nucname.openmc("Am-242m") == "Am242_m1"
    assert nucname.openmc("U-236m") == "U236_m1"

    assert nucname.openmc(10010) == "H1"
    assert nucname.openmc(2420951) == "Am242_m1"
    assert nucname.openmc(2420950) == "Am242"
    assert nucname.openmc(2360921) == "U236_m1"

    # natural elements
    assert nucname.openmc(20000000) == "He0"
    assert nucname.openmc(920000000) == "U0"
    assert nucname.openmc(930000000) == "Np0"


@pytest.mark.parametrize("val, id", set(zip([
        "He4",
        "He4",
        "Cm244",
        "Pu239",
        "Am242_m1",
        "He4",
        "Am242",
        "Am242_m1",
        "U236_m1",
        "Am242_m4",
        "Am242_m1",
        "He0",
        "U0",
        "Np0",
        "He4",
        "Cm244",
        "Pu239",
        "Am242",
        "He4",
        "Cm244_m1",
        "Pu239",
        "Am242",
        "U0",
    ], caseids)))
def test_openmc_to_id(val, id):
    if val is not None:
        check_cases(nucname.openmc_to_id, val, id)

    # tests for invalid inputs
    pytest.raises(Exception, nucname.openmc_to_id, 92)


def test_fluka():
    assert nucname.fluka(40000000) == "BERYLLIU"
    assert nucname.fluka(640000000) == "GADOLINI"
    assert nucname.fluka(280000000) == "NICKEL"
    assert nucname.fluka(1140000000) == "UNUNQUAD"
    assert nucname.fluka(1140000000) != "UNUNQUA"


def test_fluka_to_id():
    assert nucname.fluka_to_id("BERYLLIU") == 40000000
    assert nucname.fluka_to_id("NICKEL") == 280000000
    assert nucname.fluka_to_id("LITHIU-7") == 30070000


def test_serpent():
    assert nucname.serpent(942390) == "Pu-239"
    assert nucname.serpent(952421) == "Am-242m"

    assert nucname.serpent("Pu-239") == "Pu-239"

    assert nucname.serpent(94239) == "Pu-239"
    assert nucname.serpent(95642) == "Am-242"
    assert nucname.serpent(95242) == "Am-242m"
    assert nucname.serpent(92636) == "U-236m"

    assert nucname.serpent(2390940) == "Pu-239"
    assert nucname.serpent(2420951) == "Am-242m"


@pytest.mark.parametrize("val, id",set(zip([
        "He-4",
        "He-4",
        "Cm-244",
        "Pu-239",
        "Am-242m",
        "He-4",
        "Am-242",
        "Am-242m",
        "U-236m",
        None,
        "Am-242m",
        "He-nat",
        "U-nat",
        "Np-nat",
        "He-4",
        "Cm-244",
        "Pu-239",
        "Am-242",
        "He-4",
        "Cm-244m",
        "Pu-239",
        "Am-242",
        "U-nat",
    ], caseids)))
def test_serpent_to_id(val, id):
    if val is not None:
        check_cases(nucname.serpent_to_id, val, id)


def test_nist():
    assert nucname.nist(942390) == "239Pu"
    assert nucname.nist(952421) == "242Am"

    assert nucname.nist("Pu-239") == "239Pu"

    assert nucname.nist(94239) == "239Pu"
    assert nucname.nist(95642) == "242Am"

    assert nucname.nist("10000") == "H"
    assert nucname.nist("920000") == "U"
    assert nucname.nist("940000") == "Pu"

    assert nucname.nist(2390940) == "239Pu"
    assert nucname.nist(2420951) == "242Am"


@pytest.mark.parametrize("val, id",set(zip([
        "4He",
        "4He",
        "244Cm",
        "239Pu",
        None,
        "4He",
        "242Am",
        None,
        None,
        None,
        None,
        "He",
        "U",
        "Np",
        "4He",
        "244Cm",
        "239Pu",
        "242Am",
        "4He",
        None,
        "239Pu",
        "242Am",
        "U",
    ], caseids)))
def test_nist_to_id(val,id):
    if val is not None:
        check_cases(nucname.nist_to_id, val, id)


def test_cinder():
    assert nucname.cinder(10020) == 20010
    assert nucname.cinder(952421) == 2420951

    assert nucname.cinder("H2") == 20010
    assert nucname.cinder("AM242M") == 2420951

    assert nucname.cinder(1002) == 20010
    assert nucname.cinder(95642) == 2420950
    assert nucname.cinder(95242) == 2420951
    assert nucname.cinder(92636) == 2360921

    assert nucname.cinder("Am-242m") == 2420951

    assert nucname.cinder(20010) == 20010
    assert nucname.cinder(2420951) == 2420951


@pytest.mark.parametrize("val, id",set(zip([
        40020,
        40020,
        2440960,
        2390940,
        2420951,
        40020,
        2420950,
        2420951,
        2360921,
        2420954,
        2420951,
        20,
        920,
        930,
        40020,
        2440960,
        2390940,
        2420950,
        40020,
        2440961,
        2390940,
        2420950,
        920,
    ], caseids)))
def test_cinder_to_id(val, id):
    if val is not None:
        check_cases(nucname.cinder_to_id, val, id)


def test_alara():
    assert nucname.alara(942390) == "pu:239"
    assert nucname.alara(952421) == "am:242"

    assert nucname.alara("PU239") == "pu:239"

    assert nucname.alara(94239) == "pu:239"
    assert nucname.alara(95242) == "am:242"
    assert nucname.alara(95642) == "am:242"
    assert nucname.alara(92636) == "u:236"
    assert nucname.alara(2000) == "he"

    assert nucname.alara("Am-242m") == "am:242"

    assert nucname.alara(40020) == "he:4"
    assert nucname.alara(20000) == "he"
    assert nucname.alara(2440961) == "cm:244"
    assert nucname.alara(2390940) == "pu:239"
    assert nucname.alara(2420950) == "am:242"


@pytest.mark.parametrize("val, id",set(zip([
        "he:4",
        "he:4",
        "cm:244",
        "pu:239",
        None,
        "he:4",
        "am:242",
        None,
        None,
        None,
        None,
        "he",
        "u",
        "np",
        "he:4",
        "cm:244",
        "pu:239",
        "am:242",
        "he:4",
        None,
        "pu:239",
        "am:242",
        "u",
    ], caseids)))
def test_alara_to_id(val, id):
    if val is not None:
        check_cases(nucname.alara_to_id, val, id)


def test_sza():
    assert nucname.sza(20040) == 2004

    assert nucname.sza("he4") == 2004
    assert nucname.sza("Cm-244") == 96244
    assert nucname.sza("PU239") == 94239
    assert nucname.sza("AM242M") == 1095242

    assert nucname.sza(2004) == 2004
    assert nucname.sza(95642) == 95242
    assert nucname.sza(95242) == 1095242
    assert nucname.sza(92636) == 1092236
    assert nucname.sza(95942) == 4095242

    assert nucname.sza("Am-242m") == 1095242

    assert nucname.sza("he") == 2000
    assert nucname.sza("U") == 92000
    assert nucname.sza("Np") == 93000

    assert nucname.sza("4he") == 2004
    assert nucname.sza("244CM") == 96244
    assert nucname.sza("239Pu") == 94239
    assert nucname.sza("242AM") == 95242

    assert nucname.sza(40020) == 2004
    assert nucname.sza(2440961) == 1096244
    assert nucname.sza(2390940) == 94239
    assert nucname.sza(2420950) == 95242


def test_groundstate():

    assert nucname.groundstate("he4") == 20040000
    assert nucname.groundstate("Cm-244") == 962440000
    assert nucname.groundstate("PU239") == 942390000
    assert nucname.groundstate("AM242M") == 952420000

    assert nucname.groundstate(2004) == 20040000
    assert nucname.groundstate(95642) == 952420000
    assert nucname.groundstate(95242) == 952420000
    assert nucname.groundstate(92636) == 922360000
    assert nucname.groundstate(95942) == 952420000

    assert nucname.groundstate("he") == 20000000
    assert nucname.groundstate("U") == 920000000
    assert nucname.groundstate("Np") == 930000000
    assert nucname.groundstate("Cl") == 170000000

    assert nucname.groundstate("4he") == 20040000
    assert nucname.groundstate("244CM") == 962440000
    assert nucname.groundstate("239Pu") == 942390000
    assert nucname.groundstate("242AM") == 952420000

    assert nucname.groundstate(40020) == 20040000
    assert nucname.groundstate(2440961) == 962440000
    assert nucname.groundstate(2390940) == 942390000
    assert nucname.groundstate(2420950) == 952420000
    assert nucname.groundstate(92) == 920000000

    assert nucname.groundstate("94-Pu-239") == 942390000
    assert nucname.groundstate("95-Am-242m") == 952420000
    assert nucname.groundstate("94-Pu-239") == 942390000
    assert nucname.groundstate("95-Am-242") == 952420000


@pytest.mark.parametrize("val, id",set(zip([
        2004,
        2004,
        96244,
        94239,
        1095242,
        2004,
        95242,
        1095242,
        1092236,
        4095242,
        1095242,
        2000,
        92000,
        93000,
        2004,
        96244,
        94239,
        95242,
        2004,
        1096244,
        94239,
        95242,
        92000,
    ], caseids)))
def test_sza_to_id(val, id):
    if val is not None:
        check_cases(nucname.sza_to_id, val, id)


@pytest.mark.parametrize("nuc", [922350, "U235"])
def test_isnuclide(nuc):
    assert nucname.isnuclide(nuc)


@pytest.mark.parametrize("nuc", ["U3", -30060000])
def test_isnotnuclide(nuc):
    assert not nucname.isnuclide(nuc)


@pytest.mark.parametrize("nuc", [92, "U"])
def test_iselement_U235(nuc):
    assert nucname.iselement(nuc)


@pytest.mark.parametrize("nuc", [922350, "U235"])
def test_isnotelement_U235(nuc):
    assert not nucname.iselement(nuc)


@pytest.mark.parametrize("nuc", [1, "H"])
def test_iselement_H1(nuc):
    assert nucname.iselement(nuc)


@pytest.mark.parametrize("nuc", [1001, "H1"])
def test_isnotelement_H1(nuc):
    assert not nucname.iselement(nuc)


def test_state_id_to_id():
    assert nucname.state_id_to_id(190380015) == 190380002
    assert nucname.state_id_to_id(270590001) == None  # Co59M doesn't exist


def test_id_to_state_id():
    assert nucname.id_to_state_id(190380002) == 190380015
    assert nucname.id_to_state_id(270590001) == None  # Co59M doesn't exist


def test_ensdf_to_id():
    assert nucname.ensdf_to_id("16O") == 80160000
    assert nucname.ensdf_to_id("28614") == 1142860000
    assert nucname.ensdf_to_id("269Hs") == 1082690000

