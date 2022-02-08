import os
import io
import warnings
from math import e
from hashlib import md5

import nose
from nose.tools import assert_equal
from nose import SkipTest

import numpy as np
from numpy.testing import assert_array_equal, assert_allclose, assert_array_almost_equal

from pyne.utils import QAWarning

warnings.simplefilter("ignore", QAWarning)

from pyne.endf import Library, Evaluation
from pyne.utils import endftod
from pyne.rxdata import DoubleSpinDict
from pyne.xs.data_source import ENDFDataSource
from pyne import nucname

from utils import download_file


def try_download(url, target, sha):
    try:
        download_file(url, target, sha)
    except:
        raise SkipTest(url + " not available")


def ignore_future_warnings(func):
    """This is a decorator which can be used to ignore FutureWarnings
    occurring in a function."""

    def new_func(*args, **kwargs):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning)
            return func(*args, **kwargs)

    new_func.__name__ = func.__name__
    new_func.__doc__ = func.__doc__
    new_func.__dict__.update(func.__dict__)
    return new_func


tape1path = "files_test_endf/sampletape1"
library = Library(tape1path)
nuc1002, nuc1003, nuc40000 = nucname.id(1002), nucname.id(1003), nucname.id(40000)
library._read_res(nuc1002)
library._read_res(nuc1003)
library._read_res(nuc40000)
nuc40040 = nucname.id(40040)


def array_from_ENDF(fh):
    "Convert a chunk of ENDF, stripped of metadata, into a numpy array."
    return np.genfromtxt(
        fh,
        dtype="float",
        delimiter=11,
        converters={
            0: endftod,
            1: endftod,
            2: endftod,
            3: endftod,
            4: endftod,
            5: endftod,
        },
    )


def test_endftod():

    from pyne._utils import endftod

    obs = [
        endftod(" 3.28559+12"),
        endftod(" 2.328559+4"),
        endftod(" 3.28559-12"),
        endftod(" 2.328559-2"),
        endftod("-3.28559+12"),
        endftod("-2.328559+2"),
        endftod("-3.28559-12"),
        endftod("-2.328559-2"),
        endftod("        121"),
        endftod("       -121"),
    ]
    exp = [
        3.28559e12,
        2.328559e4,
        3.28559e-12,
        2.328559e-2,
        -3.28559e12,
        -2.328559e2,
        -3.28559e-12,
        -2.328559e-2,
        121.0,
        -121.0,
    ]
    obs = np.array(obs)
    exp = np.array(exp)
    assert_allclose(obs, exp, rtol=1e-8)


def test_loadtape():
    try_download(
        "http://www.nndc.bnl.gov/endf/b6.8/tapes/tape.100",
        "endftape.100",
        "b56dd0aee38bd006c58181e473232776",
    )

    testlib = Library("endftape.100")
    exp_nuclides = set(
        [10010000, 30060000, 40090000, 50110000, 30070000, 60000000, 50100000]
    )
    obs_nuclides = set(map(int, testlib.structure.keys()))
    assert_equal(exp_nuclides, obs_nuclides)


def test_get():
    obs = library.get_rx(nuc40000, 4, 2)
    exp = [
        4.898421e3,
        6.768123e0,
        0,
        1,
        0,
        0,
        2.123124e6,
        8.123142e-6,
        2.123212e6,
        8.231231e-6,
        -2.231211e6,
        8.123421e-6,
    ]
    try:
        # try to get a bad key
        library.get_rx(111, 1, 1)
        assert False
    except ValueError:
        assert True
    assert_array_equal(exp, obs)


def test_unresolved_resonances_a():
    # Case A (ENDF Manual p.70)
    obs = library.structure[nuc1003]["data"][nuc1003]["unresolved"]
    obs_LIST = obs[1][2][2, 2]

    exp = array_from_ENDF(
        io.BytesIO(
            b""" 1.801000+3          0 1.100000+0 3.078520-1 1.000000-2 0.000000+0
 2.101000+3          1 2.100000+0 7.088000-1 2.000000-2 0.000000+0
 3.101000+3          2 3.100000+0 2.120000-1 3.000000-2 0.000000+0"""
        )
    )
    exp_LIST = dict(zip(("D", "AJ", "AMUN", "GN0", "GG"), exp.transpose()))

    for key in exp_LIST:
        assert_array_equal(exp_LIST[key], obs_LIST[key])


def test_unresolved_resonances_b():
    # Case B (ENDF Manual p. 70)
    obs = library.structure[nuc40000]["data"][nuc40040]["unresolved"]
    # For the spin=4.5, L=3, J=4 section in the first isotope
    obs_1 = obs[0][2][4.5, 3, 4]
    exp_1_a = array_from_ENDF(
        io.BytesIO(
            b""" 0.000000+0 0.000000+0 3.000000+0          3         12          0
 2.804009-5 4.000000+0 3.181707+3 3.885315-9-3.382438+3 0.000000+0
 2.376630+2 7.198625-2-5.887887-8-4.380016-5 1.747888-6-4.104291-9"""
        )
    )
    exp_1 = dict(
        zip(
            (0, 0, "L", "MUF", "NE+6", 0, "D", "AJ", "AMUN", "GN0", "GG"),
            exp_1_a[:2].flat,
        )
    )
    exp_1["GF"] = exp_1_a[2]
    del exp_1[0]

    for key in exp_1:
        assert_array_equal(obs_1[key], exp_1[key])
    # For the spin=3.5, L=4, J=5 section in the second isotope
    obs_2 = obs[1][2][3.5, 4, 5]
    exp_2_a = array_from_ENDF(
        io.BytesIO(
            b""" 0.000000+0 0.000000+0 4.000000+0          4         13          0
-9.824193-5 5.000000+0 4.676826-4-4.336597+0-9.045122+2 0.000000+0
 3.699655-9-3.919000+5 8.467144-3-3.737007+9-5.750577+7-9.588021+8
-3.280571+7                                                       """
        )
    )
    exp_2 = dict(
        zip(
            (0, 0, "L", "MUF", "NE+6", 0, "D", "AJ", "AMUN", "GN0", "GG"),
            exp_2_a[:2].flat,
        )
    )
    num_e = int(exp_2["NE+6"]) - 6
    exp_2["GF"] = exp_2_a[2:].flat[:num_e]
    del exp_2[0]
    for key in exp_2:
        assert_array_equal(obs_2[key], exp_2[key])

    # Check the energy values.
    obs_ES = obs[1][2]["ES"]
    exp_ES_a = array_from_ENDF(
        io.BytesIO(
            b"""-2.723837-2-8.755303-2 2.245337-2-9.034520+2 2.252098+5 2.666587+2
 3.747872-3                                                       """
        )
    )
    exp_ES = exp_ES_a.flat[:num_e]
    assert_array_equal(obs_ES, exp_ES)


def test_unresolved_resonances_c():
    # Case C (ENDF Manual p. 70)
    obs = library.structure[nuc40000]["data"][nuc40040]["unresolved"][2][2][0.5, 6, 9]

    exp_a = array_from_ENDF(
        io.BytesIO(
            b""" 9.000000+0 0.000000+0          2          0         18          2
 0.000000+0 0.000000+0-4.253833-1-2.269388+0 0.000000+0 4.732644-4
-5.873521-3-4.808214+9 5.089619+5 4.836683+0 2.772702-3-4.865151-8
-2.659480-9 1.044275+8-1.393749+2-4.189996-6-9.596467-4 3.942829+9"""
        )
    )
    exp = dict(zip(("ES", "D", "GX", "GN0", "GG", "GF"), exp_a[2:].transpose()))
    exp.update(
        dict(
            zip(
                (
                    "AJ",
                    0,
                    "INT",
                    0,
                    "6*NE+6",
                    "NE",
                    0,
                    0,
                    "AMUX",
                    "AMUN",
                    "AMUG",
                    "AMUF",
                ),
                exp_a[:2].flat,
            )
        )
    )
    exp["AWRI"] = -5.702860e9
    for key in obs:
        assert_array_equal(obs[key], exp[key])


def test_DoubleSpinDict():
    subject = DoubleSpinDict(
        {
            (3.499999999998, 2, 1): {"A": "a", "B": "b"},
            (2.000000000012, 3, 4): {"C": "c", "D": "d"},
        }
    )
    subject.update({(3.500000000011, 8, 9): {"E": "e", "F": "f"}})

    obs = subject[(3.48, 8, 9)]
    exp = {"E": "e", "F": "f"}
    assert_equal(exp, obs)


def test_resolved_breitwigner():
    """The section looks like this:
            EL         EH        LRU        LRF        NRO       NAPS 419 2151    3
    0.000000+0 0.000000+0          0          0         NR         NP 419 2151    4
           SPI         AP          0          0        NLS          0 419 2151    5
          AWRI         QX          L        LRX      6*NRS        NRS 419 2151    6
            ER         AJ         GT         GN         GG         GF 419 2151    7
            ER         AJ         GT         GN         GG         GF 419 2151    8
          AWRI         QX          L        LRX      6*NRS        NRS 419 2151    6
            ER         AJ         GT         GN         GG         GF 419 2151    7
            ER         AJ         GT         GN         GG         GF 419 2151    8"""

    data = library.structure[nuc1002]["data"][nuc1002]["resolved"]
    # Check to see if NRO is reading from the right place.
    # NRO = 0 case
    range_nro_0 = data[2]
    assert_equal(range_nro_0[-1]["NRO"], 0)
    # NRO = 1 case
    # Check to see if NAPS is reading from the right place
    assert_equal(range_nro_0[-1]["NAPS"], 1)
    # Check to see if SPI, NLS are reading from the right place
    assert_equal(range_nro_0[-1]["SPI"], 0.5)
    assert_equal(range_nro_0[-1]["NLS"], 1)
    # Check to see if the data is alright...
    expected = {
        "ER": [350000.0, 4500000.0],
        "AJ": [1.0, 2.0],
        "GT": [10.0, 20.0],
        "GN": [2.0, 3.0],
        "GG": [1.1, 1.2],
        "GF": [3.1, 3.2],
    }

    for key in range_nro_0[2][0.5, 1]:
        assert_allclose(range_nro_0[2][0.5, 1][key], expected[key], rtol=1e-8)


def test_resolved_reichmoore():
    """The section looks like this:
            EL         EH        LRU        LRF        NRO       NAPS 419 2151    3
    0.000000+0 0.000000+0          0          0         NR         NP 419 2151    4
           SPI         AP        LAD          0        NLS       NLSC 419 2151    5
          AWRI        APL          L          0      6*NRS        NRS 419 2151    6
            ER         AJ         GN         GG        GFA        GFB 419 2151    7"""

    subsection = library.structure[nuc1002]["data"][nuc1002]["resolved"][1]
    assert_array_equal(subsection[2]["int"]["intpoints"], [3, 7])
    assert_array_equal(subsection[2]["int"]["intschemes"], [4, 1])

    obs_data = subsection[2][2.5, 2]
    exp_data = {
        "ER": 4.127773e3,
        "AJ": -3.956950e-7,
        "GN": 3.739684e-5,
        "GG": -3.872199e7,
        "GFA": 2.259559e5,
        "GFB": -3.096948e-8,
    }
    for key in subsection[2][2.5, 2]:
        assert_allclose(obs_data[key], exp_data[key], rtol=1e-8)


def test_resolved_adleradler():
    """The section looks like this:
    [MAT, 2,151/  0.0,  0.0,    0,    0,    NR,   NP/ points / AP(E) ] TAB1
    [MAT, 2,151/  SPI,   AP,    0,    0,   NLS,    0 ] CONT
    [MAT, 2,151/ AWRI,  0.0,   LI,    0,  6*NX,   NX
                  AT1,  AT2,  AT3,  AT4,   BT1,  BT2,
                  AF1, -----------------------,  BF2,
                  AC1, -----------------------,  BC2 ] LIST
    [MAT, 2,151/  0.0,  0.0,    L,    0,   NJS,    0] CONT(l)
    [MAT, 2,151/   AJ,  0.0,    0,    0,12*NLJ,  NLJ/
                 DET1, DWT1, GRT1, GIT1,  DEF1, DWF1,
                 GRF1, GIF1, DEC1, DWC1,  GRC1, GIC1,
                 DET2, DWT2, GIC2,---- --------------
                 DET3,-------------------------------
                 ------------------------------------
                 ------------------------, GICN LJ ] LIST"""
    subsection = library.structure[nuc1002]["data"][nuc1002]["resolved"][0]

    obs_LIST = subsection[2][3.5, 5, 3]
    obs_bg = subsection[2]["bg"]

    # Test to see if the LIST records read in right
    exp_LIST = {
        "DET": -3.211014e1,
        "DWT": -2.011165e3,
        "GRT": 4.178337e5,
        "GIT": 1.640997e2,
        "DEF": 1.122313e-5,
        "DWF": 1.537114e-2,
        "GRF": 4.634918e-8,
        "GIF": -3.884155e-4,
        "DEC": 2.384144e-9,
        "DWC": -3.745465e-7,
        "GRC": -1.646941e-2,
        "GIC": -2.894650e-8,
    }

    for key in exp_LIST:
        assert_allclose(exp_LIST[key], obs_LIST[key], rtol=1e-8)

    exp_bg_string = io.BytesIO(
        b""" 9.143204-3 1.601509+1-3.312654-7-3.460776+8-3.947879-5-1.646877-5
 1.009363-5-1.861342-7-1.613360+7-7.549728+5-3.064120+9-2.961641+0
-4.390193-5 2.303605+0-4.630212+5-3.237353+1-4.037885+4-3.231304+0"""
    )
    exp_bg = dict(
        zip(
            ("A1", "A2", "A3", "A4", "B1", "B2"),
            array_from_ENDF(exp_bg_string).transpose(),
        )
    )

    for key in exp_bg:
        assert_array_equal(exp_bg[key], obs_bg[key])


def test_resolved_r_matrix_kbk_kps():

    obs_3 = library.structure[nuc1002]["data"][nuc1002]["resolved"][-1][2][3.0]
    obs_4 = library.structure[nuc1002]["data"][nuc1002]["resolved"][-1][2][-4.0]

    exp_3 = array_from_ENDF(
        io.BytesIO(
            b""" 1.960831+3-1.053619+4          0          0          3          1
 3.941056-6-0.524089+0-2.023965-8 0.000000+0 0.000000+0 0.000000+0
 0.000000+0 0.000000+0          0          0          0          1
 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0
 0.000000+0 0.000000+0          0          0          0          1
 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0
 0.000000+0 0.000000+0          0          0          1          1
 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0
 0.000000+0 0.000000+0          0          0          1          3
          3          2                                            """
            + b"""
 4.016335+2 2.076736-3 4.668090-5 9.776415+2-3.940740+2-2.296483+8
 0.000000+0 0.000000+0          0          0          3         10
          3          1          6          2         10          3
-4.803282+6-1.114539-5 9.465304+2-1.436769-9 7.889727+2 4.824983+9
 4.020763+6 2.308443-6-4.188441-2 1.778263+8-3.408683+7 2.845463+7
 3.371147+1 2.054714+3-2.746606-3-9.635977-6-1.387257-2 7.042637+0
 6.917628+9-2.912896-7                                            """
        )
    )

    ch0_obs = obs_3["ch0"]

    # lbk = 3: [MAT,2,151/ ED, EU, 0, 0, LBK, 1/ R0, SO, GA, 0.0, 0.0, 0.0]LIST
    ch0_exp = dict(
        zip(("ED", "EU", 0, 0, "LBK", 0, "R0", "SO", "GA"), (exp_3[:2].flat))
    )
    ch0_exp["LPS"] = 0.0
    del ch0_exp[0]

    ch1_obs = obs_3["ch1"]

    # lbk = 0: no data
    # lps = 1: [MAT,2,151/ 0.0, 0.0,   0,   0, LPS,   1/
    #                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0]LIST
    #          [MAT,2,151/ 0.0,0.0,0,0, NR, NP/points/PSR(E)]TAB1
    #                      0.0,0.0,0,0, NR, NP/points/PSI(E)]TAB1

    ch1_exp = {
        "PSI": {
            "intpoints": [3.0, 6.0, 10.0],
            "intschemes": [1.0, 2.0, 3.0],
            "e_int": exp_3[13:17].flat[:-4:2],
            "PSI": exp_3[13:17].flat[1:-4:2],
        },
        "PSR": {
            "intpoints": 3.0,
            "intschemes": 2.0,
            "e_int": exp_3[10].flat[::2],
            "PSR": exp_3[10].flat[1::2],
        },
        "LBK": 0.0,
        "LPS": 1.0,
    }
    for key in ch1_exp:
        if isinstance(ch1_exp[key], dict):
            for intkey in ch1_exp[key]:
                assert_array_equal(ch1_exp[key][intkey], ch1_obs[key][intkey])
        else:
            assert_array_equal(ch1_obs[key], ch1_exp[key])
    assert_equal(ch0_obs, ch0_exp)

    lbk1_obs = obs_4["ch0"]
    lbk2_obs = obs_4["ch1"]
    lbk_exp = array_from_ENDF(
        io.BytesIO(
            b""" 0.000000+0 0.000000+0          0          0          1          1
 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0
 0.000000+0 0.000000+0          0          0          4          5
          1          2          2          1          3          3
          5          4                                            """
            + b"""
 0.000000+0 0.000000+0          0          0          2          6
          4          1          6          2                      """
            + b"""
 1.960831+3-1.053619+4          0          0          2          1
 1.298772+5-2.834965+3 8.381641-6 3.338353+5-1.012675-7 0.000000+0
"""
        )
    )

    # lbk = 1: [MAT,2,151/ 0.0, 0.0,   0,   0, LBK,   1/
    #                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0]LIST
    #          [MAT,2,151/ 0.0,0.0,0,0, NR, NP/points/RBR(E)]TAB1
    #                      0.0,0.0,0,0, NR, NP/points/RBI(E)]TAB1
    lbk1_exp = {
        "LBK": 1,
        "RBR": {"intpoints": [1, 2, 3, 5], "intschemes": [2, 1, 3, 4]},
        "RBI": {"intpoints": [4, 6], "intschemes": [1, 2]},
    }

    for key in lbk1_exp:
        if isinstance(lbk1_exp[key], dict):
            for intkey in lbk1_exp[key]:
                assert_array_equal(lbk1_exp[key][intkey], lbk1_obs[key][intkey])
        else:
            assert_array_equal(lbk1_obs[key], lbk1_exp[key])

    # lbk = 2: [MAT,2,151/ ED, EU,  0,  0, LBK,   1/
    #                      R0, R1, R2, S0,  S1, 0.0]LIST

    lbk2_exp = dict(
        zip(
            ("ED", "EU", 0, 0, "LBK", 0, "R0", "R1", "R2", "S0", "S1"),
            lbk_exp[7:9].flat,
        )
    )
    lbk2_exp["LBK"] = 2
    del lbk2_exp[0]

    assert_equal(lbk2_exp, lbk2_obs)


def test_resolved_r_matrix():
    pp_exp_a = array_from_ENDF(
        io.BytesIO(
            b"""1.685738+5 1.659888-5          1          7          0          0
 0.000000+0 0.000000+0          1          3          2          1
 0.000000+0 0.000000+0          2          0         24          4
 7.916271+6-3.532347-6 4.469905+7-2.134022+4-3.500000+0-3.500000+0
 5.307428-7          0          1          7          1         -1
 2.807643-8-4.478596+0 3.274758+3-2.760395+9 1.356440+3 3.447654+4
 4.790839-8          1         -1        800          1         -1"""
        )
    )
    pp_obs = library.structure[nuc1002]["data"][nuc1002]["resolved"][-1][3]
    pp_exp = dict(
        zip(
            ("MA", "MB", "ZA", "ZB", "IA", "IB", "Q", "PNT", "SHF", "MT", "PA", "PB"),
            pp_exp_a[3:7].reshape(2, 12).transpose(),
        )
    )
    pp_exp.update(dict(zip((0, 0, "IFG", "KRM", "NJS", "KRL"), pp_exp_a[1])))
    del pp_exp[0]

    ch_exp_a = array_from_ENDF(
        io.BytesIO(
            b"""-4.000000+0          0 1.914541-3-4.683290+5         12          2
 1.981937+3 9.740279-7-2.450194+5-1.304844+4 1.856158-9-1.218463-9
-4.097837+1-2.765873-9-0.913351+1 1.591290+5-2.379063+0 2.066455-6
 3.000000+0          0-2.924403-5-4.840218-1         12          2
 0.951092+5 7.932944-5-3.716218-8 4.007761-3 1.277498-6 2.041832+6
-7.375896+8-4.822942+4 4.491725+9 3.018430+8 2.238307-5-3.591395+9"""
        )
    )
    ch_exp = {}
    ch_exp[3.0] = dict(
        zip(("IPP", "L", "SCH", "BND", "APE", "APT"), ch_exp_a[4:6].transpose())
    )
    ch_exp[-4.0] = dict(
        zip(("IPP", "L", "SCH", "BND", "APE", "APT"), ch_exp_a[1:3].transpose())
    )
    ch_obs = library.structure[nuc1002]["data"][nuc1002]["resolved"][-1][2]

    gam_4_a = array_from_ENDF(
        io.BytesIO(
            b""" 2.949030-1 1.156625+7 7.255199-6          0          0          0
 4.453964+1 5.062864-5-1.110875-3          0          0          0
 2.208407-7 9.942677-6-3.134503-8          0          0          0"""
        )
    )
    gam_4_a = gam_4_a.transpose()
    gam_4_exp = {"ER": gam_4_a[0], "GAM": gam_4_a[1:3].transpose()}
    ch_exp[-4.0].update(gam_4_exp)

    gam_3_a = array_from_ENDF(
        io.BytesIO(
            b""" 5.088175-6-2.282938+0-4.236786-6          0          0          0
 8.930267-9-3.115607+8-2.521300-4          0          0          0
 3.978418+5 4.821547-6 3.110373-3          0          0          0"""
        )
    )
    gam_3_a = gam_3_a
    gam_3_a = gam_3_a.transpose()
    gam_3_exp = {"ER": gam_3_a[0], "GAM": gam_3_a[1:3].transpose()}

    ch_exp[3.0].update(gam_3_exp)

    for key in pp_exp:
        assert_array_equal(pp_obs[key], pp_exp[key])
    for spin_group in ch_exp:
        spin_group_obs = ch_obs[spin_group]
        spin_group_exp = ch_exp[spin_group]
        for key in spin_group_exp:
            assert_array_equal(spin_group_obs[key], spin_group_exp[key])


def test_xs():
    # Read in the data
    nuc_i = nucname.id(40192)
    library._read_xs(nuc40000, 2, nuc_i)
    library._read_xs(nuc40000, 600, nuc40040)

    # Manually find where the data should be reading from and check if it is
    # consistent with what the program is doing.
    exp_2_str = io.BytesIO(
        b""" 4.284918+3 6.292347+0          0          0          0          0
 4.047593+5-4.245658-8          0-4.651348+3          7         20
          6          4          9          2         12          1
         13          5         15          3         17          4
         20          1                                            """
    )
    exp_2_a = array_from_ENDF(exp_2_str)
    exp_2 = dict(
        zip(
            ("intpoints", "intschemes"),
            (exp_2_a[2:].flat[:14:2], exp_2_a[2:].flat[1:14:2]),
        )
    )
    obs_2 = library.structure[nuc40000]["data"][nuc_i]["xs"][2][0]

    exp_600_a = array_from_ENDF(
        io.BytesIO(
            b""" 4.193742+3 6.287192+0          0          0          0          0
 3.863437-5-7.373532-7          0 8.675483-1          5         20
          4          1          8          2         12          3
         16          4         20          5                     """
        )
    )

    exp_600 = dict(
        zip(
            ("intpoints", "intschemes"),
            (exp_600_a[2:].flat[:-2:2], exp_600_a[2:].flat[1:-1:2]),
        )
    )
    obs_600 = library.structure[nuc40000]["data"][nuc40040]["xs"][600][0]

    for key in exp_2:
        assert_array_equal(obs_2[key], exp_2[key])
        assert_array_equal(obs_600[key], exp_600[key])

    # Heck, why not check the flags too?
    obs_600_flags = library.structure[nuc40000]["data"][nuc40040]["xs"][600][1]
    exp_600_flags = dict(zip(("QM", "QI", 0, "LM", "NR", "NP"), exp_600_a[1]))  #
    exp_600_flags.update({"ZA": 4.004e3, "AWR": 6.287192})
    del exp_600_flags[0]
    assert_equal(obs_600_flags, exp_600_flags)


def test_xs_data_without_res():
    nuc8017 = nucname.id(8017)
    library._read_res(nuc8017)
    library._read_xs(nuc8017, 4, nuc8017)


def test_isomeric():
    nuc61148m = nucname.id("Pm148m")
    library._read_res(nuc61148m)
    assert library.structure[nuc61148m]["matflags"]["LIS0"] == 1
    assert nuc61148m in library.structure


def test_u235():
    try_download(
        "http://t2.lanl.gov/nis/data/data/ENDFB-VII.1-neutron/U/235",
        "U235.txt",
        "1b71da3769d8b1e675c3c579ba5cb2d3",
    )

    u235 = Library("U235.txt")
    nuc = 922350000
    u235._read_res(nuc)
    u235._read_xs(nuc, 37)
    exp_a = array_from_ENDF(
        io.BytesIO(
            b""" 9.223500+4 2.330248+2          0          0          0          0
-1.788560+7-1.788560+7          0          0          1          6
          6          2                                            """
            + b"""
 1.796240+7 5.05980-10 1.800000+7 3.810030-7 1.850000+7 8.441785-5
 1.900000+7 2.387410-4 1.950000+7 1.348763-3 2.000000+7 4.785594-3
"""
        )
    )
    obs = u235.structure[nuc]["data"][nuc]["xs"][37][0]
    exp = {
        "intpoints": 6,
        "intschemes": 2,
        "e_int": exp_a[3:5].flat[::2],
        "xs": exp_a[3:5].flat[1::2],
    }

    for key in obs:
        assert_array_equal(obs[key], exp[key])


# Test ENDF Data Source
@ignore_future_warnings
def test_int_hist():
    exp_Eint = np.array([1, 4, 10, 20])
    exp_xs = np.array([15, 12, -7, 10])
    obs = library.integrate_tab_range(1, exp_Eint, exp_xs)
    exp = (3 * 15 + 6 * 12 + 10 * -7) / 19.0
    assert_allclose(exp, obs, rtol=1e-12)


def test_int_hist_interpolation():
    exp_Eint = np.array([1, 4, 10, 20])
    exp_xs = np.array([15, 12, -7, 10])
    obs = library.integrate_tab_range(1, exp_Eint, exp_xs, low=2, high=15)
    exp = (2 * 15 + 6 * 12 + 5 * -7) / 13.0
    assert_allclose(exp, obs, rtol=1e-12)


def test_int_hist_only_interpolate_one_endpoint():
    endfds = ENDFDataSource(tape1path)
    obs = endfds.integrate_dst_group(
        (1, 5),
        np.array([2, 5, 8]),
        {2: 1, 5: 1, 8: 1},
        np.array([0.0, 2, 4, 6, 8]),
        np.array([0.0, 1, 0, 0, 0]),
    )
    exp = 0.5
    assert_equal(exp, obs)


def test_int_linlin():
    exp_Eint = np.array([1, 4, 10, 20])
    exp_xs = np.array([15, 12, -7, 10])
    obs = library.integrate_tab_range(2, exp_Eint, exp_xs)
    exp = (3 * 13.5 + 6 * 2.5 + 10 * 1.5) / 19.0
    assert_allclose(exp, obs, rtol=1e-12)
    return exp


def test_int_linlin_interpolation():
    exp_Eint = np.array([1, 4, 10, 20.0])
    exp_xs = np.array([15, 12, -7, 10.0])
    obs = library.integrate_tab_range(2, exp_Eint, exp_xs, low=2, high=15)
    exp = (2 * 13 + 6 * 2.5 + 5 * (-7 + 1.5) / 2) / 13.0
    assert_allclose(exp, obs, rtol=1e-12)


def test_int_linlin_interpolation_2():
    endfds = ENDFDataSource(tape1path)
    obs = endfds.integrate_dst_group(
        (1, 5),
        np.array([3, 5, 8]),
        {3: 2, 5: 2, 8: 1},
        np.array([0.0, 2, 4, 6, 8]),
        np.array([0.0, 1, 0, 0, 0]),
    )
    exp = (0.75 + 1) / 4
    assert_equal(exp, obs)


def test_int_linlin_only_interpolate_one_endpoint():
    exp_Eint = np.array([1, 4, 10, 20.0])
    exp_xs = np.array([15, 12, -7, 10.0])
    obs = library.integrate_tab_range(2, exp_Eint, exp_xs, low=1, high=15)
    exp = (3 * 13.5 + 6 * 2.5 + 5 * (-7 + 1.5) / 2) / 14.0
    assert_allclose(exp, obs, rtol=1e-12)


def test_int_linlog():
    exp_Eint = np.array([1, e, e, e**2])
    exp_xs = np.array([1, 3, 3, 0])
    obs = library.integrate_tab_range(3, exp_Eint, exp_xs)
    exp = (e + 1 + 3 * e**2 - 6 * e) / (e**2 - 1)
    assert_allclose(exp, obs, rtol=1e-12)


def test_int_linlog_interpolation():
    exp_Eint = np.array([1, e**2, e**4, e**6])
    exp_xs = np.array([0, 2, 4, 6.0])
    obs = library.integrate_tab_range(3, exp_Eint, exp_xs, low=e, high=e**5)
    exp = 4 * e**5 / (e**5 - e)
    assert_allclose(exp, obs, rtol=1e-12)


def test_int_linlog_only_interpolate_one_endpoint():
    exp_Eint = np.array([1, e**2, e**4, e**6])
    exp_xs = np.array([0, 2, 4, 6.0])
    obs = library.integrate_tab_range(3, exp_Eint, exp_xs, low=1, high=e**5)
    exp = (1 + 4 * e**5) / (e**5 - 1)
    assert_allclose(exp, obs, rtol=1e-12)


def test_int_loglin():
    exp_Eint = np.array([1.0, 2.0, 2.0, 3.0])
    exp_xs = np.array([1, e, e**2, e])
    obs = library.integrate_tab_range(4, exp_Eint, exp_xs)
    exp = (e**2 - 1) / 2
    assert_allclose(exp, obs, rtol=1e-12)


def test_int_loglin_interpolation():
    exp_Eint = np.array([0, 2, 4, 6])
    exp_xs = np.array([1, e**2, e**4, e**6])
    obs = library.integrate_tab_range(4, exp_Eint, exp_xs, low=1, high=5)
    exp = (e**5 - e) / (5 - 1)
    assert_allclose(exp, obs, rtol=1e-12)


def test_int_loglin_only_interpolate_one_endpoint():
    exp_Eint = np.array([0, 2, 4, 6], dtype="float64")
    exp_xs = np.array([1, e**2, e**4, e**6])
    obs = library.integrate_tab_range(4, exp_Eint, exp_xs, low=2, high=5)
    exp = (e**5 - e**2) / (5 - 2)
    assert_allclose(exp, obs, rtol=1e-12)


def test_int_loglog():
    exp_Eint = np.array([1.0, 2.0, 2.0, 3.0])
    exp_xs = np.array([1 / e, 4 / e, (e**2) / 4, (e**2) / 9])
    obs = library.integrate_tab_range(5, exp_Eint, exp_xs)
    exp = (7 / (3 * e) + (e**2) / 6) / 2.0
    assert_allclose(exp, obs, rtol=1e-12)


def test_int_loglog_interpolation():
    # ln y = 2 ln x + 1
    # y = e * x ** 2
    # integral of y = e/3 * x**3
    exp_Eint = np.array([1, 3, 5, 7], dtype="float64")
    exp_xs = np.array([e, 9 * e, 25 * e, 49 * e])
    obs = library.integrate_tab_range(5, exp_Eint, exp_xs, low=2, high=6)
    exp = e / 3 * (6**3 - 2**3) / (6 - 2)
    assert_allclose(exp, obs, rtol=1e-12)


def test_int_loglog_only_interpolate_one_endpoint():
    # ln y = 2 ln x + 1
    # y = e * x ** 2
    # integral of y = e/3 * x**3
    exp_Eint = np.array([1, 3, 5, 7], dtype="float64")
    exp_xs = np.array([e, 9 * e, 25 * e, 49 * e])
    obs = library.integrate_tab_range(5, exp_Eint, exp_xs, low=2, high=5)
    exp = e / 3 * (5**3 - 2**3) / (5 - 2)
    assert_allclose(exp, obs, rtol=1e-12)


def test_discretize():
    from os.path import isfile

    try:
        import urllib.request as urllib
    except ImportError:
        import urllib

    try_download(
        "http://t2.lanl.gov/nis/data/data/ENDFB-VII.1-neutron/Ni/59", "Ni59.txt", ""
    )

    endfds = ENDFDataSource("Ni59.txt")
    nonelastic_rx = endfds.reaction("Ni59", "nonelastic")
    nonelastic_rx["dst_group_struct"] = np.logspace(7, -5, 33)
    nonelastic_c = endfds.discretize("Ni59", "nonelastic")
    exp = [
        0.54334609294912528,
        0.21206255570566626,
        0.079089998725708668,
        0.039061531003500925,
        0.056193960028285306,
        0.062581135526972767,
        0.086088778452663009,
        0.1519375415918513,
        0.015156525895127398,
        0.18503957567677801,
        0.0039443417078627837,
        0.082573739674287688,
        17.523219940338304,
        0.97176481236488554,
        0.60307330340022303,
        0.71684581122716162,
        0.99386518962022252,
        1.4726882603418707,
        2.2391970686479672,
        3.405589441800994,
        5.2453926977834398,
        8.0731410528834182,
        12.384026334168054,
        19.175694435799141,
        29.334824378652982,
        45.254982026071197,
        74.217617672501689,
        162.26091389706099,
        218.90153743636509,
        312.62178192130619,
        590.40136068709603,
        724.64216445611373,
    ]
    assert_array_almost_equal(nonelastic_c, exp)
    os.remove(Ni59.txt)


def test_photoatomic():
    try_download(
        "https://www-nds.iaea.org/fendl30/data/atom/endf/ph_3000_30-Zn.txt",
        "Zn.txt",
        "e6bda2fe6125ad9a8d51a6405ae9cc2a",
    )
    photondata = Library("Zn.txt")
    xs_data = photondata.get_xs(300000000, 501)[0]
    Eints, sigmas = xs_data["e_int"], xs_data["xs"]
    assert_equal(len(Eints), 1864)
    assert_equal(len(sigmas), 1864)
    assert_array_equal(Eints[0:5], [1.0, 1.109887, 1.217224, 1.2589, 1.334942])
    assert_array_equal(
        Eints[-5:],
        [6.30960000e10, 7.94328000e10, 7.94330000e10, 8.00000000e10, 1.00000000e11],
    )
    assert_array_almost_equal(
        sigmas[0:5], [0.00460498, 0.00710582, 0.01047864, 0.01210534, 0.01556538]
    )
    assert_array_almost_equal(
        sigmas[-5:], [6.71525, 6.716781, 6.716781, 6.716829, 6.717811]
    )


def test_evaluation_neutron():
    try_download(
        "http://t2.lanl.gov/nis/data/data/ENDFB-VII.1-neutron/U/235",
        "U235.txt",
        "1b71da3769d8b1e675c3c579ba5cb2d3",
    )
    u235 = Evaluation("U235.txt", verbose=False)
    u235.read()

    # Check descriptive data
    assert hasattr(u235, "info")
    assert u235.info["library"] == ("ENDF/B", 7, 1)
    assert u235.info["sublibrary"] == 10
    assert u235.info["format"] == 6
    assert u235.material == 9228
    assert u235.projectile["mass"] == 1.0
    assert u235.target["fissionable"]
    assert u235.target["mass"] == 233.0248
    assert u235.target["zsymam"] == " 92-U -235 "
    assert not u235.target["stable"]
    assert u235.target["isomeric_state"] == 0

    # Check components of fission energy release
    delayed_betas = np.array([[6.5e06, -7.5e-02, 0.0], [5.0e04, 7.5e-03, 0.0]])
    assert_array_almost_equal(
        u235.fission["energy_release"]["delayed_betas"], delayed_betas
    )

    # Check prompt, delayed, and total nu
    E = np.logspace(-4, 6, 10)
    nu_t = np.array(
        [
            2.4367,
            2.4367,
            2.4367,
            2.4367,
            2.4367,
            2.43586422,
            2.4338,
            2.4338,
            2.43001824,
            2.532706,
        ]
    )
    nu_p = np.array(
        [
            2.42085,
            2.42085,
            2.42085,
            2.42085,
            2.42085,
            2.42001364,
            2.41794193,
            2.41784842,
            2.41331824,
            2.516006,
        ]
    )
    nu_d = np.array(
        [
            0.01585,
            0.01585,
            0.01585,
            0.01585,
            0.01585005,
            0.01585061,
            0.01585789,
            0.01595191,
            0.0167,
            0.0167,
        ]
    )
    assert_array_almost_equal(u235.fission["nu"]["total"](E), nu_t)
    assert_array_almost_equal(u235.fission["nu"]["prompt"](E), nu_p)
    assert_array_almost_equal(u235.fission["nu"]["delayed"]["values"](E), nu_d)

    # Reactions
    assert 37 in u235.reactions
    assert 38 in u235.reactions
    assert 89 in u235.reactions

    # Check reaction cross section
    r = u235.reactions[80]  # (n,n30)
    assert r.xs.x[0] == 2309870.0
    E = np.linspace(r.xs.x[0], r.xs.x[-1], 10)
    sigma = np.array(
        [
            0.0,
            0.01433554,
            0.01672903,
            0.01649244,
            0.01556801,
            0.01457872,
            0.01341327,
            0.01215489,
            0.01094848,
            0.00987096,
        ]
    )
    assert_array_almost_equal(r.xs(E), sigma)

    # Check reaction angular distribution
    assert r.angular_distribution.center_of_mass
    assert r.angular_distribution.type == "legendre"
    assert len(r.angular_distribution.energy) == 14
    mu = np.linspace(-1.0, 1.0, 10)
    p = np.array(
        [
            0.22002253,
            0.40547169,
            0.45029205,
            0.59018363,
            0.6086658,
            0.48871545,
            0.45452016,
            0.55752484,
            0.55576971,
            0.51885337,
        ]
    )
    assert_array_almost_equal(r.angular_distribution.probability[5](mu), p)


def test_evaluation_decay():
    try_download(
        "http://t2.lanl.gov/nis/data/endf/decayVII.1/092_U_233",
        "U233.txt",
        "3db23dc650bae28eabb92942dd7d0de5",
    )
    u233 = Evaluation("U233.txt", verbose=False)
    u233.read()

    assert hasattr(u233, "info")
    assert u233.info["library"] == ("ENDF/B", 7, 1)
    assert u233.info["sublibrary"] == 4
    assert u233.material == 3513
    assert u233.target["mass"] == 231.0377
    assert u233.target["zsymam"] == " 92-U -233 "
    assert not u233.target["stable"]

    assert u233.decay["half_life"] == (5023970000000.0, 6311520000.0)
    assert u233.decay["energies"] == [
        (5043.237, 536.3191),
        (1110.218, 107.6781),
        (4888351.0, 28967.6),
    ]
    assert len(u233.decay["modes"]) == 1
    assert u233.decay["modes"][0]["branching_ratio"] == (1.0, 0.0)
    assert u233.decay["modes"][0]["energy"] == (4908500.0, 1200.0)
    assert u233.decay["modes"][0]["type"] == "alpha"

    for s in ["e-", "alpha", "xray", "gamma"]:
        assert s in u233.decay["spectra"]
    assert u233.decay["spectra"]["e-"]["energy_average"] == (5043.237, 536.3191)
    alpha_spectrum = u233.decay["spectra"]["alpha"]["discrete"]
    assert alpha_spectrum[-1]["energy"] == (4824200.0, 1200.0)
    assert alpha_spectrum[-1]["intensity"] == (0.843, 0.006)


def test_evaluation_photoatomic():
    try_download(
        "https://www-nds.iaea.org/fendl30/data/atom/endf/ph_3000_30-Zn.txt",
        "Zn.txt",
        "e6bda2fe6125ad9a8d51a6405ae9cc2a",
    )
    zn = Evaluation("Zn.txt", verbose=False)
    zn.read()

    assert zn.info["library"] == ("ENDF/B", 6, 0)
    assert zn.info["sublibrary"] == 3
    assert zn.info["format"] == 6
    assert zn.projectile["mass"] == 0.0

    assert 501 in zn.reactions
    assert 515 in zn.reactions
    assert 540 in zn.reactions

    # Check cross section
    coherent = zn.reactions[502]
    assert len(coherent.xs.x) == 788
    E = np.logspace(0, 11, 10)
    xs = np.array(
        [
            4.60408600e-03,
            2.87607904e00,
            3.00145792e02,
            3.43420255e02,
            8.18532007e00,
            3.58892209e-02,
            1.29382964e-04,
            4.64984378e-07,
            1.67106649e-09,
            6.00550000e-12,
        ]
    )
    assert_array_almost_equal(coherent.xs(E), xs)

    # Check scattering factor
    ff = coherent.scattering_factor
    x = np.logspace(-3, 9, 10)
    y = np.array(
        [
            3.00000000e01,
            2.98789879e01,
            1.54855832e01,
            4.24200000e-01,
            1.40737115e-05,
            7.34268836e-10,
            8.06893048e-14,
            1.06881338e-17,
            1.27450979e-21,
            1.40780000e-25,
        ]
    )
    assert_array_almost_equal(ff(x), y)


def test_evaluation_electroatomic():
    try:
        try_download(
            "http://t2.lanl.gov/nis/data/data/ENDFB-VII-electroatomic/Mo/nat",
            "Mo.txt",
            "2139a23258c517ae3bfa5f2cc346da4c",
        )
    except:
        raise SkipTest(
            "http://t2.lanl.gov/nis/data/data/ENDFB-VII.1-neutron not available"
        )

    mo = Evaluation("Mo.txt", verbose=False)
    mo.read()

    assert mo.info["laboratory"].startswith("LLNL")
    assert "Cullen" in mo.info["author"]
    assert mo.info["library"] == ("ENDF/B", 6, 0)
    assert mo.info["sublibrary"] == 113

    assert 526 in mo.reactions
    assert 527 in mo.reactions
    assert 543 in mo.reactions

    # Check bremsstrahlung cross section
    brem = mo.reactions[527]
    E = np.logspace(1, 11, 10)
    xs = np.array(
        [
            4276.9,
            6329.94427327,
            5350.43318579,
            1818.71320602,
            467.65208672,
            376.17541997,
            434.15352828,
            484.99421137,
            535.39666445,
            591.138,
        ]
    )
    assert_array_almost_equal(brem.xs(E), xs)

    # Check bremsstrahlung secondary distributions
    assert brem.products[0]["za"] == 11  # electrons
    E = np.logspace(1, 11)
    assert np.all(brem.products[0]["yield"](E) == 1.0)
    eout = np.array(
        [
            0.1,
            0.133352,
            0.165482,
            0.228757,
            0.316228,
            0.421697,
            0.523299,
            0.723394,
            1.0,
            1.41421,
            2.0,
            2.82843,
            4.0,
            5.65685,
            8.0,
            9.9,
            10.0,
        ]
    )
    assert_array_almost_equal(brem.products[0]["energy_out"][0], eout)


def test_evaluation_relaxation():
    try_download(
        "http://t2.lanl.gov/nis/data/endf/relaxation/Xe-relax",
        "Xe.txt",
        "40ecb69da6a45f992918a98da4d98ba0",
    )
    xe = Evaluation("Xe.txt", verbose=False)
    xe.read()

    assert xe.info["laboratory"].startswith("LLNL")
    assert "Cullen" in xe.info["author"]
    assert xe.info["library"] == ("ENDF/B", 6, 0)
    assert xe.info["sublibrary"] == 6

    # Levels
    data = xe.atomic_relaxation
    assert len(data) == 17
    assert "K" in data
    assert "O3" in data

    assert data["K"]["binding_energy"] == 34556.0
    assert data["M5"]["number_electrons"] == 6.0
    assert data["L1"]["transitions"][6] == ("N2", None, 5256.21, 0.00299258)
    assert data["N2"]["transitions"][12] == ("N5", "O3", 81.95, 0.00224012)
    assert data["O3"]["transitions"] == []


def test_contents_regexp():
    testInput = """A line like this will never happen in any ENDF-6 formatted file!!!
This line looks like a (MF,MT)=(1,451) line but NOT!              012  1451 1 34
This line should be recognized as a valid (MF,MT)=(1,451) doc line 123 1451  456
It is now checked whether an integer number starts with '0' or not0123 1451  457
    (perhaps not a valid FORTRAN77 INTEGER!)                       123 1451  458
The following two lines are valid CONT Records in ENDF-6 format    123 1451  459
but are NOT recognized as CONTENTS lines in (MF,MT)=(1,451).       123 1451  460
12345689901 2345678901  345678901   45678901    5678901     6789012345 1451  461
          1          2          3          4          5          62345 1451  462
The following two lines should be recognized as the CONTENTS lines.
   (Neither the ranges of the numbers (MF:1-99, MT:1-999)
    nor the position of the CONTENTS lines are checked)
                        345678901   45678901    5678901     6789012345 1451  463
                                3          4          5          62345 1451  464
The following three lines are DATA lines but now it is not searched for.
 9.87654+32-6.54321098  345678901   45678901    5678901     6789012345 1451  465
-1.234567+1+1.23456-22          3          4          5          62345 1451  466
 2.34567890-3.456789+1          3          4          5          62345 1451  467
LRP can be NEGATIVE!
 2.34567890-3.456789+1         -3          4          5          62345 1451  468
How about the following three???
2.34567890 -3.456789+1          3          4          5          62345 1451  469
 2.34567890-3.456789+1         3          4          5          6 2345 1451  470
 2.34567890-3.456789+1         3          4          5          6 2345 1451 471 
"""
    from pyne.endf import FILE1_R, SPACEINT11_R

    for line in testInput.splitlines():
        # print("'{}'".format(line))
        if FILE1_R.match(line):
            parts = [line[i : i + 11] for i in range(0, 66, 11)]
            # CONTENTS line
            if (
                parts[0] + parts[1] == " " * 22
                and SPACEINT11_R.match(parts[2])
                and SPACEINT11_R.match(parts[3])
                and SPACEINT11_R.match(parts[4])
                and SPACEINT11_R.match(parts[5])
            ):
                print("1451CONT: '{}'".format(line))
            # DOCUMENTS line
            else:
                print("1451DOCS: '{}'".format(line))
        else:
            print("OTHER   : '{}'".format(line))


if __name__ == "__main__":
    nose.runmodule()
