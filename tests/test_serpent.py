import os
import warnings

import numpy as np
from nose.tools import assert_equal, assert_true
from numpy.testing import assert_array_equal

from pyne.utils import QAWarning

warnings.simplefilter("ignore", QAWarning)
from pyne import serpent
from pyne import nucname
from pyne.material import Material


def test_parse_res1():
    res = serpent.parse_res("sample_res.m")
    rank0 = res["IDX"]
    assert_equal(res["idx"] + 1, rank0)
    for key in res:
        if isinstance(res[key], np.ndarray):
            assert_equal(res[key].shape[0], rank0)

    # Check values
    assert_array_equal(res["SIX_FF_ETA"][1], [1.16446e00, 0.00186])
    assert_array_equal(res["PEAKF10"][rank0 - 1], [12, 11, 1.09824e00, 0.01768])


def test_parse_res2():
    res = serpent.parse_res("serp2_res.m")
    rank0 = res["IDX"]
    assert_equal(res["idx"] + 1, rank0)
    for key in res:
        if isinstance(res[key], np.ndarray):
            assert_equal(res[key].shape[0], rank0)

    # Check values
    assert_array_equal(res["MEAN_POP_SIZE"][0], [1.00140e02, 0.00359])


def test_parse_dep1():
    """
    Tests the parse_dep function with a serpent 1 _dep.m
    file.  In this, the _VOLUME variable is a float.
    """
    dep = serpent.parse_dep("sample1_dep.m")
    nuc_set = set([nucname.id(int(nuc)) for nuc in dep["ZAI"][:-2]])
    shape = (len(dep["ZAI"]), len(dep["DAYS"]))
    for key in dep:
        if key.endswith("_VOLUME"):
            assert_true(isinstance(dep[key], float))
        elif key.endswith("_MATERIAL"):
            assert_equal(len(dep[key]), shape[1])
            for n in range(shape[1]):
                assert_equal(nuc_set, set(dep[key][n].comp.keys()))
        elif key.startswith("MAT_") and key.endswith("_BURNUP"):
            assert_equal(dep[key].shape, (1, shape[1]))
        elif key.startswith("MAT_") or key.startswith("TOT_"):
            assert_equal(dep[key].shape, shape)

    # Check values
    assert_array_equal(dep["BU"], [0.00000e00, 8.40000e01, 1.68000e02])
    assert_equal(dep["i952421"], 123)
    assert_array_equal(dep["MAT_fuelp1r2_H"][3], [0.00000e00, 5.56191e-11, 3.22483e-10])


def test_parse_dep2():
    """
    Tests the parse_dep function with a simple _dep.m file,
    sample2_dep.m, which was generated using the 2d pin burnup
    example from the serpent wiki, changed to have fewer depsteps.
    In this, the _VOLUME variable is a numpy array.
    http://serpent.vtt.fi/mediawiki/index.php/2D_PWR_pin-cell_burnup_example
    """
    dep = serpent.parse_dep("sample2_dep.m")
    nuc_set = set([nucname.id(int(nuc)) for nuc in dep["ZAI"][:-2]])
    shape = (len(dep["ZAI"]), len(dep["DAYS"]))
    for key in dep:
        if key.endswith("_VOLUME"):
            assert_equal(len(dep[key]), shape[1])
        elif key.endswith("_vol"):
            assert_equal(len(dep[key]), shape[1])
        elif key.endswith("_MATERIAL"):
            assert_equal(len(dep[key]), shape[1])
            for n in range(shape[1]):
                assert_equal(nuc_set, set(dep[key][n].comp.keys()))
        elif key.endswith("_BURNUP"):
            assert_equal(len(dep[key]), shape[1])
        elif key.startswith("MAT_") or key.startswith("TOT_"):
            assert_equal(len(dep[key]), shape[0])
            assert_equal(len(dep[key][0]), shape[1])

    # Check values
    assert_array_equal(dep["BU"], [0.00000e00, 1.00000e-01, 1.00000e00])
    assert_equal(dep["i621510"], 22)
    assert_array_equal(dep["MAT_fuel_A"][6], [0.00000e00, 1.73176e05, 4.96485e06])


def test_parse_det1():
    det = serpent.parse_det("sample_det.m")
    for key in det:
        if "_" in key:
            assert_true(isinstance(det[key], int))
        elif isinstance(det[key], np.ndarray):
            assert_true(det[key].shape[1] in [3, 13])

    # Check values
    assert_array_equal(
        det["DETphi"][6], [7, 7, 1, 1, 1, 1, 1, 1, 1, 1, 2.92709e-02, 0.00857, 16768]
    )
    assert_array_equal(det["DETphiE"][-3], [1.49182e01, 1.69046e01, 1.49182e01])


def test_parse_det2():
    det = serpent.parse_det("serp2_det.m")
    for key in det:
        if "_" in key:
            assert_true(isinstance(det[key], int))
        elif isinstance(det[key], np.ndarray):
            assert_true(det[key].shape[1] in [3, 13])

    # Check values
    assert_array_equal(
        det["DET1"][4], [5, 1, 5, 1, 1, 1, 1, 1, 1, 1, 1, 5.11865e05, 0.00417]
    )
    assert_array_equal(det["DET1E"][-3], [5.25306e-05, 3.80731e-03, 1.92992e-03])
