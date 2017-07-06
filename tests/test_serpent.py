import warnings

import numpy as np
from nose.tools import assert_equal, assert_true
from numpy.testing import assert_array_equal

from pyne.utils import QAWarning
warnings.simplefilter("ignore", QAWarning)
from pyne import serpent
from pyne import nucname


def test_parse_res1():
    res = serpent.parse_res('sample_res.m')
    rank0 = res['IDX']
    assert_equal(res['idx']+1, rank0)
    for key in res:
        if isinstance(res[key], np.ndarray):
            assert_equal(res[key].shape[0], rank0)

    # Check values
    assert_array_equal(res['SIX_FF_ETA'][1],  [1.16446E+00, 0.00186])
    assert_array_equal(res['PEAKF10'][rank0-1], [12, 11, 1.09824E+00, 0.01768])


def test_parse_res2():
    res = serpent.parse_res('serp2_res.m')
    rank0 = res['IDX']
    assert_equal(res['idx']+1, rank0)
    for key in res:
        if isinstance(res[key], np.ndarray):
            assert_equal(res[key].shape[0], rank0)

    # Check values
    assert_array_equal(res['MEAN_POP_SIZE'][0], [1.00140E+02, 0.00359])


def test_parse_dep1():
    dep = serpent.parse_dep('sample_dep.m')
    nuc_set = set([nucname.id(int(nuc)) for nuc in dep['ZAI'][:-2]])
    shape = (len(dep['ZAI']), len(dep['DAYS'])) 
    for key in dep:
        if key.endswith('_VOLUME'):
            assert_true(isinstance(dep[key], float))
        elif key.endswith('_MATERIAL'):
            assert_equal(len(dep[key]), shape[1])
            for n in range(shape[1]):
                assert_equal(nuc_set, set(dep[key][n].comp.keys()))
        elif key.startswith('MAT_') and key.endswith('_BURNUP'):
            assert_equal(dep[key].shape, (1, shape[1]))
        elif key.startswith('MAT_') or key.startswith('TOT_'):
            assert_equal(dep[key].shape, shape)

    # Check values
    assert_array_equal(dep['BU'], [0.00000E+00, 8.40000E+01, 1.68000E+02])
    assert_equal(dep['i952421'], 123)
    assert_array_equal(dep['MAT_fuelp1r2_H'][3], [0.00000E+00, 5.56191E-11, 3.22483E-10])


def test_parse_det1():
    det = serpent.parse_det('sample_det.m')
    for key in det:
        if '_' in key:
            assert_true(isinstance(det[key], int))
        elif isinstance(det[key], np.ndarray):
            assert_true(det[key].shape[1] in [3, 13])

    # Check values
    assert_array_equal(det['DETphi'][6], 
                       [7, 7, 1, 1, 1, 1, 1, 1, 1, 1, 2.92709E-02, 0.00857, 16768])
    assert_array_equal(det['DETphiE'][-3], [1.49182E+01, 1.69046E+01, 1.49182E+01])

def test_parse_det2():
    det = serpent.parse_det('serp2_det.m')
    for key in det:
        if '_' in key:
            assert_true(isinstance(det[key], int))
        elif isinstance(det[key], np.ndarray):
            assert_true(det[key].shape[1] in [3, 13])

    # Check values
    assert_array_equal(det['DET1'][4], 
        [5, 1, 5, 1, 1, 1, 1, 1, 1, 1, 1, 5.11865E+05, 0.00417])
    assert_array_equal(det['DET1E'][-3], [5.25306E-05, 3.80731E-03, 1.92992E-03])

def test_parse_coe1():
    coe = serpent.parse_coe('serp2.coe')

    # test some vales. Not all, of course.
    # check BU step 1, universe 2, fuel0+rod0 branch
    # P0 scattering matrix for 4 group MSRE
    assert_array_equal(coe[1]["2"]["fuel0"]["rod0"]["INF_S0"],
                       [3.19826E-01, 5.83480E-03, 6.21411E-09, 0.00000E+00,
                        0.00000E+00, 2.73182E-01, 5.47793E-03, 0.00000E+00, 
                        0.00000E+00, 0.00000E+00, 2.59886E-01, 1.05679E-02,
                        0.00000E+00, 0.00000E+00, 7.47597E-04, 2.73214E-01])
    assert coe[1]["2"]["fuel0"]["rod0"]["INF_KINF"] == 1.61559E+00

    # check that all branches are present
    # in all universes
    fuelBranches = ["fuel"+str(i) for i in range(15)]
    rodBranches = ["rod"+str(i) for i in range(5)]
    # only fuel=2 and moder=3 universes have GCs here
    for uni in ["2", "3"]:
        for fuelBr in fuelBranches:
            assert fuelBr in coe[1][uni]
        for rodBr in rodBranches:
            for thisFuel in fuelBranches:
                assert rodBr in coe[1][uni][thisFuel]
