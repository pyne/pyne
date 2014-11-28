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
