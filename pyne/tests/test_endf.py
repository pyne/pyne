import warnings

import numpy as np
from numpy.testing import assert_array_equal # , assert_array_almost_equal

from pyne.endf import Library
from pyne.rxdata import DoubleSpinDict

import nose 
from nose.tools import assert_equal # , assert_not_equal, assert_raises, raises, assert_in


library = Library('endftest_small.txt')   


def test_mats():
    for mat_id in library.mats:
        assert_in(mat_id, [128, 131, 419])


def test_get():
    obs = library.get(419, 4, 2)

    exp = np.array([4.898421e+3, 6.768123e+0, 0, 
                    1, 0, 0, 2.123124e+6, 8.123142e-6, 
                    2.123212e+6, 8.231231e-6,
                    -2.231211e+6, 8.123421e-6])    
    badkey = library.get(111, 1, 1)
    assert_array_equal(exp, obs)
    assert_equal(badkey, False)


def test_unresolved_resonances_a():
    # Case A (ENDF Manual p.70)
    
    obs = library.mat131['rxdata']['unresolved']
    # print obs
    # assert(1==2)

    obs_D = obs[0][2][3.5,1,10]['D']
    exp_D = 2.110000e3

    obs_GNO = obs[0][2][3.5,1,10]['GNO']
    exp_GNO = 8.497500e-1
    
    new_obs = obs[1][2][2,2,1]

    new_obs_D = new_obs['D']
    new_exp_D = 2.101000e3

    new_obs_GNO = new_obs['GNO']
    new_exp_GNO = 7.088000e-1

    assert_array_equal(exp_D, obs_D)
    assert_array_equal(exp_GNO, obs_GNO)

    assert_array_equal(new_exp_D, new_obs_D)
    assert_array_equal(new_exp_GNO, new_obs_GNO)


def test_unresolved_resonances_b():
    # Case B (ENDF Manual p. 70)
    obs = library.mat419['rxdata']['unresolved'][0][2]

    obs_ES = obs[3.5,0,419]['ES']
    exp_ES = 100 * np.array([10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]) 

    obs_419_GF = obs[3.5,0,419]['GF']
    exp_419_GF = np.array([2.000000e3, 2.100000e3, 2.200000e3, 2.300000e3, 2.400000e3, 2.500000e3,
                        2.600000e3, 2.700000e3, 2.800000e3, 2.900000e3, 3.000000e3])
    
    assert_array_equal(exp_ES, obs_ES)
    assert_array_equal(exp_419_GF, obs_419_GF)


def test_unresolved_resonances_c():
    # Case C (ENDF Manual p. 70)

    obs = library.mat128['rxdata']['unresolved'][0][2][3.5,1,4]

    obs_ES = obs['ES']
    exp_ES = np.array([1.74e3, 2.04e3, 3.04e3])

    obs_D = obs['D']
    exp_D = np.array([7.762320e3, 6.766400e3, 2.780300e3])

    assert_array_equal(exp_ES, obs_ES)
    assert_array_equal(exp_D, obs_D)

    
def test_DoubleSpinDict():
    subject = DoubleSpinDict({(3.499999999998, 2, 1):{'A':'a', 'B':'b'},
                              (2.000000000012, 3, 4):{'C':'c', 'D':'d'}})
    subject.update({(3.500000000011,8,9):{'E':'e', 'F':'f'}})
    
    obs = subject[(3.48, 8, 9)]
    exp = {'E':'e', 'F':'f'}
    assert_equal(exp, obs)


if __name__ == "__main__":
    nose.main()
