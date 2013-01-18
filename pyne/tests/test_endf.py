import numpy as np
from pyne.endf import Library
import nose 
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, assert_in
import warnings


library = Library('endftest_small.txt')   

def test_mats():
    for mat_id in library.mats:
        assert_in(mat_id, [128, 131, 419])

def test_get():
    obs = library.get(419, 4, 2)
    # print result
    # print library.mat128
    exp = np.array([4.898421e+3, 6.768123e+0, 0, 
                    1, 0, 0, 2.123124e+6, 8.123142e-6, 
                    2.123212e+6, 8.231231e-6,
                    -2.231211e+6, 8.123421e-6])    
    badkey = library.get(111, 1, 1)
    assert(np.array_equal(exp, obs))
    assert_equal(badkey, False)

def test_unresolved_resonances_a():
    # Case A (ENDF Manual p.70)
    
    obs = library.mat131['RxData']['Unresolved'][0][2][3.5][1.0]

    obs_D = obs['D']
    exp_D = np.array([1.810000e3,
                      2.110000e3,
                      3.110000e3])

    obs_GNO = obs['GNO']
    exp_GNO = np.array([4.489400e-1,
                        8.497500e-1,
                        9.524900e-1])

    assert(np.array_equal(exp_D, obs_D))
    assert(np.array_equal(exp_GNO, obs_GNO))

def test_unresolved_resonances_b():
    # Case B (ENDF Manual p. 70)
    obs = library.mat419['RxData']['Unresolved'][-1][2][3.5]
    obs_ES = obs['ES']
    exp_ES = 100 * np.array([10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])

    assert(np.array_equal(exp_ES, obs_ES))

def test_unresolved_resonances_c():
    # Case C (ENDF Manual p. 70)
    obs = library.mat128['RxData']['Unresolved'][0][2][3.5][1.0][4.0]

    obs_ES = obs['ES']
    exp_ES = np.array([1.74e3, 2.04e3, 3.04e3])

    obs_D = obs['D']
    exp_D = np.array([7.762320e3, 6.766400e3, 2.780300e3])

    assert(np.array_equal(exp_ES, obs_ES))
    assert(np.array_equal(exp_D, obs_D))

    

# def do_not_test_write():
    # write library to test.txt
    # library.write('test.txt', 'endf')
    # result = library.write('test.txt', 'txt')
    # written = Library('test')
    # result = written.get(419, 4, 2)
    # expected = library.get(419, 4, 2)
    # expected = 'test.txt'
    # assert_equal(result, expected)

if __name__ == "__main__":
    nose.main()
