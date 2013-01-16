import numpy as np
from pyne.endf import Library
import nose 
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, assert_in
import warnings


library = Library('endftest_small.txt')   
cs137 = Library('137.txt') 

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
    
    obs = library.debug
    print obs
    assert (True == False)

def test_unresolved_resonances_c():
    # Case C (ENDF Manual p. 70)
    print library.structure
    obs = library.mat128['RxData']['Unresolved'][0]['data'][3.5][1.0][4.0]
    obs_ES = obs['ES']
    exp_ES = np.array([1.7e3, 2.0e3, 3.0e3])
    obs_D = obs['D']
    exp_D = np.array([7.762320e3, 6.766400e3, 2.780300e3])

    assert(np.array_equal(exp_ES, obs_ES))
    assert(np.array_equal(exp_D, obs_D))

    

# def test_unresolved_resonances():
#     resonances = library.make_resonances()
    

def test_flags():
    print library.mat128['flags']
    assert(True == True)
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
