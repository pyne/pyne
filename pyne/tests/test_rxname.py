"""rxname tests"""
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_in

from pyne import rxname

def _hash(s):
    MVAL = 2**32  # integer size
    h = 32
    for c in s:
        c = ord(c)
        h = ((h << 5) + h) ^ c
    h = h%MVAL
    return h

def test_hash():
    rxs = ["a", "hello", "total", "wowza", "z_2na", "z_3na", "absorption", "np", 
           "n2a", "z_2n2a", "nd", "nt", "nHe3", "nd3a", "nt2a", "z_4n", 
           "fission_fourth", "z_2np", "z_3np", "n2p", "npa", "n_0", "n_1", "n_2", 
           "n_3", "n_4", "n_5", "n_6", "n_7", "n_8",]
    for rx in rxs:
        yield assert_equal, rxname.hash(rx), _hash(rx)



if __name__ == "__main__":
    nose.main()

