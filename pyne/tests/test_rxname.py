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


def test_name_names():
    assert_equal(rxname.name("a"), "a")
    assert_equal(rxname.name("total"), "total")

def test_name_alts():
    assert_equal(rxname.name("alpha"), "a")
    assert_equal(rxname.name("tot"), "total")

def test_name_ids():
    assert_equal(rxname.name(_hash("a")), "a")
    assert_equal(rxname.name(_hash("total")), "total")    

    assert_equal(rxname.name(long(_hash("a"))), "a")
    assert_equal(rxname.name(long(_hash("total"))), "total")    

    assert_equal(rxname.name(str(_hash("a"))), "a")
    assert_equal(rxname.name(str(_hash("total"))), "total")    

def test_name_mts():
    assert_equal(rxname.name(107), "a")
    assert_equal(rxname.name(1), "total")

    assert_equal(rxname.name(107L), "a")
    assert_equal(rxname.name(1L), "total")

    assert_equal(rxname.name("107"), "a")
    assert_equal(rxname.name("1"), "total")

def test_name_nucdelta():
    assert_equal(rxname.name("U235", "U236"), "absorption")
    assert_equal(rxname.name("U235", "Np236", "p"), "absorption")
    assert_equal(rxname.name(922350, 912350), "p")

def test_name_not():
    assert_raises(RuntimeError, rxname.name, "Waka waka")
    assert_raises(RuntimeError, rxname.name, 0)


if __name__ == "__main__":
    nose.main()

