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

def test_id_names():
    assert_equal(rxname.id("a"), _hash("a"))
    assert_equal(rxname.id("total"), _hash("total"))

def test_id_alts():
    assert_equal(rxname.id("alpha"), _hash("a"))
    assert_equal(rxname.id("tot"), _hash("total"))

def test_id_ids():
    assert_equal(rxname.id(_hash("a")), _hash("a"))
    assert_equal(rxname.id(_hash("total")), _hash("total"))

    assert_equal(rxname.id(long(_hash("a"))), _hash("a"))
    assert_equal(rxname.id(long(_hash("total"))), _hash("total"))    

    assert_equal(rxname.id(str(_hash("a"))), _hash("a"))
    assert_equal(rxname.id(str(_hash("total"))), _hash("total"))    

def test_id_mts():
    assert_equal(rxname.id(107), _hash("a"))
    assert_equal(rxname.id(1), _hash("total"))

    assert_equal(rxname.id(107L), _hash("a"))
    assert_equal(rxname.id(1L), _hash("total"))

    assert_equal(rxname.id("107"), _hash("a"))
    assert_equal(rxname.id("1"), _hash("total"))

def test_id_nucdelta():
    assert_equal(rxname.id("U235", "U236"), _hash("absorption"))
    assert_equal(rxname.id("U235", "Np236", "p"), _hash("absorption"))
    assert_equal(rxname.id(922350, 912350), _hash("p"))

def test_id_not():
    assert_raises(RuntimeError, rxname.id, "Waka waka")
    assert_raises(RuntimeError, rxname.id, 0)

def test_mt_names():
    assert_equal(rxname.mt("a"), 107)
    assert_equal(rxname.mt("total"), 1)

def test_mt_alts():
    assert_equal(rxname.mt("alpha"), 107)
    assert_equal(rxname.mt("tot"), 1)

def test_mt_ids():
    assert_equal(rxname.mt(_hash("a")), 107)
    assert_equal(rxname.mt(_hash("total")), 1)

    assert_equal(rxname.mt(long(_hash("a"))), 107)
    assert_equal(rxname.mt(long(_hash("total"))), 1)    

    assert_equal(rxname.mt(str(_hash("a"))), 107)
    assert_equal(rxname.mt(str(_hash("total"))), 1)    

def test_mt_mts():
    assert_equal(rxname.mt(107), 107)
    assert_equal(rxname.mt(1), 1)

    assert_equal(rxname.mt(107L), 107)
    assert_equal(rxname.mt(1L), 1)

    assert_equal(rxname.mt("107"), 107)
    assert_equal(rxname.mt("1"), 1)

def test_mt_nucdelta():
    assert_equal(rxname.mt("U235", "U236"), 27)
    assert_equal(rxname.mt("U235", "Np236", "p"), 27)
    assert_equal(rxname.mt(922350, 912350), 103)

def test_mt_not():
    assert_raises(RuntimeError, rxname.mt, "Waka waka")
    assert_raises(RuntimeError, rxname.mt, 0)



if __name__ == "__main__":
    nose.main()

