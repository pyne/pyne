"""rxname tests"""
from __future__ import unicode_literals, division
import sys
import warnings

import nose
from nose.tools import (
    assert_equal,
    assert_not_equal,
    assert_raises,
    raises,
    assert_in,
    assert_greater_equal,
)

from pyne.utils import QAWarning

warnings.simplefilter("ignore", QAWarning)
from pyne import rxname

if sys.version_info[0] > 2:
    long = int


def _hash(s):
    MVAL = 2**32  # integer size
    h = 32
    for c in s:
        c = ord(c)
        h = ((h << 5) + h) ^ c
    h = h % MVAL
    return h


def test_hash():
    rxs = [
        "a",
        "hello",
        "total",
        "wowza",
        "z_2na",
        "z_3na",
        "absorption",
        "np",
        "n2a",
        "z_2n2a",
        "nd",
        "nt",
        "nHe3",
        "nd3a",
        "nt2a",
        "z_4n",
        "fission_fourth",
        "z_2np",
        "z_3np",
        "n2p",
        "npa",
        "n_0",
        "n_1",
        "n_2",
        "n_3",
        "n_4",
        "n_5",
        "n_6",
        "n_7",
        "n_8",
    ]
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

    assert_equal(rxname.name(long(107)), "a")
    assert_equal(rxname.name(long(1)), "total")

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

    assert_equal(rxname.id(long(107)), _hash("a"))
    assert_equal(rxname.id(long(1)), _hash("total"))

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

    assert_equal(rxname.mt(long(107)), 107)
    assert_equal(rxname.mt(long(1)), 1)

    assert_equal(rxname.mt("107"), 107)
    assert_equal(rxname.mt("1"), 1)


def test_mt_nucdelta():
    assert_equal(rxname.mt("U235", "U236"), 27)
    assert_equal(rxname.mt("U235", "Np236", "p"), 27)
    assert_equal(rxname.mt(922350, 912350), 103)


def test_mt_not():
    assert_raises(RuntimeError, rxname.mt, "Waka waka")
    assert_raises(RuntimeError, rxname.mt, 0)


def test_child():
    assert_equal(rxname.child("U235", "absorption"), 922360000)
    assert_equal(rxname.child(922350000, "absorption"), 922360000)
    assert_equal(rxname.child("Co58M", "gamma"), 270590000)
    assert_equal(rxname.child("Co58M", "gamma_1"), 270590001)
    assert_equal(rxname.child("Co58M", "gamma_1", "n"), 270590001)
    assert_equal(rxname.child("Co58M", "gamma_1", b"n"), 270590001)


def test_parent():
    assert_equal(rxname.parent("U235", "absorption"), 922340000)
    assert_equal(rxname.parent(922350000, "absorption"), 922340000)
    assert_equal(rxname.parent(922350000, "absorption", "n"), 922340000)
    assert_equal(rxname.parent(922350000, "absorption", b"n"), 922340000)


alabel = "(z,a)"
plabel = "(z,p)"
abslabel = "(z,abs) Absorption"
totlabel = "(z,total)"


def test_label_names():
    assert_equal(rxname.label("a"), alabel)
    assert_equal(rxname.label("total"), totlabel)


def test_label_alts():
    assert_equal(rxname.label("alpha"), alabel)
    assert_equal(rxname.label("tot"), totlabel)


def test_label_ids():
    assert_equal(rxname.label(_hash("a")), alabel)
    assert_equal(rxname.label(_hash("total")), totlabel)

    assert_equal(rxname.label(long(_hash("a"))), alabel)
    assert_equal(rxname.label(long(_hash("total"))), totlabel)

    assert_equal(rxname.label(str(_hash("a"))), alabel)
    assert_equal(rxname.label(str(_hash("total"))), totlabel)


def test_label_mts():
    assert_equal(rxname.label(107), alabel)
    assert_equal(rxname.label(1), totlabel)

    assert_equal(rxname.label(long(107)), alabel)
    assert_equal(rxname.label(long(1)), totlabel)

    assert_equal(rxname.label("107"), alabel)
    assert_equal(rxname.label("1"), totlabel)


def test_label_nucdelta():
    assert_equal(rxname.label("U235", "U236"), abslabel)
    assert_equal(rxname.label("U235", "Np236", "p"), abslabel)
    assert_equal(rxname.label(922350, 912350), plabel)


def test_label_not():
    assert_raises(RuntimeError, rxname.label, "Waka waka")
    assert_raises(RuntimeError, rxname.label, 0)


adoc = "(z,a) Production of alpha"
pdoc = "(z,p) Production of p"
absdoc = "(n,abs) Absorption"
totdoc = "(n,total) Neutron total"


def test_doc_names():
    assert_equal(rxname.doc("a"), adoc)
    assert_equal(rxname.doc("total"), totdoc)


def test_doc_alts():
    assert_equal(rxname.doc("alpha"), adoc)
    assert_equal(rxname.doc("tot"), totdoc)


def test_doc_ids():
    assert_equal(rxname.doc(_hash("a")), adoc)
    assert_equal(rxname.doc(_hash("total")), totdoc)

    assert_equal(rxname.doc(long(_hash("a"))), adoc)
    assert_equal(rxname.doc(long(_hash("total"))), totdoc)

    assert_equal(rxname.doc(str(_hash("a"))), adoc)
    assert_equal(rxname.doc(str(_hash("total"))), totdoc)


def test_doc_mts():
    assert_equal(rxname.doc(107), adoc)
    assert_equal(rxname.doc(1), totdoc)

    assert_equal(rxname.doc(long(107)), adoc)
    assert_equal(rxname.doc(long(1)), totdoc)

    assert_equal(rxname.doc("107"), adoc)
    assert_equal(rxname.doc("1"), totdoc)


def test_doc_nucdelta():
    assert_equal(rxname.doc("U235", "U236"), absdoc)
    assert_equal(rxname.doc("U235", "Np236", "p"), absdoc)
    assert_equal(rxname.doc(922350, 912350), pdoc)


def test_doc_not():
    assert_raises(RuntimeError, rxname.doc, "Waka waka")
    assert_raises(RuntimeError, rxname.doc, 0)


def test_unique_ids():
    assert_equal(len(rxname.id_name), len(rxname.name_id))


def test_no_id_mt_clash():
    badids = [rxid for rxid in rxname.id_name if rxid < 1000]
    assert_equal(0, len(badids))


if __name__ == "__main__":
    nose.runmodule()
