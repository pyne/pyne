"""rxname tests"""
from __future__ import unicode_literals, division
import sys
import warnings

import pytest

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


@pytest.mark.parametrize("rx",[
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
    ])
def test_hash(rx):
    assert rxname.hash(rx) == _hash(rx)


def test_name_names():
    assert rxname.name("a") == "a"
    assert rxname.name("total") == "total"


def test_name_alts():
    assert rxname.name("alpha") == "a"
    assert rxname.name("tot") == "total"


def test_name_ids():
    assert rxname.name(_hash("a")) == "a"
    assert rxname.name(_hash("total")) == "total"

    assert rxname.name(long(_hash("a"))) == "a"
    assert rxname.name(long(_hash("total"))) == "total"

    assert rxname.name(str(_hash("a"))) == "a"
    assert rxname.name(str(_hash("total"))) == "total"


def test_name_mts():
    assert rxname.name(107) == "a"
    assert rxname.name(1) == "total"

    assert rxname.name(long(107)) == "a"
    assert rxname.name(long(1)) == "total"

    assert rxname.name("107") == "a"
    assert rxname.name("1") == "total"


def test_name_nucdelta():
    assert rxname.name("U235", "U236") == "absorption"
    assert rxname.name("U235", "Np236", "p") == "absorption"
    assert rxname.name(922350, 912350) == "p"


def test_name_not():
    pytest.raises(RuntimeError, rxname.name, "Waka waka")
    pytest.raises(RuntimeError, rxname.name, 0)


def test_id_names():
    assert rxname.id("a") == _hash("a")
    assert rxname.id("total") == _hash("total")


def test_id_alts():
    assert rxname.id("alpha") == _hash("a")
    assert rxname.id("tot") == _hash("total")


def test_id_ids():
    assert rxname.id(_hash("a")) == _hash("a")
    assert rxname.id(_hash("total")) == _hash("total")

    assert rxname.id(long(_hash("a"))) == _hash("a")
    assert rxname.id(long(_hash("total"))) == _hash("total")

    assert rxname.id(str(_hash("a"))) == _hash("a")
    assert rxname.id(str(_hash("total"))) == _hash("total")


def test_id_mts():
    assert rxname.id(107) == _hash("a")
    assert rxname.id(1) == _hash("total")

    assert rxname.id(long(107)) == _hash("a")
    assert rxname.id(long(1)) == _hash("total")

    assert rxname.id("107") == _hash("a")
    assert rxname.id("1") == _hash("total")


def test_id_nucdelta():
    assert rxname.id("U235", "U236") == _hash("absorption")
    assert rxname.id("U235", "Np236", "p") == _hash("absorption")
    assert rxname.id(922350, 912350) == _hash("p")


def test_id_not():
    pytest.raises(RuntimeError, rxname.id, "Waka waka")
    pytest.raises(RuntimeError, rxname.id, 0)


def test_mt_names():
    assert rxname.mt("a") == 107
    assert rxname.mt("total") == 1


def test_mt_alts():
    assert rxname.mt("alpha") == 107
    assert rxname.mt("tot") == 1


def test_mt_ids():
    assert rxname.mt(_hash("a")) == 107
    assert rxname.mt(_hash("total")) == 1

    assert rxname.mt(long(_hash("a"))) == 107
    assert rxname.mt(long(_hash("total"))) == 1

    assert rxname.mt(str(_hash("a"))) == 107
    assert rxname.mt(str(_hash("total"))) == 1


def test_mt_mts():
    assert rxname.mt(107) == 107
    assert rxname.mt(1) == 1

    assert rxname.mt(long(107)) == 107
    assert rxname.mt(long(1)) == 1

    assert rxname.mt("107") == 107
    assert rxname.mt("1") == 1


def test_mt_nucdelta():
    assert rxname.mt("U235", "U236") == 27
    assert rxname.mt("U235", "Np236", "p") == 27
    assert rxname.mt(922350, 912350) == 103


def test_mt_not():
    pytest.raises(RuntimeError, rxname.mt, "Waka waka")
    pytest.raises(RuntimeError, rxname.mt, 0)


def test_child():
    assert rxname.child("U235", "absorption") == 922360000
    assert rxname.child(922350000, "absorption") == 922360000
    assert rxname.child("Co58M", "gamma") == 270590000
    assert rxname.child("Co58M", "gamma_1") == 270590001
    assert rxname.child("Co58M", "gamma_1", "n") == 270590001
    assert rxname.child("Co58M", "gamma_1", b"n") == 270590001


def test_parent():
    assert rxname.parent("U235", "absorption") == 922340000
    assert rxname.parent(922350000, "absorption") == 922340000
    assert rxname.parent(922350000, "absorption", "n") == 922340000
    assert rxname.parent(922350000, "absorption", b"n") == 922340000


alabel = "(z,a)"
plabel = "(z,p)"
abslabel = "(z,abs) Absorption"
totlabel = "(z,total)"


def test_label_names():
    assert rxname.label("a") == alabel
    assert rxname.label("total") == totlabel


def test_label_alts():
    assert rxname.label("alpha") == alabel
    assert rxname.label("tot") == totlabel


def test_label_ids():
    assert rxname.label(_hash("a")) == alabel
    assert rxname.label(_hash("total")) == totlabel

    assert rxname.label(long(_hash("a"))) == alabel
    assert rxname.label(long(_hash("total"))) == totlabel

    assert rxname.label(str(_hash("a"))) == alabel
    assert rxname.label(str(_hash("total"))) == totlabel


def test_label_mts():
    assert rxname.label(107) == alabel
    assert rxname.label(1) == totlabel

    assert rxname.label(long(107)) == alabel
    assert rxname.label(long(1)) == totlabel

    assert rxname.label("107") == alabel
    assert rxname.label("1") == totlabel


def test_label_nucdelta():
    assert rxname.label("U235", "U236") == abslabel
    assert rxname.label("U235", "Np236", "p") == abslabel
    assert rxname.label(922350, 912350) == plabel


def test_label_not():
    pytest.raises(RuntimeError, rxname.label, "Waka waka")
    pytest.raises(RuntimeError, rxname.label, 0)


adoc = "(z,a) Production of alpha"
pdoc = "(z,p) Production of p"
absdoc = "(n,abs) Absorption"
totdoc = "(n,total) Neutron total"


def test_doc_names():
    assert rxname.doc("a") == adoc
    assert rxname.doc("total") == totdoc


def test_doc_alts():
    assert rxname.doc("alpha") == adoc
    assert rxname.doc("tot") == totdoc


def test_doc_ids():
    assert rxname.doc(_hash("a")) == adoc
    assert rxname.doc(_hash("total")) == totdoc

    assert rxname.doc(long(_hash("a"))) == adoc
    assert rxname.doc(long(_hash("total"))) == totdoc

    assert rxname.doc(str(_hash("a"))) == adoc
    assert rxname.doc(str(_hash("total"))) == totdoc


def test_doc_mts():
    assert rxname.doc(107) == adoc
    assert rxname.doc(1) == totdoc

    assert rxname.doc(long(107)) == adoc
    assert rxname.doc(long(1)) == totdoc

    assert rxname.doc("107") == adoc
    assert rxname.doc("1") == totdoc


def test_doc_nucdelta():
    assert rxname.doc("U235", "U236") == absdoc
    assert rxname.doc("U235", "Np236", "p") == absdoc
    assert rxname.doc(922350, 912350) == pdoc


def test_doc_not():
    pytest.raises(RuntimeError, rxname.doc, "Waka waka")
    pytest.raises(RuntimeError, rxname.doc, 0)


def test_unique_ids():
    assert len(rxname.id_name) == len(rxname.name_id)


def test_no_id_mt_clash():
    badids = [rxid for rxid in rxname.id_name if rxid < 1000]
    assert 0 == len(badids)

