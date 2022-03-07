#!/usr/bin/env python
from __future__ import unicode_literals
import os

from nose.tools import assert_equal, assert_in, assert_almost_equal

import pyne.ace


def setup():
    try:
        import urllib.request as urllib
    except ImportError:
        import urllib

    if not os.path.isfile("C012-n.ace"):
        urllib.urlretrieve(
            "https://raw.githubusercontent.com/pyne/data/master/C012-n.ace",
            "C012-n.ace",
        )
    with open("C012-n.ace") as f:
        lines = f.readlines()

    lines.insert(0, "The next lines contain the old-style header\n")
    lines.insert(0, "11.896900 2.5263E-08 12/06/13     3\n")
    lines.insert(0, "2.0.0" + 5 * " " + " 6000.000nc" + 13 * " " + "TENDL\n")

    with open("C012-n-2p0.ace", "w") as f:
        f.writelines(lines)


def test_convert_c12():
    pyne.ace.ascii_to_binary("C012-n.ace", "C12-binary.ace")


def test_read_c12_ascii():
    c12 = pyne.ace.Library("C012-n.ace")
    c12.read()

    assert_in("6000.00c", c12.tables)
    table = c12.tables["6000.00c"]

    assert_equal(table.nxs[1], 38937)
    assert_equal(table.nxs[2], 6000)
    assert_equal(table.nxs[3], 1513)
    assert_equal(table.jxs[1], 1)

    assert_in(2, table.reactions)
    assert_in(107, table.reactions)
    assert_in(204, table.reactions)
    assert_in(444, table.reactions)

    assert_almost_equal(table.energy[0], 1.0e-11)

    assert_equal(table.reactions[2].sigma[0], 78.04874)
    assert_equal(table.reactions[2].sigma[-1], 1.00772)


def test_read_c12_2p0_ascii():
    c12 = pyne.ace.Library("C012-n-2p0.ace")
    c12.read()

    assert_in("6000.000nc", c12.tables)
    table = c12.tables["6000.000nc"]

    assert_equal(table.nxs[1], 38937)
    assert_equal(table.nxs[2], 6000)
    assert_equal(table.nxs[3], 1513)
    assert_equal(table.jxs[1], 1)

    assert_in(2, table.reactions)
    assert_in(107, table.reactions)
    assert_in(204, table.reactions)
    assert_in(444, table.reactions)

    assert_almost_equal(table.energy[0], 1.0e-11)

    assert_equal(table.reactions[2].sigma[0], 78.04874)
    assert_equal(table.reactions[2].sigma[-1], 1.00772)


def test_read_c12_binary():
    c12 = pyne.ace.Library("C12-binary.ace")
    c12.read()

    assert_in("6000.00c", c12.tables)
    table = c12.tables["6000.00c"]

    assert_equal(table.nxs[1], 38937)
    assert_equal(table.nxs[2], 6000)
    assert_equal(table.nxs[3], 1513)
    assert_equal(table.jxs[1], 1)

    assert_in(2, table.reactions)
    assert_in(107, table.reactions)
    assert_in(204, table.reactions)
    assert_in(444, table.reactions)

    assert_almost_equal(table.energy[0], 1.0e-11)

    assert_equal(table.reactions[2].sigma[0], 78.04874)
    assert_equal(table.reactions[2].sigma[-1], 1.00772)


def teardown():
    if os.path.exists("C12-binary.ace"):
        os.remove("C12-binary.ace")
