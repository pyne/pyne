#!/usr/bin/env python

import os

from nose.tools import assert_equal, assert_in

import pyne.ace

def test_convert_h1():
    pyne.ace.ascii_to_binary('H_001_293.6K.ace', 'H1-binary.ace')

def test_read_h1_ascii():
    h1 = pyne.ace.Library('H_001_293.6K.ace')
    h1.read()

    assert_in('1001.71c', h1.tables)
    table = h1.tables['1001.71c']

    assert_equal(table.nxs[1], 8177)
    assert_equal(table.nxs[2], 1001)
    assert_equal(table.nxs[3], 590)
    assert_equal(table.jxs[1], 1)

    assert_in(2, table.reactions)
    assert_in(102, table.reactions)
    assert_in(204, table.reactions)
    assert_in(444, table.reactions)

    assert_equal(table.energy[0], 1.0e-11)

    assert_equal(table.reactions[2].sigma[0], 1160.546)
    assert_equal(table.reactions[2].sigma[-1], 0.4827462)

def test_read_h1_binary():
    h1 = pyne.ace.Library('H1-binary.ace')
    h1.read()

    assert_in('1001.71c', h1.tables)
    table = h1.tables['1001.71c']

    assert_equal(table.nxs[1], 8177)
    assert_equal(table.nxs[2], 1001)
    assert_equal(table.nxs[3], 590)
    assert_equal(table.jxs[1], 1)

    assert_in(2, table.reactions)
    assert_in(102, table.reactions)
    assert_in(204, table.reactions)
    assert_in(444, table.reactions)

    assert_equal(table.energy[0], 1.0e-11)

    assert_equal(table.reactions[2].sigma[0], 1160.546)
    assert_equal(table.reactions[2].sigma[-1], 0.4827462)

def teardown():
    if os.path.exists('H1-binary.ace'):
        os.remove('H1-binary.ace')
