#!/usr/bin/env python

"""Fortran reader and writer tests"""
import os
import unittest
import nose
import struct
import warnings
import tables
import nose.tools

from nose.tools import (
    assert_almost_equal,
    assert_equal,
    assert_true,
    assert_not_equal,
    assert_false,
    assert_raises,
)
from nose.plugins.skip import SkipTest
from pyne.utils import QAWarning
from pyne.fortranformat import FortranRecordReader
from pyne.fortranformat import FortranRecordWriter


def read_line(line, ffs):
    """Read a line with the given fortran format specification"""

    ff = FortranRecordReader(ffs)
    words = ff.read(line)
    return words


def write_line(words, ffs):
    """Write a line with the given fortran format specification"""

    ff = FortranRecordWriter(ffs)
    line = ff.write(words)
    return line


##################
#
# the following set of tests, test reading and writing of different formats
#
##################

# test reading simple integers
def test_read_int():

    set_integer = "    2    4"
    set_integer1 = "   60"
    word = read_line(set_integer, "(2I5)")
    word1 = read_line(set_integer1, "(I5)")
    assert_equal(word, [2, 4])
    assert_equal(word1, [60])


# test writing simple integers
def test_write_int():

    set_integer = [2, 4]
    set_integer1 = [60]
    word = write_line(set_integer, "(2I5)")
    word1 = write_line(set_integer1, "(I5)")
    assert_equal(word, "    2    4")
    assert_equal(word1, "   60")


# Tests reading a mix of integers and strings
def test_read_mix_int_string():

    set_int_string = "ntal     2"
    set_int_string1 = "tally    4   -1    0    0"
    set_int_string2 = "   2  tin  he  40"
    word = read_line(set_int_string, "(A4,I6)")
    word1 = read_line(set_int_string1, "(A5,4I5)")
    word2 = read_line(set_int_string2, "(I4,A5,A4,I4)")

    assert_equal(word, ["ntal", 2])
    assert_equal(word1, ["tally", 4, -1, 0, 0])
    assert_equal(word2, [2, "  tin", "  he", 40])


# tests writing integers and strings mixed
def test_write_mix_int_string():
    set_int_string = ["ntal", 2]
    set_int_string1 = ["tally", 4, -1, 0, 0]
    set_int_string2 = [2, "tin", "he", 40]
    word = write_line(set_int_string, "(A4,I6,1X,A5,I6)")
    word1 = write_line(set_int_string1, "(A5,4I5)")
    word2 = write_line(set_int_string2, "(I4,A5,A4,I4)")
    assert_equal(word, "ntal     2")
    assert_equal(word1, "tally    4   -1    0    0")
    assert_equal(word2, "   2  tin  he  40")


# test reading strings
def test_read_string():

    set_string = " Sample Problem Input Deck"
    set_string1 = "Hello World!"
    word = read_line(set_string, "(1x,A30)")
    word1 = read_line(set_string1, "(A13)")
    assert_equal(word, ["Sample Problem Input Deck     "])
    assert_equal(word1, ["Hello World! "])


# tests writing strings
def test_write_string():

    string_set = ["Sample Problem Input Deck"]
    set_string1 = ["Hello World! "]
    word = write_line(string_set, "(1x,A26)")
    word1 = write_line(set_string1, "(A13)")
    assert_equal(word, "  Sample Problem Input Deck")
    assert_equal(word1, "Hello World! ")


# tests reading float numbers
def test_read_floating():

    floating_set = "  2.06784E-05 0.7239  1.80389E-05 0.7299  2.62211E-05 0.7411  0.00000E+00 0.0000"
    word = read_line(floating_set, "(4(1PE13.5,0PF7.4))")
    assert_equal(
        word,
        [
            2.06784e-05,
            0.7239,
            1.80389e-05,
            0.7299,
            2.62211e-05,
            0.7411,
            0.00000e00,
            0.000,
        ],
    )


# tests writing float numbers
def test_write_floating():

    floating_set = [
        2.06784e-05,
        0.7239,
        1.80389e-05,
        0.7299,
        2.62211e-05,
        0.7411,
        0.00000e00,
        0.0000,
    ]
    word = write_line(floating_set, "(4(1PE13.5,0PF7.4))")
    assert_equal(
        word,
        "  2.06784E-05 0.7239  1.80389E-05 0.7299  2.62211E-05 0.7411  0.00000E+00 0.0000",
    )


# tests reading integers and floats mixed
def test_read_mix_int_float():

    int_float_set = "       8000  2.15273E-04  9.99937E-01"
    word = read_line(int_float_set, "(I11,1P2E13.5)")
    assert_equal(word, [8000, 2.15273e-04, 9.99937e-01])


# tests writing integers and floads mixed
def test_write_mix_int_float():

    int_float_set = [8000, 2.15273e-04, 9.99937e-01]
    word = write_line(int_float_set, "(I11,1P3E13.5)")
    assert_equal(word, "       8000  2.15273E-04  9.99937E-01")


# tests a  mctal file
def test_mctal_file():

    filename = "EXAMPLE.INPm"
    f = open(filename, "r")
    configs = [
        (
            "(2A8,A19,I5,2I16)",
            "                                       2          100000          821020",
        ),
        ("(1x,A25)", " Sample Problem Input Deck"),
        ("(A4,I6,1X,A5,I6)", "ntal     2      "),
        ("(16I5)", "    2    4"),
        ("(A5,I5,I21,2I5)", "tally    2                   -1    0    0"),
        (
            "(40I2)",
            " 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
        ),
        ("(A2,I8)", "f        1"),
        ("(11I7)", "      8"),
        ("(A2,I8)", "d        1"),
        ("(A2,I8)", "u        0"),
        ("(A2,I8)", "s        0"),
        ("(A2,I8)", "m        0"),
        ("(A2,I8,I4)", "c        0"),
        ("(A2,I8,I4)", "et      15"),
        (
            "(1P6E13.5)",
            "  1.00000E+00  2.00000E+00  3.00000E+00  4.00000E+00  5.00000E+00  6.00000E+00",
        ),
        (
            "(1P6E13.5)",
            "  7.00000E+00  8.00000E+00  9.00000E+00  1.00000E+01  1.10000E+01  1.20000E+01",
        ),
        ("(1P6E13.5)", "  1.30000E+01  1.40000E+01"),
        ("(A2,I8,I4)", "t        0"),
        ("(A4)", "vals"),
        (
            "(4(1PE13.5,0PF7.4))",
            "  2.06784E-05 0.7239  1.80389E-05 0.7299  2.62211E-05 0.7411  0.00000E+00 0.0000",
        ),
        (
            "(4(1PE13.5,0PF7.4))",
            "  0.00000E+00 0.0000  0.00000E+00 0.0000  6.15160E-06 1.0000  6.96320E-06 1.0000",
        ),
        (
            "(4(1PE13.5,0PF7.4))",
            "  5.85347E-06 1.0000  0.00000E+00 0.0000  0.00000E+00 0.0000  0.00000E+00 0.0000",
        ),
        (
            "(4(1PE13.5,0PF7.4))",
            "  0.00000E+00 0.0000  5.17611E-04 0.3562  6.01518E-04 0.3105",
        ),
        (
            "(A3,I5,8I8)",
            "tfc   13       1       1       1       1       1       1      15       1",
        ),
        ("(I11,1P3E13.5)", "       8000  2.15273E-04  9.99937E-01"),
        ("(I11,1P3E13.5)", "      16000  2.72525E-04  5.32249E-01"),
        ("(I11,1P3E13.5)", "      24000  4.37049E-04  3.75945E-01"),
        ("(I11,1P3E13.5)", "      32000  3.95198E-04  3.34352E-01"),
        ("(I11,1P3E13.5)", "      40000  3.63456E-04  3.06290E-01"),
        ("(I11,1P3E13.5)", "      48000  6.46706E-04  4.38910E-01"),
        ("(I11,1P3E13.5)", "      56000  5.87226E-04  4.16305E-01"),
        ("(I11,1P3E13.5)", "      64000  5.39816E-04  3.97730E-01"),
        ("(I11,1P3E13.5)", "      72000  7.10568E-04  3.61496E-01"),
        ("(I11,1P3E13.5)", "      80000  6.47902E-04  3.57051E-01"),
        ("(I11,1P3E13.5)", "      88000  6.15867E-04  3.42404E-01"),
        ("(I11,1P3E13.5)", "      96000  5.89509E-04  3.28822E-01"),
        ("(I11,1P3E13.5)", "     100000  6.01518E-04  3.10516E-01"),
    ]

    for fs, fs1 in configs:
        line = f.readline()
        word = read_line(line, fs)
        word1 = [item for item in word if item is not None]
        l = write_line(word1, fs)
    assert_equal(l, fs1)
