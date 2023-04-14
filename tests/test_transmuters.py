"""Transmuter tests"""
import pytest
from pyne import data
from pyne import cram
from pyne import nucname
from pyne import transmuters


def test_transmuters_cram():
    n0 = {"H3": 1.0}
    A = -cram.DECAY_MATRIX * data.half_life("H3")
    n1 = transmuters.cram(A, n0, order=16)
    assert 2 == len(n1)
    assert 0.5 == pytest.approx(n1[nucname.id("H3")])
    assert 0.5 == pytest.approx(n1[nucname.id("He3")])

