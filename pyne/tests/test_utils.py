"""PyNE utility tests"""
import os

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, assert_in

from pyne import utils
import numpy as np


def test_to_sec():
    assert_equal(120.0, utils.to_sec(2, 'M'))


def test_to_barns():
    assert_equal(3E3, utils.to_barns(3, 'KB'))


def test_message():
    if os.name is 'posix':
        assert_equal('\033[1;32mHello\033[0m', utils.message('Hello'))
    else:
        assert_equal('*** MESSAGE ***: Hello', utils.message('Hello'))


def test_failure():
    if os.name is 'posix':
        assert_equal('\033[1;31mWorld\033[0m', utils.failure('World'))
    else:
        assert_equal('*** FAILURE ***: World', utils.failure('World'))


if __name__ == "__main__":
    nose.main()

