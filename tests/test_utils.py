"""PyNE utility tests"""
import os

import nose

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_in, assert_true

from numpy.testing import assert_array_equal, assert_array_almost_equal

from pyne import utils
import numpy as np


def test_to_sec():
    # as
    assert_equal(2e-18, utils.to_sec(2, 'AS'))
    assert_equal(2e-18, utils.to_sec(2, 'attosec'))
    assert_equal(2e-18, utils.to_sec(2, 'attosecond'))
    assert_equal(2e-18, utils.to_sec(2, 'attoseconds'))
    # fs
    assert_equal(2e-15, utils.to_sec(2, 'FS'))
    assert_equal(2e-15, utils.to_sec(2, 'femtosec'))
    assert_equal(2e-15, utils.to_sec(2, 'femtosecond'))
    assert_equal(2e-15, utils.to_sec(2, 'femtoseconds'))
    # ps
    assert_equal(2e-12, utils.to_sec(2, 'PS'))
    assert_equal(2e-12, utils.to_sec(2, 'picosec'))
    assert_equal(2e-12, utils.to_sec(2, 'picosecond'))
    assert_equal(2e-12, utils.to_sec(2, 'picoseconds'))
    # ns
    assert_equal(2e-9, utils.to_sec(2, 'NS'))
    assert_equal(2e-9, utils.to_sec(2, 'nanosec'))
    assert_equal(2e-9, utils.to_sec(2, 'nanosecond'))
    assert_equal(2e-9, utils.to_sec(2, 'nanoseconds'))
    # us
    assert_equal(2e-6, utils.to_sec(2, 'US'))
    assert_equal(2e-6, utils.to_sec(2, 'microsec'))
    assert_equal(2e-6, utils.to_sec(2, 'microsecond'))
    assert_equal(2e-6, utils.to_sec(2, 'microseconds'))
    # ms
    assert_equal(2e-3, utils.to_sec(2, 'MS'))
    assert_equal(2e-3, utils.to_sec(2, 'millisec'))
    assert_equal(2e-3, utils.to_sec(2, 'millisecond'))
    assert_equal(2e-3, utils.to_sec(2, 'milliseconds'))
    # s
    assert_equal(2.0, utils.to_sec(2, 'S'))
    assert_equal(2.0, utils.to_sec(2, 'sec'))
    assert_equal(2.0, utils.to_sec(2, 'second'))
    assert_equal(2.0, utils.to_sec(2, 'seconds'))
    # m
    assert_equal(120.0, utils.to_sec(2, 'M'))
    assert_equal(120.0, utils.to_sec(2, 'min'))
    assert_equal(120.0, utils.to_sec(2, 'minute'))
    assert_equal(120.0, utils.to_sec(2, 'minutes'))
    # h
    assert_equal(7200.0, utils.to_sec(2, 'H'))
    assert_equal(7200.0, utils.to_sec(2, 'hour'))
    assert_equal(7200.0, utils.to_sec(2, 'hours'))
    # d
    assert_equal(172800.0, utils.to_sec(2, 'D'))
    assert_equal(172800.0, utils.to_sec(2, 'day'))
    assert_equal(172800.0, utils.to_sec(2, 'days'))
    # w
    assert_equal(1209600.0, utils.to_sec(2, 'W'))
    assert_equal(1209600.0, utils.to_sec(2, 'week'))
    assert_equal(1209600.0, utils.to_sec(2, 'weeks'))
    # y
    assert_equal(63115200.0, utils.to_sec(2, 'Y'))
    assert_equal(63115200.0, utils.to_sec(2, 'year'))
    assert_equal(63115200.0, utils.to_sec(2, 'years'))
    # c
    assert_equal(6311520000.0, utils.to_sec(2, 'C'))
    assert_equal(6311520000.0, utils.to_sec(2, 'century'))
    assert_equal(6311520000.0, utils.to_sec(2, 'centuries'))
    # undifined unit trigs ValueError
    assert_raises(ValueError, utils.to_sec, 2, 'month')


def test_to_barns():
    assert_equal(3E3, utils.to_barns(3, 'KB'))


def test_from_barns():
    assert_equal(3E-3, utils.from_barns(3, 'KB'))


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


def test_remove():
    # test file removal
    assert_true('hello.txt' not in os.listdir('.'))
    with open('hello.txt', 'w') as f:
        f.write("Hello PyNE!")
    assert_true('hello.txt' in os.listdir('.'))

    utils.remove('hello.txt')
    assert_true('hello.txt' not in os.listdir('.'))

    # test recursive dir removal
    assert_true('rawr' not in os.listdir('.'))
    os.mkdir('rawr')

    assert_true('hello.txt' not in os.listdir('rawr'))
    with open('rawr/hello.txt', 'w') as f:
        f.write("Hello PyNE!")
    assert_true('hello.txt' in os.listdir('rawr'))

    utils.remove('rawr')
    assert_true('rawr' not in os.listdir('.'))

    # Test pass through on non-existant
    assert_true('Medic!' not in os.listdir('.'))
    utils.remove('Medic!')


def test_fromstring_split1():
    s = "1 2.3 4"
    obs = utils.fromstring_split(s)
    exp = np.fromstring(s, sep=" ")
    assert_array_equal(obs, exp)


def test_fromstring_split2():
    s = "1-2.3-4"
    obs = utils.fromstring_split(s, sep='-')
    exp = np.fromstring(s, sep="-")
    assert_array_equal(obs, exp)


def test_fromstring_split3():
    s = "1,2.3,4"
    obs = utils.fromstring_split(s, sep=',')
    exp = np.fromstring(s, sep=",")
    assert_array_equal(obs, exp)


def test_fromstring_split4():
    s = "1\n 23 \t 4"
    obs = utils.fromstring_split(s, dtype=int)
    exp = np.array([1, 23, 4])
    assert_array_equal(obs, exp)


def test_fromstring_token1():
    s = "1 2.3 4"
    obs = utils.fromstring_token(s)
    exp = np.fromstring(s, sep=" ")
    assert_array_equal(obs, exp)


def test_fromstring_token2():
    s = "1-2.3-4"
    obs = utils.fromstring_token(s, sep='-')
    exp = np.fromstring(s, sep="-")
    assert_array_equal(obs, exp)


def test_fromstring_token3():
    s = "1,  2.3, 4"
    obs = utils.fromstring_token(s, sep=' ,')
    exp = np.fromstring(s, sep=", ")
    assert_array_equal(obs, exp)


def test_fromstring_token4():
    s = "1, 2.3 ,4"
    obs = utils.fromstring_token(s, sep=' ,', inplace=True)
    exp = np.array([1.0, 2.3, 4.0])
    assert_array_equal(obs, exp)


def test_use_warnings():
    first = utils.use_warnings()
    second = utils.use_warnings()
    assert_equal(first, second)


def test_toggle_warnings():
    state = utils.use_warnings()
    observed = utils.toggle_warnings()
    assert_equal(state, not observed)


if __name__ == "__main__":
    nose.runmodule()
