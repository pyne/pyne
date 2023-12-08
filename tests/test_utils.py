"""PyNE utility tests"""
import os

import pytest


from numpy.testing import assert_array_equal, assert_array_almost_equal

from pyne import utils
import numpy as np


def test_to_sec():
    # as
    assert 2e-18 == utils.to_sec(2, "AS")
    assert 2e-18 == utils.to_sec(2, "attosec")
    assert 2e-18 == utils.to_sec(2, "attosecond")
    assert 2e-18 == utils.to_sec(2, "attoseconds")
    # fs
    assert 2e-15 == utils.to_sec(2, "FS")
    assert 2e-15 == utils.to_sec(2, "femtosec")
    assert 2e-15 == utils.to_sec(2, "femtosecond")
    assert 2e-15 == utils.to_sec(2, "femtoseconds")
    # ps
    assert 2e-12 == utils.to_sec(2, "PS")
    assert 2e-12 == utils.to_sec(2, "picosec")
    assert 2e-12 == utils.to_sec(2, "picosecond")
    assert 2e-12 == utils.to_sec(2, "picoseconds")
    # ns
    assert 2e-9 == utils.to_sec(2, "NS")
    assert 2e-9 == utils.to_sec(2, "nanosec")
    assert 2e-9 == utils.to_sec(2, "nanosecond")
    assert 2e-9 == utils.to_sec(2, "nanoseconds")
    # us
    assert 2e-6 == utils.to_sec(2, "US")
    assert 2e-6 == utils.to_sec(2, "microsec")
    assert 2e-6 == utils.to_sec(2, "microsecond")
    assert 2e-6 == utils.to_sec(2, "microseconds")
    # ms
    assert 2e-3 == utils.to_sec(2, "MS")
    assert 2e-3 == utils.to_sec(2, "millisec")
    assert 2e-3 == utils.to_sec(2, "millisecond")
    assert 2e-3 == utils.to_sec(2, "milliseconds")
    # s
    assert 2.0 == utils.to_sec(2, "S")
    assert 2.0 == utils.to_sec(2, "sec")
    assert 2.0 == utils.to_sec(2, "second")
    assert 2.0 == utils.to_sec(2, "seconds")
    # m
    assert 120.0 == utils.to_sec(2, "M")
    assert 120.0 == utils.to_sec(2, "min")
    assert 120.0 == utils.to_sec(2, "minute")
    assert 120.0 == utils.to_sec(2, "minutes")
    # h
    assert 7200.0 == utils.to_sec(2, "H")
    assert 7200.0 == utils.to_sec(2, "hour")
    assert 7200.0 == utils.to_sec(2, "hours")
    # d
    assert 172800.0 == utils.to_sec(2, "D")
    assert 172800.0 == utils.to_sec(2, "day")
    assert 172800.0 == utils.to_sec(2, "days")
    # w
    assert 1209600.0 == utils.to_sec(2, "W")
    assert 1209600.0 == utils.to_sec(2, "week")
    assert 1209600.0 == utils.to_sec(2, "weeks")
    # y
    assert 63115200.0 == utils.to_sec(2, "Y")
    assert 63115200.0 == utils.to_sec(2, "year")
    assert 63115200.0 == utils.to_sec(2, "years")
    # c
    assert 6311520000.0 == utils.to_sec(2, "C")
    assert 6311520000.0 == utils.to_sec(2, "century")
    assert 6311520000.0 == utils.to_sec(2, "centuries")
    # undifined unit trigs ValueError
    pytest.raises(ValueError, utils.to_sec, 2, "month")


def test_to_barns():
    assert 3e3 == utils.to_barns(3, "KB")


def test_from_barns():
    assert 3e-3 == utils.from_barns(3, "KB")


def test_message():
    if os.name == "posix":
        assert "\033[1;32mHello\033[0m" == utils.message("Hello")
    else:
        assert "*** MESSAGE ***: Hello" == utils.message("Hello")


def test_failure():
    if os.name == "posix":
        assert "\033[1;31mWorld\033[0m" == utils.failure("World")
    else:
        assert "*** FAILURE ***: World" == utils.failure("World")


def test_remove():
    # test file removal
    assert "hello.txt" not in os.listdir(".")
    with open("hello.txt", "w") as f:
        f.write("Hello PyNE!")
    assert "hello.txt" in os.listdir(".")

    utils.remove("hello.txt")
    assert "hello.txt" not in os.listdir(".")

    # test recursive dir removal
    assert "rawr" not in os.listdir(".")
    os.mkdir("rawr")

    assert "hello.txt" not in os.listdir("rawr")
    with open("rawr/hello.txt", "w") as f:
        f.write("Hello PyNE!")
    assert "hello.txt" in os.listdir("rawr")

    utils.remove("rawr")
    assert "rawr" not in os.listdir(".")

    # Test pass through on non-existant
    assert "Medic!" not in os.listdir(".")
    utils.remove("Medic!")


def test_fromstring_split1():
    s = "1 2.3 4"
    obs = utils.fromstring_split(s)
    exp = np.fromstring(s, sep=" ")
    assert_array_equal(obs, exp)


def test_fromstring_split2():
    s = "1-2.3-4"
    obs = utils.fromstring_split(s, sep="-")
    exp = np.fromstring(s, sep="-")
    assert_array_equal(obs, exp)


def test_fromstring_split3():
    s = "1,2.3,4"
    obs = utils.fromstring_split(s, sep=",")
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
    obs = utils.fromstring_token(s, sep="-")
    exp = np.fromstring(s, sep="-")
    assert_array_equal(obs, exp)


def test_fromstring_token3():
    s = "1,  2.3, 4"
    obs = utils.fromstring_token(s, sep=" ,")
    exp = np.fromstring(s, sep=", ")
    assert_array_equal(obs, exp)


def test_fromstring_token4():
    s = "1, 2.3 ,4"
    obs = utils.fromstring_token(s, sep=" ,", inplace=True)
    exp = np.array([1.0, 2.3, 4.0])
    assert_array_equal(obs, exp)


def test_use_warnings():
    first = utils.use_warnings()
    second = utils.use_warnings()
    assert first == second


def test_toggle_warnings():
    state = utils.use_warnings()
    observed = utils.toggle_warnings()
    assert state == (not observed)

def test_str_to_unicode():
    """
    Convert binary str to unicode str.
    """
    exp_answer = "test"
    
    # default str
    s = "test"
    assert utils.str_to_unicode(s) == exp_answer
    
    # binary str
    s = b"test"
    assert utils.str_to_unicode(s) == exp_answer
    
    # unicode str
    s = "test"
    assert utils.str_to_unicode(s) == exp_answer

    # list of str
    s = ["test1", "test2", b"test3"]
    exp_answer = ["test1", "test2", "test3"]
    assert_array_equal(utils.str_to_unicode(s), exp_answer)
    
    # set of str
    s = {"test1", "test2", b"test3"}
    exp_answer = {"test1", "test2", "test3"}
    assert utils.str_to_unicode(s) == exp_answer

    # tuple of str
    s = ("test1", "test2", b"test3")
    exp_answer = ("test1", "test2", "test3")
    assert utils.str_to_unicode(s) == exp_answer


def test_str_almost_same():
    """
    Test utils.str_almost_same.
    """
    # exactly the same
    s1 = "test1"
    s2 = "test1"
    assert utils.str_almost_same(s1, s2)

    # exactly the same, can be converted to floats
    s1 = "1.2e-3"
    s2 = "1.2e-3"
    assert utils.str_almost_same(s1, s2)

    # almost the same, with default rel_tol=1e-9
    s1 = "9.5"
    s2 = "9.500000000001"
    assert utils.str_almost_same(s1, s2)

    # almost the same, with big rel_tol=1e-6
    s1 = "1.2e-3"
    s2 = "1.20000005e-3"
    assert utils.str_almost_same(s1, s2, rel_tol=1e-6)

    # different str, can be converted too floats
    s1 = "1.2e-3"
    s2 = "1.20000005e-3"
    assert utils.str_almost_same(s1, s2) == False

    # different str, can not be converted to floats
    s1 = "test1"
    s2 = "test2"
    assert utils.str_almost_same(s1, s2) == False


def test_line_almost_same():
    """
    Test utils.line_almost_same.
    """
    # exactly the same lines, w/o numbers
    l1 = "test strings"
    l2 = "test strings"
    assert utils.line_almost_same(l1, l2)

    # almost same lines, with numbers
    l1 = "test data 9.5"
    l2 = "test data 9.50000000001"
    assert utils.line_almost_same(l1, l2)

    # different lines, w/o numbers
    l1 = "test1 strings"
    l2 = "test2 strings"
    assert utils.line_almost_same(l1, l2) == False

    # different lines, with numbers
    l1 = "test data 9.5"
    l2 = "test data 9.5001"
    assert utils.line_almost_same(l1, l2) == False


def test_file_almost_same():
    """
    Test utils.file_almost_same.
    """
    # exactly the same, w/o numbers
    f1 = """l1\nl2 string"""
    f2 = """l1\nl2 string"""
    assert utils.file_almost_same(f1, f2)

    # almost the same, with numbers
    f1 = """l1\nl2 data 9.5"""
    f2 = """l1\nl2 data 9.5000000000001"""
    assert utils.file_almost_same(f1, f2)

    # different contents
    f1 = """l1\nl2 string1"""
    f2 = """l1\nl2 string2"""
    assert utils.file_almost_same(f1, f2) == False


def test_block_in_blocks():
    # exactly in
    block1 = """test1"""
    blocks2 = ["""test1""", """test2"""]
    assert utils.block_in_blocks(block1, blocks2)

    # with float number tolerable difference
    block1 = """test data  9.5"""
    blocks2 = ["""test data 9.500000000001""", """test2"""]
    assert utils.block_in_blocks(block1, blocks2)

    # block not in blocks
    block1 = """test1"""
    blocks2 = ["""test2""", """test3"""]
    assert utils.block_in_blocks(block1, blocks2) == False


def test_file_block_almost_same():
    # exactly same file
    f1 = """block1\n\nblock2"""
    f2 = """block1\n\nblock2"""
    assert utils.file_block_almost_same(f1, f2)

    # same block, different sequence
    f1 = """block1\n\nblock2"""
    f2 = """block2\n\nblock1"""
    assert utils.file_block_almost_same(f1, f2)

    # almost same block, same sequence
    f1 = """test data 9.5\n\nblock2"""
    f2 = """test data 9.500000000001\n\nblock2"""
    assert utils.file_block_almost_same(f1, f2)

    # almost same block, different sequence
    f1 = """test data 9.5\n\nblock2"""
    f2 = """block2\n\ntest data 9.500000000001"""
    assert utils.file_block_almost_same(f1, f2)

    # different block
    f1 = """block1\n\nblock2"""
    f2 = """block1\n\nblock3"""
    assert utils.file_block_almost_same(f1, f2) == False


def test_check_iterable():
    # list
    obj = ["a", 1, 1.0]
    assert utils.check_iterable(obj)
    # tuple
    obj = ("a", 1, 1.0)
    assert utils.check_iterable(obj)
    # dict
    obj = {1: "a", 2: 1, 3: 1.0}
    assert utils.check_iterable(obj)
    # set
    obj = set(["a", 1, 1.0])
    assert utils.check_iterable(obj)


def test_ifbar():
    loops = 100
    bar = utils.IfBar(
        "if bar print", max=loops, suffix="%(percent).1f%% - %(eta)ds", show=True
    )
    for i in range(loops):
        bar.next()
    bar.finish()
    bar = utils.IfBar(
        "if bar print", max=loops, suffix="%(percent).1f%% - %(eta)ds", show=False
    )
    for i in range(loops):
        bar.next()
    bar.finish()

