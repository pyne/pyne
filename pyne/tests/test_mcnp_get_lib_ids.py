from pyne.mcnp import Xsdir
from nose.tools import assert_equal

xsFile= Xsdir("sample_xsdir")

def test_get_lib_ids1():
    assert_equal(xsFile.get_lib_ids('922350',"None" ), {'922350': '67y'})

def test_get_lib_ids2():
    assert_equal(xsFile.get_lib_ids('922380',"None" ), {'922380': '22c'})

def test_get_lib_ids3():
    assert_equal(xsFile.get_lib_ids('2003',"None" ), {'2003': '21c'})

def test_get_lib_ids4():
    assert_equal(xsFile.get_lib_ids('922350',"21c 22c 67y" ), {'922350': '21c'})

def test_get_lib_ids5():
    assert_equal(xsFile.get_lib_ids('922380',"21c 22c 67y" ), {'922380': '21c'})

def test_get_lib_ids6():
    assert_equal(xsFile.get_lib_ids('2003',"24y 22c 67y" ), {'2003': '24y'})

def test_get_lib_ids7():
    assert_equal(xsFile.get_lib_ids('922350',"" ), 0)

def test_get_lib_ids8():
    assert_equal(xsFile.get_lib_ids('992233',"None" ), 0)
