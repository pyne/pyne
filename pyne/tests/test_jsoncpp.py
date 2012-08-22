"""JsonCpp tests"""
import os
try:
    import simplejson as json
except ImportError:
    import json

import nose 
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, assert_in, \
                       assert_true, assert_false

from pyne.jsoncpp import Value


#def test_reader():
#    pass

def test_strvalue():
    s = "No one expects the Spanish Inquisition!!!!"
    v = Value(s)
    assert_equal(len(v), 42)
    assert_equal(str(v), s)
    
def test_fltvalue():
    d = 42.0
    v = Value(d)
    assert_equal(float(v), d)
    assert_equal(int(v), 42)

def test_intvalue():
    i = 42
    v = Value(i)
    assert_equal(float(v), 42.0)
    assert_equal(int(v), i)
    
def test_truevalue():
    b = True
    v = Value(b)
    assert_true(v)

def test_falsevalue():
    b = False
    v = Value(b)
    assert_false(v)
    

#r = jsoncpp.Reader()
#root = r.parse({'a': 10, 'b': 'Hello', 'c': {'d': [1, 2, 10.0]}})

"""\
print root
print root['a']
print root['b']
print root['c']['d'][1]

q = root['c']['d']
print q
r = root['c']['d']
print r

print "len root = ", len(root)
root['ORLY'] = 'YARLY'
print "len root = ", len(root)
print len(root['ORLY'])
print root['ORLY']


t = root['b']
print t
t = "Ho!"
print root['b']
"""


if __name__ == "__main__":
    nose.runmodule()

