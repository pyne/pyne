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
    assert_true(v.isstring())
    assert_equal(v.type(), 4)
    assert_equal(v.type_name(), 'string')
    
def test_fltvalue():
    d = 42.0
    v = Value(d)
    assert_equal(float(v), d)
    assert_equal(int(v), 42)
    assert_true(v.isfloat())
    assert_equal(v.type(), 3)
    assert_equal(v.type_name(), 'real')

def test_intvalue():
    i = 42
    v = Value(i)
    assert_equal(float(v), 42.0)
    assert_equal(int(v), i)
    assert_true(v.isint())
    assert_equal(v.type(), 1)
    assert_equal(v.type_name(), 'int')
    
def test_truevalue():
    b = True
    v = Value(b)
    assert_true(v)
    assert_true(v.isbool())
    assert_equal(v.type(), 5)
    assert_equal(v.type_name(), 'boolean')

def test_falsevalue():
    b = False
    v = Value(b)
    assert_false(v)
    assert_true(v.isbool())
    assert_equal(v.type(), 5)
    assert_equal(v.type_name(), 'boolean')
    
def test_nonevalue():
    n = None
    v = Value(n)
    assert_true(v.isnull())
    assert_true(v.isarray())
    assert_true(v.isobject())
    assert_equal(v.type(), 0)
    assert_equal(v.type_name(), 'null')

def test_arrvalue():
    a = [1, 2, 5, 3]
    v = Value(a)
    assert_equal(len(v), len(a))
    assert_equal([v[i] for i in range(len(a))], a)
    assert_true(v.isarray())
    assert_equal(v.type(), 6)
    assert_equal(v.type_name(), 'array')

def test_tuplevalue():
    a = (1, 2, 5, 3)
    v = Value(a)
    assert_equal(len(v), len(a))
    assert_equal(tuple([v[i] for i in range(len(a))]), a)
    assert_true(v.isarray())
    assert_equal(v.type(), 6)
    assert_equal(v.type_name(), 'array')

def test_setvalue():
    a = set([1, 2, 5, 3])
    v = Value(a)
    assert_equal(len(v), len(a))
    assert_equal(set([v[i] for i in range(len(a))]), a)
    assert_true(v.isarray())
    assert_equal(v.type(), 6)
    assert_equal(v.type_name(), 'array')

def test_mapvalue():
    m = {'name': 'Terry Jones', 'age': 42.0}
    v = Value(m)
    assert_equal(len(v), len(m))
    assert_equal(dict([(k, v[k]) for k in m.keys()]), m)
    assert_true(v.isobject())
    assert_equal(v.type(), 7)
    assert_equal(v.type_name(), 'object')

def test_badmapvalue():
    m = {'name': 'Terry Jones', 42: 42.0}
    assert_raises(KeyError, Value, m)

def test_nestedvalue():
    lofl = [[1, 2, 3, True, 5, 6], '1', {'a': {'b': 42.0}}]
    lofl = [[1, 2, 3, False, 5, 6], '1', {'a': {'b': 42.0}}]
    lofl = [[1, 2, 3, None, 5, 6], '1', {'a': {'b': 42.0}}]
    v = Value(lofl)
    assert_equal(len(v), len(lofl))
    assert_equal(len(v[0]), len(lofl[0]))
    assert_equal(v[0][1], lofl[0][1])
    assert_equal(v[-1]['a']['b'], lofl[-1]['a']['b'])

def test_arrsetitem():
    l = ['Terry Jones', 1, None, 42.0]
    v = Value([None] * len(l))
    for i, value in enumerate(l):
        v[i] = value
    assert_equal([v[i] for i in range(len(l))], l)

def test_mapsetitem():
    m = {'name': 'Terry Jones', 'age': 42.0}
    v = Value({})
    for key, value in m.items():
        v[key] = value
    assert_equal(dict([(k, v[k]) for k in m.keys()]), m)

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

