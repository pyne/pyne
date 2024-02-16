"""JsonCpp tests"""
import os

try:
    import simplejson as json
except ImportError:
    import json
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import pytest

from pyne.jsoncpp import Value, Reader, FastWriter, StyledWriter, CustomWriter


def test_strvalue():
    s = "No one expects the Spanish Inquisition!!!!"
    v = Value(s)
    assert len(v) == 42
    assert repr(v) == s
    assert v.isstring()
    assert v.type() == 4
    assert v.type_name() == "string"


def test_fltvalue():
    d = 42.0
    v = Value(d)
    assert float(v) == d
    assert int(v) == 42
    assert v.isfloat()
    assert v.type() == 3
    assert v.type_name() == "real"


def test_intvalue():
    i = 42
    v = Value(i)
    assert float(v) == 42.0
    assert int(v) == i
    assert v.isint()
    assert v.type() == 1
    assert v.type_name() == "int"


def test_truevalue():
    b = True
    v = Value(b)
    assert v
    assert v.isbool()
    assert v.type() == 5
    assert v.type_name() == "boolean"


def test_falsevalue():
    b = False
    v = Value(b)
    assert not v
    assert v.isbool()
    assert v.type() == 5
    assert v.type_name() == "boolean"


def test_nonevalue():
    n = None
    v = Value(n)
    assert v.isnull()
    assert v.isarray()
    assert v.isobject()
    assert v.type() == 0
    assert v.type_name() == "null"


def test_arrvalue():
    a = [1, 2, 5, 3]
    v = Value(a)
    assert len(v) == len(a)
    assert [v[i] for i in range(len(a))] == a
    assert v.isarray()
    assert v.type() == 6
    assert v.type_name() == "array"


def test_tuplevalue():
    a = (1, 2, 5, 3)
    v = Value(a)
    assert len(v) == len(a)
    assert tuple([v[i] for i in range(len(a))]) == a
    assert v.isarray()
    assert v.type() == 6
    assert v.type_name() == "array"


def test_setvalue():
    a = set([1, 2, 5, 3])
    v = Value(a)
    assert len(v) == len(a)
    assert set([v[i] for i in range(len(a))]) == a
    assert v.isarray()
    assert v.type() == 6
    assert v.type_name() == "array"


def test_mapvalue():
    m = {"name": "Terry Jones", "age": 42.0}
    v = Value(m)
    assert len(v) == len(m)
    assert dict([(k, v[k]) for k in m]) == m
    assert v.isobject()
    assert v.type() == 7
    assert v.type_name() == "object"


def test_badmapvalue():
    m = {"name": "Terry Jones", 42: 42.0}
    pytest.raises(KeyError, Value, m)


def test_nestedvalue():
    lofl = [[1, 2, 3, True, 5, 6], "1", {"a": {"b": 42.0}}]
    lofl = [[1, 2, 3, False, 5, 6], "1", {"a": {"b": 42.0}}]
    lofl = [[1, 2, 3, None, 5, 6], "1", {"a": {"b": 42.0}}]
    v = Value(lofl)
    assert len(v) == len(lofl)
    assert len(v[0]) == len(lofl[0])
    assert v[0][1] == lofl[0][1]
    assert v[-1]["a"]["b"] == lofl[-1]["a"]["b"]


def test_arrsetitem():
    l = ["Terry Jones", 1, None, 42.0]
    v = Value([None] * len(l))
    for i, value in enumerate(l):
        v[i] = value
    assert [v[i] for i in range(len(l))] == l


def test_mapsetitem():
    m = {"name": "Terry Jones", "age": 42.0}
    v = Value({})
    for key, value in m.items():
        v[key] = value
    assert dict([(k, v[k]) for k in m]) == m


def test_getslice():
    a = [1, 2, 5, 3]
    v = Value(a)

    t = v[1:-1]
    obs = [t[i] for i in range(2)]
    exp = [2, 5]
    assert obs == exp

    t = v[::-1]
    obs = [t[i] for i in range(4)]
    exp = a[::-1]
    assert obs == exp

    t = v[-3::-1]
    obs = [t[i] for i in range(2)]
    exp = a[-3::-1]
    assert obs == exp


def test_setslice():
    a = [1, 2, 5, 3]
    v = Value(a)

    v[1:-1] = [42, 65]
    obs = [v[i] for i in range(len(a))]
    exp = [1, 42, 65, 3]
    assert obs == exp

    v = Value(a)
    v[::-1] = "abcd"
    obs = [v[i] for i in range(4)]
    exp = ["d", "c", "b", "a"]
    assert obs == exp

    v = Value(a)
    v[-3::-1] = [65, 42]
    obs = [v[i] for i in range(4)]
    exp = [42, 65, 5, 3]
    assert obs == exp


def test_setvalue():
    a = Value({"i": 10, "j": "rawr"})
    b = Value(65.0)
    a["counter"] = b
    assert a["i"] == 10
    assert a["j"] == "rawr"
    assert float(b) == 65.0
    assert a["counter"] == 65.0

    a = Value({"i": 10, "j": "rawr"})
    b = Value("burninating")
    a["counter"] = b
    assert a["i"] == 10
    assert a["j"] == "rawr"
    assert str(b) == "burninating"
    assert a["counter"] == "burninating"

    a = Value({"i": 10, "j": "rawr"})
    b = Value([1, 2, 5, 3])
    a["counter"] = b
    assert a["i"] == 10
    assert a["j"] == "rawr"
    assert [b[i] for i in range(4)] == [1, 2, 5, 3]
    assert [a["counter"][i] for i in range(4)] == [b[i] for i in range(4)]

    a = Value([1, 2, 5, 3])
    b = Value([42, 65])
    a[1:-1] = b
    assert [b[i] for i in range(2)] == [42, 65]
    assert [a[i] for i in range(4)] == [1, 42, 65, 3]


def test_delitem_contains():
    a = Value({"i": 10, "j": "rawr"})
    assert "i" in a
    del a["i"]
    assert not "i" in a

    a = Value([1, 2, 5, 3])
    assert len(a) == 4
    assert 2 in a
    del a[1]
    assert not 2 in a
    assert len(a) == 3

    a = Value([1, 2, 5, 3])
    assert len(a) == 4
    assert 2 in a
    assert 5 in a
    del a[1:-1]
    assert not 2 in a
    assert not 5 in a
    assert len(a) == 2

    a = Value([1, 2, 5, 3])
    assert len(a) == 4
    assert 1 in a
    assert 2 in a
    del a[-3::-1]
    assert not 1 in a
    assert not 2 in a
    assert len(a) == 2


def test_keys():
    a = Value({"i": 10, "j": "rawr"})
    assert a.keys() == ["i", "j"]


def test_vals():
    a = Value({"i": 10, "j": "rawr"})
    assert a.values() == [10, "rawr"]


def test_items():
    a = Value({"i": 10, "j": "rawr"})
    assert a.items() == [("i", 10), ("j", "rawr")]


def test_iter():
    a = Value({"i": 10, "j": "rawr"})
    assert a.keys() == [k for k in a]

    a = Value([1, 2, 5, 3])
    assert [i for i in a] == [1, 2, 5, 3]

    a = Value("rawr")
    assert [i for i in a] == ["r", "a", "w", "r"]


def test_get():
    a = Value({"i": 10, "j": "rawr"})
    assert a.get("i") == 10
    assert a.get("wahhh") == None
    assert a.get("wahhh", 42.0) == 42.0


def test_cmp():
    a = Value({"i": 10, "j": "rawr"})
    assert a == {"i": 10, "j": "rawr"}
    assert a != {"i": 10}

    a = Value(10)
    assert a == 10
    assert a != 11
    assert a < 11
    assert a <= 11
    assert a > 9
    assert a >= 9


def test_mutablemap():
    a = Value({"i": 10, "j": "rawr"})
    assert a.pop("i") == 10
    assert a.popitem("j") == ("j", "rawr")
    a.setdefault("z", "man")
    assert a == {"z": "man"}
    a.update({"i": 10, "j": "rawr"})
    assert a == {"i": 10, "j": "rawr", "z": "man"}
    a.clear()
    assert len(a) == 0


def test_mutableseq():
    pya = [1, 2, 5, 3]
    a = Value(pya)
    assert [i for i in reversed(a)] == pya[::-1]
    assert a.index(5) == 2
    assert a.index(5, 2) == 2
    assert a.index(2, 1, -1) == 1

    pya = [1, 2, 5, 3, 1, 1, 6]
    a = Value(pya)
    assert a.count(1) == pya.count(1)
    assert a.count(5) == pya.count(5)

    assert len(a) == len(pya)
    a.append(42)
    assert len(a) == len(pya) + 1

    pya = [1, 2, 5, 3]
    a = Value(pya)
    a.insert(2, 42)
    pya.insert(2, 42)
    assert [i for i in a] == pya
    a.insert(-3, 65)
    pya.insert(-3, 65)
    assert [i for i in a] == pya

    pya = [1, 2, 5, 3]
    a = Value(pya)
    a.reverse()
    assert [i for i in a] == pya[::-1]
    pya = [1, 2, 42, 5, 3]
    a = Value(pya)
    a.reverse()
    assert [i for i in a] == pya[::-1]

    pya = [1, 2, 5, 3]
    a = Value(pya)
    a.extend([42, 65])
    pya.extend([42, 65])
    assert [i for i in a] == pya

    pya = [1, 2, 5, 3]
    a = Value(pya)
    a.pop(-2)
    pya.pop(-2)
    assert [i for i in a] == pya

    pya = [1, 2, 5, 3]
    a = Value(pya)
    a.remove(5)
    pya.remove(5)
    assert [i for i in a] == pya

    pya = [1, 2, 5, 3]
    a = Value(pya)
    a += [42, 65]
    pya += [42, 65]
    assert [i for i in a] == pya


def test_str_repr():
    assert repr(Value({"hello": 1})) == '{"hello":1}'
    assert str(Value({"hello": 1})) == '{\n   "hello" : 1\n}'
    assert repr(Value("hello")) == "hello"
    assert str(Value("hello")) == "hello"


def test_reader():
    r = Reader()
    value = Value({"hello": 1})
    valstr = repr(value)
    obs = r.parse(valstr)
    assert value == obs

    strio = StringIO(valstr)
    obs = r.parse(valstr)
    assert value == obs


def test_fastwriter():
    fw = FastWriter()
    assert fw.write(Value({"hello": 1})) == '{"hello":1}\n'
    assert fw.write(Value("hello")) == '"hello"\n'


def test_styledwriter():
    sw = StyledWriter()
    assert sw.write(Value({"hello": 1})) == '{\n   "hello" : 1\n}\n'
    assert sw.write(Value("hello")) == '"hello"\n'


def test_customwriter():
    cw = CustomWriter(colon=": ", closecurly="\n}")
    assert cw.write(Value({"hello": 1})) == '{"hello": 1\n}\n'
    assert cw.write(Value("hello")) == '"hello"\n'

