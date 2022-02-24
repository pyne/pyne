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

import nose
from nose.tools import (
    assert_equal,
    assert_not_equal,
    assert_raises,
    raises,
    assert_in,
    assert_true,
    assert_false,
)

from pyne.jsoncpp import Value, Reader, FastWriter, StyledWriter, CustomWriter


def test_strvalue():
    s = "No one expects the Spanish Inquisition!!!!"
    v = Value(s)
    assert_equal(len(v), 42)
    assert_equal(repr(v), s)
    assert_true(v.isstring())
    assert_equal(v.type(), 4)
    assert_equal(v.type_name(), "string")


def test_fltvalue():
    d = 42.0
    v = Value(d)
    assert_equal(float(v), d)
    assert_equal(int(v), 42)
    assert_true(v.isfloat())
    assert_equal(v.type(), 3)
    assert_equal(v.type_name(), "real")


def test_intvalue():
    i = 42
    v = Value(i)
    assert_equal(float(v), 42.0)
    assert_equal(int(v), i)
    assert_true(v.isint())
    assert_equal(v.type(), 1)
    assert_equal(v.type_name(), "int")


def test_truevalue():
    b = True
    v = Value(b)
    assert_true(v)
    assert_true(v.isbool())
    assert_equal(v.type(), 5)
    assert_equal(v.type_name(), "boolean")


def test_falsevalue():
    b = False
    v = Value(b)
    assert_false(v)
    assert_true(v.isbool())
    assert_equal(v.type(), 5)
    assert_equal(v.type_name(), "boolean")


def test_nonevalue():
    n = None
    v = Value(n)
    assert_true(v.isnull())
    assert_true(v.isarray())
    assert_true(v.isobject())
    assert_equal(v.type(), 0)
    assert_equal(v.type_name(), "null")


def test_arrvalue():
    a = [1, 2, 5, 3]
    v = Value(a)
    assert_equal(len(v), len(a))
    assert_equal([v[i] for i in range(len(a))], a)
    assert_true(v.isarray())
    assert_equal(v.type(), 6)
    assert_equal(v.type_name(), "array")


def test_tuplevalue():
    a = (1, 2, 5, 3)
    v = Value(a)
    assert_equal(len(v), len(a))
    assert_equal(tuple([v[i] for i in range(len(a))]), a)
    assert_true(v.isarray())
    assert_equal(v.type(), 6)
    assert_equal(v.type_name(), "array")


def test_setvalue():
    a = set([1, 2, 5, 3])
    v = Value(a)
    assert_equal(len(v), len(a))
    assert_equal(set([v[i] for i in range(len(a))]), a)
    assert_true(v.isarray())
    assert_equal(v.type(), 6)
    assert_equal(v.type_name(), "array")


def test_mapvalue():
    m = {"name": "Terry Jones", "age": 42.0}
    v = Value(m)
    assert_equal(len(v), len(m))
    assert_equal(dict([(k, v[k]) for k in m]), m)
    assert_true(v.isobject())
    assert_equal(v.type(), 7)
    assert_equal(v.type_name(), "object")


def test_badmapvalue():
    m = {"name": "Terry Jones", 42: 42.0}
    assert_raises(KeyError, Value, m)


def test_nestedvalue():
    lofl = [[1, 2, 3, True, 5, 6], "1", {"a": {"b": 42.0}}]
    lofl = [[1, 2, 3, False, 5, 6], "1", {"a": {"b": 42.0}}]
    lofl = [[1, 2, 3, None, 5, 6], "1", {"a": {"b": 42.0}}]
    v = Value(lofl)
    assert_equal(len(v), len(lofl))
    assert_equal(len(v[0]), len(lofl[0]))
    assert_equal(v[0][1], lofl[0][1])
    assert_equal(v[-1]["a"]["b"], lofl[-1]["a"]["b"])


def test_arrsetitem():
    l = ["Terry Jones", 1, None, 42.0]
    v = Value([None] * len(l))
    for i, value in enumerate(l):
        v[i] = value
    assert_equal([v[i] for i in range(len(l))], l)


def test_mapsetitem():
    m = {"name": "Terry Jones", "age": 42.0}
    v = Value({})
    for key, value in m.items():
        v[key] = value
    assert_equal(dict([(k, v[k]) for k in m]), m)


def test_getslice():
    a = [1, 2, 5, 3]
    v = Value(a)

    t = v[1:-1]
    obs = [t[i] for i in range(2)]
    exp = [2, 5]
    assert_equal(obs, exp)

    t = v[::-1]
    obs = [t[i] for i in range(4)]
    exp = a[::-1]
    assert_equal(obs, exp)

    t = v[-3::-1]
    obs = [t[i] for i in range(2)]
    exp = a[-3::-1]
    assert_equal(obs, exp)


def test_setslice():
    a = [1, 2, 5, 3]
    v = Value(a)

    v[1:-1] = [42, 65]
    obs = [v[i] for i in range(len(a))]
    exp = [1, 42, 65, 3]
    assert_equal(obs, exp)

    v = Value(a)
    v[::-1] = "abcd"
    obs = [v[i] for i in range(4)]
    exp = ["d", "c", "b", "a"]
    assert_equal(obs, exp)

    v = Value(a)
    v[-3::-1] = [65, 42]
    obs = [v[i] for i in range(4)]
    exp = [42, 65, 5, 3]
    assert_equal(obs, exp)


def test_setvalue():
    a = Value({"i": 10, "j": "rawr"})
    b = Value(65.0)
    a["counter"] = b
    assert_equal(a["i"], 10)
    assert_equal(a["j"], "rawr")
    assert_equal(float(b), 65.0)
    assert_equal(a["counter"], 65.0)

    a = Value({"i": 10, "j": "rawr"})
    b = Value("burninating")
    a["counter"] = b
    assert_equal(a["i"], 10)
    assert_equal(a["j"], "rawr")
    assert_equal(str(b), "burninating")
    assert_equal(a["counter"], "burninating")

    a = Value({"i": 10, "j": "rawr"})
    b = Value([1, 2, 5, 3])
    a["counter"] = b
    assert_equal(a["i"], 10)
    assert_equal(a["j"], "rawr")
    assert_equal([b[i] for i in range(4)], [1, 2, 5, 3])
    assert_equal([a["counter"][i] for i in range(4)], [b[i] for i in range(4)])

    a = Value([1, 2, 5, 3])
    b = Value([42, 65])
    a[1:-1] = b
    assert_equal([b[i] for i in range(2)], [42, 65])
    assert_equal([a[i] for i in range(4)], [1, 42, 65, 3])


def test_delitem_contains():
    a = Value({"i": 10, "j": "rawr"})
    assert_true("i" in a)
    del a["i"]
    assert_false("i" in a)

    a = Value([1, 2, 5, 3])
    assert_equal(len(a), 4)
    assert_true(2 in a)
    del a[1]
    assert_false(2 in a)
    assert_equal(len(a), 3)

    a = Value([1, 2, 5, 3])
    assert_equal(len(a), 4)
    assert_true(2 in a)
    assert_true(5 in a)
    del a[1:-1]
    assert_false(2 in a)
    assert_false(5 in a)
    assert_equal(len(a), 2)

    a = Value([1, 2, 5, 3])
    assert_equal(len(a), 4)
    assert_true(1 in a)
    assert_true(2 in a)
    del a[-3::-1]
    assert_false(1 in a)
    assert_false(2 in a)
    assert_equal(len(a), 2)


def test_keys():
    a = Value({"i": 10, "j": "rawr"})
    assert_equal(a.keys(), ["i", "j"])


def test_vals():
    a = Value({"i": 10, "j": "rawr"})
    assert_equal(a.values(), [10, "rawr"])


def test_items():
    a = Value({"i": 10, "j": "rawr"})
    assert_equal(a.items(), [("i", 10), ("j", "rawr")])


def test_iter():
    a = Value({"i": 10, "j": "rawr"})
    assert_equal(a.keys(), [k for k in a])

    a = Value([1, 2, 5, 3])
    assert_equal([i for i in a], [1, 2, 5, 3])

    a = Value("rawr")
    assert_equal([i for i in a], ["r", "a", "w", "r"])


def test_get():
    a = Value({"i": 10, "j": "rawr"})
    assert_equal(a.get("i"), 10)
    assert_equal(a.get("wahhh"), None)
    assert_equal(a.get("wahhh", 42.0), 42.0)


def test_cmp():
    a = Value({"i": 10, "j": "rawr"})
    assert_true(a == {"i": 10, "j": "rawr"})
    assert_true(a != {"i": 10})

    a = Value(10)
    assert_true(a == 10)
    assert_true(a != 11)
    assert_true(a < 11)
    assert_true(a <= 11)
    assert_true(a > 9)
    assert_true(a >= 9)


def test_mutablemap():
    a = Value({"i": 10, "j": "rawr"})
    assert_equal(a.pop("i"), 10)
    assert_equal(a.popitem("j"), ("j", "rawr"))
    a.setdefault("z", "man")
    assert_true(a == {"z": "man"})
    a.update({"i": 10, "j": "rawr"})
    assert_true(a == {"i": 10, "j": "rawr", "z": "man"})
    a.clear()
    assert_equal(len(a), 0)


def test_mutableseq():
    pya = [1, 2, 5, 3]
    a = Value(pya)
    assert_equal([i for i in reversed(a)], pya[::-1])
    assert_equal(a.index(5), 2)
    assert_equal(a.index(5, 2), 2)
    assert_equal(a.index(2, 1, -1), 1)

    pya = [1, 2, 5, 3, 1, 1, 6]
    a = Value(pya)
    assert_equal(a.count(1), pya.count(1))
    assert_equal(a.count(5), pya.count(5))

    assert_equal(len(a), len(pya))
    a.append(42)
    assert_equal(len(a), len(pya) + 1)

    pya = [1, 2, 5, 3]
    a = Value(pya)
    a.insert(2, 42)
    pya.insert(2, 42)
    assert_equal([i for i in a], pya)
    a.insert(-3, 65)
    pya.insert(-3, 65)
    assert_equal([i for i in a], pya)

    pya = [1, 2, 5, 3]
    a = Value(pya)
    a.reverse()
    assert_equal([i for i in a], pya[::-1])
    pya = [1, 2, 42, 5, 3]
    a = Value(pya)
    a.reverse()
    assert_equal([i for i in a], pya[::-1])

    pya = [1, 2, 5, 3]
    a = Value(pya)
    a.extend([42, 65])
    pya.extend([42, 65])
    assert_equal([i for i in a], pya)

    pya = [1, 2, 5, 3]
    a = Value(pya)
    a.pop(-2)
    pya.pop(-2)
    assert_equal([i for i in a], pya)

    pya = [1, 2, 5, 3]
    a = Value(pya)
    a.remove(5)
    pya.remove(5)
    assert_equal([i for i in a], pya)

    pya = [1, 2, 5, 3]
    a = Value(pya)
    a += [42, 65]
    pya += [42, 65]
    assert_equal([i for i in a], pya)


def test_str_repr():
    assert_equal(repr(Value({"hello": 1})), '{"hello":1}')
    assert_equal(str(Value({"hello": 1})), '{\n   "hello" : 1\n}')
    assert_equal(repr(Value("hello")), "hello")
    assert_equal(str(Value("hello")), "hello")


def test_reader():
    r = Reader()
    value = Value({"hello": 1})
    valstr = repr(value)
    obs = r.parse(valstr)
    assert_true(value == obs)

    strio = StringIO(valstr)
    obs = r.parse(valstr)
    assert_true(value == obs)


def test_fastwriter():
    fw = FastWriter()
    assert_equal(fw.write(Value({"hello": 1})), '{"hello":1}\n')
    assert_equal(fw.write(Value("hello")), '"hello"\n')


def test_styledwriter():
    sw = StyledWriter()
    assert_equal(sw.write(Value({"hello": 1})), '{\n   "hello" : 1\n}\n')
    assert_equal(sw.write(Value("hello")), '"hello"\n')


def test_customwriter():
    cw = CustomWriter(colon=": ", closecurly="\n}")
    assert_equal(cw.write(Value({"hello": 1})), '{"hello": 1\n}\n')
    assert_equal(cw.write(Value("hello")), '"hello"\n')


if __name__ == "__main__":
    nose.runmodule()
