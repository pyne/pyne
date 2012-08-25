"""Python wrapper for jsoncpp."""
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.stdlib cimport malloc, free
from cython cimport pointer
from libc.string cimport const_char, memcpy
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector

# Python imports
import collections

# local imports
cimport cpp_jsoncpp

try:
    import simplejson as json
except ImportError:
    import json

cdef cpp_jsoncpp.Value * toboolval(bint b):
    # NOTE: This is a little hack-y but has to be done since
    # Cython bints are not actually C++ bools
    cdef cpp_jsoncpp.Value * cval = \
            new cpp_jsoncpp.Value(<cpp_jsoncpp.ValueType> cpp_jsoncpp.booleanValue)
    cdef cpp_jsoncpp.Reader reader= cpp_jsoncpp.Reader()
    if b:
        reader.parse('true', deref(cval), 0)
    return cval


cdef cpp_jsoncpp.Value * tocppval(object doc) except NULL:
    cdef cpp_jsoncpp.Value * cval = NULL
    if isinstance(doc, Value):
        cval = new cpp_jsoncpp.Value(<cpp_jsoncpp.Value &> (<Value> doc)._inst[0])
    elif isinstance(doc, collections.Mapping):
        cval = new cpp_jsoncpp.Value(<cpp_jsoncpp.ValueType> cpp_jsoncpp.objectValue)
        for k, v in doc.items():
            if not isinstance(k, basestring):
                raise KeyError('object keys must be strings, got {0}'.format(k))
            cval[0][<const_char *> k].swap(deref(tocppval(v)))
    elif isinstance(doc, basestring):
        # string must come before other sequences
        cval = new cpp_jsoncpp.Value(<char *> doc)
    elif isinstance(doc, collections.Sequence) or isinstance(doc, collections.Set):
        cval = new cpp_jsoncpp.Value(<cpp_jsoncpp.ValueType> cpp_jsoncpp.arrayValue)
        cval.resize(<int> len(doc))
        for i, d in enumerate(doc):
            cval[0][<int> i].swap(deref(tocppval(d)))
    elif isinstance(doc, float):
        cval = new cpp_jsoncpp.Value(<double> doc)
    elif isinstance(doc, bool):
        # NOTE: bool must go before int!
        # Python bools are ints, but ints are not bools.
        cval = toboolval(<bint> doc)
    elif isinstance(doc, int):
        cval = new cpp_jsoncpp.Value(<int> doc)
    elif doc is None:
        cval = new cpp_jsoncpp.Value(<cpp_jsoncpp.ValueType> cpp_jsoncpp.nullValue)
        #cval = <cpp_jsoncpp.Value *> &((new cpp_jsoncpp.Value()).null)
    else:
        raise ValueError("{0} not of known type".format(doc))
    return cval


cdef int toposindex(int i, int I) except -1:
    cdef int valid_i = i 
    if valid_i < 0:
        valid_i = I + valid_i
    if (I <= valid_i) or (valid_i < 0):
        raise IndexError
    return valid_i

cdef class Value(object):

    _value_type_names = ['null', 'int', 'uint', 'real', 'string', 'boolean',
                         'array', 'object']

    def __cinit__(self, document=None, bint view=False):
        """Value C++ constuctor."""
        self._view = view
        if view:
            self._inst = NULL
        elif document is not None:
            self._inst = tocppval(document)
        else:
            self._inst = new cpp_jsoncpp.Value()

    def __dealloc__(self):
        """Value C++ destructor."""
        if not self._view:
            del self._inst

    def __getitem__(self, pykey):
        cdef cpp_jsoncpp.Value * cvalue
        cdef Value pyvalue = Value(view=True)

        # convert key and get value
        if isinstance(pykey, basestring):
            cvalue = &self._inst[0][<const_char *> pykey]
        elif isinstance(pykey, int) and (self._inst.type() == cpp_jsoncpp.arrayValue):
            pykey = toposindex(pykey, self._inst[0].size())
            cvalue = &self._inst[0][<int> pykey]
        elif isinstance(pykey, slice) and (self._inst.type() == cpp_jsoncpp.arrayValue):
            pyvalue._view = False
            cvalue = new cpp_jsoncpp.Value(<cpp_jsoncpp.ValueType> cpp_jsoncpp.arrayValue)
            curr_size = self._inst.size()
            r = range(*pykey.indices(curr_size))
            cvalue.resize(<int> len(r))
            for i, j in enumerate(r):
                j = toposindex(j, curr_size)
                cvalue[0][<int> i].swap(cpp_jsoncpp.Value(self._inst[0][<int> j]))
        else:
            if (isinstance(pykey, int) or isinstance(pykey, slice)) and not \
               (self._inst.type() == cpp_jsoncpp.arrayValue):
                raise KeyError('key is int or slice but object is not an array')
            else:            
                raise KeyError('key not of appropriate type, got {0}'.format(type(pykey)))

        # convert value
        if (cvalue.type() == cpp_jsoncpp.objectValue) or \
           (cvalue.type() == cpp_jsoncpp.arrayValue):
            pyvalue._inst = cvalue
            return pyvalue
        elif cvalue.isString():
            return <char *> cvalue.asCString()
        elif cvalue.isDouble():
            return cvalue.asDouble()
        elif cvalue.isBool():
            return cvalue.asBool()
        elif cvalue.isIntegral():
            return cvalue.asInt()
        elif cvalue.isNull():
            return None
        else:
            raise ValueError("{0} not of known type".format(pykey))

    def __setitem__(self, key, value):
        cdef cpp_jsoncpp.Value * ckey = NULL
        if isinstance(key, basestring):
            ckey = &self._inst[0][<const_char *> key]
            ckey.swap(deref(tocppval(value)))
        elif isinstance(key, int):
            curr_size = self._inst[0].size()
            key = toposindex(key, curr_size)
            ckey = &self._inst[0][<int> key]
            ckey.swap(deref(tocppval(value)))
        elif isinstance(key, slice):
            curr_size = self._inst[0].size()
            r = range(*key.indices(curr_size))
            for i, v in zip(r, value):
                i = toposindex(i, curr_size)
                ckey = &self._inst[0][<int> i]
                ckey.swap(deref(tocppval(v)))
        else:
            raise KeyError('key not of appropriate type, got {0}'.format(type(key)))

    def __delitem__(self, key):
        cdef int i, ikey, curr_size, end_size
        cdef cpp_jsoncpp.Value ctemp
        if isinstance(key, basestring) and (self._inst.type() == cpp_jsoncpp.objectValue):
            self._inst.removeMember(<const_char *> key)
        elif isinstance(key, int) and (self._inst.type() == cpp_jsoncpp.arrayValue):
            curr_size = self._inst[0].size()
            ikey = key
            ikey = toposindex(ikey, curr_size)
            for i in range(ikey+1, curr_size):
                self._inst[0][i-1].swap(self._inst[0][i])
            self._inst.resize(curr_size-1)
        elif isinstance(key, slice) and (self._inst.type() == cpp_jsoncpp.arrayValue):
            curr_size = self._inst[0].size()
            r = range(curr_size)
            del r[key]
            end_size = len(r)
            ctemp = cpp_jsoncpp.Value(cpp_jsoncpp.arrayValue)
            ctemp.resize(end_size)
            for i, r_i in enumerate(r):
                ctemp[i].swap(self._inst[0][<int> r_i])
            for i in range(end_size):
                self._inst[0][i].swap(ctemp[i])
            self._inst.resize(end_size)
        else:
            raise KeyError('key or object not of appropriate type')

    def __contains__(self, item):
        cdef int i, curr_size
        if isinstance(item, basestring) and (self._inst.type() == cpp_jsoncpp.objectValue):
            return self._inst.isMember(<const_char *> item)
        elif (self._inst.type() == cpp_jsoncpp.arrayValue):
            i = 0
            curr_size = self._inst[0].size()
            for i in range(curr_size):
                if item == self[i]:
                    return True
            return False
        else:
            raise KeyError('key or object not of appropriate type')

    def __iter__(self):
        if (self._inst.type() == cpp_jsoncpp.objectValue):
            return iter(self.keys())
        elif (self._inst.type() == cpp_jsoncpp.arrayValue):
            return iter([self[i] for i in range(len(self))])
        elif (self._inst.type() == cpp_jsoncpp.stringValue):
            return iter(str(self))
        else:
            raise TypeError('JSON Value not of appropriate type')
        

    def __len__(self):
        if self._inst.isObject() or self._inst.isArray():
            return self._inst.size()
        elif self._inst.isString():
            return len(str(self))
        else:
            raise TypeError("JSON Value has no length")

    def __str__(self):
        cdef const_char * s = NULL
        if self._inst.isString():
            s = self._inst.asCString()
        else:
            # FIXME: add writer here
            s = "Nothing here yet!"
        return s

    def __float__(self):
        if self._inst.isNumeric():
            return self._inst.asDouble()
        else:
            return NotImplemented

    def __int__(self):
        if self._inst.isNumeric():
            return self._inst.asInt()
        else:
            return NotImplemented

    def __nonzero__(self):
        if self._inst.isBool():
            return self._inst.asBool()
        else:
            return NotImplemented

    def isnull(self):
        """True if JSON null, False otherwise."""
        return self._inst.isNull()

    def isbool(self):
        """True if JSON boolean, False otherwise."""
        return self._inst.isBool()

    def isint(self):
        """True if is any JSON integer type, False otherwise."""
        return self._inst.isIntegral()

    def isfloat(self):
        """True if is any JSON float or double type, False otherwise."""
        return self._inst.isDouble()

    def isstring(self):
        """True if JSON string, False otherwise."""
        return self._inst.isString()

    def isarray(self):
        """True if JSON array or null, False otherwise."""
        return self._inst.isArray()

    def isobject(self):
        """True if JSON object or null, False otherwise."""
        return self._inst.isObject()

    def type(self):
        """The type number of this JSON value."""
        return self._inst.type()

    def type_name(self):
        """The type name of this JSON value."""
        return self._value_type_names[self._inst.type()]

    def keys(self):
        """Returns a list of keys in JSON object."""
        cdef std_vector[std_string] ckeys
        if (self._inst.type() != cpp_jsoncpp.objectValue):
            raise TypeError("no keys, not JSON object.")
        ckeys = self._inst.getMemberNames()
        return ckeys

    def values(self):
        """Returns a list of values in JSON object."""
        cdef std_vector[std_string] ckeys
        if (self._inst.type() != cpp_jsoncpp.objectValue):
            raise TypeError("no values, not JSON object.")
        ckeys = self._inst.getMemberNames()
        vals = [self[k] for k in ckeys]
        return vals

    def items(self):
        """Returns a list of items in JSON object."""
        cdef std_vector[std_string] ckeys
        if (self._inst.type() != cpp_jsoncpp.objectValue):
            raise TypeError("no values, not JSON object.")
        ckeys = self._inst.getMemberNames()
        its = [(k, self[k]) for k in ckeys]
        return its

    def get(self, key, default=None):
        """Returns key if present, or default otherwise."""
        if (self._inst.type() != cpp_jsoncpp.objectValue):
            raise TypeError("no keys, not JSON object.")
        if self._inst.isMember(<const_char *> key):
            return self[key]
        else:
            return default

    def __richcmp__(self, other, int op):
        cdef Value cother = Value(other)
        if op == 0:
            return ((<Value> self)._inst[0] < cother._inst[0])
        elif op == 1:
            return ((<Value> self)._inst[0] <= cother._inst[0])
        elif op == 2:
            return ((<Value> self)._inst[0] == cother._inst[0])
        elif op == 3:
            return ((<Value> self)._inst[0] != cother._inst[0])
        elif op == 4:
            return ((<Value> self)._inst[0] > cother._inst[0])
        elif op == 5:
            return ((<Value> self)._inst[0] >= cother._inst[0])

    def __cmp__(self, other):
        cdef Value cother = Value(other)
        return self._inst.compare((<Value> other)._inst[0])

cdef class Reader:
    def __cinit__(self):
        """Reader C++ constuctor."""
        self._inst = new cpp_jsoncpp.Reader()

    def __dealloc__(self):
        """Value C++ destructor."""
        del self._inst

    def parse(self, document, bint collect_comments=True):
        """Read a Value from a JSON document.

        Parameters
        ----------
        document : string, dict, file-object, list, etc
            Anything sufficiently JSON-izable
        collect_comments : bool, optional 
            True to collect comment and allow writing them back during
            serialization, and False to discard comments.

        Returns
        -------
        root : Value 
            The root value of the document if it was successfully parsed.

        """
        cdef Value root = Value()
        cdef std_string cdoc
        pydoc = json.JSONEncoder(separators=(',', ':')).encode(document)
        cdoc = pydoc
        self._inst.parse(cdoc, root._inst[0], collect_comments)
        return root
