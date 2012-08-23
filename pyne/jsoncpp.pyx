"""Python wrapper for jsoncpp."""
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.stdlib cimport malloc, free
from cython cimport pointer
from libc.string cimport const_char

# Python imports
import collections

# local imports
cimport std
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


cdef cpp_jsoncpp.Value * tocppval(object doc):
    cdef cpp_jsoncpp.Value * cval = NULL
#    if cvalue.isObject() or cvalue.isArray():
#            #pyvalue._inst = &cvalue
#            pyvalue._inst[0] = cvalue
#            return pyvalue
    if False:
        pass
    elif isinstance(doc, basestring):
        cval = new cpp_jsoncpp.Value(<char *> doc)
    elif isinstance(doc, float):
        cval = new cpp_jsoncpp.Value(<double> doc)
    elif isinstance(doc, bool):
        # NOTE: bool must go before int!
        # Python bools are ints, but ints are not bools.
        cval = toboolval(<bint> doc)
    elif isinstance(doc, int):
        cval = new cpp_jsoncpp.Value(<int> doc)
#        elif cvalue.isNull():
#            return None
    else:
        raise ValueError("{0} not of know type".format(doc))
    return cval



cdef class Value(object):
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
        cdef cpp_jsoncpp.Value cvalue
        cdef Value pyvalue = Value()
        #cdef Value pyvalue = Value(view=True)
        cdef std.string cstrvalue
        cdef int cintkey
        cdef std.string cstrkey

        # convert key and get value
        if isinstance(pykey, basestring):
            cstrkey = std.string(pykey)
            cvalue = self._inst.get(cstrkey, self._inst.null)
        elif isinstance(pykey, int) and self._inst.isArray():
            cintkey = pykey
            cvalue = self._inst.get(cintkey, self._inst.null)
        elif isinstance(pykey, slice) and self._inst.isArray():
            N = self._inst.size()
            #for 
        else:
            if (isinstance(pykey, int) or isinstance(pykey, slice)) and not self._inst.isArray():
                raise KeyError('key is int but object is not an array')
            else:            
                raise KeyError('key not of appropriate type, got {0}'.format(type(pykey)))

        # convert value
        if cvalue.isObject() or cvalue.isArray():
            #pyvalue._inst = &cvalue
            pyvalue._inst[0] = cvalue
            return pyvalue
        elif cvalue.isString():
            cstrvalue = cvalue.asString()
            return cstrvalue.c_str()
        elif cvalue.isDouble():
            return cvalue.asDouble()
        elif cvalue.isIntegral():
            return cvalue.asInt()
        elif cvalue.isBool():
            return cvalue.asBool()
        elif cvalue.isNull():
            return None
        else:
            raise ValueError("{0} not of know type".format(pykey))

    def __setitem__(self, key, value):
        cdef std.string cstrkey 
        cdef std.string cstrval
        cdef cpp_jsoncpp.Value * ckey 
        cdef cpp_jsoncpp.Value * cvalue

        cstrkey = std.string(key)
        ckey = &self._inst[0][cstrkey]
        
        cstrval = std.string(value)
        cvalue = new cpp_jsoncpp.Value(cstrval)
        ckey.swap(deref(cvalue))

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
            pass
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
            print "Here"
            return self._inst.asBool()
        else:
            print "Not Here"
            return NotImplemented


cdef class Reader:
    def __cinit__(self):
        """Reade C++ constuctor."""
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
        cdef std.string cdoc
        pydoc = json.JSONEncoder(separators=(',', ':')).encode(document)
        cdoc = std.string(pydoc)
        self._inst.parse(cdoc, root._inst[0], collect_comments)
        return root
