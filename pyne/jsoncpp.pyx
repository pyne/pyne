"""Python wrapper for jsoncpp."""
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.stdlib cimport malloc, free
from cython cimport pointer

# Python imports
import collections

# local imports
cimport std
cimport cpp_jsoncpp

try:
    import simplejson as json
except ImportError:
    import json


cdef class Value:
    def __cinit__(self):
        """Value C++ constuctor."""
        self._inst = new cpp_jsoncpp.Value()

    def __dealloc__(self):
        """Value C++ destructor."""
        del self._inst

    def __getitem__(self, pykey):
        cdef cpp_jsoncpp.Value cvalue
        #cdef cpp_jsoncpp.Value * cvalueref
        cdef Value pyvalue = Value()
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
            for 
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
            #return const_cast<char *> cvalue.asCString()
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

    def __len__(self):
        return self._inst.size()


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
