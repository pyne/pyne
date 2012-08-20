"""C++ wrapper for jsoncpp."""
cimport std

cdef extern from "json/json.h" namespace "Json":
    cdef cppclass Value:
        Value null

        Value() except +

        char * asCString() except +
        std.string asString() except +
        int asInt() except +
        #uint asUInt() except +
        #double int asInt64() except +
        #double uint asUInt64() except +
        float asFloat() except +
        double asDouble() except +
        bint asBool() except +

        bint isNull() except +
        bint isBool() except +
        bint isInt() except +
        bint isUInt() except +
        bint isIntegral() except +
        bint isDouble() except +
        bint isNumeric() except +
        bint isString() except +
        bint isArray() except +
        bint isObject() except +

        Value get(int, Value) except +
        Value get(std.string, Value) except +

        int size() except +


    cdef cppclass Reader:
        Reader() except +
        bint parse(std.string, Value) except +
        bint parse(std.string, Value, bint) except +
