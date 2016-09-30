"""C++ wrapper for jsoncpp."""
from libc.string cimport const_char
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector

cdef extern from "json.h" namespace "Json":

    cdef enum ValueType:
        nullValue,
        intValue,
        uintValue,
        realValue,
        stringValue,
        booleanValue,
        arrayValue,
        objectValue

    cdef cppclass Value:
        Value null

        Value() except +
        Value(ValueType) except +
        Value(char *) except +
        Value(std_string) except +
        Value(double) except +
        Value(int) except +
        Value(bint) except +
        Value(bint, enum) except +
        Value(Value &) except +

        const_char * asCString() except +
        std_string asString() except +
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
        ValueType type() except +

        Value get(int, Value) except +
        Value get(std_string, Value) except +
        Value & operator[](int)
        Value & operator[](std_string)
        Value & operator[](const_char *)
        void swap(Value &) except +
        #Value & operator=(Value &) except +
        Value removeMember(std_string) except +
        Value removeMember(const_char *) except +

        bint isMember(std_string) except +
        bint isMember(const_char *) except +
        bint operator<(Value &) except +
        bint operator<=(Value &) except +
        bint operator==(Value &) except +
        bint operator!=(Value &) except +
        bint operator>(Value &) except +
        bint operator>=(Value &) except +
        int compare(Value &) except +

        std_vector[std_string] getMemberNames() except +

        int size() except +
        void resize(int) except +
        void clear() except +
        Value & append(Value &) except +

    cdef cppclass Reader:
        Reader() except +
        bint parse(std_string, Value) except +
        bint parse(std_string, Value, bint) except +

    cdef cppclass FastWriter:
        FastWriter() except +
        void enableYAMLCompatibility() except +
        std_string write(Value &)

    cdef cppclass StyledWriter:
        StyledWriter() except +
        std_string write(Value &)

cdef extern from "jsoncustomwriter.h" namespace "Json":

    cdef cppclass CustomWriter:
        CustomWriter() except +
        CustomWriter(std_string) except +
        CustomWriter(std_string, std_string) except +
        CustomWriter(std_string, std_string, std_string) except +
        CustomWriter(std_string, std_string, std_string, std_string) except +
        CustomWriter(std_string, std_string, std_string, std_string, std_string) except +
        CustomWriter(std_string, std_string, std_string, std_string, std_string,
                     std_string) except +
        CustomWriter(std_string, std_string, std_string, std_string, std_string,
                     std_string, std_string) except +
        CustomWriter(std_string, std_string, std_string, std_string, std_string,
                     std_string, std_string, int) except +
        std_string write(Value &)

