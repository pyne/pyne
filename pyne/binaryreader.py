#!/usr/bin/env python

"""Module for manipulating binary files, especially those with Fortran formatted
records.

"""

import struct

class _FortranRecord(object):
    """A single Fortran formatted record.

    Parameters
    ----------
    data : str
        A string of binary data.
    numBytes : int
        Total number of bytes in record.

    """

    def __init__(self, data, numBytes):
        """
        Initialize instance of Record object.
        """

        self.data       = data
        self.numBytes   = numBytes

        self.reset()
        self.intSize    = struct.calcsize('i')
        self.floatSize  = struct.calcsize('f')
        self.doubleSize = struct.calcsize('d')
    
    def get_int(self, n=1):
        """
        Returns one or more integers at the current position within
        the data list. If more than one integer is read, the integers
        are returned in a list.
        """

        if self.pos >= self.numBytes:
            raise BinaryReaderError("Already read all data from record")
        
        values = struct.unpack('{0}i'.format(n), self.data[self.pos:self.pos+self.intSize*n])
        self.pos += self.intSize * n
        if n == 1:
            return values[0]
        else:
            return list(values)
                             
    def get_float(self, n=1):
        """
        Returns one or more floats at the current position within the
        data list. If more than one float is read, the floats are
        returned in a list.
        """

        if self.pos >= self.numBytes:
            raise BinaryReaderError("Already read all data from record")
        
        values = struct.unpack('{0}f'.format(n), self.data[self.pos:self.pos+self.floatSize*n])
        self.pos += self.floatSize * n
        if n == 1:
            return values[0]
        else:
            return list(values)
    
    def get_double(self, n=1):
        """
        Returns one or more doubles at the current position within the
        data list. If more than one double is read, the doubles are
        returned in a list.
        """

        if self.pos >= self.numBytes:
            raise BinaryReaderError("Already read all data from record")
        
        values = struct.unpack('{0}d'.format(n),self.data[self.pos:self.pos+self.doubleSize*n])
        self.pos += self.doubleSize * n
        if n == 1:
            return values[0]
        else:
            return list(values)
    
    def get_string(self, length, n=1):
        """
        Returns a string of a specified length starting at the current
        position in the data list.
        """

        if self.pos >= self.numBytes:
            raise BinaryReaderError("Already read all data from record")
        
        relevantData = self.data[self.pos:self.pos+length*n]
        (s,) = struct.unpack('{0}s'.format(length*n), relevantData)
        self.pos += length*n
        if n == 1:
            return s
        else:
            return [s[i*length:(i+1)*length] for i in range(n)]

    def reset(self):
        self.pos = 0

    def __repr__(self):
        return "<Record: {0} bytes>".format(self.numBytes)


class _BinaryReader(object):
    """Reads a binary file according to CCCC standards. This was created
    following Prof. James Paul Holloway's (hagar@umich.edu) alpha release of
    ccccutils written in C++ from 2001.
    """
    
    def __init__(self, filename):
        self.intSize = struct.calcsize('i')
        self.floatSize = struct.calcsize('d')
        self.f = open(filename, 'rb')
        
    def get_int(self):
        (i,) = struct.unpack('i',self.f.read(self.intSize))
        return i
                             
    def get_float(self):
        (f,) =  struct.unpack('d',self.f.read(self.floatSize))
        return f

    def get_fortran_record(self):
        """Fortran formatted records start with an integer and end with the same
        integer that represents the number of bytes that the record is.
        """

        numBytes = self.get_int()
        if numBytes % self.intSize:
            raise InvalidFortranRecordError("Number of bytes in Fortran formatted record is not a multiple of the integer size")

        # Read numBytes from the record
        data = self.f.read(numBytes)
        
        # now read end of record
        numBytes2 = self.get_int()
        if numBytes2 != numBytes:
            raise InvalidFortranRecordError("Starting and matching integers in Fortran formatted record do not match.")
        
        return _FortranRecord(data, numBytes)


class BinaryReaderError(Exception):
    """Case class for all binary reader errors."""

    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return repr(self.msg)


class InvalidFortranRecordError(BinaryReaderError):
    pass
