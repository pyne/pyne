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
        self.longSize   = struct.calcsize('q')
        self.floatSize  = struct.calcsize('f')
        self.doubleSize = struct.calcsize('d')

    def get_data(self, n, typeCode, itemSize):
        """
        Returns one or more items of a specified type at the current
        position within the data list. If more than one item is read,
        the items are returned in a list.
        """
        if self.pos >= self.numBytes:
            raise BinaryReaderError("Already read all data from record")
        
        values = struct.unpack('{0}{1}'.format(n,typeCode), self.data[self.pos:self.pos+itemSize*n])
        self.pos += itemSize * n
        return list(values)

    def get_int(self, n=1):
        """
        Returns one or more 4-byte integers.
        """
        return self.get_data(n,'i',self.intSize)
                             
    def get_long(self, n=1):
        """
        Returns one or more 8-byte integers.
        """
        return self.get_data(n,'q',self.longSize)
                             
    def get_float(self, n=1):
        """
        Returns one or more floats.
        """
        return self.get_data(n,'f',self.floatSize)
                             
    def get_double(self, n=1):
        """
        Returns one or more double
        """
        return self.get_data(n,'d',self.doubleSize)
                             
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
        return [s[i*length:(i+1)*length] for i in range(n)]

    def put_data(self, newdata, format, itemSize):
        """
        Packs a list of data objects at the current position with a
        specified format and data size.
        """
        try: 
            iter(newdata) 
        except TypeError: 
            newdata = [newdata]

        for i in range(len(newdata)):
            self.data += struct.pack(format,newdata[i])
            self.pos += itemSize
            self.numBytes += itemSize

    def put_int(self, data):
        """
        Pack a list of 4-byte integers.
        """
        self.put_data(data,'1i',self.intSize)

    def put_long(self,data):
        """
        Pack a list of 8-byte integers.
        """
        self.put_data(data,'1q',self.longSize)

    def put_float(self,data):
        """
        Pack a list of floats
        """
        self.put_data(data,'1f',self.floatSize)
        
    def put_double(self,data):
        """
        Pack a list of doubles
        """
        self.put_data(data,'1d',self.doubleSize)
        
    
    def put_string(self, data, length, n=1):
        """
        Packs a list of one or more double at the current
        position within the data list.
        """
        self.put_data(data,'{0}s'.format(length),length)
    
    def reset(self):
        self.pos = 0

    def __repr__(self):
        return "<Record: {0} bytes>".format(self.numBytes)


class _BinaryReader(object):
    """Reads/writes a binary file according to CCCC standards. This
    was created following Prof. James Paul Holloway's
    (hagar@umich.edu) alpha release of ccccutils written in C++ from
    2001.
    """
    
    def __init__(self,filename,mode='rb'):
        self.intSize = struct.calcsize('i')
        self.longSize   = struct.calcsize('q')
        self.f = open(filename, mode)
        
    def close(self):
        self.f.close()

    def get_int(self):
        (i,) = struct.unpack('i',self.f.read(self.intSize))
        return i

    def put_int(self,data):
        self.f.write(struct.pack('i',data))
                             
    def put_fortran_record(self,record):
        """Fortran formatted records start with an integer and end with the same
        integer that represents the number of bytes that the record is.
        """
        self.put_int(record.numBytes)
        self.f.write(record.data)
        self.put_int(record.numBytes)

        return 1

    def get_fortran_record(self):
        """Fortran formatted records start with an integer and end with the same
        integer that represents the number of bytes that the record is.
        """

        numBytes = self.get_int()

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

class FortranRecordError(BinaryReaderError):
    pass
