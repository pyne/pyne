#!/usr/bin/env python

"""Module for manipulating binary files, especially those with
Fortran formatted records.

"""
import struct

try:
    from collections.abc import Iterable
except ImportError:
    from collections import Iterable

from pyne.utils import QA_warn

QA_warn(__name__)


class _FortranRecord(object):
    """A single Fortran formatted record.

    Parameters
    ----------
    data : str
        A string of binary data.
    num_bytes : int
        Total number of bytes in record.

    """

    def __init__(self, data, num_bytes):
        """Initialize instance of Record object."""
        if isinstance(data, str):
            data = data.encode()
        self.data = data
        self.num_bytes = num_bytes

        self.reset()
        self.int_size = struct.calcsize("i")
        self.long_size = struct.calcsize("q")
        self.float_size = struct.calcsize("f")
        self.double_size = struct.calcsize("d")

    def get_data(self, n, typeCode, item_size):
        """
        Returns one or more items of a specified type at the current
        position within the data list. If more than one item is read,
        the items are returned in a list.
        """
        if self.pos >= self.num_bytes:
            raise ValueError(
                "All data read from record, pos="
                + str(self.pos)
                + " >= num_bytes="
                + str(self.num_bytes)
            )

        values = struct.unpack(
            "{0}{1}".format(n, typeCode), self.data[self.pos : self.pos + item_size * n]
        )
        self.pos += item_size * n
        return list(values)

    def get_int(self, n=1):
        """
        Returns one or more 4-byte integers.
        """
        return self.get_data(n, "i", self.int_size)

    def get_long(self, n=1):
        """
        Returns one or more 8-byte integers.
        """
        return self.get_data(n, "q", self.long_size)

    def get_float(self, n=1):
        """Returns one or more floats."""
        return self.get_data(n, "f", self.float_size)

    def get_double(self, n=1):
        """
        Returns one or more double
        """
        return self.get_data(n, "d", self.double_size)

    def get_string(self, length, n=1):
        """Returns a string of a specified length starting at the current
        position in the data list.
        """

        if self.pos >= self.num_bytes:
            raise ValueError(
                "All data read from record, pos="
                + str(self.pos)
                + " >= num_bytes="
                + str(self.num_bytes)
            )

        relevantData = self.data[self.pos : self.pos + length * n]
        (s,) = struct.unpack("{0}s".format(length * n), relevantData)
        self.pos += length * n
        return [s[i * length : (i + 1) * length].decode() for i in range(n)]

    def put_data(self, newdata, format, item_size):
        """
        Packs a list of data objects at the current position with a
        specified format and data size.
        """
        if not isinstance(newdata, Iterable):
            newdata = [newdata]

        for i in range(len(newdata)):
            nd = newdata[i]
            if isinstance(nd, str):
                nd = nd.encode()
            self.data += struct.pack(format, nd)
            self.pos += item_size
            self.num_bytes += item_size

    def put_int(self, data):
        """
        Pack a list of 4-byte integers.
        """
        self.put_data(data, "1i", self.int_size)

    def put_long(self, data):
        """
        Pack a list of 8-byte integers.
        """
        self.put_data(data, "1q", self.long_size)

    def put_float(self, data):
        """Pack a list of floats"""
        self.put_data(data, "1f", self.float_size)

    def put_double(self, data):
        """Pack a list of doubles."""
        self.put_data(data, "1d", self.double_size)

    def put_string(self, data, length, n=1):
        """Packs a list of one or more double at the current
        position within the data list.
        """
        self.put_data(data, "{0}s".format(length), length)

    def reset(self):
        self.pos = 0

    def __repr__(self):
        return "<Record: {0} bytes>".format(self.num_bytes)


class _BinaryReader(object):
    """Reads/writes a binary file according to CCCC standards. This
    was created following Prof. James Paul Holloway's
    (hagar@umich.edu) alpha release of ccccutils written in C++ from
    2001.
    """

    def __init__(self, filename, mode="rb"):
        self.int_size = struct.calcsize("i")
        self.long_size = struct.calcsize("q")
        self.f = open(filename, mode)

    def close(self):
        self.f.close()

    def get_int(self):
        (i,) = struct.unpack("i", self.f.read(self.int_size))
        return i

    def put_int(self, data):
        self.f.write(struct.pack("i", data))

    def put_fortran_record(self, record):
        """Fortran formatted records start with an integer and end with
        the same integer that represents the number of bytes
        that the record is.
        """
        self.put_int(record.num_bytes)
        self.f.write(record.data)
        self.put_int(record.num_bytes)
        return 1

    def get_fortran_record(self):
        """Fortran formatted records start with an integer and
        end with the same integer that represents the number of bytes
        that the record is.
        """

        num_bytes = self.get_int()

        # Read num_bytes from the record
        data = self.f.read(num_bytes)
        if isinstance(data, str):
            data = bytearray(data)

        # now read end of record
        num_bytes2 = self.get_int()
        if num_bytes2 != num_bytes:
            raise ValueError(
                "Fortran formatted record Mismatch"
                + " in starting and matching integers, "
                + str(self.num_bytes2)
                + " != "
                + str(self.num_bytes)
            )

        return _FortranRecord(data, num_bytes)
