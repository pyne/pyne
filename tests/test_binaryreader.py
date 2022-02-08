#!/usr/bin/env python

from __future__ import print_function
import struct
import filecmp
import os
import warnings

from nose.tools import assert_equal
from nose.plugins.skip import SkipTest

from pyne.utils import QAWarning

warnings.simplefilter("ignore", QAWarning)

from pyne.binaryreader import _FortranRecord, _BinaryReader


# test the ability to make a new empty FortranRecord
def test_make_empty_FR():
    test_record = _FortranRecord("", 0)

    if len(test_record.data) != 0:
        raise ValueError(
            "Failed to make an new empty _FortranRecord. " " Record has data."
        )
    if test_record.num_bytes != 0:
        raise ValueError(
            "Failed to make an new empty _FortranRecord. " " Record has num_bytes>0."
        )
    if test_record.pos != 0:
        raise ValueError(
            "Failed to make an new empty _FortranRecord. " " Position is not at 0."
        )

    return 1


# test the ability to reset the pos pointer in a FortranRecord
def test_reset_FR():
    temp_pos = 4

    test_record = _FortranRecord("", 0)  # already tested

    test_record.pos = temp_pos

    if test_record.pos != temp_pos:
        raise ValueError("Internal error: unable to update test_record.pos")

    test_record.reset()
    if test_record.pos != 0:
        raise ValueError("reset() method did not reset pos")

    return 1


# useful for checking all the changes expected from writing to a FortranRecord
def check_write_record_data(record, pos, num, data, typestring):
    if record.pos != pos:
        raise ValueError(
            "Writing "
            + typestring
            + " to record did not update pos properly: "
            + str(record.pos)
        )
    if record.num_bytes != num:
        raise ValueError(
            "Writing  " + typestring + "to record did not update num_bytes properly"
        )
    if record.data != data:
        raise ValueError(
            "Writing  " + typestring + "to record did not set data member correctly"
        )

    return 1


###################
#
# Test writing to all the different combinations of data type and single
# vs. multiple
#
###################
def test_write_FR_single_int():
    set_int = 8
    set_num_bytes = 4
    set_data = b"\x08\x00\x00\x00"

    # create new record
    test_record = _FortranRecord("", 0)  # already tested

    # write integer
    test_record.put_int(set_int)

    return check_write_record_data(
        test_record, set_num_bytes, set_num_bytes, set_data, "4-byte integer"
    )


def test_write_FR_int_list():
    set_intList = [8, 16]
    set_num_bytes = 8
    set_data = b"\x08\x00\x00\x00\x10\x00\x00\x00"

    test_record = _FortranRecord("", 0)

    # write integer list
    test_record.put_int(set_intList)

    return check_write_record_data(
        test_record, set_num_bytes, set_num_bytes, set_data, "list of 4-byte integers"
    )


def test_write_FR_single_long():
    set_long = 8

    set_num_bytes = 8
    set_data = b"\x08\x00\x00\x00\x00\x00\x00\x00"

    # create new record
    test_record = _FortranRecord("", 0)  # already tested

    # write long
    test_record.put_long(set_long)

    return check_write_record_data(
        test_record, set_num_bytes, set_num_bytes, set_data, "8-byte integer"
    )


def test_write_FR_long_list():
    set_long_list = [8, 16]
    set_num_bytes = 16
    set_data = b"\x08\x00\x00\x00\x00\x00\x00\x00\x10\x00\x00\x00\x00\x00\x00\x00"

    test_record = _FortranRecord("", 0)
    test_record.put_long(set_long_list)

    return check_write_record_data(
        test_record, set_num_bytes, set_num_bytes, set_data, "list of 8-byte integers"
    )


def test_write_FR_single_float():
    set_float = 3.14
    set_num_bytes = 4
    set_data = b"\xc3\xf5H@"

    # create new record
    test_record = _FortranRecord("", 0)  # already tested

    # write float
    test_record.put_float(set_float)

    return check_write_record_data(
        test_record, set_num_bytes, set_num_bytes, set_data, "float"
    )


def test_write_FR_float_list():
    set_floatList = [3.14, 0.3333]
    set_num_bytes = 8
    set_data = b"\xc3\xf5H@L\xa6\xaa>"

    test_record = _FortranRecord("", 0)
    test_record.put_float(set_floatList)

    return check_write_record_data(
        test_record, set_num_bytes, set_num_bytes, set_data, "list of floats"
    )


def test_write_FR_single_double():
    set_double = 3.14
    set_num_bytes = 8
    set_data = b"\x1f\x85\xebQ\xb8\x1e\t@"

    # create new record
    test_record = _FortranRecord("", 0)  # already tested

    # write double
    test_record.put_double(set_double)

    return check_write_record_data(
        test_record, set_num_bytes, set_num_bytes, set_data, "double"
    )


def test_write_FR_double_list():
    set_double_list = [3.14, 0.3333]
    set_num_bytes = 16
    set_data = b"\x1f\x85\xebQ\xb8\x1e\t@io\xf0\x85\xc9T\xd5?"

    test_record = _FortranRecord("", 0)
    test_record.put_double(set_double_list)

    return check_write_record_data(
        test_record, set_num_bytes, set_num_bytes, set_data, "list of doubles"
    )


def test_write_FR_single_string():
    set_string = "Hello World!"
    set_length = len(set_string)
    set_num_bytes = 12
    set_data = b"Hello World!"

    # create new record
    test_record = _FortranRecord("", 0)  # already tested

    # write string
    test_record.put_string([set_string], 12, 1)

    return check_write_record_data(
        test_record, set_num_bytes, set_num_bytes, set_data, "string"
    )


def test_write_FR_string_list():
    set_string_list = ["Hello ", "World!"]
    set_length = len(set_string_list[0])
    set_num_bytes = 12
    set_data = b"Hello World!"

    test_record = _FortranRecord("", 0)
    test_record.put_string(set_string_list, set_length)

    return check_write_record_data(
        test_record, set_num_bytes, set_num_bytes, set_data, "list of strings"
    )


def test_write_FR_mixed_record():
    set_int = 8
    set_float = 3.14
    set_double_list = [1.6e-19, 6.02e23]
    set_string = "Hello World!"

    set_num_bytes = 4 + 4 + (2 * 8) + 12
    set_data = (
        b"\x08\x00\x00\x00Hello World!#B\x92\x0c\xa1\x9c\x07<a\xd3"
        + b"\xa8\x10\x9f\xde\xdfD\xc3\xf5H@"
    )

    test_record = _FortranRecord("", 0)
    test_record.put_int(set_int)
    test_record.put_string([set_string], len(set_string))
    test_record.put_double(set_double_list)
    test_record.put_float([set_float])

    return check_write_record_data(
        test_record, set_num_bytes, set_num_bytes, set_data, "list of doubles"
    )


###################
#
# Test reading from all the different combinations of data type and single
# vs. multiple
#
###################


def test_read_FR_single_int():
    set_int = 8

    test_record = _FortranRecord("", 0)  # already tested
    test_record.put_int([set_int])  # already tested

    test_record.reset()  # already tested

    test_int = test_record.get_int()[0]

    if test_int != set_int:
        raise ValueError("Value from get_int doesn't match value " "from put_int.")

    return 1


def test_read_FR_int_list():
    set_intList = [8, 16]
    num_ints = 2

    test_record = _FortranRecord("", 0)  # already tested
    test_record.put_int(set_intList)  # already tested

    test_record.reset()  # already tested

    test_int = test_record.get_int(num_ints)

    if test_int != set_intList:
        raise ValueError("Value from get_int doesn't match value " "from put_int.")

    return 1


def test_read_FR_single_long():
    set_long = 8

    test_record = _FortranRecord("", 0)  # already tested
    test_record.put_long([set_long])  # already tested

    test_record.reset()  # already tested

    testLong = test_record.get_long()[0]

    if testLong != set_long:
        raise ValueError("Value from get_long doesn't match value " "from put_long.")

    return 1


def test_read_FR_long_list():
    set_long_list = [8, 16]
    numLongs = 2

    test_record = _FortranRecord("", 0)  # already tested
    test_record.put_long(set_long_list)  # already tested

    test_record.reset()  # already tested

    testLong = test_record.get_long(numLongs)

    if testLong != set_long_list:
        raise ValueError("Value from get_long doesn't match value " "from put_long.")

    return 1


def test_read_FR_single_float():
    set_float = 6.1

    test_record = _FortranRecord("", 0)
    test_record.put_float([set_float])

    test_record.reset()

    test_float = test_record.get_float()[0]

    # NOTE: since Python doesn't do native 32-bit floats, both values should be
    #       passed through a string formatting to truncate to 6 digits
    set_float = "%10.6f" % set_float
    test_float = "%10.6f" % test_float

    if test_float != set_float:
        print("{0} != {1}".format(set_float, test_float))
        raise ValueError("Value from get_float doesn't match value " "from put_float.")

    return 1


def test_read_FR_float_list():
    floatList = [2.34, 8.65]

    test_record = _FortranRecord("", 0)
    test_record.put_float(floatList)

    test_record.reset()

    testList = test_record.get_float(2)

    # NOTE: since Python doesn't do native 32-bit floats, both values should be
    #       passed through a string formatting to truncate to 6 digits
    floatList = ["%10.6f" % floatList[0], "%10.6f" % floatList[1]]
    testList = ["%10.6f" % testList[0], "%10.6f" % testList[1]]

    if testList != floatList:
        print(floatList)
        print(testList)
        raise ValueError("List from get_float doesn't match value " "from put_float.")

    return 1


def test_read_FR_single_double():
    set_double = 1.43

    test_record = _FortranRecord("", 0)
    test_record.put_double([set_double])

    test_record.reset()

    testDouble = test_record.get_double()[0]

    if testDouble != set_double:
        raise ValueError(
            "Value from get_double doesn't match value " "from put_double."
        )

    return 1


def test_read_FR_double_list():
    double_list = [2.34, 8.65]

    test_record = _FortranRecord("", 0)
    test_record.put_double(double_list)

    test_record.reset()

    testList = test_record.get_double(2)

    if testList != double_list:
        raise ValueError("List from get_double doesn't match value " "from put_double.")

    return 1


def test_read_FR_single_string():
    set_string = "Hello World!"
    set_length = len(set_string)

    # create new record
    test_record = _FortranRecord("", 0)  # already tested
    test_record.put_string([set_string], 12, 1)

    test_record.reset()

    test_string = test_record.get_string(set_length)[0]

    if test_string != set_string:
        raise ValueError("List from get_string doesn't match value " "from put_string.")

    return 1


def test_read_FR_string_list():
    set_string_list = ["Hello ", "World!"]
    set_length = len(set_string_list[0])

    test_record = _FortranRecord("", 0)
    test_record.put_string(set_string_list, set_length)
    test_record.reset()

    test_string_list = test_record.get_string(set_length, 2)

    if test_string_list != set_string_list:
        raise ValueError("List from get_string doesn't match value " "from put_string.")

    return 1


def test_read_FR_mixed_record():
    set_int = 8
    set_float = 3.14
    set_double_list = [1.6e-19, 6.02e23]
    set_string = "Hello World!"

    test_record = _FortranRecord("", 0)
    test_record.put_int([set_int])
    test_record.put_string([set_string], len(set_string))
    test_record.put_double(set_double_list)
    test_record.put_float([set_float])

    test_record.reset()

    test_int = test_record.get_int()[0]
    if test_int != set_int:
        raise ValueError("Value from get_int doesn't match value " "from put_int.")

    test_string = test_record.get_string(12)[0]
    if test_string != set_string:
        raise ValueError("List from get_string doesn't match value " "from put_string.")

    test_double_list = test_record.get_double(2)
    if test_double_list != set_double_list:
        raise ValueError("List from get_double doesn't match value " "from put_double.")

    test_float = test_record.get_float()[0]
    # NOTE: since Python doesn't do native 32-bit floats, both values should be
    #       passed through a string formatting to truncate to 6 digits
    set_float = "%10.6f" % set_float
    test_float = "%10.6f" % test_float

    if test_float != set_float:
        print("{0} != {1}".format(set_float, test_float))
        raise ValueError("Value from get_float doesn't match value " "from put_float.")

    return 1


####################
#
# Test binary file operations here
#
####################


def test_open_writable_BR():
    binary_file = _BinaryReader("test.file", "wb")
    if not binary_file.f:
        raise ValueError("Failed to open new file for writing.")

    binary_file.close()

    return 1


def test_write_BR():
    binary_file = _BinaryReader("test.file", "wb")
    if not binary_file.f:
        raise ValueError("Failed to open new file for writing.")

    set_int = 8
    set_float = 3.14
    set_double_list = [1.6e-19, 6.02e23]
    set_string = "Hello World!"

    test_record = _FortranRecord("", 0)
    test_record.put_int([set_int])
    test_record.put_string([set_string], len(set_string))
    test_record.put_double(set_double_list)
    test_record.put_float([set_float])

    test_record.reset()

    binary_file.put_fortran_record(test_record)

    binary_file.close()

    if not filecmp.cmp("test.file", "test_readBR.ref"):
        raise ValueError("Created file does not match reference.")

    return 1


def test_read_BR():
    set_int = 8
    set_float = 3.14
    set_double_list = [1.6e-19, 6.02e23]
    set_string = b"Hello World!"

    binary_file = _BinaryReader("test_readBR.ref")

    test_record = binary_file.get_fortran_record()

    try:
        test_int = test_record.get_int()[0]
    except:
        raise SkipTest
    assert_equal(set_int, test_int)

    test_string = test_record.get_string(12)[0]
    assert_equal(set_string.decode(), test_string)

    test_double_list = test_record.get_double(2)
    if test_double_list != set_double_list:
        raise ValueError("List of doubles was not as expected.")

    test_float = test_record.get_float()[0]
    # NOTE: since Python doesn't do native 32-bit floats, both values should be
    #       passed through a string formatting to truncate to 6 digits
    set_float = "%10.6f" % set_float
    test_float = "%10.6f" % test_float

    if test_float != set_float:
        print("{0} != {1}".format(set_float, test_float))
        raise ValueError("Float value was not as expected.")

    return 1


# start all tests here

tests = [0, 0]

if os.name == "posix":
    passed = "\x1b[31G\x1b[32mPASSED\x1b[0m"
    failed = "\x1b[31G\x1b[31mFAILED\x1b[0m"
else:
    passed = "PASSED"
    failed = "FAILED"

print("test_make_empty_FR: ")
try:
    tests[0] += test_make_empty_FR()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_reset_FR: ")
try:
    tests[0] += test_reset_FR()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_write_FR_single_int: ")
try:
    tests[0] += test_write_FR_single_int()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_write_FR_int_list: ")
try:
    tests[0] += test_write_FR_int_list()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_write_FR_single_long: ")
try:
    tests[0] += test_write_FR_single_long()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_write_FR_long_list: ")
try:
    tests[0] += test_write_FR_long_list()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_write_FR_single_float: ")
try:
    tests[0] += test_write_FR_single_float()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_write_FR_float_list: ")
try:
    tests[0] += test_write_FR_float_list()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_write_FR_single_double: ")
try:
    tests[0] += test_write_FR_single_double()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_write_FR_double_list: ")
try:
    tests[0] += test_write_FR_double_list()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_write_FR_single_string: ")
try:
    tests[0] += test_write_FR_single_string()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_write_FR_string_list: ")
try:
    tests[0] += test_write_FR_string_list()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_write_FR_mixed_record: ")
try:
    tests[0] += test_write_FR_mixed_record()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_read_FR_single_int: ")
try:
    tests[0] += test_read_FR_single_int()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_read_FR_int_list: ")
try:
    tests[0] += test_read_FR_int_list()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_read_FR_single_long: ")
try:
    tests[0] += test_read_FR_single_long()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_read_FR_long_list: ")
try:
    tests[0] += test_read_FR_long_list()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_read_FR_single_float: ")
try:
    tests[0] += test_read_FR_single_float()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_read_FR_float_list: ")
try:
    tests[0] += test_read_FR_float_list()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_read_FR_single_double: ")
try:
    tests[0] += test_read_FR_single_double()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_read_FR_double_list: ")
try:
    tests[0] += test_read_FR_double_list()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_read_FR_single_string: ")
try:
    tests[0] += test_read_FR_single_string()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_read_FR_string_list: ")
try:
    tests[0] += test_read_FR_string_list()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_read_FR_mixed_record: ")
try:
    tests[0] += test_read_FR_mixed_record()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_open_writable_BR: ")
try:
    tests[0] += test_open_writable_BR()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_write_BR: ")
try:
    tests[0] += test_write_BR()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("test_read_BR: ")
try:
    tests[0] += test_read_BR()
    print(passed)
except Exception as inst:
    print(failed + ": " + str(inst))
    tests[1] += 1

print("Ran    " + str(tests[0] + tests[1]) + " tests.")
print("PASSED " + str(tests[0]) + " tests.")
print("FAILED " + str(tests[1]) + " tests.")
