#!/usr/bin/env python

import struct
import filecmp
import os

from pyne.binaryreader import (_FortranRecord, FortranRecordError,
                               _BinaryReader, BinaryReaderError)

# test the ability to make a new empty FortranRecord
def test_make_empty_FR():
    testRecord = _FortranRecord('',0)

    if testRecord.data != '':
        raise FortranRecordError("Failed to make an new empty _FortranRecord.  "
                                 "Record has data.")
    if testRecord.numBytes != 0:
        raise FortranRecordError("Failed to make an new empty _FortranRecord.  "
                                 "Record has numBytes>0.")
    if testRecord.pos != 0:
        raise FortranRecordError("Failed to make an new empty _FortranRecord.  "
                                 "Position is not at 0.")

    return 1

# test the ability to reset the pos pointer in a FortranRecord
def test_reset_FR():

    tempPos = 4

    testRecord = _FortranRecord('',0)  # already tested

    testRecord.pos = tempPos

    if testRecord.pos != tempPos:
        raise FortranRecordError("Internal error: unable to update "
                                 "testRecord.pos")

    testRecord.reset()
    if testRecord.pos != 0:
        raise FortranRecordError("reset() method did not reset pos")
        
    return 1
    

# useful for checking all the changes expected from writing to a FortranRecord
def check_write_record_data(record, pos, num, data, typeString):

    if record.pos != pos:
        raise FortranRecordError("Writing " + typeString + 
                                 " to record did not update pos properly: "
                                 + str(record.pos))
    if record.numBytes != num:
        raise FortranRecordError("Writing " + typeString +
                                 " to record did not update numBytes properly")
    if record.data != data:
        raise FortranRecordError("Writing " + typeString +
                                 " to record did not set data member correctly")
    
    return 1

###################
#
# Test writing to all the different combinations of data type and single
# vs. multiple
#
###################
def test_write_FR_single_int():

    setInt = 8
    setNumBytes = 4
    setData = '\x08\x00\x00\x00'
    
    # create new record
    testRecord = _FortranRecord('',0)  # already tested
    
    # write integer
    testRecord.put_int(setInt)

    return check_write_record_data(testRecord, setNumBytes, setNumBytes,
                                   setData,"4-byte integer")

def test_write_FR_int_list():

    setIntList = [8,16]
    setNumBytes = 8
    setData = '\x08\x00\x00\x00\x10\x00\x00\x00'

    testRecord = _FortranRecord('',0)

    # write integer list
    testRecord.put_int(setIntList)
    
    return check_write_record_data(testRecord, setNumBytes, setNumBytes,
                                   setData, "list of 4-byte integers")

def test_write_FR_single_long():

    setLong = 8
    setNumBytes = 8
    setData = '\x08\x00\x00\x00\x00\x00\x00\x00'
    
    # create new record
    testRecord = _FortranRecord('',0)  # already tested
    
    # write long
    testRecord.put_long(setLong)

    return check_write_record_data(testRecord, setNumBytes, setNumBytes,
                                   setData,"8-byte integer")

def test_write_FR_long_list():

    setLongList = [8,16]
    setNumBytes  = 16
    setData = '\x08\x00\x00\x00\x00\x00\x00\x00\x10\x00\x00\x00\x00\x00\x00\x00'

    testRecord = _FortranRecord('',0)
    testRecord.put_long(setLongList)
    
    return check_write_record_data(testRecord, setNumBytes, setNumBytes,
                                   setData,"list of 8-byte integers")

def test_write_FR_single_float():

    setFloat = 3.14
    setNumBytes = 4
    setData = '\xc3\xf5H@'
    
    # create new record
    testRecord = _FortranRecord('',0)  # already tested
    
    # write float
    testRecord.put_float(setFloat)

    return check_write_record_data(testRecord, setNumBytes, setNumBytes,
                                   setData,"float")

def test_write_FR_float_list():

    setFloatList = [3.14,0.3333]
    setNumBytes  = 8
    setData = '\xc3\xf5H@L\xa6\xaa>'

    testRecord = _FortranRecord('',0)
    testRecord.put_float(setFloatList)
    
    return check_write_record_data(testRecord, setNumBytes, setNumBytes,
                                   setData,"list of floats")

def test_write_FR_single_double():

    setDouble = 3.14
    setNumBytes = 8
    setData = '\x1f\x85\xebQ\xb8\x1e\t@'
    
    # create new record
    testRecord = _FortranRecord('',0)  # already tested
    
    # write double
    testRecord.put_double(setDouble)

    return check_write_record_data(testRecord, setNumBytes, setNumBytes,
                                   setData,"double")

def test_write_FR_double_list():

    setDoubleList = [3.14,0.3333]
    setNumBytes  = 16
    setData = '\x1f\x85\xebQ\xb8\x1e\t@io\xf0\x85\xc9T\xd5?'

    testRecord = _FortranRecord('',0)
    testRecord.put_double(setDoubleList)
    
    return check_write_record_data(testRecord, setNumBytes, setNumBytes,
                                   setData,"list of doubles")

def test_write_FR_single_string():

    setString = "Hello World!"
    setLength = len(setString)
    setNumBytes = 12
    setData = 'Hello World!'
    
    # create new record
    testRecord = _FortranRecord('',0)  # already tested
    
    # write string
    testRecord.put_string([setString],12,1)
    
    return check_write_record_data(testRecord, setNumBytes, setNumBytes,
                                   setData,"string")

def test_write_FR_string_list():

    setStringList = ["Hello ", "World!"]
    setLength = len(setStringList[0])
    setNumBytes  = 12
    setData = 'Hello World!'

    testRecord = _FortranRecord('',0)
    testRecord.put_string(setStringList,setLength)
    
    return check_write_record_data(testRecord, setNumBytes, setNumBytes,
                                   setData,"list of strings")

def test_write_FR_mixed_record():

    setInt = 8
    setFloat = 3.14
    setDoubleList = [1.6e-19,6.02e23]
    setString = "Hello World!"
    
    setNumBytes = 4+4+(2*8)+12
    setData = '\x08\x00\x00\x00Hello World!#B\x92\x0c\xa1\x9c\x07<a\xd3' + \
        '\xa8\x10\x9f\xde\xdfD\xc3\xf5H@'

    testRecord = _FortranRecord('',0)
    testRecord.put_int(setInt)
    testRecord.put_string([setString],len(setString))
    testRecord.put_double(setDoubleList)
    testRecord.put_float([setFloat])

    return check_write_record_data(testRecord, setNumBytes, setNumBytes,
                                   setData,"list of doubles")

###################
#
# Test reading from all the different combinations of data type and single
# vs. multiple
#
###################

def test_read_FR_single_int():

    setInt = 8

    testRecord = _FortranRecord('',0)  # already tested
    testRecord.put_int([setInt])       # already tested

    testRecord.reset()                 # already tested

    testInt = testRecord.get_int()[0]

    if testInt != setInt:
        raise FortranRecordError("Value from get_int doesn't match value "
                                 "from put_int.")
        
    return 1

def test_read_FR_int_list():

    setIntList = [8,16]
    numInts = 2

    testRecord = _FortranRecord('',0)  # already tested
    testRecord.put_int(setIntList)       # already tested

    testRecord.reset()                 # already tested

    testInt = testRecord.get_int(numInts)

    if testInt != setIntList:
        raise FortranRecordError("Value from get_int doesn't match value "
                                 "from put_int.")
        
    return 1

def test_read_FR_single_long():

    setLong = 8

    testRecord = _FortranRecord('',0)  # already tested
    testRecord.put_long([setLong])       # already tested

    testRecord.reset()                 # already tested

    testLong = testRecord.get_long()[0]

    if testLong != setLong:
        raise FortranRecordError("Value from get_long doesn't match value "
                                 "from put_long.")
        
    return 1

def test_read_FR_long_list():

    setLongList = [8,16]
    numLongs = 2

    testRecord = _FortranRecord('',0)  # already tested
    testRecord.put_long(setLongList)       # already tested

    testRecord.reset()                 # already tested

    testLong = testRecord.get_long(numLongs)

    if testLong != setLongList:
        raise FortranRecordError("Value from get_long doesn't match value "
                                 "from put_long.")
        
    return 1

def test_read_FR_single_float():

    setFloat =  6.1

    testRecord = _FortranRecord('',0)
    testRecord.put_float([setFloat])

    testRecord.reset()
    
    testFloat = testRecord.get_float()[0]

    # NOTE: since Python doesn't do native 32-bit floats, both values should be
    #       passed through a string formatting to truncate to 6 digits
    setFloat  = '%10.6f' % setFloat
    testFloat = '%10.6f' % testFloat

    if testFloat != setFloat:
        print str(setFloat) + " != " + str(testFloat)
        raise FortranRecordError("Value from get_float doesn't match value "
                                 "from put_float.")
        
    return 1

def test_read_FR_float_list():

    floatList = [2.34,8.65]
    
    testRecord = _FortranRecord('',0)
    testRecord.put_float(floatList)
    
    testRecord.reset()

    testList = testRecord.get_float(2)

    # NOTE: since Python doesn't do native 32-bit floats, both values should be
    #       passed through a string formatting to truncate to 6 digits
    floatList  = ['%10.6f' % floatList[0], '%10.6f' % floatList[1]]
    testList = ['%10.6f' % testList[0], '%10.6f' % testList[1]]

    if testList != floatList:
        print floatList
        print testList
        raise FortranRecordError("List from get_float doesn't match value "
                                 "from put_float.")

    return 1

def test_read_FR_single_double():

    setDouble = 1.43

    testRecord = _FortranRecord('',0)
    testRecord.put_double([setDouble])

    testRecord.reset()
    
    testDouble = testRecord.get_double()[0]

    if testDouble != setDouble:
        raise FortranRecordError("Value from get_double doesn't match value "
                                 "from put_double.")
        
    return 1

def test_read_FR_double_list():

    doubleList = [2.34,8.65]
    
    testRecord = _FortranRecord('',0)
    testRecord.put_double(doubleList)
    
    testRecord.reset()

    testList = testRecord.get_double(2)

    if testList != doubleList:
        raise FortranRecordError("List from get_double doesn't match value "
                                 "from put_double.")

    return 1

def test_read_FR_single_string():

    setString = "Hello World!"
    setLength = len(setString)
    
    # create new record
    testRecord = _FortranRecord('',0)  # already tested
    testRecord.put_string([setString],12,1)

    testRecord.reset()

    testString = testRecord.get_string(setLength)[0]

    if testString != setString:
        raise FortranRecordError("List from get_string doesn't match value "
                                 "from put_string.")
        
    return 1

def test_read_FR_string_list():

    setStringList = ["Hello ", "World!"]
    setLength = len(setStringList[0])

    testRecord = _FortranRecord('',0)
    testRecord.put_string(setStringList,setLength)
    testRecord.reset()

    testStringList = testRecord.get_string(setLength,2)

    if testStringList != setStringList:
        raise FortranRecordError("List from get_string doesn't match value "
                                 "from put_string.")
        
    return 1

def test_read_FR_mixed_record():

    setInt = 8
    setFloat = 3.14
    setDoubleList = [1.6e-19,6.02e23]
    setString = "Hello World!"
    
    testRecord = _FortranRecord('',0)
    testRecord.put_int([setInt])
    testRecord.put_string([setString],len(setString))
    testRecord.put_double(setDoubleList)
    testRecord.put_float([setFloat])

    testRecord.reset()

    testInt = testRecord.get_int()[0]
    if testInt != setInt:
        raise FortranRecordError("Value from get_int doesn't match value "
                                 "from put_int.")

    testString = testRecord.get_string(12)[0]
    if testString != setString:
        raise FortranRecordError("List from get_string doesn't match value "
                                 "from put_string.")

    testDoubleList = testRecord.get_double(2)
    if testDoubleList != setDoubleList:
        raise FortranRecordError("List from get_double doesn't match value "
                                 "from put_double.")

    testFloat = testRecord.get_float()[0]
    # NOTE: since Python doesn't do native 32-bit floats, both values should be
    #       passed through a string formatting to truncate to 6 digits
    setFloat  = '%10.6f' % setFloat
    testFloat = '%10.6f' % testFloat

    if testFloat != setFloat:
        print str(setFloat) + " != " + str(testFloat)
        raise FortranRecordError("Value from get_float doesn't match value "
                                 "from put_float.")
    
    return 1


####################
#
# Test binary file operations here
#
####################

def test_open_writable_BR():
    
    binaryFile = _BinaryReader('test.file','wb')
    if not binaryFile.f:
        raise BinaryReaderError("Failed to open new file for writing.")
    
    binaryFile.close()
    
    return 1

def test_write_BR():

    binaryFile = _BinaryReader('test.file','wb')
    if not binaryFile.f:
        raise BinaryReaderError("Failed to open new file for writing.")

    setInt = 8
    setFloat = 3.14
    setDoubleList = [1.6e-19,6.02e23]
    setString = "Hello World!"
    
    testRecord = _FortranRecord('',0)
    testRecord.put_int([setInt])
    testRecord.put_string([setString],len(setString))
    testRecord.put_double(setDoubleList)
    testRecord.put_float([setFloat])

    testRecord.reset()

    binaryFile.put_fortran_record(testRecord)

    binaryFile.close()

    if not filecmp.cmp('test.file','test_readBR.ref'):
        raise BinaryReaderError('Created file does not match reference.')

    return 1

def test_read_BR():

    setInt = 8
    setFloat = 3.14
    setDoubleList = [1.6e-19,6.02e23]
    setString = "Hello World!"
    
    binaryFile = _BinaryReader('test_readBR.ref')
    
    testRecord = binaryFile.get_fortran_record()

    testInt = testRecord.get_int()[0]
    if testInt != setInt:
        raise BinaryReaderError("Integer value was not as expected.")

    testString = testRecord.get_string(12)[0]
    if testString != setString:
        raise BinaryReaderError("String was not as expected.")

    testDoubleList = testRecord.get_double(2)
    if testDoubleList != setDoubleList:
        raise BinaryReaderError("List of doubles was not as expected.")

    testFloat = testRecord.get_float()[0]
    # NOTE: since Python doesn't do native 32-bit floats, both values should be
    #       passed through a string formatting to truncate to 6 digits
    setFloat  = '%10.6f' % setFloat
    testFloat = '%10.6f' % testFloat

    if testFloat != setFloat:
        print str(setFloat) + " != " + str(testFloat)
        raise BinaryReaderError("Float value was not as expected.")
    
    return 1
    

# start all tests here

tests = [0,0]

if os.name == 'posix':
    passed = '\x1b[31G\x1b[32mPASSED\x1b[0m'
    failed = '\x1b[31G\x1b[31mFAILED\x1b[0m'
else:
    passed = 'PASSED'
    failed = 'FAILED'

print "test_make_empty_FR: ",
try:
    tests[0] += test_make_empty_FR()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_reset_FR: ",
try:
    tests[0] += test_reset_FR()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_write_FR_single_int: ",
try:
    tests[0] += test_write_FR_single_int()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_write_FR_int_list: ",
try:
    tests[0] += test_write_FR_int_list()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_write_FR_single_long: ",
try:
    tests[0] += test_write_FR_single_long()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_write_FR_long_list: ",
try:
    tests[0] += test_write_FR_long_list()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_write_FR_single_float: ",
try:
    tests[0] += test_write_FR_single_float()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_write_FR_float_list: ",
try:
    tests[0] += test_write_FR_float_list()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_write_FR_single_double: ",
try:
    tests[0] += test_write_FR_single_double()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_write_FR_double_list: ",
try:
    tests[0] += test_write_FR_double_list()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_write_FR_single_string: ",
try:
    tests[0] += test_write_FR_single_string()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_write_FR_string_list: ",
try:
    tests[0] += test_write_FR_string_list()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_write_FR_mixed_record: ",
try:
    tests[0] += test_write_FR_mixed_record()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_read_FR_single_int: ",
try:
    tests[0] += test_read_FR_single_int()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_read_FR_int_list: ",
try:
    tests[0] += test_read_FR_int_list()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_read_FR_single_long: ",
try:
    tests[0] += test_read_FR_single_long()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_read_FR_long_list: ",
try:
    tests[0] += test_read_FR_long_list()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_read_FR_single_float: ",
try:
    tests[0] += test_read_FR_single_float()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_read_FR_float_list: ",
try:
    tests[0] += test_read_FR_float_list()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_read_FR_single_double: ",
try:
    tests[0] += test_read_FR_single_double()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_read_FR_double_list: ",
try:
    tests[0] += test_read_FR_double_list()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_read_FR_single_string: ",
try:
    tests[0] += test_read_FR_single_string()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_read_FR_string_list: ",
try:
    tests[0] += test_read_FR_string_list()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_read_FR_mixed_record: ",
try:
    tests[0] += test_read_FR_mixed_record()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_open_writable_BR: ",
try:
    tests[0] += test_open_writable_BR()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_write_BR: ",
try:
    tests[0] += test_write_BR()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1

print "test_read_BR: ",
try:
    tests[0] += test_read_BR()
    print passed
except Exception as inst:
    print failed + ": " + str(inst)
    tests[1] += 1



print "Ran    " + str(tests[0]+tests[1]) + " tests."
print "PASSED " + str(tests[0]) + " tests."
print "FAILED " + str(tests[1]) + " tests."
