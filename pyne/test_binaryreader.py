#!/usr/bin/env python

import struct
import filecmp
from binaryreader import _FortranRecord, FortranRecordError, _BinaryReader, BinaryReaderError

# test the ability to make a new empty FortranRecord
def test_makeEmptyFR():
    testRecord = _FortranRecord('',0)

    if testRecord.data != '':
        raise FortranRecordError("Failed to make an new empty _FortranRecord.  Record has data.")
    if testRecord.numBytes != 0:
        raise FortranRecordError("Failed to make an new empty _FortranRecord.  Record has numBytes>0.")
    if testRecord.pos != 0:
        raise FortranRecordError("Failed to make an new empty _FortranRecord.  Position is not at 0.")

    return 1

# test the ability to reset the pos pointer in a FortranRecord
def test_resetFR():

    tempPos = 4

    testRecord = _FortranRecord('',0)  # already tested

    testRecord.pos = tempPos

    if testRecord.pos != tempPos:
        raise FortranRecordError("Internal error: unable to update testRecord.pos")

    testRecord.reset()
    if testRecord.pos != 0:
        raise FortranRecordError("reset() method did not reset pos")
        
    return 1
    

# useful for checking all the changes expected from writing to a FortranRecord
def checkWriteRecordData(record,pos,num,data,typeString):

    if record.pos != pos:
        raise FortranRecordError("Writing " + typeString + " to record did not update pos properly: " + str(record.pos))
    if record.numBytes != num:
        raise FortranRecordError("Writing " + typeString + " to record did not update numBytes properly")
    if record.data != data:
        raise FortranRecordError("Writing " + typeString + " to record did not set data member correctly")
    
    return 1

###################
#
# Test writing to all the different combinations of data type and single vs. multiple
#
###################
def test_writeFR_singleInt():

    setInt = 8
    setNumBytes = 4
    setData = '\x08\x00\x00\x00'
    
    # create new record
    testRecord = _FortranRecord('',0)  # already tested
    
    # write integer
    testRecord.put_int([setInt])

    return checkWriteRecordData(testRecord,setNumBytes,setNumBytes,setData,"4-byte integer")

def test_writeFR_intList():

    setIntList = [8,16]
    setNumBytes = 8
    setData = '\x08\x00\x00\x00\x10\x00\x00\x00'

    testRecord = _FortranRecord('',0)

    # write integer list
    testRecord.put_int(setIntList)
    
    return checkWriteRecordData(testRecord,setNumBytes,setNumBytes,setData,"list of 4-byte integers")

def test_writeFR_singleLong():

    setLong = 8
    setNumBytes = 8
    setData = '\x08\x00\x00\x00\x00\x00\x00\x00'
    
    # create new record
    testRecord = _FortranRecord('',0)  # already tested
    
    # write long
    testRecord.put_long([setLong])

    return checkWriteRecordData(testRecord,setNumBytes,setNumBytes,setData,"8-byte integer")

def test_writeFR_longList():

    setLongList = [8,16]
    setNumBytes  = 16
    setData = '\x08\x00\x00\x00\x00\x00\x00\x00\x10\x00\x00\x00\x00\x00\x00\x00'

    testRecord = _FortranRecord('',0)
    testRecord.put_long(setLongList)
    
    return checkWriteRecordData(testRecord,setNumBytes,setNumBytes,setData,"list of 8-byte integers")

def test_writeFR_singleFloat():

    setFloat = 3.14
    setNumBytes = 4
    setData = '\xc3\xf5H@'
    
    # create new record
    testRecord = _FortranRecord('',0)  # already tested
    
    # write float
    testRecord.put_float([setFloat])

    return checkWriteRecordData(testRecord,setNumBytes,setNumBytes,setData,"float")

def test_writeFR_floatList():

    setFloatList = [3.14,0.3333]
    setNumBytes  = 8
    setData = '\xc3\xf5H@L\xa6\xaa>'

    testRecord = _FortranRecord('',0)
    testRecord.put_float(setFloatList)
    
    return checkWriteRecordData(testRecord,setNumBytes,setNumBytes,setData,"list of floats")

def test_writeFR_singleDouble():

    setDouble = 3.14
    setNumBytes = 8
    setData = '\x1f\x85\xebQ\xb8\x1e\t@'
    
    # create new record
    testRecord = _FortranRecord('',0)  # already tested
    
    # write double
    testRecord.put_double([setDouble])

    return checkWriteRecordData(testRecord,setNumBytes,setNumBytes,setData,"double")

def test_writeFR_doubleList():

    setDoubleList = [3.14,0.3333]
    setNumBytes  = 16
    setData = '\x1f\x85\xebQ\xb8\x1e\t@io\xf0\x85\xc9T\xd5?'

    testRecord = _FortranRecord('',0)
    testRecord.put_double(setDoubleList)
    
    return checkWriteRecordData(testRecord,setNumBytes,setNumBytes,setData,"list of doubles")

def test_writeFR_singleString():

    setString = "Hello World!"
    setLength = len(setString)
    setNumBytes = 12
    setData = 'Hello World!'
    
    # create new record
    testRecord = _FortranRecord('',0)  # already tested
    
    # write string
    testRecord.put_string([setString],12,1)
    
    return checkWriteRecordData(testRecord,setNumBytes,setNumBytes,setData,"string")

def test_writeFR_stringList():

    setStringList = ["Hello ", "World!"]
    setLength = len(setStringList[0])
    setNumBytes  = 12
    setData = 'Hello World!'

    testRecord = _FortranRecord('',0)
    testRecord.put_string(setStringList,setLength)
    
    return checkWriteRecordData(testRecord,setNumBytes,setNumBytes,setData,"list of strings")

def test_writeFR_mixedRecord():

    setInt = 8
    setFloat = 3.14
    setDoubleList = [1.6e-19,6.02e23]
    setString = "Hello World!"
    
    setNumBytes = 4+4+(2*8)+12
    setData = '\x08\x00\x00\x00Hello World!#B\x92\x0c\xa1\x9c\x07<a\xd3\xa8\x10\x9f\xde\xdfD\xc3\xf5H@'

    testRecord = _FortranRecord('',0)
    testRecord.put_int([setInt])
    testRecord.put_string([setString],len(setString))
    testRecord.put_double(setDoubleList)
    testRecord.put_float([setFloat])

    return checkWriteRecordData(testRecord,setNumBytes,setNumBytes,setData,"list of doubles")

###################
#
# Test reading from all the different combinations of data type and single vs. multiple
#
###################

def test_readFR_singleInt():

    setInt = 8

    testRecord = _FortranRecord('',0)  # already tested
    testRecord.put_int([setInt])       # already tested

    testRecord.reset()                 # already tested

    testInt = testRecord.get_int()

    if testInt != setInt:
        raise FortranRecordError("Value from get_int doesn't match value from put_int.")
        
    return 1

def test_readFR_intList():

    setIntList = [8,16]
    numInts = 2

    testRecord = _FortranRecord('',0)  # already tested
    testRecord.put_int(setIntList)       # already tested

    testRecord.reset()                 # already tested

    testInt = testRecord.get_int(numInts)

    if testInt != setIntList:
        raise FortranRecordError("Value from get_int doesn't match value from put_int.")
        
    return 1

def test_readFR_singleLong():

    setLong = 8

    testRecord = _FortranRecord('',0)  # already tested
    testRecord.put_long([setLong])       # already tested

    testRecord.reset()                 # already tested

    testLong = testRecord.get_long()

    if testLong != setLong:
        raise FortranRecordError("Value from get_long doesn't match value from put_long.")
        
    return 1

def test_readFR_longList():

    setLongList = [8,16]
    numLongs = 2

    testRecord = _FortranRecord('',0)  # already tested
    testRecord.put_long(setLongList)       # already tested

    testRecord.reset()                 # already tested

    testLong = testRecord.get_long(numLongs)

    if testLong != setLongList:
        raise FortranRecordError("Value from get_long doesn't match value from put_long.")
        
    return 1

def test_readFR_singleFloat():

    setFloat =  6.1

    testRecord = _FortranRecord('',0)
    testRecord.put_float([setFloat])

    testRecord.reset()
    
    testFloat = testRecord.get_float()

    # NOTE: since Python doesn't do native 32-bit floats, both values should be passed through
    #       a string formatting to truncate to 6 digits
    setFloat  = '%10.6f' % setFloat
    testFloat = '%10.6f' % testFloat

    if testFloat != setFloat:
        print str(setFloat) + " != " + str(testFloat)
        raise FortranRecordError("Value from get_float doesn't match value from put_float.")
        
    return 1

def test_readFR_floatList():

    floatList = [2.34,8.65]
    
    testRecord = _FortranRecord('',0)
    testRecord.put_float(floatList)
    
    testRecord.reset()

    testList = testRecord.get_float(2)

    # NOTE: since Python doesn't do native 32-bit floats, both values should be passed through
    #       a string formatting to truncate to 6 digits
    floatList  = ['%10.6f' % floatList[0], '%10.6f' % floatList[1]]
    testList = ['%10.6f' % testList[0], '%10.6f' % testList[1]]

    if testList != floatList:
        print floatList
        print testList
        raise FortranRecordError("List from get_float doesn't match value from put_float.")

    return 1

def test_readFR_singleDouble():

    setDouble = 1.43

    testRecord = _FortranRecord('',0)
    testRecord.put_double([setDouble])

    testRecord.reset()
    
    testDouble = testRecord.get_double()

    if testDouble != setDouble:
        raise FortranRecordError("Value from get_double doesn't match value from put_double.")
        
    return 1

def test_readFR_doubleList():

    doubleList = [2.34,8.65]
    
    testRecord = _FortranRecord('',0)
    testRecord.put_double(doubleList)
    
    testRecord.reset()

    testList = testRecord.get_double(2)

    if testList != doubleList:
        raise FortranRecordError("List from get_double doesn't match value from put_double.")

    return 1

def test_readFR_singleString():

    setString = "Hello World!"
    setLength = len(setString)
    
    # create new record
    testRecord = _FortranRecord('',0)  # already tested
    testRecord.put_string([setString],12,1)

    testRecord.reset()

    testString = testRecord.get_string(setLength)

    if testString != setString:
        raise FortranRecordError("List from get_string doesn't match value from put_string.")
        
    return 1

def test_readFR_stringList():

    setStringList = ["Hello ", "World!"]
    setLength = len(setStringList[0])

    testRecord = _FortranRecord('',0)
    testRecord.put_string(setStringList,setLength)
    testRecord.reset()

    testStringList = testRecord.get_string(setLength,2)

    if testStringList != setStringList:
        raise FortranRecordError("List from get_string doesn't match value from put_string.")
        
    return 1

def test_readFR_mixedRecord():

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

    testInt = testRecord.get_int()
    if testInt != setInt:
        raise FortranRecordError("Value from get_int doesn't match value from put_int.")

    testString = testRecord.get_string(12)
    if testString != setString:
        raise FortranRecordError("List from get_string doesn't match value from put_string.")

    testDoubleList = testRecord.get_double(2)
    if testDoubleList != setDoubleList:
        raise FortranRecordError("List from get_double doesn't match value from put_double.")

    testFloat = testRecord.get_float()
    # NOTE: since Python doesn't do native 32-bit floats, both values should be passed through
    #       a string formatting to truncate to 6 digits
    setFloat  = '%10.6f' % setFloat
    testFloat = '%10.6f' % testFloat

    if testFloat != setFloat:
        print str(setFloat) + " != " + str(testFloat)
        raise FortranRecordError("Value from get_float doesn't match value from put_float.")
    
    return 1


####################
#
# Test binary file operations here
#
####################

def test_openWritableBR():
    
    binaryFile = _BinaryReader('test.file','wb')
    if not binaryFile.f:
        raise BinaryReaderError("Failed to open new file for writing.")
    
    binaryFile.close()
    
    return 1

def test_writeBR():

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

def test_readBR():

    setInt = 8
    setFloat = 3.14
    setDoubleList = [1.6e-19,6.02e23]
    setString = "Hello World!"
    
    binaryFile = _BinaryReader('test_readBR.ref')
    
    testRecord = binaryFile.get_fortran_record()

    testInt = testRecord.get_int()
    if testInt != setInt:
        raise BinaryReaderError("Integer value was not as expected.")

    testString = testRecord.get_string(12)
    if testString != setString:
        raise BinaryReaderError("String was not as expected.")

    testDoubleList = testRecord.get_double(2)
    if testDoubleList != setDoubleList:
        raise BinaryReaderError("List of doubles was not as expected.")

    testFloat = testRecord.get_float()
    # NOTE: since Python doesn't do native 32-bit floats, both values should be passed through
    #       a string formatting to truncate to 6 digits
    setFloat  = '%10.6f' % setFloat
    testFloat = '%10.6f' % testFloat

    if testFloat != setFloat:
        print str(setFloat) + " != " + str(testFloat)
        raise BinaryReaderError("Float value was not as expected.")
    
    return 1
    

# start all tests here

tests = [0,0]

print "test_makeEmptyFR: ",
try:
    tests[0] += test_makeEmptyFR()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_resetFR: ",
try:
    tests[0] += test_resetFR()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_writeFR_singleInt: ",
try:
    tests[0] += test_writeFR_singleInt()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_writeFR_intList: ",
try:
    tests[0] += test_writeFR_intList()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_writeFR_singleLong: ",
try:
    tests[0] += test_writeFR_singleLong()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_writeFR_longList: ",
try:
    tests[0] += test_writeFR_longList()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_writeFR_singleFloat: ",
try:
    tests[0] += test_writeFR_singleFloat()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_writeFR_floatList: ",
try:
    tests[0] += test_writeFR_floatList()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_writeFR_singleDouble: ",
try:
    tests[0] += test_writeFR_singleDouble()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_writeFR_doubleList: ",
try:
    tests[0] += test_writeFR_doubleList()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_writeFR_singleString: ",
try:
    tests[0] += test_writeFR_singleString()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_writeFR_stringList: ",
try:
    tests[0] += test_writeFR_stringList()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_writeFR_mixedRecord: ",
try:
    tests[0] += test_writeFR_mixedRecord()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_readFR_singleInt: ",
try:
    tests[0] += test_readFR_singleInt()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_readFR_intList: ",
try:
    tests[0] += test_readFR_intList()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_readFR_singleLong: ",
try:
    tests[0] += test_readFR_singleLong()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_readFR_longList: ",
try:
    tests[0] += test_readFR_longList()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_readFR_singleFloat: ",
try:
    tests[0] += test_readFR_singleFloat()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_readFR_floatList: ",
try:
    tests[0] += test_readFR_floatList()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_readFR_singleDouble: ",
try:
    tests[0] += test_readFR_singleDouble()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_readFR_doubleList: ",
try:
    tests[0] += test_readFR_doubleList()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_readFR_singleString: ",
try:
    tests[0] += test_readFR_singleString()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_readFR_stringList: ",
try:
    tests[0] += test_readFR_stringList()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_readFR_mixedRecord: ",
try:
    tests[0] += test_readFR_mixedRecord()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_openWritableBR: ",
try:
    tests[0] += test_openWritableBR()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_writeBR: ",
try:
    tests[0] += test_writeBR()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1

print "test_readBR: ",
try:
    tests[0] += test_readBR()
    print "PASSED"
except Exception as inst:
    print "FAILED: " + str(inst)
    tests[1] += 1



print "Ran    " + str(tests[0]+tests[1]) + " tests."
print "Passed " + str(tests[0]) + " tests."
print "FAILED " + str(tests[1]) + " tests."
