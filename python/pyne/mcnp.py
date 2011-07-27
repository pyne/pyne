#!/usr/bin/env python

"""
Module for parsing MCNP output data. MCNP is a general-purpose Monte Carlo
N-Particle code developed at Los Alamos National Laboratory that can be used for
neutron, photon, electron, or coupled neutron/photon/electron transport. Further
information on MCNP can be obtained from http://mcnp-green.lanl.gov/

Mctal and Runtpe classes still need work. Also should add Meshtal and Outp
classes.

"""

import struct

class Record(object):
    """
    Stores data from a single fortran record.

    Attributes:
      data       = list of bytes read in from a binary file
      pos        = position in data list
      intSize    = size of integer on running machine
      floatSize  = size of float on running machine
      doubleSize = size of double on running machine
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
    
    def getInt(self, n = 1):
        """
        Returns one or more integers at the current position within
        the data list. If more than one integer is read, the integers
        are returned in a list.
        """

        if self.pos >= self.numBytes:
            raise Exception("Already read all data from record")
        
        values = struct.unpack('{0}i'.format(n), self.data[self.pos:self.pos+self.intSize*n])
        self.pos += self.intSize * n
        if n == 1:
            return values[0]
        else:
            return list(values)
                             
    def getFloat(self, n = 1):
        """
        Returns one or more floats at the current position within the
        data list. If more than one float is read, the floats are
        returned in a list.
        """

        if self.pos >= self.numBytes:
            raise Exception("Already read all data from record")
        
        values = struct.unpack('{0}f'.format(n), self.data[self.pos:self.pos+self.floatSize*n])
        self.pos += self.floatSize * n
        if n == 1:
            return values[0]
        else:
            return list(values)
    
    def getDouble(self, n = 1):
        """
        Returns one or more doubles at the current position within the
        data list. If more than one double is read, the doubles are
        returned in a list.
        """

        if self.pos >= self.numBytes:
            raise Exception("Already read all data from record")
        
        values = struct.unpack('{0}d'.format(n),self.data[self.pos:self.pos+self.doubleSize*n])
        self.pos += self.doubleSize * n
        if n == 1:
            return values[0]
        else:
            return list(values)
    
    def getString(self, length):
        """
        Returns a string of a specified length starting at the current
        position in the data list.
        """

        if self.pos >= self.numBytes:
            raise Exception("Already read all data from record")
        
        relevantData = self.data[self.pos:self.pos+length]
        (s,) = struct.unpack('%ds'%length,relevantData)
        self.pos += length
        return s

    def reset(self):
        self.pos = 0

    def __repr__(self):
        return "<Record: {0} bytes>".format(self.numBytes)


class Mctal(object):
    def __init__(self):
        pass

    def read(self, filename):
        """
        Parses a 'mctal' tally output file from MCNP. Currently this
        only supports reading the kcode data- the remaining tally data
        will not be read.
        """

        # open file
        self.f = open(filename, 'r')

        # get code name, version, date/time, etc
        words = self.f.readline().split()
        self.codeName = words[0]
        self.codeVersion = words[1]
        self.codeDate = words[2]
        self.codeTime = words[3]
        self.n_dump = words[4]
        self.n_histories = int(words[5])
        self.n_prn       = int(words[6])

        # comment line of input file
        self.comment = self.f.readline().strip()

        # read tally line
        words = self.f.readline().split()
        self.n_tallies = words[1]
        if len(words) > 2:
            # perturbation tallies present
            pass

        # read tally numbers
        tally_nums = [int(i) for i in self.f.readline().split()]

        # read tallies
        for i_tal in tally_nums:
            pass

        # read kcode information
        words = self.f.readline().split()
        self.n_cycles = int(words[1])
        self.n_inactive = int(words[2])
        vars_per_cycle = int(words[3])

        self.k_col = []
        self.k_abs = []
        self.k_path = []
        self.prompt_life_col = []
        self.prompt_life_path = []
        self.avg_k_col = []
        self.avg_k_abs = []
        self.avg_k_path = []
        self.avg_k_combined = []
        self.avg_k_combined_active = []
        self.prompt_life_combined = []
        self.cycle_histories = []
        self.avg_k_combined_FOM = []

        for cycle in range(self.n_cycles):
            # read keff and prompt neutron lifetimes
            if vars_per_cycle == 0 or vars_per_cycle == 5:
                values = [float(i) for i in getWords(self.f, lines = 1)]
            elif vars_per_cycle == 19:
                values = [float(i) for i in getWords(self.f, lines = 4)]
            
            self.k_col.append(values[0])
            self.k_abs.append(values[1])
            self.k_path.append(values[2])
            self.prompt_life_col.append(values[3])
            self.prompt_life_path.append(values[4])

            if vars_per_cycle <= 5:
                continue

            avg, stdev = (values[5],values[6])
            self.avg_k_col.append((avg,stdev))
            avg, stdev = (values[7],values[8])
            self.avg_k_abs.append((avg,stdev))
            avg, stdev = (values[9],values[10])
            self.avg_k_path.append((avg,stdev))
            avg, stdev = (values[11],values[12])
            self.avg_k_combined.append((avg,stdev))
            avg, stdev = (values[13],values[14])
            self.avg_k_combined_active.append((avg,stdev))
            avg, stdev = (values[15],values[16])
            self.prompt_life_combined.append((avg,stdev))
            self.cycle_histories.append(values[17])
            self.avg_k_combined_FOM.append(values[18])
            
def getWords(f, lines = 1):
    words = []
    for i in range(lines):
        local_words = f.readline().split()
        words += local_words
    return words

class Srctp(object):

    def __init__(self):
        pass

    def read(self, filename):

        # open file
        self.f = open(filename,'rb')

        # read header block
        header = getRecord(self.f)

        # interpret header block
        k              = header.getInt() # unique code (947830)
        self.loc_next  = header.getInt() # location of next site in FSO array (ixak)
        self.n_run     = header.getInt() # source particles yet to be run (nsa)
        self.loc_store = header.getInt() # where to store next source neutron (ist)
        self.n_source  = header.getInt() # number of source points in fso

        # read source site array
        fso = getRecord(self.f)
        
        self.sites = []
        for i in range(self.n_source):
            vals = fso.getDouble(11)

            site = SourceSite()
            site.x = vals[0]
            site.y = vals[1]
            site.z = vals[2]
            site.E = vals[3]

            self.sites.append(site)

    def remainingSites(self):
        index = self.loc_next - 1
        return self.sites[index : index + self.n_run]

    def __repr__(self):
        return "<Srctp: {0}>".format(self.f.name)
        
        
class SourceSite(object):
    
    def __init__(self):
        pass

    def __repr__(self):
        return "<SourceSite: ({0.x},{0.y},{0.z})>".format(self)

class Runtpe(object):

    def __init__(self):
        pass

    def read(self, filename):

        # open file
        self.f = open(filename,'rb')

        # read identification block
        header = getRecord(self.f)

        # parse identification block
        self.codeName = header.getString(8)
        self.codeVersion = header.getString(5)
        self.codeDate = header.getString(8)
        header.getString(19) # machine designator, date and time
        self.chargeCode = header.getString(10)
        self.problemID = header.getString(19)
        self.problemIDsurf = header.getString(19)
        self.title = header.getString(80)
        header.pos += 3*6*11 # skip user file characteristics
        self.n_tables = header.getInt()

        # read cross-section tables
        self.tables = []
        for i in range(self.n_tables):
            self.tables.append(getRecord(self.f))
            

    def __repr__(self):
        return "<Runtpe: {0}>".format(self.f.name)


def getRecord(f):
    """
    Fortran unformatted records start with an int and end with the
    same int. This int represents the number of bytes that the record
    is. That makes it easy to read.
    """

    intSize = struct.calcsize('i')
    (numBytes,) = struct.unpack('i', f.read(intSize))

    data = f.read(numBytes)
        
    # now read end of record
    (numBytes2,) = struct.unpack('i', f.read(intSize))
    if numBytes2 != numBytes:
        raise Exception("Error while reading unformatted record.")
        
    return Record(data, numBytes)
