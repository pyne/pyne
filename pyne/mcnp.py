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
import math 

from binaryreader import _BinaryReader

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
                values = [float(i) for i in get_words(self.f, lines = 1)]
            elif vars_per_cycle == 19:
                values = [float(i) for i in get_words(self.f, lines = 4)]
            
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
            
def get_words(f, lines = 1):
    words = []
    for i in range(lines):
        local_words = f.readline().split()
        words += local_words
    return words

class SourceSurf(object):
    def __init__(self):
        pass

class TrackData(object):
    def __init__(self):
        pass

class SurfSrc(_BinaryReader):

    def __init__(self, filename, mode="rb"):
        super(SurfSrc, self).__init__(filename,mode)

    def __str__(self):
        return self.print_header()

    def print_header(self):
        headerString  = "Code: {0} (version: {1}) [{2}]\n".format(self.kod, self.ver, self.loddat)
        headerString += "Problem info: ({0}) {1}\n{2}\n".format(self.idtm, self.probid, self.aid)
        headerString += "Showing dump #{0}\n".format(self.knod)
        headerString += "{0} histories, {1} tracks, {2} record size, {3} surfaces, {4} histories\n".format(self.np1, self.nrss, self.ncrd, self.njsw, self.niss)
        headerString += "{0} cells, source particle: {1}, macrobody facet flag: {2}\n".format(self.niwr, self.mipts, self.kjaq)
        for i in self.surflist:
            headerString += "Surface {0}: facet {1}, type {2} with {3} parameters: (".format(i.id, i.facetId, i.type, i.numParams)
            if i.numParams > 1:
                for j in i.surfParams:
                    headerString += " {0}".format(j)
            else:
                headerString += " {0}".format(i.surfParams)
            headerString += ")\n"
        headerString += "Summary Table: " + str(self.summaryTable)

        return headerString

    def print_tracklist(self):
        trackData = "Track Data\n"
  #                       1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
        trackData +=     "       nps   BITARRAY        WGT        ERG        TME             X             Y             Z          U          V     COSINE  |       W\n"
        for j in self.tracklist:
            trackData += "%10d %10g %10.5g %10.5g %10.5g %13.5e %13.5e %13.5e %10.5f %10.5f %10.5f  | %10.5f " % (j.nps, j.bitarray, j.wgt, j.erg, j.tme, j.x, j.y, j.z, j.u, j.v, j.cs, j.w) + "\n"

        return trackData

    def compare(self,other):
        if other.kod != self.kod:
            print "kod does not match"
            return False
        if other.ver != self.ver:
            print "ver does not match"
            return False
        if other.loddat != self.loddat:
            print "loddat does not match"
            return False

        if other.ncrd != self.ncrd:
            print "ncrd does not match"
            return False
        if other.njsw != self.njsw:
            print "njsw does not match"
            return False

        if other.niwr  != self.niwr:
            print "niwr does not match"
            return False
        if other.mipts != self.mipts:
            print "mipts does not match"
            return False
        if other.kjaq  != self.kjaq:
            print "kjaq does not match"
            return False

        for surf in range(len(self.surflist)):
            if other.surflist[surf].id         != self.surflist[surf].id:
                print "surf " + str(surf) + " ID doesn't match"
                return False
            if other.surflist[surf].facetId    != self.surflist[surf].facetId:
                print "surf " + str(surf) + " facetId doesn't match"
                return False
            if other.surflist[surf].type       != self.surflist[surf].id:
                print "surf " + str(surf) + " type doesn't match"
                return False
            if other.surflist[surf].numParams  != self.surflist[surf].numParams:
                print "surf " + str(surf) + " numParams doesn't match"
                return False
            if other.surflist[surf].surfParams != self.surflist[surf].surfParams:
                print "surf " + str(surf) + " surfParams doesn't match"
                return False

        return True

    def read_header(self):
        # read header block
        header = self.get_fortran_record()

        # interpret header block
        self.kod    = header.get_string(8)  # code identifier
        self.ver    = header.get_string(5)  # code version identifier
        self.loddat = header.get_string(8)  # code version date
        self.idtm   = header.get_string(19) # current date and time
        self.probid = header.get_string(19) # problem identification string
        self.aid    = header.get_string(80) # title card of initial run
        self.knod   = header.get_int()      # dump number

        # read various counts and sizes
        tablelengths = self.get_fortran_record()

        # interpret table lengths
        self.np1  = tablelengths.get_long()  # number of histories used to generate source
        self.nrss = tablelengths.get_long()  # number of tracks written to surface source
        self.ncrd = tablelengths.get_int()  # number of values in surface source record
                                            # 6 for a spherical source
                                            # 11 otherwise
        self.njsw = tablelengths.get_int()  # number of surfaces 
        self.niss = tablelengths.get_long()  # number of histories written to surface source

        if self.np1 < 0:
            # read more size info
            tablelengths = self.get_fortran_record()

            self.niwr  = tablelengths.get_int()  # number of cells in surface source card
            self.mipts = tablelengths.get_int()  # source particle type
            self.kjaq  = tablelengths.get_int()  # macrobody facet flag

        self.np1 = abs(self.np1)

        # get info for each surface
        self.surflist = []
        for j in range(self.njsw):
            surfaceinfo = self.get_fortran_record()
            
            surfinfo = SourceSurf()
            surfinfo.id = surfaceinfo.get_int()            # surface ID
            if self.kjaq == 1:
                surfinfo.facetId = surfaceinfo.get_int()   # facet ID
            else:
                surfinfo.facetId = -1                      # dummy facet ID

            surfinfo.type = surfaceinfo.get_int()                      # surface type
            surfinfo.numParams = surfaceinfo.get_int()                 # number of surface parameters
            surfinfo.surfParams = surfaceinfo.get_double(surfinfo.numParams)

            self.surflist.append(surfinfo)                  


        for j in range(self.njsw,self.njsw+self.niwr):
            self.get_fortran_record()
            print j

        summaryInfo = self.get_fortran_record()            # summary table
        self.summaryTable = summaryInfo.get_int((2+4*self.mipts)*(self.njsw+self.niwr)+1)
        self.summaryTable = self.summaryTable[1:]
        
    def read_tracklist(self):
        self.tracklist = []
        for j in range(self.nrss):
            trackInfo = self.get_fortran_record()
            trackData = TrackData()
            trackData.record   = trackInfo.get_double(self.ncrd)
            trackData.nps      = trackData.record[0]
            trackData.bitarray = trackData.record[1]
            trackData.wgt      = trackData.record[2]
            trackData.erg      = trackData.record[3]
            trackData.tme      = trackData.record[4]
            trackData.x        = trackData.record[5]
            trackData.y        = trackData.record[6]
            trackData.z        = trackData.record[7]
            trackData.u        = trackData.record[8]
            trackData.v        = trackData.record[9]
            trackData.cs       = trackData.record[10]
            trackData.w        = math.copysign(math.sqrt(1 - trackData.u*trackData.u - trackData.v*trackData.v),trackData.bitarray)
            # trackData.bitarray = abs(trackData.bitarray)
            
            self.tracklist.append(trackData)

class Srctp(_BinaryReader):

    def __init__(self, filename):
        super(Srctp, self).__init__(filename)

    def read(self):
        # read header block
        header = self.get_fortran_record()

        # interpret header block
        k = header.get_int() # unique code (947830)
        self.loc_next = header.get_int() # location of next site in FSO array (ixak)
        self.n_run = header.get_int() # source particles yet to be run (nsa)
        self.loc_store = header.get_int() # where to store next source neutron (ist)
        self.n_source = header.get_int() # number of source points in fso

        # read source site array
        fso = self.get_fortran_record()
        
        self.sites = []
        for i in range(self.n_source):
            vals = fso.get_double(11)

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


class Runtpe(_BinaryReader):

    def __init__(self, filename):
        super(Runtpe, self).__init__(filename)

    def read(self, filename):
        # read identification block
        header = self.get_fortran_record()

        # parse identification block
        self.codeName = header.get_string(8)
        self.codeVersion = header.get_string(5)
        self.codeDate = header.get_string(8)
        header.get_string(19) # machine designator, date and time
        self.chargeCode = header.get_string(10)
        self.problemID = header.get_string(19)
        self.problemIDsurf = header.get_string(19)
        self.title = header.get_string(80)
        header.pos += 3*6*11 # skip user file characteristics
        self.n_tables = header.get_int()

        # read cross-section tables
        self.tables = []
        for i in range(self.n_tables):
            self.tables.append(self.get_fortran_record())
            

    def __repr__(self):
        return "<Runtpe: {0}>".format(self.f.name)
