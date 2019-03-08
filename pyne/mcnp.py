#!/usr/bin/env python
"""Module for parsing MCNP output data. MCNP is a general-purpose Monte Carlo
N-Particle code developed at Los Alamos National Laboratory that can be used
for neutron, photon, electron, or coupled neutron/photon/electron transport.
Further information on MCNP can be obtained from http://mcnp.lanl.gov/

Mctal and Runtpe classes still need work. Also should add Meshtal and Outp
classes.

If PyMOAB is not installed, then Wwinp, Meshtal, and Meshtally will not be
available to use.

"""
from __future__ import print_function, division
from pyne.mesh import Mesh, StatMesh, HAVE_PYMOAB
import sys
import struct
import math
import os
import linecache
import datetime
from warnings import warn

import numpy as np
import tables

from pyne.utils import QAWarning
from pyne.material import Material
from pyne.material import MultiMaterial
from pyne import nucname
from pyne.binaryreader import _BinaryReader, _FortranRecord

warn(__name__ + " is not yet QA compliant.", QAWarning)

# Mesh specific imports

if HAVE_PYMOAB:
    from pyne.mesh import NativeMeshTag
else:
    warn("The PyMOAB optional dependency could not be imported. "
         "Some aspects of the mcnp module may be incomplete.",
         QAWarning)

if sys.version_info[0] > 2:
    def cmp(a, b):
        return (a > b) - (a < b)


class Mctal(object):
    def __init__(self):
        pass

    def read(self, filename):
        """Parses a 'mctal' tally output file from MCNP. Currently this
        only supports reading the kcode data- the remaining tally data
        will not be read.
        """

        # open file
        self.f = open(filename, 'r')

        # get code name, version, date/time, etc
        words = self.f.readline().split()
        self.code_name = words[0]
        self.code_version = words[1]
        self.code_date = words[2]
        self.code_time = words[3]
        self.n_dump = words[4]
        self.n_histories = int(words[5])
        self.n_prn = int(words[6])

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
                values = [float(i) for i in get_words(self.f, lines=1)]
            elif vars_per_cycle == 19:
                values = [float(i) for i in get_words(self.f, lines=4)]

            self.k_col.append(values[0])
            self.k_abs.append(values[1])
            self.k_path.append(values[2])
            self.prompt_life_col.append(values[3])
            self.prompt_life_path.append(values[4])

            if vars_per_cycle <= 5:
                continue

            avg, stdev = (values[5], values[6])
            self.avg_k_col.append((avg, stdev))
            avg, stdev = (values[7], values[8])
            self.avg_k_abs.append((avg, stdev))
            avg, stdev = (values[9], values[10])
            self.avg_k_path.append((avg, stdev))
            avg, stdev = (values[11], values[12])
            self.avg_k_combined.append((avg, stdev))
            avg, stdev = (values[13], values[14])
            self.avg_k_combined_active.append((avg, stdev))
            avg, stdev = (values[15], values[16])
            self.prompt_life_combined.append((avg, stdev))
            self.cycle_histories.append(values[17])
            self.avg_k_combined_FOM.append(values[18])


def get_words(f, lines=1):
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
    """Enables manipulating both the header and tracklists in surface source
    files.

    Example use cases include adding source particles from other codes, and
    combining multiple files together. Note that typically additional code
    will be needed to supplement this class in order to modify the header or
    track information in a way suitable to the use case.

    Parameters
    ----------
    filename : str
        Path to surface source file being read or written.
    mode : str, optional
        String indicating file opening mode to be used (defaults to 'rb').

    """

    def __init__(self, filename, mode="rb"):
        super(SurfSrc, self).__init__(filename, mode)

    def __str__(self):
        return self.print_header()

    def print_header(self):
        """Returns contents of SurfSrc's header as an informative string.

        Returns
        -------
        header_string : str
            A line-by-line listing of the contents of the SurfSrc's header.

        """
        header_string = "Code: {0} (version: {1}) [{2}]\n".format(
            self.kod, self.ver, self.loddat)
        header_string += "Problem info: ({0}) {1}\n{2}\n".format(
            self.idtm, self.probid, self.aid)
        header_string += "Showing dump #{0}\n".format(self.knod)
        header_string += (
            "{0} histories, {1} tracks, {2} record size, "
            "{0} surfaces, {1} histories\n").format(
            self.np1, self.nrss, self.ncrd,
            self.njsw, self.niss)
        header_string += (
            "{0} cells, source particle: {1},"
            " macrobody facet flag: {2}\n").format(
            self.niwr, self.mipts, self.kjaq)
        for i in self.surflist:
            header_string += (
                "Surface {0}: facet {1},"
                " type {2} with {3} parameters: (").format(
                i.id, i.facet_id, i.type, i.num_params)
            if i.num_params > 1:
                for j in i.surf_params:
                    header_string += " {0}".format(j)
            else:
                header_string += " {0}".format(i.surf_params)
            header_string += ")\n"
        header_string += "Summary Table: " + str(self.summary_table)

        return header_string

    def print_tracklist(self, max_tracks=None):
        """Returns tracklists in SurfSrc as a string.

        Parameters
        ----------
        max_tracks : int, optional
            Maximum number of tracks to print. Defaults to all tracks.

        Returns
        -------
        track_data : str
            Single string with data for one track on each line.
        """

        if max_tracks is None:
            max_tracks = self.nrss

        track_data = "Track Data\n"
        track_data += \
            "       nps   BITARRAY        WGT        ERG        TME" \
            "             X             Y             Z" \
            "          U          V     COSINE  |       W\n"
        for cnt, j in enumerate(self.tracklist):

            format_string = "%10d %10g %10.5g %10.5g %10.5g" \
                            " %13.5e %13.5e %13.5e" \
                            " %10.5f %10.5f %10.5f  | %10.5f "
            track_data += format_string % (
                j.nps, j.bitarray, j.wgt, j.erg, j.tme,
                j.x, j.y, j.z, j.u, j.v, j.cs, j.w) + "\n"
            if cnt > max_tracks:
                break

        return track_data

    def __eq__(self, other):
        rtn = self.__cmp__(other)
        return rtn == 0

    def __cmp__(self, other):
        """ Comparison is not completely robust. Tracklists are not compared!!!
        """

        if other.kod != self.kod:
            # kod does not match
            return cmp(other.kod, self.kod)
        if other.ver != self.ver:
            # ver does not match
            return cmp(other.ver, self.ver)
        if other.loddat != self.loddat:
            # loddat does not match
            return cmp(other.loddat, self.loddat)

        if other.ncrd != self.ncrd:
            # ncrd does not match
            return cmp(other.nrcd, self.nrcd)
        if other.njsw != self.njsw:
            # njsw does not match
            return cmp(other.njsw, self.njsw)

        if other.np1 != self.np1:
            # np1 does not match
            return cmp(other.np1, self.np1)
        if other.nrss != self.nrss:
            # nrss does not match
            return cmp(other.nrss, self.nrss)
        if other.niss != self.niss:
            # nrss does not match
            return cmp(other.niss, self.niss)

        if other.niwr != self.niwr:
            # niwr does not match
            return cmp(other.niwr, self.niwr)
        if other.mipts != self.mipts:
            # mipts does not match
            return cmp(other.mipts, self.mipts)
        if other.kjaq != self.kjaq:
            # kjaq does not match
            return cmp(other.kjaq, self.kjaq)

        for surf in range(len(self.surflist)):
            if other.surflist[surf].id != self.surflist[surf].id:
                # ID doesn't match
                return cmp(other.surflist[surf].id, self.surflist[surf].id)
            if other.surflist[surf].facet_id != self.surflist[surf].facet_id:
                # facet_id doesn't match
                return cmp(other.surflist[surf].facet_id,
                           self.surflist[surf].facet_id)
            if other.surflist[surf].type != self.surflist[surf].type:
                # type doesn't match
                return cmp(other.surflist[surf].type,
                           self.surflist[surf].type)
            if other.surflist[surf].num_params != \
                    self.surflist[surf].num_params:
                # num_params ddoesn't match
                return cmp(other.surflist[surf].num_params,
                           self.surflist[surf].num_params)
            if other.surflist[surf].surf_params != \
                    self.surflist[surf].surf_params:
                # surf_params doesn't match
                return cmp(other.surflist[surf].surf_params,
                           self.surflist[surf].surf_params)

        return 0

    def read_header(self):
        """Read in the header block data. This block comprises 4 fortran
        records which we refer to as: header, table1, table2, summary.
        """
        # read header record
        header = self.get_fortran_record()

        # interpret header
        self.kod = header.get_string(8)[0]  # code identifier

        if 'SF_00001' not in self.kod:
            self.ver = header.get_string(5)[0]  # code version identifier

            if '2.6.0' in self.ver:
                self.loddat = header.get_string(28)[0]  # code version date
            elif '5    ' in self.ver:
                self.loddat = header.get_string(8)[0]  # code version date
            else:
                raise NotImplementedError("MCNP5/X Version:" +
                                          self.ver.rstrip() + " not supported")

            self.idtm = header.get_string(19)[0]    # current date and time
            self.probid = header.get_string(19)[0]  # problem id string
            self.aid = header.get_string(80)[0]     # title card of initial run
            self.knod = header.get_int()[0]         # dump number

            # read table 1 record; various counts and sizes
            tablelengths = self.get_fortran_record()

            # interpret table lengths
            if '2.6.0' in self.ver:
                self.np1 = tablelengths.get_int()[0]    # hist used to gen. src
                self.nrss = tablelengths.get_int()[0]   # #tracks to surf src
            else:
                self.np1 = tablelengths.get_long()[0]   # hist used to gen. src
                self.nrss = tablelengths.get_long()[0]  # #tracks to surf src

	    # values in surf src record
            # 6 for a spherical source
            # 11 otherwise
            self.ncrd = tablelengths.get_int()[0]
            self.njsw = tablelengths.get_int()[0]  # number of surfaces
            self.niss = tablelengths.get_int()[0]  # #histories to surf src
            self.table1extra = list()
            while tablelengths.num_bytes > tablelengths.pos:
                self.table1extra += tablelengths.get_int()

        elif 'SF_00001' in self.kod:
            header = self.get_fortran_record()
            self.ver = header.get_string(12)[0]     # code version identifier
            self.loddat = header.get_string(9)[0]   # code version date
            self.idtm = header.get_string(19)[0]    # current date and time
            self.probid = header.get_string(19)[0]  # problem id string
            self.aid = header.get_string(80)[0]     # title card of initial run
            self.knod = header.get_int()[0]         # dump number

            # read table 1 record; various counts and sizes
            tablelengths = self.get_fortran_record()

            # interpret table lengths
            self.np1 = tablelengths.get_int()[0]     # hist used to gen.source
            self.notsure0 = tablelengths.get_int()[0]  # vals in surf src rec.
            self.nrss = tablelengths.get_int()[0]    # tracks writ. to surf.src
            self.notsure1 = tablelengths.get_int()[0]  # number of surfaces
            self.ncrd = tablelengths.get_int()[0]      # histories to surf.src
            self.njsw = tablelengths.get_int()[0]      # number of surfaces
            self.niss = tablelengths.get_int()[0]      # histories to surf.src
            self.table1extra = list()
            while tablelengths.num_bytes > tablelengths.pos:
                self.table1extra += tablelengths.get_int()

        if self.np1 < 0:
            # read table 2 record; more size info
            tablelengths = self.get_fortran_record()

            self.niwr = tablelengths.get_int()[0]   # #cells in surf.src card
            self.mipts = tablelengths.get_int()[0]  # source particle type
            self.kjaq = tablelengths.get_int()[0]   # macrobody facet flag
            self.table2extra = list()
            while tablelengths.num_bytes > tablelengths.pos:
                self.table2extra += tablelengths.get_int()

        else:
            pass

        # Since np1 can be negative, preserve the actual np1 value while
        # taking the absolute value so that np1 can be used mathematically
        self.orignp1 = self.np1
        self.np1 = abs(self.np1)

        # get info for each surface
        self.surflist = list()
        for j in range(self.njsw):
            # read next surface info record
            self.surfaceinfo = self.get_fortran_record()

            surfinfo = SourceSurf()
            surfinfo.id = self.surfaceinfo.get_int()            # surface ID
            if self.kjaq == 1:
                surfinfo.facet_id = self.surfaceinfo.get_int()  # facet ID
            else:
                surfinfo.facet_id = -1                        # dummy facet ID

            surfinfo.type = self.surfaceinfo.get_int()           # surface type
            surfinfo.num_params = self.surfaceinfo.get_int()[0]  # #surface prm
            surfinfo.surf_params = \
                self.surfaceinfo.get_double(surfinfo.num_params)

            self.surflist.append(surfinfo)

        # we read any extra records as determined by njsw+niwr...
        # no known case of their actual utility is known currently
        for j in range(self.njsw, self.njsw+self.niwr):
            self.get_fortran_record()
            warn("Extra info in header not handled: {0}".format(j),
                 RuntimeWarning)

        # read summary table record
        summary_info = self.get_fortran_record()
        self.summary_table = summary_info.get_int(
            (2+4*self.mipts)*(self.njsw+self.niwr)+1)
        self.summary_extra = list()
        while summary_info.num_bytes > summary_info.pos:
            self.summary_extra += summary_info.get_int()

    def read_tracklist(self):
        """Reads in track records for individual particles."""
        self.tracklist = []
        for j in range(self.nrss):
            track_info = self.get_fortran_record()
            track_data = TrackData()
            track_data.record = track_info.get_double(abs(self.ncrd))
            track_data.nps = track_data.record[0]
            track_data.bitarray = track_data.record[1]
            track_data.cell = abs(track_data.bitarray) // 8 % 100000000
            track_data.wgt = track_data.record[2]
            track_data.erg = track_data.record[3]
            track_data.tme = track_data.record[4]
            track_data.x = track_data.record[5]
            track_data.y = track_data.record[6]
            track_data.z = track_data.record[7]
            track_data.u = track_data.record[8]
            track_data.v = track_data.record[9]
            track_data.cs = track_data.record[10]
            track_data.w = math.copysign(
                math.sqrt(1 - track_data.u*track_data.u -
                          track_data.v*track_data.v), track_data.bitarray)
            # track_data.bitarray = abs(track_data.bitarray)

            self.tracklist.append(track_data)
        return

    def put_header(self):
        """Write the header part of the header to the surface source file"""
        if 'SF_00001' in self.kod:
            rec = [self.kod]
            joinrec = "".join(rec)
            newrecord = _FortranRecord(joinrec, len(joinrec))
            self.put_fortran_record(newrecord)

            rec = [self.ver, self.loddat, self.idtm, self.probid, self.aid]
            joinrec = "".join(rec)
            newrecord = _FortranRecord(joinrec, len(joinrec))
            newrecord.put_int([self.knod])
            self.put_fortran_record(newrecord)
        else:
            rec = [self.kod, self.ver, self.loddat,
                   self.idtm, self.probid, self.aid]
            joinrec = "".join(rec)
            newrecord = _FortranRecord(joinrec, len(joinrec))
            newrecord.put_int([self.knod])
            self.put_fortran_record(newrecord)
        return

    def put_table_1(self):
        """Write the table1 part of the header to the surface source file"""
        newrecord = _FortranRecord("", 0)

        if '2.6.0' in self.ver:
            newrecord.put_int([self.np1])
            newrecord.put_int([self.nrss])
        else:
            newrecord.put_long([self.np1])
            newrecord.put_long([self.nrss])

        newrecord.put_int([self.ncrd])
        newrecord.put_int([self.njsw])
        newrecord.put_int([self.niss])  # MCNP needs 'int', could be 'long' ?
        newrecord.put_int(self.table1extra)
        self.put_fortran_record(newrecord)
        return

    def put_table_2(self):
        """Write the table2 part of the header to the surface source file"""
        newrecord = _FortranRecord("", 0)
        newrecord.put_int([self.niwr])
        newrecord.put_int([self.mipts])
        newrecord.put_int([self.kjaq])
        newrecord.put_int(self.table2extra)
        self.put_fortran_record(newrecord)
        return

    def put_surface_info(self):
        """Write the record for each surface to the surface source file"""

        for cnt, s in enumerate(self.surflist):
            newrecord = _FortranRecord("", 0)
            newrecord.put_int(s.id)
            if self.kjaq == 1:
                newrecord.put_int(s.facet_id)  # don't add a 'dummy facet ID'
            # else no macrobody flag byte in the record

            newrecord.put_int(s.type)
            newrecord.put_int(s.num_params)
            newrecord.put_double(s.surf_params)

            self.put_fortran_record(newrecord)
        return

    def put_summary(self):
        """Write the summary part of the header to the surface source file"""
        newrecord = _FortranRecord("", 0)
        newrecord.put_int(list(self.summary_table))
        newrecord.put_int(list(self.summary_extra))
        self.put_fortran_record(newrecord)
        return

    def write_header(self):
        """Write the first part of the MCNP surface source file. The header content
        comprises five parts shown below.
        """
        self.put_header()
        self.put_table_1()
        self.put_table_2()
        self.put_surface_info()
        self.put_summary()

    def write_tracklist(self):
        """Write track records for individual particles. Second part of the MCNP
        surface source file.  Tracklist is also known as a 'phase space'.
        """

        for j in range(self.nrss):  # nrss is the size of tracklist
            newrecord = _FortranRecord("", 0)
            # 11 records comprising particle information
            newrecord.put_double(self.tracklist[j].nps)
            newrecord.put_double(self.tracklist[j].bitarray)
            newrecord.put_double(self.tracklist[j].wgt)
            newrecord.put_double(self.tracklist[j].erg)
            newrecord.put_double(self.tracklist[j].tme)
            newrecord.put_double(self.tracklist[j].x)
            newrecord.put_double(self.tracklist[j].y)
            newrecord.put_double(self.tracklist[j].z)
            newrecord.put_double(self.tracklist[j].u)
            newrecord.put_double(self.tracklist[j].v)
            newrecord.put_double(self.tracklist[j].cs)
            self.put_fortran_record(newrecord)
        return

    def update_tracklist(self, surf_src):
        """ Update tracklist from another surface source.
        This updates the surface source in-place.
        """

        # Catch for improper non-SurfSrc type
        if type(surf_src) != SurfSrc:
            raise TypeError('Surface Source is not of type SurfSrc')

        # Because 'kod' is the first header attribute
        elif not hasattr(surf_src, 'kod'):
            raise AttributeError(
                'No header attributes for surface source argument')
        elif not hasattr(self, 'kod'):
            raise AttributeError(
                'No header attributes read for surface source')

        # Because 'tracklist' forms the non-header portion
        elif not hasattr(surf_src, 'tracklist'):
            raise AttributeError(
                'No tracklist read for surface source argument')
        elif not hasattr(self, 'tracklist'):
            raise AttributeError(
                'No tracklist read for surface source')

        # No point in updating with self
        elif self == surf_src:
            raise ValueError('Tracklist cannot be updated with itself')

        self.tracklist = surf_src.tracklist
        self.nrss = surf_src.nrss

    def __del__(self):
        """Destructor. The only thing to do is close the file."""
        self.f.close()


class Srctp(_BinaryReader):
    """This class stores source site data from a 'srctp' file written by
    MCNP. The source sites are stored in the 'fso' array in MCNP.

    Parameters
    ----------
    filename : str
        Path to Srctp file being worked with.

    """

    def __init__(self, filename):
        super(Srctp, self).__init__(filename)

    def read(self):
        # read header block
        header = self.get_fortran_record()

        # interpret header block
        # NOTUSED k = header.get_int()  # unique code (947830)
        self.loc_next = header.get_int()  # loc. of next site in FSO arr (ixak)
        self.n_run = header.get_int()  # source particles yet to be run (nsa)
        self.loc_store = header.get_int()  # where to put nxt src neutron (ist)
        self.n_source = header.get_int()  # # of source points in fso (mrl)

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

    def remaining_sites(self):
        index = self.loc_next - 1
        if (self.loc_next + self.n_run) >= self.n_source:
            return (self.sites[index:] +
                    self.sites[:self.n_run - (self.n_source - index)])
        else:
            return self.sites[index: index + self.n_run]

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
        self.code_name = header.get_string(8)
        self.code_version = header.get_string(5)
        self.code_date = header.get_string(8)
        header.get_string(19)  # machine designator, date and time
        self.charge_code = header.get_string(10)
        self.problem_ID = header.get_string(19)
        self.problem_ID_surf = header.get_string(19)
        self.title = header.get_string(80)
        header.pos += 3*6*11  # skip user file characteristics
        self.n_tables = header.get_int()

        # read cross-section tables
        self.tables = []
        for i in range(self.n_tables):
            self.tables.append(self.get_fortran_record())

    def __repr__(self):
        return "<Runtpe: {0}>".format(self.f.name)


class Xsdir(object):
    """This class stores the information contained in a single MCNP xsdir file.

    Attributes
    ----------
    f : file handle
        The xsdir file.
    filename : str
        Path to the xsdir file.
    directory : str
        Path to the directory containing the xsdir file.
    datapath : str
        The data path specified in the first line of the xsdir file, if it
        exists.
    awr : dict
        Maps material ids to their atomic weight ratios.
    tables : list
        Entries are XsdirTable objects, that appear in the same order as the
        xsdir table lines.

    Notes
    -----
    See MCNP5 User's Guide Volume 3 Appendix K for more information.
    """

    def __init__(self, filename):
        """Parameters
        ----------
        filename : str
            Path to xsdir file.
        """
        self.f = open(filename, 'r')
        self.filename = os.path.abspath(filename)
        self.directory = os.path.dirname(filename)
        self.awr = {}
        self.tables = []

        self.read()

    def read(self):
        """Populate the Xsdir object by reading the file.
        """
        # Go to beginning of file
        self.f.seek(0)

        # Read first section (DATAPATH)
        line = self.f.readline()
        words = line.split()
        if words:
            if words[0].lower().startswith('datapath'):
                index = line.index('=')
                self.datapath = line[index+1:].strip()

        # Read second section
        line = self.f.readline()
        words = line.split()
        assert len(words) == 3
        assert words[0].lower() == 'atomic'
        assert words[1].lower() == 'weight'
        assert words[2].lower() == 'ratios'

        while True:
            line = self.f.readline()
            words = line.split()

            # Check for end of second section
            if len(words) % 2 != 0 or words[0] == 'directory':
                break

            for zaid, awr in zip(words[::2], words[1::2]):
                self.awr[zaid] = awr

        # Read third section
        while words[0] != 'directory':
            words = self.f.readline().split()

        while True:
            words = self.f.readline().split()
            if not words:
                break

            # Handle continuation lines
            while words[-1] == '+':
                extraWords = self.f.readline().split()
                words = words + extraWords
            assert len(words) >= 7

            # Create XsdirTable object and add to line
            table = XsdirTable()
            self.tables.append(table)

            # All tables have at least 7 attributes
            table.name = words[0]
            table.awr = float(words[1])
            table.filename = words[2]
            table.access = words[3]
            table.filetype = int(words[4])
            table.address = int(words[5])
            table.tablelength = int(words[6])

            if len(words) > 7:
                table.recordlength = int(words[7])
            if len(words) > 8:
                table.entries = int(words[8])
            if len(words) > 9:
                table.temperature = float(words[9])
            if len(words) > 10:
                table.ptable = (words[10] == 'ptable')

    def find_table(self, name):
        """Find all tables for a given ZIAD.

        Parameters
        ----------
        name : str
            The ZIAD name.

        Returns
        -------
        tables : list
            All XsdirTable objects for a given ZIAD.
        """
        tables = []
        for table in self:
            if name in table.name:
                tables.append(table)
        return tables

    def to_xsdata(self, filename):
        """Writes a Serpent xsdata file for all continuous energy xs tables.

        Parameters
        ----------
        filename : str
            The output filename.

        """
        xsdata = open(filename, 'w')
        for table in self.tables:
            if table.serpent_type == 1:
                xsdata.write(table.to_serpent() + '\n')
        xsdata.close()

    def __iter__(self):
        for table in self.tables:
            yield table

    def nucs(self):
        """Provides a set of the valid nuclide ids for nuclides contained
        in the xsdir.

        Returns
        -------
        valid_nucs : set
            The valid nuclide ids.
        """
        valid_nucs = set(nucname.id(table.name.split('.')[0])
                         for table in self.tables if
                         nucname.isnuclide(table.name.split('.')[0]))
        return valid_nucs


class XsdirTable(object):
    """Stores all information that describes a xsdir table entry, which appears
    as a single line in xsdir file. Attribute names are based off of those
    found in the MCNP5 User's Guide Volume 3, appendix K.

    Attributes
    ----------
    name : str
        The ZAID and library identifier, delimited by a '.'.
    awr : float
        The atomic mass ratio of the nuclide.
    filename : str
        The relative path of the file containing the xs table.
    access : str
       Additional string to specify an access route, such as UNIX directory.
       This entry is typically 0.
    filetype : int
        Describes whether the file contains formated (1) or unformated (2)
        file.
    address : int
        If filetype is 1, address is the line number of the xsdir table. If
        filetype is 2, address is the record number.
    tablelength : int
        Length of the second block of a data table.
    recordlength : int
        Unused for filetype = 1. For filetype = 2, recordlength is the number
        of entires per record times the size (in bytes) of each entry.
    entries : int
        Unused for filetype = 1. For filetype = 2, it is the number of entries
        per record
    temperature : float
        Temperature in MeV for neutron data only.
    ptable : bool
        True if xs table describes continuous energy neutron data with
        unresolved resonance range probability tables.
    """

    def __init__(self):
        self.name = None
        self.awr = None
        self.filename = None
        self.access = None
        self.filetype = None
        self.address = None
        self.tablelength = None
        self.recordlength = None
        self.entries = None
        self.temperature = None
        self.ptable = False

    @property
    def alias(self):
        """Returns the name of the table entry <ZIAD>.<library id>.
        """
        return self.name

    @property
    def serpent_type(self):
        """Converts cross section table type to Serpent format:
            :1: continuous energy (c).
            :2: dosimetry table (y).
            :3: termal (t).
        """
        if self.name.endswith('c'):
            return 1
        elif self.name.endswith('y'):
            return 2
        elif self.name.endswith('t'):
            return 3
        else:
            return None

    @property
    def metastable(self):
        """Returns 1 is xsdir table nuclide is metastable. Returns zero
        otherwise.
        """
        # Only valid for neutron cross-sections
        if not self.name.endswith('c'):
            return

        # Handle special case of Am-242 and Am-242m
        if self.zaid == '95242':
            return 1
        elif self.zaid == '95642':
            return 0

        # All other cases
        A = int(self.name.split('.')[0]) % 1000
        if A > 600:
            return 1
        else:
            return 0

    @property
    def zaid(self):
        """Returns the ZIAD of the nuclide.
        """
        return self.name[:self.name.find('.')]

    def to_serpent(self, directory=''):
        """Converts table to serpent format.

        Parameters
        ----------
        directory : str
            The directory where Serpent data is to be stored.
        """
        # Adjust directory
        if directory:
            if not directory.endswith('/'):
                directory = directory.strip() + '/'

        return "{0} {0} {1} {2} {3} {4} {5:.11e} {6} {7}".format(
            self.name,
            self.serpent_type, self.zaid, 1 if self.metastable else 0,
            self.awr, self.temperature/8.6173423e-11, self.filetype - 1,
            directory + self.filename)

    def __repr__(self):
        return "<XsDirTable: {0}>".format(self.name)


class PtracEvent(tables.IsDescription):
    """This class holds one Ptrac event and serves as a table definition
    for saving Ptrac data to a HDF5 file.
    """
    event_type = tables.Int32Col()
    node = tables.Float32Col()
    nsr = tables.Float32Col()
    nsf = tables.Float32Col()
    nxs = tables.Float32Col()
    ntyn = tables.Float32Col()
    ipt = tables.Float32Col()
    ncl = tables.Float32Col()
    mat = tables.Float32Col()
    ncp = tables.Float32Col()
    xxx = tables.Float32Col()
    yyy = tables.Float32Col()
    zzz = tables.Float32Col()
    uuu = tables.Float32Col()
    vvv = tables.Float32Col()
    www = tables.Float32Col()
    erg = tables.Float32Col()
    wgt = tables.Float32Col()
    tme = tables.Float32Col()


class PtracReader(object):
    """Class to read _binary_ PTRAC files generated by MCNP.
    """

    def __init__(self, filename):
        """Construct a new Ptrac reader for a given filename, determine the
        number format and read the file's headers.
        """
        self.variable_mappings = {
            1: "nps",
            3: "ncl",
            4: "nsf",  # surface id
            8: "node",
            9: "nsr",
            10: "nxs",
            11: "ntyn",
            12: "nsf",
            16: "ipt",
            17: "ncl",
            18: "mat",
            19: "ncp",
            20: "xxx",  # position x
            21: "yyy",  # position y
            22: "zzz",  # position z
            23: "uuu",  # cos(x-direction)
            24: "vvv",  # cos(y-direction)
            25: "www",  # cos(z-direction)
            26: "erg",  # energy
            27: "wgt",  # mass
            28: "tme"
        }

        self.eightbytes = False

        self.f = open(filename, 'rb')
        self.determine_endianness()
        self.read_headers()
        self.read_variable_ids()

        self.next_event = 0

    def __del__(self):
        """Destructor. The only thing to do is close the Ptrac file.
        """
        self.f.close()

    def determine_endianness(self):
        """Determine the number format (endianness) used in the Ptrac file.
        For this, the file's first entry is used. It is always minus one
        and has a length of 4 bytes.
        """

        # read and unpack first 4 bytes
        b = self.f.read(4)
        should_be_4 = struct.unpack('<i', b)[0]
        if should_be_4 == 4:
            self.endianness = '<'
        else:
            self.endianness = '>'

        # discard the next 8 bytes (the value -1 und another 4)
        self.f.read(8)

    def read_next(self, format, number=1, auto=False, raw_format=False):
        """Helper method for reading records from the Ptrac file.
        All binary records consist of the record content's length in bytes,
        the content itself and then the length again.
        format can be one of the struct module's format characters (i.e. i
        for an int, f for a float, s for a string).
        The length of the record can either be hard-coded by setting the
        number parameter (e.g. to read 10 floats) or determined automatically
        by setting auto=True.
        Setting the parameter raw_format to True means that the format string
        will not be expanded by number, but will be used directly.
        """

        if self.eightbytes and (not raw_format) and format == 'f':
            format = 'd'
        if self.eightbytes and (not raw_format) and format == 'i':
            format = 'q'

        # how long is one field of the read values
        format_length = 1
        if format in ['h', 'H'] and not raw_format:
            format_length = 2
        elif format in ['i', 'I', 'l', 'L', 'f'] and not raw_format:
            format_length = 4
        elif format in ['d', 'q', 'Q'] and not raw_format:
            format_length = 8

        if auto and not raw_format:
            b = self.f.read(4)

            if b == b'':
                raise EOFError

            length = struct.unpack(self.endianness.encode() + b'i', b)[0]
            number = length // format_length

            b = self.f.read(length + 4)
            tmp = struct.unpack(b"".join([self.endianness.encode(),
                                          (format*number).encode(), b'i']), b)
            length2 = tmp[-1]
            tmp = tmp[:-1]
        else:
            bytes_to_read = number * format_length + 8
            b = self.f.read(bytes_to_read)
            if b == b'':
                raise EOFError

            fmt_string = self.endianness + "i"
            if raw_format:
                fmt_string += format + "i"
            else:
                fmt_string += format * number + "i"

            tmp = struct.unpack(fmt_string.encode(), b)
            length = tmp[0]
            length2 = tmp[-1]
            tmp = tmp[1:-1]

        assert length == length2

        if format == 's':
            # return just one string
            return b''.join(tmp).decode()
        elif number == 1:
            # just return the number and not a tuple containing just the number
            return tmp[0]
        else:
            # convert tuple to list
            return list(tmp)

    def read_headers(self):
        """Read and save the MCNP version and problem description from the
        Ptrac file.
        """
        # mcnp version info
        self.mcnp_version_info = self.read_next('s', auto=True)
        # problem title
        self.problem_title = self.read_next('s', auto=True).strip()

        # ptrac input data. can be omitted for now,
        # but has to be parsed, because it has variable length.
        # Also, this is the first difference between a file generated
        # with 4-byte and 8-byte numbers.
        line = self.read_next('f', auto=True)
        # if this line doesn't consist of 10 floats, then we've read them with
        # the wrong byte length and re have to re-read them (and every
        # following float) with 8 bytes length.
        if len(line) != 10:
            self.eightbytes = True
            tmp = struct.pack(self.endianness + "f"*20, *line)
            line = list(struct.unpack(self.endianness + "d"*10, tmp))

        # the first item is always 13. afterwards, there is 13 times the
        # following scheme:
        # N x_0 ... x_N,
        # where N is the number of values for the current input variable and
        # the x_i are its N values.
        num_variables = int(line[0])  # should always be 13.
        current_pos = 1
        current_variable = 1

        while current_variable <= num_variables:
            n = int(line[current_pos])
            if current_variable < num_variables and (current_pos + n + 1) >= \
                    len(line):
                line += self.read_next('f', 10)
            current_pos += n + 1
            current_variable += 1

    def read_variable_ids(self):
        """Read the list of variable IDs that each record type in the Ptrac
        file is comprised of. The variables can vary for different problems.
        Consult the MCNP manual for details.
        """

        variable_nums = dict()
        variable_ids = dict()

        if self.eightbytes:
            variable_info = self.read_next(
                "qqqqqqqqqqqiiiiiiiii", 124, raw_format=True)
        else:
            variable_info = self.read_next('i', 20)

        variable_nums["nps"] = variable_info[0]
        variable_nums["src"] = variable_info[1] + variable_info[2]
        variable_nums["bnk"] = variable_info[3] + variable_info[4]
        variable_nums["sur"] = variable_info[5] + variable_info[6]
        variable_nums["col"] = variable_info[7] + variable_info[8]
        variable_nums["ter"] = variable_info[9] + variable_info[10]

        num_vars_total = sum(variable_info[:11])

        if self.eightbytes:
            # only the NPS vars are in 8 byte, the other ones are still 4
            fmt_string = "q" * variable_info[0] + \
                "i" * sum(variable_info[1:11])
            fmt_length = 8 * variable_info[0] + 4 * sum(variable_info[1:11])
            all_var_ids = self.read_next(
                fmt_string, fmt_length, raw_format=True)
        else:
            all_var_ids = self.read_next('i', num_vars_total)

        for l in ["nps", "src", "bnk", "sur", "col", "ter"]:
            variable_ids[l] = all_var_ids[:variable_nums[l]]
            all_var_ids = all_var_ids[variable_nums[l]:]

        self.variable_nums = variable_nums
        self.variable_ids = variable_ids

    def read_nps_line(self):
        """Read an NPS record and save the type of the next event.
        """
        nps_line = self.read_next('i', self.variable_nums["nps"])
        self.next_event = nps_line[1]

    def read_event_line(self, ptrac_event):
        """Read an event record and save it to a given PtracParticle instance.
        """

        # save for current event, because this record
        # contains only the next event's type
        event_type = self.next_event

        if event_type == 1000:
            e = "src"
        elif event_type == 3000:
            e = "sur"
        elif event_type == 4000:
            e = "col"
        elif event_type == 5000:
            e = "ter"
        else:
            e = "bnk"

        evt_line = self.read_next('f', self.variable_nums[e])

        self.next_event = evt_line[0]

        for i in range(1, len(self.variable_ids[e])):
            if self.variable_ids[e][i] in self.variable_mappings:
                ptrac_event[self.variable_mappings[
                    self.variable_ids[e][i]]] = \
                    evt_line[i]
        ptrac_event["event_type"] = event_type

    def write_to_hdf5_table(self, hdf5_table, print_progress=0):
        """Writes the events contained in this Ptrac file to a given HDF5
        table. The table must already exist and have rows that match the
        PtracEvent definition.
        If desired, the number of processed events can be printed to the
        console each N events by passing the print_progress=N parameter.
        """

        ptrac_event = hdf5_table.row
        counter = 0

        while True:
            try:
                self.read_nps_line()
            except EOFError:
                break  # no more entries

            while self.next_event != 9000:
                self.read_event_line(ptrac_event)
                ptrac_event.append()

                counter += 1
                if print_progress > 0 and counter % print_progress == 0:
                    print("processing event {0}".format(counter))


def _is_cell_line(line):
    is_cell = False
    if len(line.split()) > 3:
        if line.split()[0].isdigit() and \
           line.split()[1].isdigit() and \
           not line.split()[2][0].isalpha() and \
           line[0:5] != '     ' and \
           line.split()[1] != '0':
            is_cell = True
    return is_cell


def mats_from_inp(inp):
    """This function reads an MCNP inp file and returns a mapping of material
    numbers to material objects.

    Parameters
    ----------
    inp : str
        MCNP input file

    Returns
    --------
    materials : dict
       Keys are MCNP material numbers and values are PyNE material objects (for
       single density materials) and MultiMaterial objects (for multiple density
       materials).
    """

    mat_lines = []  # line of lines that begin material cards
    densities = {}  # dictionary to be populuated with material # and densities
    mat_nums = []  # list of material numbers to be printed to stdout
    materials = {}  # to be populated with material objectes and returned

    line_count = 1
    line = linecache.getline(inp, line_count)
    # scroll through every line of the mcnp inp file
    while line != '':
        line = linecache.getline(inp, line_count)
        # check to see if line contains a cell card. If so, grab the density.
        # information is stored in a dictionary where:
        # key = material number, value = list of densities
        if _is_cell_line(line):
            mat_num = int(line.split()[1])
            den = float(line.split()[2])

            if mat_num not in densities.keys():
                densities[mat_num] = [den]

            else:
                same_bool = False
                for j in range(0, len(densities[mat_num])):
                    if abs((den - densities[mat_num][j])/den) < 1E-4:
                        same_bool = True

                if same_bool is False:
                    densities[mat_num].append(den)

        # check line to see if it contain a material card, in the form
        # m* where * is a digit. If so store the line num. and material number
        if line.split() != []:
            if line.split()[0][0] == 'm' or line.split()[0][0] == 'M':
                if line.split()[0][1].isdigit() is True:
                    mat_lines.append(line_count)
                    mat_nums.append(int(line.split()[0][1:]))

        line_count += 1
        line = linecache.getline(inp, line_count)

    for i in range(0, len(mat_nums)):
        if mat_nums[i] in densities.keys():
            materials[mat_nums[i]] = mat_from_inp_line(inp, mat_lines[i],
                                                       densities[mat_nums[i]])
        else:
            materials[mat_nums[i]] = mat_from_inp_line(inp, mat_lines[i])
    return materials


def mat_from_inp_line(filename, mat_line, densities='None'):
    """ This function reads an MCNP material card from a file and returns a
    Material or Multimaterial object for the material described by the card.
    This function is used by :func:`mats_from_inp`.

    Parameters
    ----------
    filename : str
        Name of the MCNP input file
    mat_line : int
        Line number of the material card or interest
    densities : list of floats
        The densities associated with the material

    Returns
    -------
    finished_mat : Material or MultiMaterial
        A Material object is returned if there is 1 density supplied. If
        multiple densities are supplied a MultiMaterial is returned.
    """

    data_string = linecache.getline(filename, mat_line).split('$')[0]
    # collect all material card data on one string
    line_index = 1
    line = linecache.getline(filename, mat_line + line_index)
    # people sometimes put comments in materials and then this loop breaks                                                                                       # so we need to keep reading if we encounter comments
    while len(line.split()) > 0 and (line[0:5] == '     ' or line[0].lower() == 'c'):
        # make sure element/isotope is not commented out
        if line.split()[0][0] != 'c' and line.split()[0][0] != 'C':
            data_string += line.split('$')[0]
            line_index += 1
            line = linecache.getline(filename, mat_line + line_index)
        # otherwise this not a line we care about, move on and
        # skip lines that start with c or C
        else:
            line_index += 1
            line = linecache.getline(filename, mat_line + line_index)

    # create dictionaries nucvec and table_ids
    nucvec = {}
    table_ids = {}
    for i in range(1, len(data_string.split())):
        if i & 1 == 1:
            zzzaaam = str(nucname.zzaaam(
                nucname.mcnp_to_id(data_string.split()[i].split('.')[0])))

            # this allows us to read nuclides that are repeated
            if zzzaaam in nucvec.keys():
                nucvec[zzzaaam] += float(data_string.split()[i+1])
            else:
                nucvec[zzzaaam] = float(data_string.split()[i+1])

            if len(data_string.split()[i].split('.')) > 1:
                table_ids[str(zzzaaam)] = data_string.split()[i].split('.')[1]

    # Check to see it material is definted my mass or atom fracs.
    # Do this by comparing the first non-zero fraction to the rest
    # If atom fracs, convert.
    nucvecvals = list(nucvec.values())
    n = 0
    isatom = 0 < nucvecvals[n]
    while 0 == nucvecvals[n]:
        n += 1
        isatom = 0 < nucvecvals[n]
    for value in nucvecvals[n+1:]:
        if isatom != (0 <= value):
            msg = 'Mixed atom and mass fractions not supported.'
            ' See material defined on line {0}'.format(mat_line)
            warn(msg)

    # apply all data to material object
    if isatom:
        mat = Material()
        mat.from_atom_frac(nucvec)
    else:
        # set nucvec attribute to the nucvec dict from above
        mat = Material(nucvec=nucvec)

    mat.metadata['table_ids'] = table_ids
    mat.metadata['mat_number'] = data_string.split()[0][1:]

    # collect metadata, if present
    mds = ['source', 'comments', 'name']
    line_index = 1
    mds_line = linecache.getline(filename, mat_line - line_index)
    # while reading non-empty comment lines
    while mds_line.strip() not in set('cC') \
            and mds_line.split()[0] in ['c', 'C']:
        if mds_line.split()[0] in ['c', 'C'] \
                and len(mds_line.split()) > 1:
            possible_md = mds_line.split()[1].split(':')[0].lower()
            if possible_md in mds:
                if possible_md.lower() == 'comments':
                    comments_string = str(
                        ''.join(mds_line.split(':')[1:]).split('\n')[0])
                    comment_index = 1
                    comment_line = linecache.getline(
                        filename, mat_line - line_index + comment_index)
                    while comment_line.split()[0] in ['c', 'C']:
                        if comment_line.split()[1].split(':')[0].lower() in mds:
                            break
                        comments_string += ' ' + ' '.join(
                            comment_line.split()[1:])
                        comment_index += 1
                        comment_line = \
                            linecache.getline(filename,
                                              mat_line - line_index +
                                              comment_index)
                    mat.metadata[possible_md] = comments_string
                else:
                    mat.metadata[possible_md] = ''.join(
                        mds_line.split(':')[1:]).split('\n')[0]
                    # set metadata
        line_index += 1
        mds_line = linecache.getline(filename, mat_line - line_index)

    # Check all the densities. If they are atom densities, convert them to mass
    # densities. If they are mass densities they willl be negative, so make
    # them positive.
    if densities != 'None':
        converted_densities = []
        for den in densities:
            if den <= 0:
                converted_densities.append(-1*float(den))
            else:
                converted_densities.append(mat.mass_density(float(den)*1E24))

        # check to see how many densities are associated with this material.
        # if there is more than one, create a multimaterial"""
        if len(converted_densities) == 1:
            mat.density = converted_densities[0]
            finished_mat = mat

        elif len(converted_densities) > 1:
            mat_dict = {}
            for density in converted_densities:
                mat2 = Material()
                mat2.comp = mat.comp
                mat2.atoms_per_molecule = mat.atoms_per_molecule
                mat2.mass = mat.mass
                mat2.metadata = mat.metadata
                mat2.density = density
                mat_dict[mat2] = 1
            finished_mat = MultiMaterial(mat_dict)
    else:
        finished_mat = mat

    return finished_mat


class Wwinp(Mesh):
    """A Wwinp object stores all of the information from a single MCNP WWINP
    file. Weight window lower bounds are stored on a structured mesh. Only
    Cartesian mesh WWINP files are supported. Neutron, photon, and
    simotaneous neutron and photon WWINP files are supported.

    Attributes
    ----------
    ni : number of integers on card 2.
        ni = 1 for neutron WWINPs, ni = 2 for photon WWINPs or neutron + photon WWINPs.
    nr : int
        10 for rectangular, 16 for cylindrical.
    ne : list of number of energy groups for neutrons and photons.
        If ni = 1 the list is only 1 value long,
        to represent the number of neutron energy groups
    nf : list of numbers
        of fine mesh points in the i, j, k dimensions
    nft : int
        total number of fine mesh points
    origin : list of i, j, k
        minimums.
    nc : list
        number of coarse mesh points in the i, j, k dimensions
    nwg : int
        1 for rectangular, 2 for cylindrical.
    cm : list of lists
        of coarse mesh points in the i, j, k dimensions. Note
        the origin is not considered a coarse mesh point (as in MCNP).
    fm : list of lists
        of number of fine mesh points between the coarses
        mesh points in the i, j, k dimensions.
    e : list of lists
        of energy upper bounds for neutrons, photons. If
        ni = 1, the e will look like [[]]. If ni = 2, e will look like
        [[], []].
    bounds : list of lists
        of spacial bounds in the i, j, k dimensions.
    mesh : Mesh object
        with a structured mesh containing all the neutron and/or
        photon weight window lower bounds. These tags have the form
        "ww_X" where X is n or p The mesh has rootSet tags in the form
        X_e_upper_bounds.

    Notes
    -----
    Attribute names are identical to names speficied in WWINP file
    description in the MCNP5 User's Guide Volume 3 Appendix J.

    """

    def __init__(self):
        if not HAVE_PYMOAB:
            raise RuntimeError("PyMOAB is not available, "
                               "unable to create Wwinp Mesh.")
        pass

    def read_wwinp(self, filename):
        """This method creates a Wwinp object from the WWINP file <filename>.
        """
        with open(filename, 'r') as f:
            self._read_block1(f)
            self._read_block2(f)
            self._read_block3(f)

    def _read_block1(self, f):
        # Retrieves all of the information from block 1 of a wwinp file.

        line_1 = f.readline()
        self.ni = int(line_1.split()[2])
        self.nr = int(line_1.split()[3])

        line_2 = f.readline()
        self.ne = [int(x) for x in line_2.split()]

        if self.nr == 10:  # Cartesian
            line_3 = f.readline()
            self.nf = [int(float(x)) for x in line_3.split()[0:3]]
            self.nft = self.nf[0]*self.nf[1]*self.nf[2]
            self.origin = [float(x) for x in line_3.split()[3:6]]

            line_4 = f.readline()
            self.nc = [int(float(x)) for x in line_4.split()[0:3]]
            self.nwg = int(float(line_4.split()[3]))

        if self.nr == 16:  # Cylindrical
            raise ValueError('Cylindrical WWINP not currently supported')

    def _read_block2(self, f):
        # Retrieves all of the information from block 2 of a wwinp file.

        self.bounds = [[], [], []]
        self.cm = [[], [], []]
        self.fm = [[], [], []]

        for i in [0, 1, 2]:
            # Create a list of raw block 2 values.
            raw = []
            while len(raw) < 3*self.nc[i] + 1:
                raw += [float(x) for x in f.readline().split()]

            # Remove all the rx(i), ry(i), rz(i) values that
            # contaminated the raw list.
            removed_values = [raw[0]]
            for j in range(1, len(raw)):
                if j % 3 != 0:
                    removed_values.append(raw[j])

            # Expaned out nfx/nfy/nfx values to get structured mesh bounds.
            for j in range(0, len(removed_values)):
                if j % 2 == 0:
                    self.bounds[i].append(removed_values[j])
                    if j != 0:
                        self.cm[i].append(removed_values[j])

                else:
                    self.fm[i].append(removed_values[j])
                    for k in range(1, int(removed_values[j])):
                        self.bounds[i].append(
                            (removed_values[j+1] - removed_values[j-1])
                            * k / removed_values[j] + removed_values[j-1])

    def _read_block3(self, f):
        # Retrives all the information of the block 3 of a wwinp file.

        self.e = [[]]
        if self.ne[0] != 0:
            while len(self.e[0]) < self.ne[0]:
                self.e[0] += [float(x) for x in f.readline().split()]

            self._read_wwlb('n', f)

        if len(self.ne) == 2 and self.ne[1] != 0:
            self.e.append([])
            while len(self.e[-1]) < self.ne[1]:
                self.e[-1] += [float(x) for x in f.readline().split()]

            self._read_wwlb('p', f)

    def _read_wwlb(self, particle, f):
        # Reads the weight window lower bounds from block 3 and returns a
        # mesh.

        # If this is the first time this method is called then created a mesh,
        # otherwise (in the case of n and p in the same WWINP) add to the
        # preexisting mesh.
        if not hasattr(self, 'mesh'):
            super(Wwinp, self).__init__(structured_coords=[self.bounds[0],
                                                           self.bounds[1], self.bounds[2]],
                                        structured=True)

        volume_elements = list(self.structured_iterate_hex('zyx'))

        if particle == 'n':
            particle_index = 0

        elif particle == 'p':
            particle_index = 1

        # read in WW data for a single particle type
        ww_data = np.empty(shape=(self.ne[particle_index], self.nft))
        for i in range(0, self.ne[particle_index]):
            count = 0
            ww_row = []
            while count < self.nft:
                ww_row += [float(x) for x in f.readline().split()]
                count += 6  # number of entries per row in WWINP

            ww_data[i] = ww_row

        # create vector tags for data
        ww_tag_name = "ww_{0}".format(particle)
        self.tag(ww_tag_name, size=self.ne[particle_index],
                 dtype=float, tagtype='nat_mesh')
        tag_ww = self.get_tag(ww_tag_name)

        # tag vector data to mesh
        for i, volume_element in enumerate(volume_elements):
            tag_ww[volume_element] = ww_data[:, i]

        # Save energy upper bounds to rootset.
        e_bounds_tag_name = '{0}_e_upper_bounds'.format(particle)
        self.tag(e_bounds_tag_name,
                 size=len(self.e[particle_index]),
                 dtype=float, tagtype='nat_mesh')
        tag_e_bounds = self.get_tag(e_bounds_tag_name)
        tag_e_bounds[self] = self.e[particle_index]

    def write_wwinp(self, filename):
        """This method writes a complete WWINP file to <filename>.
        """
        with open(filename, 'w') as f:
            self._write_block1(f)
            self._write_block2(f)
            self._write_block3(f)

    def _write_block1(self, f):
        # Writes the all block 1 data to WWINP file

        block1 = ''

        # Create a MCNP formated time string.
        now = datetime.datetime.now()
        time = '{0:02d}/{1}/{2} {3}:{4}:{5}'.format(
            now.month, now.day, str(now.year)[2:],
            now.hour, now.minute, now.second)

        # Append line 1.
        block1 += \
            "{0:10.0f}{1:10.0f}{2:10.0f}{3:10.0f}" \
            "{4:>38}\n".format(1, 1,  self.ni, self.nr, time)

        # Append line 2.
        for i in self.ne:
            block1 += '{0:10.0f}'.format(int(i))

        block1 += '\n'

        # Append line 3.
        block1 += \
            '{0:13.5E}{1:13.5E}{2:13.5E}{3:13.5E}{4:13.5E}{5:13.5E}\n'\
            .format(self.nf[0], self.nf[1], self.nf[2],
                    self.origin[0], self.origin[1], self.origin[2])

        # Append line 4.
        block1 += '{0:13.5E}{1:13.5E}{2:13.5E}{3:13.5E}\n'\
            .format(self.nc[0], self.nc[1], self.nc[2], self.nwg)

        f.write(block1)

    def _write_block2(self, f):
        # Writes the all block 2 data to WWINP file

        # Create an array of values to be print in block 2.
        block2_array = [[], [], []]
        for i in [0, 1, 2]:
            block2_array[i].append(self.origin[i])
            for j in range(0, len(self.cm[i])):
                block2_array[i] += [self.fm[i][j], self.cm[i][j], 1.0000]

        # Translate block2 vector into a string with appropriate text wrapping.
        block2 = ""
        for i in range(0, 3):
            line_count = 0  # number of entries printed to current line, max=6
            for j in range(0, len(block2_array[i])):
                block2 += '{0:13.5E}'.format(block2_array[i][j])
                line_count += 1
                if line_count == 6:
                    block2 += '\n'
                    line_count = 0

            if line_count != 0:
                block2 += '\n'

        f.write(block2)

    def _write_block3(self, f):
        # Writes the all block 3 data to WWINP file

        if self.ne[0] != 0:
            self._write_block3_single('n', f)

        if len(self.ne) == 2:
            self._write_block3_single('p', f)

    def _write_block3_single(self, particle, f):
        # Write all of block 3 a single time (e.g. for WWINP with only n or
        #   p). This function is called twice in the case of the WWINP having
        #   both n and p.

        if particle == 'n':
            particle_index = 0

        elif particle == 'p':
            particle_index = 1

        # Append energy line.
        block3 = ''
        line_count = 0

        for e_upper_bound in self.e[particle_index]:
            block3 += '{0:13.5E}'.format(e_upper_bound)
            line_count += 1
            if line_count == 6:
                block3 += '\n'
                line_count = 0

        if line_count != 0:
            block3 += '\n'

        # Get ww_data.
        ww_data = np.empty(shape=(self.nft, self.ne[particle_index]))
        volume_elements = list(self.structured_iterate_hex('zyx'))
        for i, volume_element in enumerate(volume_elements):
            ww_data[i] = self.get_tag("ww_{0}".format(particle))[
                volume_element]

        for i in range(0, self.ne[particle_index]):
            # Append ww_data to block3 string.
            line_count = 0
            for ww in ww_data[:, i]:
                block3 += '{0:13.5E}'.format(ww)
                line_count += 1

                if line_count == 6:
                    block3 += '\n'
                    line_count = 0

            if line_count != 0:
                block3 += '\n'

        f.write(block3)

    def read_mesh(self, mesh):
        """This method creates a Wwinp object from a structured mesh object.
        The mesh must have tags in the form "ww_X" where X is n
        or p. For every particle there must be a rootSet tag in the form
        X_e_upper_bounds containing a list of energy upper bounds.
        """

        super(Wwinp, self).__init__(mesh=mesh, structured=True)

        # Set geometry related attributes.
        self.nr = 10
        self.nwg = 1

        # Set energy related attributes.
        self.e = []
        self.ne = []
        all_tags = [x.name for x in self.get_all_tags()]

        if 'n_e_upper_bounds' in all_tags:
            n_e_upper_bounds = self.n_e_upper_bounds[self]
            # In the single energy group case, the "E_upper_bounds" tag
            # returns a non-iterable float. If this is the case, put this
            # float into an array so that it can be iterated over
            if isinstance(n_e_upper_bounds, float):
                n_e_upper_bounds = [n_e_upper_bounds]

            self.e.append(n_e_upper_bounds)
            self.ne.append(int(len(n_e_upper_bounds)))

        else:
            self.e.append([])
            self.ne.append(0)

        if 'p_e_upper_bounds' in all_tags:
            p_e_upper_bounds = self.p_e_upper_bounds[self]
            if isinstance(p_e_upper_bounds, float):
                p_e_upper_bounds = [p_e_upper_bounds]

            self.e.append(p_e_upper_bounds)
            self.ne.append(int(len(p_e_upper_bounds)))

        self.ni = int(len(self.ne))

        # Set space related attributes.
        self.bounds = [self.structured_get_divisions('x'),
                       self.structured_get_divisions('y'),
                       self.structured_get_divisions('z')]

        self.origin = [self.bounds[0][0], self.bounds[1][0], self.bounds[2][0]]

        # Cycle through the rest of bounds to get the fine/coarse information.
        self.cm = [[], [], []]
        self.fm = [[], [], []]
        for i, points in enumerate(self.bounds):
            # There exists at least 1 fine mesh point.
            self.fm[i].append(1)
            # Loop through remaining points to determine which are coarse,
            # and the number of fine meshes between them.
            j = 1
            while j < len(points) - 1:
                # Floating point comparison characterizes coarse vs. fine.
                if abs((points[j] - points[j-1]) -
                       (points[j+1] - points[j])) <= 1.01E-4:
                    self.fm[i][len(self.cm[i])] += 1
                else:
                    self.cm[i].append(points[j])
                    self.fm[i].append(1)

                j += 1

            # Append last point as coarse point, as this is always the case.
            self.cm[i].append(points[-1])

        self.nc = [len(self.cm[0]), len(self.cm[1]), len(self.cm[2])]
        self.nf = [sum(self.fm[0]), sum(self.fm[1]), sum(self.fm[2])]
        self.nft = self.nf[0]*self.nf[1]*self.nf[2]


class Meshtal(object):
    """This class stores all the information from an MCNP meshtal file with
    single or multiple fmesh4 neutron or photon tallies. The "tally" attribute
    provides key/value access to invidial MeshTally objects.

    Attributes
    ----------
    filename : string
        Path to an MCNP meshtal file
    version : float
        The MCNP verison number
    ld : string
        The MCNP verison date
    title : string
        Title card from the MCNP input
    histories : int
        Number of histories from the MCNP simulation
    tally : dict
        A dictionary with MCNP fmesh4 tally numbers
        (e.g. 4, 14, 24) as keys and
        MeshTally objects as values.
    tags : dict
        Maps integer tally numbers to iterables containing four strs, the
        results tag name, the relative error tag name, the total results
        tag name, and the total relative error tag name. If tags is None
        the tags are named 'x_result', 'x_rel_error', 'x_result_total',
        'x_rel_error_total' where x is n or p for neutrons or photons.

    """

    def __init__(self, filename, tags=None, meshes_have_mats=False):
        """Parameters
        ----------
        filename : str
            MCNP meshtal file.
        tags : dict, optional
            Maps integer tally numbers to iterables containing four strs: the
            results tag name, the relative error tag name, the total results
            tag name, and the total relative error tag name. If tags is None
            the tags are named 'x_result', 'x_rel_error', 'x_result_total',
            'x_rel_error_total' where x is n or p for neutrons or photons.
        meshes_have_mats : bool
             If false, Meshtally objects will be created without PyNE material
             material objects.
        """

        if not HAVE_PYMOAB:
            raise RuntimeError("PyMOAB is not available, "
                               "unable to create Meshtal.")

        self.tally = {}
        self.tags = tags
        self._meshes_have_mats = meshes_have_mats

        with open(filename, 'r') as f:
            self._read_meshtal_head(f)
            self._read_tallies(f)

    def _read_meshtal_head(self, f):
        """Get the version, ld, title card and number of histories.
        """

        line_1 = f.readline()
        # set mcnp version
        self.version = line_1.split()[2]
        # get version date ("ld" in MCNP User's Manual)
        self.ld = line_1.split()[3][3:]

        line_2 = f.readline()
        # store title card
        self.title = line_2.strip()

        line_3 = f.readline()
        # get number of histories
        self.histories = int(float(line_3.split()[-1]))

    def _read_tallies(self, f):
        """Read in all of the mesh tallies from the meshtal file.
        """
        line = f.readline()

        while line != "":
            if line.split()[0:3] == ['Mesh', 'Tally', 'Number']:
                tally_num = int(line.split()[3])
                if self.tags is not None and tally_num in self.tags.keys():
                    self.tally[tally_num] = MeshTally(f, tally_num,
                                                      self.tags[tally_num],
                                                      mesh_has_mats=self._meshes_have_mats)
                else:
                    self.tally[tally_num] = MeshTally(f, tally_num,
                                                      mesh_has_mats=self._meshes_have_mats)

            line = f.readline()


class MeshTally(StatMesh):
    """This class stores all information from all single MCNP mesh tally that
    exists within some meshtal file. Header information is stored as attributes
    and the "mesh" attribute is a MOAB mesh with all result and relative error
    data tagged. This class inherits from StatMesh, exposing all statistical
    mesh manipulation methods.

    Attributes
    ----------
    tally_number : int
        The MCNP tally number. Must end in 4 (e.g. 4, 14, 214).
    particle : string
        Either "neutron" for a neutron mesh tally or "photon" for a photon mesh
        tally.
    dose_response : bool
        True is the tally is modified by a dose response function.
    x_bounds : list of floats
        The locations of mesh vertices in the x direction.
    y_bounds : list of floats
        The locations of mesh vertices in the y direction.
    z_bounds : list of floats
        The locations of mesh vertices in the z direction.
    e_bounds : list of floats
        The minimum and maximum bounds for energy bins
    mesh :
        An PyMOAB core instance tagged with all results and
        relative errors
    tag_names : iterable
        Four strs that specify the tag names for the results, relative errors,
        total results, and relative errors of the total results.

    Notes
    -----
    All Mesh/StatMesh attributes are also present via a super() call to
    StatMesh.__init__().

    """

    def __init__(self, f, tally_number, tag_names=None, mesh_has_mats=False):
        """Create MeshTally object from a filestream open to the second
        line of a mesh tally header (the neutron/photon line). MeshTally objects
        should be instantiated only through the Meshtal class.

        Parameters
        ----------
        f : filestream
            Open to the neutron/photon line.
        tally_number : int
            The MCNP fmesh4 tally number (e.g. 4, 14, 24).
        tag_names : iterable, optional
            Four strs that specify the tag names for the results, relative
            errors, total results and relative errors of the total results.
            This should come from the Meshtal.tags attribute dict.
        mesh_has_mats : bool
             If false, Meshtally objects will be created without PyNE material
             objects.
        """

        if not HAVE_PYMOAB:
            raise RuntimeError("PyMOAB is not available, "
                               "unable to create Meshtally Mesh.")

        self.tally_number = tally_number
        self._read_meshtally_head(f)
        self._read_column_order(f)

        if tag_names is None:
            self.tag_names = ("{0}_result".format(self.particle),
                              "{0}_result_rel_error".format(self.particle),
                              "{0}_result_total".format(self.particle),
                              "{0}_result_total_rel_error".format(self.particle))
        else:
            self.tag_names = tag_names

        self._create_mesh(f, mesh_has_mats)

    def _read_meshtally_head(self, f):
        """Get the particle type, spacial and energy bounds, and whether or
        not flux-to-dose conversion factors are being used.
        """
        line = f.readline()
        if ('neutron' in line):
            self.particle = 'neutron'
        elif ('photon' in line):
            self.particle = 'photon'

        # determine if meshtally flux-to-dose conversion factors are being used.
        line = f.readline()
        dr_str = 'This mesh tally is modified by a dose response function.'
        if line.strip() == dr_str:
            self.dose_response = True
        else:
            self.dose_response = False

        # advance the file to the line where x, y, z, bounds start
        while line.strip() != 'Tally bin boundaries:':
            line = f.readline()

        self.x_bounds = [float(x) for x in f.readline().split()[2:]]
        self.y_bounds = [float(x) for x in f.readline().split()[2:]]
        self.z_bounds = [float(x) for x in f.readline().split()[2:]]
        # "Energy bin boundaries" contain one more word than "X boundaries"
        self.e_bounds = [float(x) for x in f.readline().split()[3:]]

        self.dims = [0, 0, 0] + [len(self.x_bounds) - 1,
                                 len(self.y_bounds) - 1,
                                 len(self.z_bounds) - 1]

        # skip blank line between enery bin boundaries and table headings
        f.readline()

    def _read_column_order(self, f):
        """Create dictionary with table headings as keys and their column
        location as values. Dictionary is the private attribute _column_idx.
        """
        line = f.readline()
        column_names = line.replace('Rel ', 'Rel_').replace(
            'Rslt * ', 'Rslt_*_').strip().split()
        self._column_idx = dict(zip(column_names, range(0, len(column_names))))

    def _create_mesh(self, f, mesh_has_mats):
        """Instantiate a Mesh object and tag the PyMOAB core instance
           with results and relative errors.
        """

        mats = () if mesh_has_mats is True else None
        super(MeshTally, self).__init__(structured_coords=[self.x_bounds,
                                                           self.y_bounds, self.z_bounds],
                                        structured=True, mats=mats)

        num_vol_elements = (len(self.x_bounds)-1) * (len(self.y_bounds)-1)\
            * (len(self.z_bounds)-1)
        num_e_groups = len(self.e_bounds)-1

        # get result and relative error data from file
        result = np.empty(shape=(num_e_groups, num_vol_elements))
        rel_error = np.empty(shape=(num_e_groups, num_vol_elements))
        for i in range(0, num_e_groups):
            result_row = []
            rel_error_row = []
            for j in range(0, num_vol_elements):
                line = f.readline().split()
                result_row.append(float(line[self._column_idx["Result"]]))
                rel_error_row.append(
                    float(line[self._column_idx["Rel_Error"]]))

            result[i] = result_row
            rel_error[i] = rel_error_row

        # Tag results and error vector to mesh
        self.tag(self.tag_names[0], tagtype='nat_mesh',
                 size=num_e_groups, dtype=float)
        res_tag = self.get_tag(self.tag_names[0])
        self.tag(self.tag_names[1], tagtype='nat_mesh',
                 size=num_e_groups, dtype=float)
        rel_err_tag = self.get_tag(self.tag_names[1])

        if num_e_groups == 1:
            res_tag[:] = result[0]
            rel_err_tag[:] = rel_error[0]
        else:
            res_tag[:] = result.transpose()
            rel_err_tag[:] = rel_error.transpose()

        # If "total" data exists (i.e. if there is more than
        # 1 energy group) get it and tag it onto the mesh.
        if num_e_groups > 1:
            result = []
            rel_error = []
            for i in range(0, num_vol_elements):
                line = f.readline().split()
                result.append(float(line[self._column_idx["Result"]]))
                rel_error.append(
                    float(line[self._column_idx["Rel_Error"]]))

            self.tag(self.tag_names[2], size=1,
                     dtype=float, tagtype='nat_mesh')
            res_tot_tag = self.get_tag(self.tag_names[2])

            self.tag(self.tag_names[3], size=1,
                     dtype=float, tagtype='nat_mesh')
            rel_err_tot_tag = self.get_tag(self.tag_names[3])

            res_tot_tag[:] = result
            rel_err_tot_tag[:] = rel_error


def mesh_to_geom(mesh, frac_type='mass', title_card="Generated from PyNE Mesh"):
    """This function reads a structured Mesh object and returns the geometry
    portion of an MCNP input file (cells, surfaces, materials), prepended by a
    title card. The mesh must be axis aligned. Surfaces and cells are written
    in xyz iteration order (z changing fastest).

    Parameters
    ----------
    mesh : PyNE Mesh object
        A structured Mesh object with materials and valid densities.
    frac_type : str, optional
        Either 'mass' or 'atom'. The type of fraction to use for the material
        definition.
    title_card : str, optional
        The MCNP title card to appear at the top of the input file.

    Returns
    -------
    geom : str
        The title, cell, surface, and material cards of an MCNP input file in
        the proper order.

    """
    mesh._structured_check()
    divs = (mesh.structured_get_divisions('x'),
            mesh.structured_get_divisions('y'),
            mesh.structured_get_divisions('z'))

    cell_cards = _mesh_to_cell_cards(mesh, divs)
    surf_cards = _mesh_to_surf_cards(mesh, divs)
    mat_cards = _mesh_to_mat_cards(mesh, divs, frac_type)

    return "{0}\n{1}\n{2}\n{3}".format(title_card, cell_cards,
                                       surf_cards, mat_cards)


def _mesh_to_cell_cards(mesh, divs):
    """Prepares the cell cards for mesh_to_geom."""
    cell_cards = ""
    count = 1
    idx = mesh.iter_structured_idx('xyz')

    # Establish min and max idx values for each dimension.
    x_min = 1
    x_max = len(divs[0])
    y_min = x_max + 1
    y_max = x_max + len(divs[1])
    z_min = y_max + 1
    z_max = y_max + len(divs[2])

    for i in range(1, len(divs[0])):
        for j in range(1, len(divs[1])):
            for k in range(1, len(divs[2])):
                # Cell number, mat number, density
                cell_cards += "{0} {1} {2} ".format(count, count,
                                                    mesh.density[idx.next()])
                # x, y, and z surfaces
                cell_cards += "{0} -{1} {2} -{3} {4} -{5}\n".format(
                              i, i + 1, j + x_max, j + x_max + 1,
                              k + y_max, k + y_max + 1)
                count += 1

    # Append graveyard.
    cell_cards += "{0} 0 -{1}:{2}:-{3}:{4}:-{5}:{6}\n".format(
                  count, x_min, x_max, y_min, y_max, z_min, z_max)

    return cell_cards


def _mesh_to_surf_cards(mesh, divs):
    """Prepares the surface cards for mesh_to_geom."""
    surf_cards = ""
    count = 1
    for i, dim in enumerate("xyz"):
        for div in divs[i]:
            surf_cards += "{0} p{1} {2}\n".format(count, dim, div)
            count += 1

    return surf_cards


def _mesh_to_mat_cards(mesh, divs, frac_type):
    """Prepares the material cards for mesh_to_geom."""
    mat_cards = ""
    idx = mesh.iter_structured_idx('xyz')
    for i in idx:
        mesh.mats[i].metadata['mat_number'] = i + 1
        mat_cards += mesh.mats[i].mcnp(frac_type=frac_type)

    return mat_cards
