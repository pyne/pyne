"""Module for parsing and manipulating data from ENDF evaluations. Currently, it
only can read several MTs from File 1, but with time it will be expanded to
include the entire ENDF format.

All the classes and functions in this module are based on document
ENDF-102 titled "Data Formats and Procedures for the Evaluated Nuclear
Data File ENDF-6". The latest version from June 2009 can be found at
http://www-nds.iaea.org/ndspub/documents/endf/endf102/endf102.pdf

For more information on the Evaluation class, contact Paul Romano
<paul.k.romano@gmail.com>. For more information on the Library class, contact
John Xia <john.danger.xia@gmail.com>.
"""
from __future__ import print_function, division, unicode_literals

import re
import os
try:
    from collections.abc import OrderedDict, Iterable
except ImportError:
    from collections import OrderedDict, Iterable
from warnings import warn
from pyne.utils import QA_warn

cimport numpy as np
import numpy as np
from numpy.polynomial.polynomial import Polynomial
from numpy.polynomial.legendre import Legendre
from scipy.interpolate import interp1d
cimport cython

from pyne cimport cpp_nucname
from pyne import nucname
from pyne import rxdata
from pyne.rxname import label
from pyne.utils import fromendf_tok, endftod

np.import_array()

QA_warn(__name__)

libraries = {0: 'ENDF/B', 1: 'ENDF/A', 2: 'JEFF', 3: 'EFF',
             4: 'ENDF/B High Energy', 5: 'CENDL', 6: 'JENDL',
             31: 'INDL/V', 32: 'INDL/A', 33: 'FENDL', 34: 'IRDF',
             35: 'BROND', 36: 'INGDB-90', 37: 'FENDL/A', 41: 'BROND'}
FILE1END = r'([1-9]\d{3}| [1-9]\d{2}|  [1-9]\d|   [1-9]) 1451(?= *[1-9]\d*$)[ \d]{5}$'
FILE1_R = re.compile(r'^.{66}'+FILE1END)
SPACEINT11_R = re.compile('^((?= *-?[1-9]\d*$)[ \d-]{11}| {10}0)$')      # I11 (FORTRAN77)

def _radiation_type(value):
    p = {0: 'gamma', 1: 'beta-', 2: 'ec/beta+', 3: 'IT',
         4: 'alpha', 5: 'neutron', 6: 'sf', 7: 'proton',
         8: 'e-', 9: 'xray', 10: 'unknown'}
    if value % 1.0 == 0:
        return p[int(value)]
    else:
        return (p[int(value)], p[int(10*value % 10)])


class Library(rxdata.RxLib):
    """A class for a file which contains multiple ENDF evaluations."""
    def __init__(self, fh):
        self.mts = {}
        self.structure = {}
        self.mat_dict = {}
        self.more_files = True
        self.intdict = {1: self._histogram, 2: self._linlin, 3: self._linlog, 4:
                        self._loglin, 5: self._loglog, 6:
                        self._chargedparticles, 11: self._histogram,
                        12: self._linlin, 13: self._linlog, 14: self._loglin,
                        15: self._loglog, 21: self._histogram, 22: self._linlin,
                        23: self._linlog, 24: self._loglin, 25: self._loglog}
        self.chars_til_now = 0 # offset (byte) from the top of the file for seek()ing
        self.fh = fh
        self._set_line_length()
        # read first line (Tape ID)
        self._read_tpid()
        # read headers for all materials
        while self.more_files:
            self._read_headers()

    def _set_line_length(self):
        opened_here = False
        if isinstance(self.fh, str):
            fh = open(self.fh, 'r')
            opened_here = True
        else:
            fh = self.fh

        # Make sure the newlines attribute is set: read a couple of lines
        # and rewind to the beginning.
        # Yes, we need to read two lines to make sure that fh.newlines gets
        # set. One line seems to be enough for *nix-style line terminators, but
        # two seem to be necessary for Windows-style terminators.
        fh.seek(0)
        fh.readline()
        line = fh.readline()
        self.line_length = len(line) # actual chars/line read
        # self.offset now stores the diff. between the length of a line read
        # and the length of a line in the ENDF-6 formatted file.
        self.offset = len(fh.newlines) - (self.line_length - 80)
        fh.seek(0)

        if opened_here:
            fh.close()

    def load(self):
        """load()
        Read the ENDF file into a NumPy array.

        Returns
        -------
        data : np.array, 1d, float64
            Returns a 1d float64 NumPy array.
        """
        opened_here = False
        if isinstance(self.fh, basestring):
            fh = open(self.fh, 'r')
            opened_here = True
        else:
            fh = self.fh
        fh.readline()
        data = fromendf_tok(fh.read())
        fh.seek(0)
        if opened_here:
            fh.close()
        return data

    def _read_tpid(self):
        if self.chars_til_now == 0:
            opened_here = False
            if isinstance(self.fh, basestring):
                fh = open(self.fh, 'r')
                opened_here = True
            else:
                fh = self.fh
            line = fh.readline()
            self.chars_til_now = len(line) + self.offset
        else:
            warn('TPID is the first line, has been read already', UserWarning)

    def _isContentLine(self,parts):
        """Check whether a line is consisted of 22*spaces and 4*(I11)s (FORTRAN77).

        Parameters
        -----------
        parts: list
            made by dividing 1-66 chars of an input line into 6 parts of equal length
        """
        return parts[0]+parts[1]==' '*22 and \
        SPACEINT11_R.match(parts[2]) and SPACEINT11_R.match(parts[3]) and \
        SPACEINT11_R.match(parts[4]) and SPACEINT11_R.match(parts[5])

    def _read_headers(self):
        cdef int nuc
        cdef int mat_id
        cdef double nucd
        opened_here = False
        if isinstance(self.fh, basestring):
            fh = open(self.fh, 'r')
            opened_here = True
        else:
            fh = self.fh
        # Go to current file position
        fh.seek(self.chars_til_now)
        # get mat_id
        line = fh.readline()
        mat_id = int(line[66:70].strip() or -1)
        # check for isomer (LIS0/LISO entry)
        matflagstring = line + fh.read(3*self.line_length)
        flagkeys = ['ZA', 'AWR', 'LRP', 'LFI', 'NLIB', 'NMOD', 'ELIS',
                    'STA', 'LIS', 'LIS0', 0, 'NFOR', 'AWI', 'EMAX',
                    'LREL', 0, 'NSUB', 'NVER', 'TEMP', 0, 'LDRV',
                    0, 'NWD', 'NXC']
        flags = dict(zip(flagkeys, fromendf_tok(matflagstring)))
        nuc = cpp_nucname.id(<int> (<int> flags['ZA']*10000 + flags['LIS0']))
        # Make a new dict in self.structure to contain the material data.
        if nuc not in self.structure:
            self.structure.update(
                {nuc: {'styles': '', 'docs': [], 'particles': [], 'data': {},
                       'matflags': {}}})
            self.mat_dict.update({nuc: {'end_line': [],
                                        'mfs': {}}})
        # Parse header (all lines with 1451)
        mf = 1
        start = self.chars_til_now//(self.line_length+self.offset) # present (the first) line number
        stop = start  # if no 451 can be found
        line = fh.readline() # get the next line; start parsing from the 6th line
        while FILE1_R.search(line):
            # divide 1-66 chars of the line into six 11-char parts
            lineparts = [line[i:i+11] for i in range(0, 66, 11)]
            # parse contents section
            if self._isContentLine(lineparts):
                # When MF and MT change, add offset due to SEND/FEND records.
                old_mf = mf
                mf, mt = int(lineparts[2]), int(lineparts[3])
                if old_mf != mf:
                    start += 1
                mt_length = int(lineparts[4])
                stop = start + mt_length
                # The first number in the tuple is the offset in the file to seek(),
                # whereas the second stands for the number of characters to be read().
                self.mat_dict[nuc]['mfs'][mf, mt] = ((self.line_length+self.offset)*start,
                                                     self.line_length*stop)
                start = stop + 1
                line = fh.readline()
            # parse comment
            else:
                self.structure[nuc]['docs'].append(line[0:66])
                line = fh.readline()
        # Find where the end of the material is and then jump to it.
        # The end is 3 lines after the last mf,mt
        # combination (SEND, FEND, MEND)
        self.chars_til_now = (stop + 3)*(self.line_length+self.offset) # at the end of a MAT
        fh.seek(self.chars_til_now)
        nextline = fh.readline()
        self.more_files = (nextline != '' and nextline[68:70] != '-1')
        # Update materials dict
        if mat_id != -1:
            self.mat_dict[nuc]['end_line'] = \
                (self.chars_til_now+self.offset)//self.line_length
            setattr(self, 'mat{0}'.format(nuc), self.structure[nuc])
        self._read_mat_flags(nuc)
        fh.seek(0)
        if opened_here:
            fh.close()

    def _read_mat_flags(self, nuc):
        """Reads the global flags for a certain material.

        Parameters
        -----------
        nuc: int
            ZZAAAM of material.
        """
        mf1 = self.get_rx(nuc, 1, 451, lines=4)
        flagkeys = ['ZA', 'AWR', 'LRP', 'LFI', 'NLIB', 'NMOD', 'ELIS',
                    'STA', 'LIS', 'LIS0', 0, 'NFOR', 'AWI', 'EMAX',
                    'LREL', 0, 'NSUB', 'NVER', 'TEMP', 0, 'LDRV',
                    0, 'NWD', 'NXC']
        flags = dict(zip(flagkeys, mf1[:12]))
        del flags[0]
        self.structure[nuc]['matflags'] = flags

    def _get_cont(self, keys, line):
        """Read one line of the array, treating it as a CONT record.

        Parameters
        -----------
        keys: iterable
            An iterable containing the labels for each field in the CONT record.
            For empty/unassigned fields, use 0.
        line: array-like
            The line to be read.

        Returns
        --------
        cont : dict
            Contains labels and values mapped to each other.
        """
        cont = dict(zip(keys, line.flat[:6]))
        if 0 in cont:
            del cont[0]
        return cont

    def _get_head(self, keys, line):
        """Read one line of the array, treating it as a HEAD record.

        Parameters
        -----------
        keys: iterable
            An iterable containing the labels for each field in the HEAD record.
            For empty/unassigned fields, use 0.
        line: array-like
            The line to be read.

        Returns
        --------
        cont : dict
            Contains labels and values mapped to each other.
        """
        # Just calls self._get_cont because HEAD is just a special case of CONT
        if (keys[0] == 'ZA' and keys[1] == 'AWR'):
            return self._get_cont(keys, line)
        else:
            raise ValueError('This is not a HEAD record: {}'.format(
                dict(zip(keys, line))))

    def _get_list(self, headkeys, itemkeys, lines):
        """Read some lines of the array, treating it as a LIST record.

        Parameters
        -----------
        headkeys: iterable
            An iterable containing the labels for each field in the first
            record. For empty/unassigned fields, use 0.
        itemkeys: iterable
            An iterable containing the labels for each field in the next
            records. For empty/unassigned fields, use 0. If itemkeys has length
            1, the array is flattened and assigned to that key.
        lines: two-dimensional array-like
            The lines to be read. Each line should have 6 elements. The first
            line should be the first line of the LIST record; since we don't
            know the length of the LIST record, the last line should be the last
            line it is plausible for the LIST record to end.

        Returns
        --------
        head: dict
            Contains elements of the first line paired with their labels.
        items: dict
            Contains columns of the LIST array paired with their labels, unless
            itemkeys has length 1, in which case items contains the flattened
            LIST array paired with its label.
        total_lines: int
            The number of lines the LIST record takes up.
        """

        head = dict(zip(headkeys, lines[0:].flat[:len(headkeys)]))
        if 0 in head:
            del head[0]
        npl = int(lines[0][4])
        headlines = (len(headkeys)-1)//6 + 1
        arraylines = (npl-1)//6 + 1
        if len(itemkeys) == 1:
            array_len = npl - (headlines-1) * 6
            items = {itemkeys[0]: lines[headlines:].flat[:array_len]}
        else:
            array_width = ((len(itemkeys)-1)//6 + 1)*6
            items_transposed = np.transpose(
                lines[headlines:headlines+arraylines].reshape(-1,
                                                              array_width))
            items = dict(zip(itemkeys, items_transposed))
        if 0 in items:
            del items[0]

        total_lines = 1+arraylines
        return head, items, total_lines

    def _get_tab1(self, headkeys, xykeys, lines):
        """Read some lines of the array, treating it as a TAB1 record.

        Parameters
        -----------
        headkeys: iterable, length 6
            An iterable containing the labels for each field in the first
            line. For empty/unassigned fields, use 0.
        xykeys: iterable, length 2
            An iterable containing the labels for the interpolation data. The
            first key should be xint, the second should be y(x).
        lines: two-dimensional array-like
            The lines to be read. Each line should have 6 elements. The first
            line should be the first line of the TAB1 record; since we don't
            know the length of the TAB1 record, the last line should be the last
            line it is plausible for the TAB1 record to end.

        Returns
        --------
        head: dict
            Contains elements of the first card paired with their labels.
        intdata: dict
            Contains the interpolation data.
        total_lines: int
            The number of lines the TAB1 record takes up.
        """
        head = dict(zip(headkeys, lines[0]))
        if 0 in head:
            del head[0]
        nr, np_ = int(lines[0][4]), int(lines[0][5])
        meta_len = (nr*2-1)//6 + 1
        data_len = (np_*2-1)//6 + 1
        intmeta = dict(zip(('intpoints', 'intschemes'),
                           (np.asarray(lines[1:1+meta_len].flat[:nr*2:2], dtype=int),
                            lines[1:1+meta_len].flat[1:nr*2:2])))
        intdata = dict(zip(xykeys,
            (lines[1+meta_len:1+meta_len+data_len].flat[:np_*2:2],
             lines[1+meta_len:1+meta_len+data_len].flat[1:np_*2:2])))
        intdata.update(intmeta)
        total_lines = 1 + meta_len + data_len
        return head, intdata, total_lines

    def _histogram(self, e_int, xs, low, high):
        if low in e_int:
            # truncate at lower bound
            xs = xs[e_int >= low]
            e_int = e_int[e_int >= low]
        elif low is not None and low > e_int[0]:
            # truncate at lower bound and prepend interpolated endpoint
            low_xs = xs[e_int < low][-1]
            xs = np.insert(xs[e_int > low], 0, low_xs)
            e_int = np.insert(e_int[e_int > low], 0, low)
        if high in e_int:
            # truncate at higher bound
            xs = xs[e_int <= high]
            e_int = e_int[e_int <= high]
        elif high is not None:
            # truncate at higher bound and prepend interpolated endpoint
            high_xs = xs[e_int < high][-1]
            xs = np.append(xs[e_int < high], high_xs)
            e_int = np.append(e_int[e_int < high], high)
        de_int = float(e_int[-1]-e_int[0])
        return np.nansum((e_int[1:]-e_int[:-1]) * xs[:-1]/de_int)

    def _linlin(self, e_int, xs, low, high):
        if low is not None or high is not None:
            interp = interp1d(e_int, xs)
            if low in e_int:
                xs = xs[e_int >= low]
                e_int = e_int[e_int >= low]
            elif low is not None and low > e_int[0]:
                low_xs = interp(low)
                xs = np.insert(xs[e_int > low], 0, low_xs)
                e_int = np.insert(e_int[e_int > low], 0, low)
            if high in e_int:
                xs = xs[e_int <= high]
                e_int = e_int[e_int <= high]
            elif high is not None:
                high_xs = interp(high)
                xs = np.append(xs[e_int < high], high_xs)
                e_int = np.append(e_int[e_int < high], high)
        de_int = float(e_int[-1]-e_int[0])
        return np.nansum((e_int[1:]-e_int[:-1]) * (xs[1:]+xs[:-1])/2./de_int)

    def _linlog(self, e_int, xs, low, high):
        if low is not None or high is not None:
            interp = interp1d(np.log(e_int), xs)
            if low in e_int:
                xs = xs[e_int >= low]
                e_int = e_int[e_int >= low]
            elif low is not None and low > e_int[0]:
                low_xs = interp(np.log(low))
                xs = np.insert(xs[e_int > low], 0, low_xs)
                e_int = np.insert(e_int[e_int > low], 0, low)
            if high in e_int:
                xs = xs[e_int <= high]
                e_int = e_int[e_int <= high]
            elif high is not None:
                high_xs = interp(np.log(high))
                xs = np.append(xs[e_int < high], high_xs)
                e_int = np.append(e_int[e_int < high], high)

        de_int = float(e_int[-1]-e_int[0])
        x1 = e_int[:-1]
        x2 = e_int[1:]
        y1 = xs[:-1]
        y2 = xs[1:]
        A = (y1-y2)/(np.log(x1/x2))
        B = y1-A*np.log(x1)
        return np.nansum(A*(x2*np.log(x2) -
                            x1*np.log(x1)-x2+x1) + B*(x2-x1))/de_int

    def _loglin(self, e_int, xs, low, high):
        if low is not None or high is not None:
            interp = interp1d(e_int, np.log(xs))
            if low in e_int:
                xs = xs[e_int >= low]
                e_int = e_int[e_int >= low]
            elif low is not None and low > e_int[0]:
                low_xs = np.e ** interp(low)
                xs = np.insert(xs[e_int > low], 0, low_xs)
                e_int = np.insert(e_int[e_int > low], 0, low)
            if high in e_int:
                xs = xs[e_int <= high]
                e_int = e_int[e_int <= high]
            elif high is not None:
                high_xs = np.e ** interp(high)
                xs = np.append(xs[e_int < high], high_xs)
                e_int = np.append(e_int[e_int < high], high)

        de_int = float(e_int[-1]-e_int[0])
        x1 = e_int[:-1]
        x2 = e_int[1:]
        y1 = xs[:-1]
        y2 = xs[1:]
        A = (np.log(y1)-np.log(y2))/(x1-x2)
        B = np.log(y1) - A*x1
        return np.nansum((y2-y1)/A)/de_int

    def _loglog(self, e_int, xs, low, high):
        if low is not None or high is not None:
            interp = interp1d(np.log(e_int), np.log(xs))
            if low in e_int:
                xs = xs[e_int >= low]
                e_int = e_int[e_int >= low]
            elif low is not None and low > e_int[0]:
                low_xs = np.e ** interp(np.log(low))
                xs = np.insert(xs[e_int > low], 0, low_xs)
                e_int = np.insert(e_int[e_int > low], 0, low)
            if high in e_int:
                xs = xs[e_int <= high]
                e_int = e_int[e_int <= high]
            elif high is not None:
                high_xs = np.e ** interp(np.log(high))
                xs = np.append(xs[e_int < high], high_xs)
                e_int = np.append(e_int[e_int < high], high)

        de_int = float(e_int[-1]-e_int[0])
        x1 = e_int[:-1]
        x2 = e_int[1:]
        y1 = xs[:-1]
        y2 = xs[1:]
        A = - np.log(y2/y1)/np.log(x1/x2)
        B = - (np.log(y1)*np.log(x2) - np.log(y2)*np.log(x1))/np.log(x1/x2)
        return np.nansum(np.e**B / (A+1) * (x2**(A+1) - x1**(A+1))/de_int)

    def _chargedparticles(self, e_int, xs, flags=None):
        q = flags['Q']
        if q > 0:
            T = 0
        else:
            T = q
        de_int = float(e_int[-1]-e_int[0])
        x1 = e_int[:-1]
        x2 = e_int[1:]
        y1 = xs[:-1]
        y2 = xs[1:]
        B = np.log(y2*x2/(x1*y1)) / (1/(x1-T)**0.5 - 1/(x2-T)**0.5)
        A = np.e**(B/(x1-T)**0.5)*y1*x1
        # FIXME
        raise NotImplementedError('see docs for more details.')

    def integrate_tab_range(self, intscheme, e_int, xs, low=None, high=None):
        """integrate_tab_range(intscheme, e_int, xs, low=None, high=None)
        Integrates across one tabulation range.

        Parameters
        ----------
        intscheme : int or float
            The interpolation scheme used in this range.
        e_int : array
            The energies at which we have xs data.
        xs : array
            The xs data corresponding to e_int.
        low, high : float
            Lower and upper bounds within the tabulation range to start/stop at.

        Returns
        -------
        sigma_g : float
            The group xs.
        """
        with np.errstate(divide='ignore', invalid='ignore'):
            # each of these functions returns a normalized integration
            # over the range
            return self.intdict[intscheme](e_int, xs, low, high)

    def _cont_and_update(self, flags, keys, data, total_lines):
        flags.update(self._get_cont(keys, data[total_lines]))
        return flags, total_lines+1

    def _nls_njs_loop(self, L_keys, j_keys, itemkeys, data, total_lines,
                      range_flags, subsection_dict):
        nls = int(range_flags['NLS'])
        for nls_iter in range(nls):
            if j_keys is None:
                L_flags, items, lines = self._get_list(
                    L_keys, itemkeys, data[total_lines:])
                total_lines += lines
                spi, L = range_flags['SPI'], L_flags['L']
                subsection_dict[spi, L] = items
            else:
                L_flags = self._get_cont(L_keys, data[total_lines])
                total_lines += 1
                njs = int(L_flags['NJS'])
                for njs_iter in range(njs):
                    j_flags, items, lines = self._get_list(
                        j_keys, itemkeys, data[total_lines:])
                    total_lines += lines
                    items.update(j_flags)
                    spi, L, aj = range_flags['SPI'], L_flags['L'], j_flags['AJ']
                    subsection_dict[(spi, L, aj)] = items
        return total_lines

    def _read_res(self, mat_id):
        """_read_res(mat_id)
        Read the resonance data from one material in the library and updates
        self.structure.

        Parameters
        -----------
        mat_id: int
            Material id .
        """
        lrp = self.structure[mat_id]['matflags']['LRP']
        if (lrp == -1 or mat_id in (-1, 0)):
            # If the LRP flag for the material is -1,
            # there's no resonance data.
            # Also if the mat id is invalid.
            #
            # However other methods expects _read_res to set
            # structur[nuc]['data'], so fill it with a single
            # entry for mat_id and empty values:
            self.structure[mat_id]['data'].update(
                {mat_id: {'resolved': [],
                          'unresolved': [],
                          'datadocs': [],
                          'xs': {},
                          'output': {'channel1': [],
                                     'channel2': []},
                          'isotope_flags': {}}})
            pass
        else:
            # Load the resonance data.
            mf2 = self.get_rx(mat_id, 2, 151).reshape(-1, 6)

            self.structure[mat_id]['matflags'].update(
                self._get_head(['ZA', 'AWR', 0, 0, 'NIS', 0], mf2[0]))
            total_lines = 1
            for isotope_num in range(
                    int(self.structure[mat_id]['matflags']['NIS'])):
                total_lines += self._read_nis(mf2[total_lines:], lrp, mat_id)
        for isotope in self.structure[mat_id]['data'].values():
            isotope['resolved'].sort()
            isotope['unresolved'].sort()

    def _read_nis(self, isotope_data, lrp, mat_id):
        """_read_nis(isotope_data, lrp, mat_id)
        Read resonance data for a specific isotope.

        Parameters
        -----------
        isotope_data: 2D array
            The section of the resonance data to read. The isotope starts at the
            top of this.
        lrp: int
            A flag denoting the type of data in the isotope. Exact meaning of
            this flag can be found in ENDF Manual pp.50-51.
        mat_id: int
            Material ZZAAAM.

        Returns
        --------
        total_lines: int
            The number of lines the isotope takes up.
        """
        isotope_flags = self._get_cont(['ZAI', 'ABN', 0, 'LFW', 'NER', 0],
                                       isotope_data[0])
        # according to endf manual, there is no specification
        # for metastable states in ZAI
        # if we have a LIS0 != 0 we add the state to all isotopes
        if self.structure[mat_id]['matflags']['LIS0'] == 0:
            nuc_i = nucname.id(int(isotope_flags['ZAI']*10))
        else:
            nuc_i = nucname.id(int(isotope_flags['ZAI']*10 +
                               self.structure[mat_id]['matflags']['LIS0']))

        self.structure[mat_id]['data'].update(
            {nuc_i: {'resolved': [],
                     'unresolved': [],
                     'datadocs': [],
                     'xs': {},
                     'output': {'channel1': [],
                                'channel2': []},
                     'isotope_flags': isotope_flags}})
        total_lines = 1
        for er in range(int(isotope_flags['NER'])):
            total_lines += self._read_subsection(isotope_data[total_lines:],
                                                 isotope_flags,
                                                 mat_id,
                                                 nuc_i)

        return total_lines

    def _read_subsection(self, subsection, isotope_flags, mat_id, nuc_i):
        """Read resonance data for a specific energy range subsection.

        Parameters
        -----------
        subsection: 2D array
            The section of the resonance data to read. The energy range
            subsection starts at the top of this.
        range_flags: dict
            Dictionary of flags inherited from the range.
        isotope_flags: dict
            Dictionary of flags inherited from the isotope.
        mat_id: int
            Material ZZAAAM.
        nuc_i: int
            Isotope ZZAAAM.

        Returns
        --------
        total_lines: int
            The number of lines the energy range subsection takes up.
        """
        range_flags = self._get_cont(('EL', 'EH', 'LRU', 'LRF', 'NRO', 'NAPS'),
                                     subsection[0])
        total_lines = 1
        lru = int(round(range_flags['LRU']))
        lru_list = [self._read_ap_only, self._read_resolved,
                    self._read_unresolved]
        total_lines += lru_list[lru](subsection[1:],
                                     range_flags,
                                     isotope_flags,
                                     mat_id,
                                     nuc_i)
        return total_lines

    def _read_resolved(self, subsection, range_flags, isotope_flags, mat_id,
                       nuc_i):
        """Read the subsection for a resolved energy range.

        Parameters
        -----------
        subsection: 2D array
            The section of the resonance data to read. The energy range
            subsection starts at the top of this.
        range_flags: dict
            Dictionary of flags inherited from the range.
        isotope_flags: dict
            Dictionary of flags inherited from the isotope.
        mat_id: int
            ZZAAAM of the material.
        nuc_i: int
            ZZAAAM of the isotope.

        Returns
        --------
        total_lines: int
            The number of lines taken up by the subsection.
        """
        def read_kbks(nch, subsection, aj_data, total_lines):
            for ch in range(nch):
                lbk = int(subsection[total_lines][4])
                lbk_list_keys = {2: ('R0', 'R1', 'R2', 'S0', 'S1', 0),
                                 3: ('R0', 'SO', 'GA', 0, 0, 0)}
                aj_data['ch{}'.format(ch)] = {'LBK': lbk}
                ch_data = aj_data['ch{}'.format(ch)]
                if lbk == 0:
                    total_lines += 2
                elif lbk == 1:
                    total_lines += 2
                    rbr, rbr_size = self._get_tab1(
                        (0, 0, 0, 0, 'NR', 'NP'), ('e_int', 'RBR'),
                        subsection[total_lines:])[1:3]
                    total_lines += rbr_size
                    ch_data['RBR'] = rbr
                    rbi, rbi_size = self._get_tab1(
                        (0, 0, 0, 0, 'NR', 'NP'), ('e_int', 'RBI'),
                        (subsection[total_lines:]))[1:3]
                    total_lines += rbi_size
                    ch_data['RBI'] = rbi
                else:
                    ch_data, total_lines = self._cont_and_update(
                        ch_data, ('ED', 'EU', 0, 0, 'LBK', 0), subsection,
                        total_lines)
                    ch_data, total_lines = self._cont_and_update(
                        ch_data, lbk_list_keys[lbk], subsection,
                        total_lines)
            return total_lines

        def read_kpss(nch, subsection, aj_data, total_lines):
            for ch in range(nch):
                ch_data = aj_data['ch{}'.format(ch)]
                lps = subsection[total_lines][4]
                ch_data['LPS'] = lps
                total_lines += 2
                if lps == 1:
                    psr, psr_size = self._get_tab1(
                        (0, 0, 0, 0, 'NR', 'NP'), ('e_int', 'PSR'),
                        subsection[total_lines:])[1:3]
                    total_lines += psr_size
                    ch_data['PSR'] = psr
                    psi, psi_size = self._get_tab1(
                        (0, 0, 0, 0, 'NR', 'NP'), ('e_int', 'PSI'),
                        (subsection[total_lines:]))[1:3]
                    total_lines += psi_size
                    ch_data['PSI'] = psi
                    total_lines += psi_size
            return total_lines

        lrf = int(range_flags['LRF'])
        subsection_dict = rxdata.DoubleSpinDict({})
        headers = [None,
                   ('SPI', 'AP', 0, 0, 'NLS', 0),
                   ('SPI', 'AP', 0, 0, 'NLS', 0),
                   ('SPI', 'AP', 'LAD', 0, 'NLS', 'NLSC'),
                   ('SPI', 'AP', 0, 0, 'NLS', 0),
                   None,
                   None,
                   (0, 0, 'IFG', 'KRM', 'NJS', 'KRL')]
        if range_flags['NRO'] > 0:
            intdata, total_lines = self._get_tab1((0, 0, 0, 0, 'NR', 'NP'),
                                                  ('E', 'AP'),
                                                  subsection)[1:3]
            subsection_dict['int'] = intdata
        else:
            total_lines = 0
        range_flags, total_lines = self._cont_and_update(
            range_flags, headers[lrf], subsection, total_lines)

        lrf_L_keys = [None,
                      ('AWRI', 'QX', 'L', 'LRX', '6*NRS', 'NRS'),
                      ('AWRI', 'QX', 'L', 'LRX', '6*NRS', 'NRS'),
                      ('AWRI', 'APL', 'L', 0, '6*NRS', 'NRS'),
                      (0, 0, 'L', 0, 'NJS', 0)]
        lrf_J_keys = [None, None, None, None, ('AJ', 0, 0, 0, '12*NLJ', 'NLJ')]
        lrf_itemkeys = [None,
                        ('ER', 'AJ', 'GT', 'GN', 'GG', 'GF'),
                        ('ER', 'AJ', 'GT', 'GN', 'GG', 'GF'),
                        ('ER', 'AJ', 'GN', 'GG', 'GFA', 'GFB'),
                        ('DET', 'DWT', 'GRT', 'GIT', 'DEF', 'DWF',
                         'GRF', 'GIF', 'DEC', 'DWC', 'GRC', 'GIC')]
        if lrf == 4:
            # Adler-Adler
            bg_flags, bg, bg_size = self._get_list(
                ('AWRI', 0, 'LI', 0, '6*NX', 'NX'),
                ('A1', 'A2', 'A3', 'A4', 'B1', 'B2'),
                subsection[total_lines:])
            total_lines += bg_size
            subsection_dict['bg'] = bg

        if lrf < 5:
            total_lines = self._nls_njs_loop(lrf_L_keys[lrf],
                                             lrf_J_keys[lrf],
                                             lrf_itemkeys[lrf],
                                             subsection,
                                             total_lines,
                                             range_flags,
                                             subsection_dict)
        if lrf == 7:
            # R-Matrix Limited Format (ENDF Manual pp. 62-67)
            # Particle pair descriptions for the whole range
            particle_pair_data, pp_size = self._get_list(
                (0, 0, 'NPP', 0, '12*NPP', '2*NPP'),
                ('MA', 'MB', 'ZA', 'ZB', 'IA', 'IB',
                 'Q', 'PNT', 'SHF', 'MT', 'PA', 'PB'),
                subsection[total_lines:])[1:3]
            total_lines += pp_size
            range_flags.update(particle_pair_data)
            for aj_section in range(int(range_flags['NJS'])):
                # Read first LIST record, with channel descriptions
                aj_flags, ch_items, ch_size = self._get_list(
                    ('AJ', 'PJ', 'KBK', 'KPS', '6*NCH', 'NCH'),
                    ('IPP', 'L', 'SCH', 'BND', 'APE', 'APT'),
                    subsection[total_lines:])
                total_lines += ch_size
                # Second LIST record, with resonance energies and widths.
                er_flags, er_data, er_size = self._get_list(
                    (0, 0, 0, 'NRS', '6*NX', 'NX'), ('ER',),
                    subsection[total_lines:])
                total_lines += er_size
                nch = int(aj_flags['NCH'])
                er_array_width = (nch//6+1)*6
                er_data = er_data['ER'].reshape(-1, er_array_width).transpose()
                aj_data = {'ER': er_data[0],
                           'GAM': er_data[1:1+nch].transpose()}
                aj_data.update(ch_items)
                aj = aj_flags['AJ']
                # Additional records
                if aj_flags['KBK'] > 0:
                    lbk_list_keys = ((), (),  # ('ED','EU',0,0,'LBK',0),
                                     ('R0', 'R1', 'R2', 'S0', 'S1', 0),
                                     ('R0', 'SO', 'GA', 0, 0, 0))
                    total_lines = read_kbks(nch, subsection,
                                            aj_data, total_lines)
                if aj_flags['KPS'] > 0:
                    total_lines = read_kpss(nch, subsection,
                                            aj_data, total_lines)
                subsection_dict[aj] = aj_data

        el, eh = range_flags['EL'], range_flags['EH']
        subsection_data = (el, eh, subsection_dict, range_flags)
        isotope_dict = self.structure[mat_id]['data'][nuc_i]
        isotope_dict['resolved'].append(subsection_data)
        return total_lines

    def _read_unresolved(self, subsection, range_flags, isotope_flags, mat_id,
                         nuc_i):
        """Read unresolved resonances of an energy subsection.

        Parameters
        -----------
        subsection: array
            Contains data for energy subsection.
        range_flags: dict
            Contains metadata flags for energy range.
        isotope_flags: dict
            Contiains flags for isotope.
        mat_id: int
            Material ZZAAAM.
        nuc_i: int
            Isotope ZZAAAM.

        Returns
        --------
        total_lines: int
        """
        head_cont = ('SPI', 'AP', 'LSSF', 0, 'NLS', 0)
        has_head_cont = {(0, 1): True, (1, 1): False,
                         (0, 2): True, (1, 2): True}
        L_keys = {(0, 1): ('AWRI', 0, 'L', 0, '6*NJS', 'NJS'),
                  (1, 1): ('AWRI', 0, 'L', 0, 'NJS', 0),
                  (0, 2): ('AWRI', 0, 'L', 0, 'NJS', 0),
                  (1, 2): ('AWRI', 0, 'L', 0, 'NJS', 0)}
        j_keys = {(0, 1): None,
                  (1, 1): (0, 0, 'L', 'MUF', 'NE+6', 0,
                           'D', 'AJ', 'AMUN', 'GN0', 'GG', 0),
                  (0, 2): ('AJ', 0, 'INT', 0, '6*NE+6', 'NE',
                           0, 0, 'AMUX', 'AMUN', 'AMUG', 'AMUF'),
                  (1, 2): ('AJ', 0, 'INT', 0, '6*NE+6', 'NE',
                           0, 0, 'AMUX', 'AMUN', 'AMUG', 'AMUF')}
        itemkeys = {(0, 1): ('D', 'AJ', 'AMUN', 'GN0', 'GG', 0),
                    (1, 1): ('GF', ),
                    (0, 2): ('ES', 'D', 'GX', 'GN0', 'GG', 'GF'),
                    (1, 2): ('ES', 'D', 'GX', 'GN0', 'GG', 'GF')}

        lfw, lrf = int(isotope_flags['LFW']), int(range_flags['LRF'])
        subsection_dict = rxdata.DoubleSpinDict({})
        if range_flags['NRO'] > 0:
            tabhead, intdata, total_lines = self._get_tab1((0, 0, 0,
                                                            0, 'NR', 'NP'),
                                                           ('E', 'AP'),
                                                           subsection)
            subsection_dict['int'] = intdata
        else:
            total_lines = 0
        if has_head_cont[(lfw, lrf)]:
            range_flags, total_lines = self._cont_and_update(
                range_flags, head_cont, subsection, total_lines)
        if (lfw, lrf) == (1, 1):
            # Case B in ENDF manual p.70
            head_flags, es_array, lines = self._get_list(
                ('SPI', 'AP', 'LSSF', 0, 'NE', 'NLS'),
                ('ES', ),
                subsection[total_lines:])
            subsection_dict['ES'] = es_array['ES']
            total_lines += lines
            range_flags.update(head_flags)
        total_lines = self._nls_njs_loop(L_keys[(lfw, lrf)],
                                         j_keys[(lfw, lrf)],
                                         itemkeys[(lfw, lrf)],
                                         subsection,
                                         total_lines,
                                         range_flags,
                                         subsection_dict)
        el, eh = range_flags['EL'], range_flags['EH']
        subsection_data = (el, eh, subsection_dict, range_flags)
        isotope_dict = self.structure[mat_id]['data'][nuc_i]
        isotope_dict['unresolved'].append(subsection_data)
        return total_lines

    def _read_ap_only(self, subsection, range_flags, isotope_flags, mat_id,
                      nuc_i):
        'Read in scattering radius when it is the only resonance data given.'
        subsection_dict = {}
        if range_flags['NRO'] > 0:
            tabhead, intdata, total_lines = self._get_tab1((0, 0, 0, 0,
                                                            'NR', 'NP'),
                                                           ('E', 'AP'),
                                                           subsection)
            subsection_dict['int'] = intdata
        else:
            total_lines = 0
        range_flags, total_lines = self._cont_and_update(
            range_flags, ('SPI', 'AP', 0, 0, 'NLS', 0), subsection, total_lines)
        return total_lines

    def _read_xs(self, nuc, mt, nuc_i=None):
        """Read in cross-section data. Read resonances with Library._read_res
        first.

        Parameters
        -----------
        nuc: int
            id of material.
        mt: int
            Reaction number to find cross-section data of.
        nuc_i: int
            Isotope to find; if None, defaults to mat_id.
        """
        nuc = nucname.id(nuc)
        if nuc_i is None:
            nuc_i = nuc
        if 600 > mt > 500:
            xsdata = self.get_rx(nuc, 23, mt).reshape(-1, 6)
        else:
            xsdata = self.get_rx(nuc, 3, mt).reshape(-1, 6)
        total_lines = 0
        head_flags = self._get_head(('ZA', 'AWR', 0, 0, 0, 0),
                                    xsdata[total_lines])
        total_lines += 1
        int_flags, int_data, int_size = self._get_tab1(
            ('QM', 'QI', 0, 'LM', 'NR', 'NP'),
            ('e_int', 'xs'),
            xsdata[total_lines:])
        int_flags.update(head_flags)
        isotope_dict = self.structure[nuc]['data'][nuc_i]
        isotope_dict['xs'].update({mt: (int_data, int_flags)})
        total_lines += int_size

    def get_xs(self, nuc, mt, nuc_i=None):
        """get_xs(nuc, mt, nuc_i=None)
        Grab cross-section data.

        Parameters
        ----------
        nuc: int
            id of nuclide to read.
        mt: int
            ENDF reaction number to read.
        nuc_i: int
            id of isotope to read. Defaults to nuc.

        Returns
        -------
        tuple
            Returns a tuple with xs data in tuple[0] and flags in tuple[1].
        """
        nuc = nucname.id(nuc)
        if not nuc_i:
            nuc_i = nuc
        else:
            nuc_i = nucname.id(nuc_i)
        if (nuc not in self.structure) or (not self.structure[nuc]['data']):
            self._read_res(nuc)
        if nuc_i not in self.structure[nuc]['data'] or \
           mt not in self.structure[nuc]['data'][nuc_i]['xs']:
            self._read_xs(nuc, mt, nuc_i)
        return self.structure[nuc]['data'][nuc_i]['xs'][mt]

    def get_rx(self, nuc, mf, mt, lines=0):
        """get_rx(nuc, mf, mt, lines=0)
        Grab the data from one reaction type.

        Parameters
        ----------
        nuc : int
            id form of material to read from.
        mf : int
            ENDF file number (MF).
        mt : int
            ENDF reaction number (MT).
        lines : int
            Number of lines to read from this reaction, starting from the top.
            Default value is 0, which reads in the entire reaction.

        Returns
        -------
        data : NumPy array
            Contains the reaction data in an Nx6 array.
        """
        nuc = nucname.id(nuc)
        if nuc in self.structure:
            return self._read_nucmfmt(nuc, mf, mt, lines)
        else:
            raise ValueError('Material {} does not exist.'.format(nuc))

    def _read_nucmfmt(self, nuc, mf, mt, lines):
        """Load in the data from one reaction into self.structure.

        Parameters
        ----------
        nuc : int
            id of nuclide.
        mf : int
            ENDF file number (MF).
        mt : int
            ENDF reaction number (MT).

        Returns
        -------
        array, 1d, float64
            1d, float64 NumPy array containing the reaction data.
        """
        opened_here = False
        if isinstance(self.fh, basestring):
            fh = open(self.fh, 'r')
            opened_here = True
        else:
            fh = self.fh
        try:
            start, stop = self.mat_dict[nuc]['mfs'][mf, mt]
        except KeyError as e:
            msg = 'MT {1} not found in File {0}.'.format(mf, mt)
            e.args = (msg,)
            raise e
        fh.readline()
        fh.seek(start)
        if lines == 0:
            s = fh.read(stop-start)
        else:
            s = fh.read(lines*self.line_length)
        if opened_here:
            fh.close
        return fromendf_tok(s)


def at_end_of_tape(f):
    """Indicate whether file is positioned at the end of an ENDF tape.

    Parameters
    ----------
    f : file_like
        File to check

    Returns
    -------
    bool
        Whether the file is at the end of the ENDF tape

    """
    position = f.tell()
    line = f.readline()
    if line == '' or line[66:70] == '  -1':
        return True
    else:
        f.seek(position)
        return False


def seek_material_end(f):
    """Position the file at the end of the ENDF material (MAT) currently being read.

    Parameters
    ----------
    f : file_like
        File to position

    """
    while True:
        line = f.readline()
        if line[66:70] == '   0':
            break


def seek_file_end(f):
    """Position the file at the end of the ENDF file (MF) currently being read.

    Parameters
    ----------
    f : file_like
        File to position

    """
    while True:
        line = f.readline()
        if line[70:72] == ' 0':
            break


def seek_section_end(f):
    """Position the file at the end of the ENDF section (MT) currently being read.

    Parameters
    ----------
    f : file_like
        File to position

    """
    while True:
        line = f.readline()
        if line[72:75] == '  0':
            break


class Evaluation(object):
    """ENDF material evaluation with multiple files/sections

    The Evaluation class provides a means to parse data from an ENDF-6 format
    file and access it as stored internal Python objects. A summary of the
    parsing capabilities is as follows:

    == === =============================================== ========
    MF MT  Description                                     Complete
    == === =============================================== ========
    1  451 Descriptive data and directory                  Yes
    1  452 Number of neutrons per fission                  Yes
    1  455 Delayed neutron data                            Yes
    1  456 Number of prompt neutrons per fission           Yes
    1  458 Components of fission energy release            Yes
    1  460 Delayed photon data                             Yes
    2  151 Resonance parameters                            Yes
    3   -  Reaction cross sections                         Yes
    4   -  Angular distributions                           Yes
    5   -  Energy distributions                            Yes
    6   -  Product energy-angle distributions              Yes
    7  2   Thermal elastic scattering                      Yes
    7  4   Thermal inelastic scattering                    Yes
    8  454 Independent fission yields                      Yes
    8  457 Radioactive decay data                          Yes
    8  459 Cumulative fission yields                       Yes
    8   -  Radioactive nuclide production                  Yes
    9   -  Multiplicities of radioactive products          Yes
    10  -  Production cross sections for radionuclides     Yes
    12  -  Photon production yield data                    Yes
    13  -  Photon production cross sections                Yes
    14  -  Photon angular distributions                    Yes
    15  -  Continuous photon energy spectra                Yes
    23  -  Photon and electron cross sections              Yes
    26  -  Secondary distributions for electro-atomic data Yes
    27  -  Atomic form factors                             Yes
    28 533 Atomic relaxation data                          Yes
    30 1   Directory and correspondance table              No
    30 2   Covariance matrix                               No
    30  -  Sensitivities                                   No
    31  -  Covariances of fission                          No
    32  -  Covariances of resonance parameters             No
    33  -  Covariances of neutron cross sections           No
    34  -  Covariances for angular distributions           No
    35  -  Covariances for energy distributions            No
    40  -  Covariances for radionuclide production         No
    == === =============================================== ========

    Attributes
    ----------
    atomic_relaxation : dict
        Dictionary containing atomic relaxation data from MF=28, MT=533. If the
        evaluation is not an atomic relaxation sublibrary, the dictionary is
        empty.
    decay : dict
        Dictionary containing decay data from MF=8. If the evaluation is not
        from a decay sublibrary, the dictionary is empty.
    fission : dict
        Dictionary containing fission-related data, such as neutrons release
        from fission (MF=1, MT=452,455,456), components of energy release (MF=1,
        MT=458), delayed photons from fission (MF=1, MT=460), and
        cumulative/independent fission yields (MF=8, MT=454,459).
    info : dict
        Miscallaneous information about the evaluation.
    target : dict
        Information about the target material, such as its mass, isomeric state,
        whether it's stable, and whether it's fissionable.
    projectile : dict
        Information about the projectile such as its mass.
    reaction_list : list of 4-tuples
        List of sections in the evaluation. The entries of the tuples are the
        file (MF), section (MT), number of records (NC), and modification
        indicator (MOD).
    reactions : collections.OrderedDict
        Dictionary whose keys are MT numbers and values are Reaction instances.
    resonances : dict
        Resolved resonance data from MF=2, MT=151.
    thermal_elastic : dict
        Coherent and/or incoherent thermal elastic data from MF=7, MT=2.
    thermal_inelastic : dict
        Incoherent thermal inelastic data from MF=7, MT=4.

    """

    def __init__(self, filename_or_handle, verbose=True):
        if hasattr(filename_or_handle, 'read'):
            self._fh = filename_or_handle
        else:
            self._fh = open(filename_or_handle, 'r')
        self._verbose = verbose
        self._veryverbose = False

        # Create public attributes
        self.atomic_relaxation = {}
        self.decay = {}
        self.fission = {'nu': {}, 'energy_release': {}, 'delayed_photon': {},
                        'yield_independent': {}, 'yield_cumulative': {}}
        self.info = {}
        self.target = {}
        self.projectile = {}
        self.reaction_list = []
        self.reactions = OrderedDict()
        self.resonances = {}
        self.thermal_elastic = {}
        self.thermal_inelastic = {}

        # Determine MAT number for this evaluation
        MF = 0
        while MF == 0:
            position = self._fh.tell()
            line = self._fh.readline()
            MF = int(line[70:72])
        self.material = int(line[66:70])

        # Save starting position for this evaluation
        self._fh.seek(position)

        # First we need to read MT=1, MT=451 which has a description of the ENDF
        # file and a list of what data exists in the file
        self._read_header()

        # Save starting position
        self._start_position = self._fh.tell()

    def read(self, reactions=None, skip_mf=[], skip_mt=[]):
        """Reads reactions from the ENDF file of the Evaluation object. If no
        arguments are provided, this method will read all the reactions in the
        file. A single reaction can be read if provided.

        Parameters
        ----------
        reactions : tuple or list of tuple, optional
            A single reaction in the following format: (MF, MT)
        skip_mf : list of int, optional
            Files (MF) which should not be read
        skip_mt : list of int, optional
            Reactions (MT) which should not be read

        """

        # Make sure file is positioned correctly
        self._fh.seek(self._start_position)

        if isinstance(reactions, tuple):
            reactions = [reactions]

        while True:
            # Find next section
            while True:
                position = self._fh.tell()
                line = self._fh.readline()
                MAT = int(line[66:70])
                MF = int(line[70:72])
                MT = int(line[72:75])
                if MT > 0 or MAT == 0:
                    self._fh.seek(position)
                    break

            # If end of material reached, exit loop
            if MAT == 0:
                break

            # If there are files/reactions requested to be skipped, check them
            if MF in skip_mf:
                seek_file_end(self._fh)
                continue
            if MT in skip_mt:
                seek_section_end(self._fh)
                continue

            # If reading is restricted to certain reactions, check here
            if reactions and (MF, MT) not in reactions:
                seek_section_end(self._fh)
                continue


            # File 1 data
            if MF == 1:
                if MT == 452:
                    # Number of total neutrons per fission
                    self._read_total_nu()
                elif MT == 455:
                    # Number of delayed neutrons per fission
                    self._read_delayed_nu()
                elif MT == 456:
                    # Number of prompt neutrons per fission
                    self._read_prompt_nu()
                elif MT == 458:
                    # Components of energy release due to fission
                    self._read_fission_energy()
                elif MT == 460:
                    self._read_delayed_photon()

            elif MF == 2:
                # Resonance parameters
                if MT == 151:
                    self._read_resonances()
                else:
                    seek_section_end(self._fh)

            elif MF == 3:
                # Reaction cross sections
                self._read_reaction_xs(MT)

            elif MF == 4:
                # Angular distributions
                self._read_angular_distribution(MT)

            elif MF == 5:
                # Energy distributions
                self._read_energy_distribution(MT)

            elif MF == 6:
                # Product energy-angle distributions
                self._read_product_energy_angle(MT)

            elif MF == 7:
                # Thermal scattering data
                if MT == 2:
                    self._read_thermal_elastic()
                if MT == 4:
                    self._read_thermal_inelastic()

            elif MF == 8:
                # decay and fission yield data
                if MT == 454:
                    self._read_independent_yield()
                elif MT == 459:
                    self._read_cumulative_yield()
                elif MT == 457:
                    self._read_decay()
                else:
                    self._read_radioactive_nuclide(MT)

            elif MF == 9:
                # multiplicities
                self._read_multiplicity(MT)

            elif MF == 10:
                # cross sections for production of radioactive nuclides
                self._read_production_xs(MT)

            elif MF == 12:
                # Photon production yield data
                self._read_photon_production_yield(MT)

            elif MF == 13:
                # Photon production cross sections
                self._read_photon_production_xs(MT)

            elif MF == 14:
                # Photon angular distributions
                self._read_photon_angular_distribution(MT)

            elif MF == 15:
                # Photon continuum energy distributions
                self._read_photon_energy_distribution(MT)

            elif MF == 23:
                # photon interaction data
                self._read_photon_interaction(MT)

            elif MF == 26:
                # secondary distributions for photon interactions
                self._read_electron_products(MT)

            elif MF == 27:
                # atomic form factors or scattering functions
                self._read_scattering_functions(MT)

            elif MF == 28:
                # atomic relaxation data
                self._read_atomic_relaxation()

            else:
                seek_file_end(self._fh)

    def _read_header(self):
        self._print_info(1, 451)

        # Information about target/projectile
        # First HEAD record
        items = self._get_head_record()
        self.target['ZA'] = items[0]
        self.target['mass'] = items[1]
        self._LRP = items[2]
        self.target['fissionable'] = (items[3] == 1)
        try:
            global libraries
            library = libraries[items[4]]
        except KeyError:
            library = 'Unknown'
        self.info['modification'] = items[5]

        # Control record 1
        items = self._get_cont_record()
        self.target['excitation_energy'] = items[0]
        self.target['stable'] = (int(items[1]) == 0)
        self.target['state'] = items[2]
        self.target['isomeric_state'] = items[3]
        self.info['format'] = items[5]
        assert self.info['format'] == 6

        # Control record 2
        items = self._get_cont_record()
        self.projectile['mass'] = items[0]
        self.info['energy_max'] = items[1]
        library_release = items[2]
        self.info['sublibrary'] = items[4]
        library_version = items[5]
        self.info['library'] = (library, library_version, library_release)

        # Control record 3
        items = self._get_cont_record()
        self.target['temperature'] = items[0]
        self.info['derived'] = (items[2] > 0)
        NWD = items[4]
        NXC = items[5]

        # Text records
        text = [self._get_text_record() for i in range(NWD)]
        if len(text) >= 5:
            self.target['zsymam'] = text[0][0:11]
            self.info['laboratory'] = text[0][11:22]
            self.info['date'] = text[0][22:32]
            self.info['author'] = text[0][32:66]
            self.info['reference'] = text[1][1:22]
            self.info['date_distribution'] = text[1][22:32]
            self.info['date_release'] = text[1][33:43]
            self.info['date_entry'] = text[1][55:63]
            self.info['identifier'] = text[2:5]
            self.info['description'] = text[5:]

        # File numbers, reaction designations, and number of records
        for i in range(NXC):
            items = self._get_cont_record(skipC=True)
            MF, MT, NC, MOD = items[2:6]
            self.reaction_list.append((MF, MT, NC, MOD))

    def _read_total_nu(self):
        self._print_info(1, 452)

        # Determine representation of total nu data
        items = self._get_head_record()
        LNU = items[3]

        # Polynomial representation
        if LNU == 1:
            coefficients = np.asarray(self._get_list_record(onlyList=True))
            self.fission['nu']['total'] = Polynomial(coefficients)

        # Tabulated representation
        elif LNU == 2:
            params, self.fission['nu']['total'] = self._get_tab1_record()

        # Skip SEND record
        self._fh.readline()

    def _read_delayed_nu(self):
        self._print_info(1, 455)

        # Create delayed nu reaction
        self.fission['nu']['delayed'] = {}

        # Determine representation of delayed nu data
        items = self._get_head_record()
        LDG = items[2]
        LNU = items[3]
        self.fission['nu']['delayed']['decay_energy_dependent'] = (LDG == 1)

        if LDG == 0:
            # Delayed-group constants energy independent
            self.fission['nu']['delayed']['decay_constants'] = np.asarray(
                self._get_list_record(onlyList=True))
        elif LDG == 1:
            # Delayed-group constants energy dependent
            raise NotImplementedError

        if LNU == 1:
            # Nu represented as polynomial
            coefficients = np.asarray(self._get_list_record(onlyList=True))
            self.fission['nu']['delayed']['values'] = Polynomial(coefficients)
        elif LNU == 2:
            # Nu represented by tabulation
            params, self.fission['nu']['delayed']['values'] = self._get_tab1_record()
        self._fh.readline()

    def _read_prompt_nu(self):
        self._print_info(1, 456)

        # Determine representation of delayed nu data
        items = self._get_head_record()
        LNU = items[3]

        if LNU == 1:
            # Polynomial representation (spontaneous fission)
            coefficients = np.asarray(self._get_list_record(onlyList=True))
            self.fission['nu']['prompt'] = Polynomial(coefficients)
        elif LNU == 2:
            # Tabulated values of nu
            params, self.fission['nu']['prompt'] = self._get_tab1_record()

        # Skip SEND record
        self._fh.readline()

    def _read_fission_energy(self):
        self._print_info(1, 458)
        er = self.fission['energy_release']

        # Skip HEAD record
        self._get_head_record()

        # Read LIST record containing components of fission energy release (or
        # coefficients)
        items, values = self._get_list_record()
        NPLY = items[3]
        er['order'] = NPLY

        values = np.asarray(values)
        values.shape = (NPLY + 1, 18)
        er['fission_products'] = np.vstack((values[:,0], values[:,1]))
        er['prompt_neutrons'] = np.vstack((values[:,2], values[:,3]))
        er['delayed_neutrons'] = np.vstack((values[:,4], values[:,5]))
        er['prompt_gammas'] = np.vstack((values[:,6], values[:,7]))
        er['delayed_gammas'] = np.vstack((values[:,8], values[:,9]))
        er['delayed_betas'] = np.vstack((values[:,10], values[:,11]))
        er['neutrinos'] = np.vstack((values[:,12], values[:,13]))
        er['total_less_neutrinos'] = np.vstack((values[:,14], values[:,15]))
        er['total'] = np.vstack((values[:,16], values[:,17]))

        # Skip SEND record
        self._fh.readline()

    def _read_reaction_xs(self, MT):
        self._print_info(3, MT)

        # Get Reaction instance
        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(3)

        # Read HEAD record with ZA and atomic mass ratio
        items = self._get_head_record()

        # Read TAB1 record with reaction cross section
        params, rx.xs = self._get_tab1_record()
        rx.Q_mass_difference = params[0]
        rx.Q_reaction = params[1]
        rx.complex_breakup_flag = params[3]

        # Skip SEND record
        self._fh.readline()

    def _read_angular_distribution(self, MT):
        # Find energy distribution
        self._print_info(4, MT)

        # Get Reaction instance
        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(4)

        adist = AngularDistribution()
        rx.angular_distribution = adist

        # Read HEAD record
        items = self._get_head_record()
        ltt = items[3]

        # Read CONT record
        items =self._get_cont_record()
        li = items[2]
        adist.center_of_mass = (items[3] == 2)

        if ltt == 0 and li == 1:
            # Purely isotropic
            adist.type = 'isotropic'

        elif ltt == 1 and li == 0:
            # Legendre polynomial coefficients
            adist.type = 'legendre'

            adist.tab2 = self._get_tab2_record()
            n_energy = adist.tab2.params[5]

            adist.energy = np.zeros(n_energy)
            adist.probability = []
            for i in range(n_energy):
                items, al = self._get_list_record()
                temperature = items[0]
                adist.energy[i] = items[1]
                coefficients = np.asarray([1.0] + al)
                for i in range(len(coefficients)):
                    coefficients[i] *= (2.*i + 1.)/2.
                adist.probability.append(Legendre(coefficients))

        elif ltt == 2 and li == 0:
            # Tabulated probability distribution
            adist.type = 'tabulated'

            adist.tab2 = self._get_tab2_record()
            n_energy = adist.tab2.params[5]

            adist.energy = np.zeros(n_energy)
            adist.probability = []
            for i in range(n_energy):
                params, f = self._get_tab1_record()
                temperature = params[0]
                adist.energy[i] = params[1]
                adist.probability.append(f)

        elif ltt == 3 and li == 0:
            # Legendre for low energies / tabulated for high energies
            adist.type = 'legendre/tabulated'

            adist.tab2_legendre = self._get_tab2_record()
            n_energy_legendre = adist.tab2_legendre.params[5]

            energy_legendre = np.zeros(n_energy_legendre)
            adist.probability = []
            for i in range(n_energy_legendre):
                items, al = self._get_list_record()
                temperature = items[0]
                energy_legendre[i] = items[1]
                coefficients = np.asarray([1.0] + al)
                for i in range(len(coefficients)):
                    coefficients[i] *= (2.*i + 1.)/2.
                adist.probability.append(Legendre(coefficients))

            adist.tab2_tabulated = self._get_tab2_record()
            n_energy_tabulated = adist.tab2_tabulated.params[5]

            energy_tabulated = np.zeros(n_energy_tabulated)
            for i in range(n_energy_tabulated):
                params, f = self._get_tab1_record()
                temperature = params[0]
                energy_tabulated[i] = params[1]
                adist.probability.append(f)

            adist.energy = np.concatenate((energy_legendre, energy_tabulated))

    def _read_energy_distribution(self, MT):
        # Find energy distribution
        self._print_info(5, MT)

        # Get Reaction instance
        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(5)

        # Read HEAD record
        items = self._get_head_record()
        nk = items[4]

        for i in range(nk):
            edist = EnergyDistribution()

            # Read TAB1 record for p(E)
            params, applicability = self._get_tab1_record()
            lf = params[3]
            if lf == 1:
                tab2 = self._get_tab2_record()
                n_energies = tab2.params[5]

                energy = np.zeros(n_energies)
                pdf = []
                for j in range(n_energies):
                    params, func = self._get_tab1_record()
                    energy[j] = params[1]
                    pdf.append(func)
                edist = ArbitraryTabulated(energy, pdf)
            elif lf == 5:
                # General evaporation spectrum
                u = params[0]
                params, theta = self._get_tab1_record()
                params, g = self._get_tab1_record()
                edist = GeneralEvaporation(theta, g, u)
            elif lf == 7:
                # Simple Maxwellian fission spectrum
                u = params[0]
                params, theta = self._get_tab1_record()
                edist = Maxwellian(theta, u)
            elif lf == 9:
                # Evaporation spectrum
                u = params[0]
                params, theta = self._get_tab1_record()
                edist = Evaporation(theta, u)
            elif lf == 11:
                # Energy-dependent Watt spectrum
                u = params[0]
                params, a = self._get_tab1_record()
                params, b = self._get_tab1_record()
                edist = Watt(a, b, u)
            elif lf == 12:
                # Energy-dependent fission neutron spectrum (Madland-Nix)
                params, tm = self._get_tab1_record()
                efl, efh = params[0:2]
                edist = MadlandNix(efl, efh, tm)

            edist.applicability = applicability
            rx.energy_distribution.append(edist)

    def _read_product_energy_angle(self, MT):
        # Find distribution
        self._print_info(6, MT)

        # Get Reaction instance
        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(6)

        # Read HEAD record
        items = self._get_head_record()
        rx.reference_frame = {1: 'laboratory', 2: 'center-of-mass',
                               3: 'light-heavy'}[items[3]]
        n_products = items[4]

        for i in range(n_products):
            product = {}
            rx.product_distribution.append(product)

            # Read TAB1 record for product yield
            params, product['yield'] = self._get_tab1_record()
            product['za'] = params[0]
            product['mass'] = params[1]
            product['lip'] = params[2]
            product['law'] = params[3]

            if product['law'] == 1:
                # Continuum energy-angle distribution
                tab2 = self._get_tab2_record()
                product['lang'] = tab2.params[2]
                product['lep'] = tab2.params[3]
                ne = tab2.params[5]
                product['energy'] = np.zeros(ne)
                product['n_discrete_energies'] = np.zeros(ne)
                product['energy_out'] = []
                product['b'] = []
                for i in range(ne):
                    items, values = self._get_list_record()
                    product['energy'][i] = items[1]
                    product['n_discrete_energies'][i] = items[2]
                    n_angle = items[3]
                    n_energy_out = items[5]
                    values = np.array(values)
                    values.shape = (n_energy_out, n_angle + 2)
                    product['energy_out'].append(values[:,0])
                    product['b'].append(values[:,1:])

            elif product['law'] == 2:
                # Discrete two-body scattering
                tab2 = self._get_tab2_record()
                ne = tab2.params[5]
                product['energy'] = np.zeros(ne)
                product['lang'] = np.zeros(ne, dtype=int)
                product['Al'] = []
                for i in range(ne):
                    items, values = self._get_list_record()
                    product['energy'][i] = items[1]
                    product['lang'][i] = items[2]
                    product['Al'].append(np.asarray(values))

            elif product['law'] == 5:
                # Charged particle elastic scattering
                tab2 = self._get_tab2_record()
                product['spin'] = tab2.params[0]
                product['identical'] = (tab2.params[2] == 1)
                ne = tab2.params[5]
                product['energies'] = np.zeros(ne)
                product['ltp'] = np.zeros(ne, dtype=int)
                product['coefficients'] = []
                for i in range(ne):
                    items, values = self._get_list_record()
                    product['energies'][i] = items[1]
                    product['lpt'][i] = items[2]
                    product['coefficients'].append(values)

            elif product['law'] == 6:
                # N-body phase-space distribution
                items = self._get_cont_record()
                product['total_mass'] = items[0]
                product['n_particles'] = items[5]

            elif product['law'] == 7:
                # Laboratory energy-angle distribution
                tab2 = self._get_tab2_record()
                ne = tab2.params[5]
                product['energies'] = np.zeros(ne)
                product['mu'] = []
                product['distribution'] = []
                for i in range(ne):
                    tab2mu = self._get_tab2_record()
                    product['energies'][i] = tab2mu.params[1]
                    nmu = tab2mu.params[5]
                    mu = np.zeros(nmu)
                    dists = []
                    for j in range(nmu):
                        params, f = self._get_tab1_record()
                        mu[j] = params[1]
                        dists.append(f)
                    product['mu'].append(mu)
                    product['distribution'].append(dists)

    def _read_delayed_photon(self):
        self._print_info(1, 460)
        dp = self.fission['delayed_photon']

        # Determine whether discrete or continuous representation
        items = self._get_head_record()
        LO = items[2]
        NG = items[4]

        # Discrete representation
        if LO == 1:
            dp['form'] = 'discrete'

            # Initialize lists for energies of photons and time dependence of
            # photon multiplicity
            dp['energy'] = np.zeros(NG)
            dp['multiplicity'] = []
            for i in range(NG):
                # Read TAB1 record with multiplicity as function of time
                params, mult = self._get_tab1_record()
                dp['multiplicity'].append(mult)

                # Determine energy
                dp['energy'][i] = params[0]

        # Continuous representation
        elif LO == 2:
            # Determine decay constant and number of precursor families
            dp['form'] = 'continuous'
            dp['decay_constant'] = self._get_list_record(onlyList=True)

    def _read_resonances(self):
        self._print_info(2, 151)
        res = self.resonances

        # Determine whether discrete or continuous representation
        items = self._get_head_record()
        NIS = items[4] # Number of isotopes
        res['isotopes'] = []

        for iso in range(NIS):
            # Create dictionary for this isotope
            isotope = {}
            res['isotopes'].append(isotope)

            items = self._get_cont_record()
            isotope['abundance'] = items[1]
            LFW = items[3] # average fission width flag
            NER = items[4] # number of resonance energy ranges

            isotope['ranges'] = []

            for j in range(NER):
                items = self._get_cont_record()
                emin, emax = items[0:2]  # min/max energies of range
                resonance_flag = items[2]  # flag for resolved (1)/unresolved (2)
                resonance_formalism = items[3]  # resonance formalism
                nro = items[4]  # flag for energy dependence of scattering radius
                naps = items[5]  # flag controlling use of channel/scattering radius

                if resonance_flag == 0 and nro == 0:
                    # Only scattering radius specified
                    erange = ScatteringRadius(emin, emax, nro, naps)
                    items = self._get_cont_record()
                    erange.spin = items[0]
                    erange.scattering_radius = items[1]

                elif resonance_flag == 1:
                    # resolved resonance region
                    erange = _formalisms[resonance_formalism](
                        emin, emax, nro, naps)
                    erange.read(self)
                    if NIS == 1:
                        res['resolved'] = erange

                elif resonance_flag == 2:
                    # unresolved resonance region
                    erange = Unresolved(emin, emax, nro, naps)
                    erange.fission_widths = (LFW == 1)
                    erange.LRF = resonance_formalism
                    erange.read(self)
                    if NIS == 1:
                        res['unresolved'] = erange

                erange.material = self
                isotope['ranges'].append(erange)


    def _read_thermal_elastic(self):
        self._print_info(7, 2)
        elast = self.thermal_elastic

        # Get head record
        items = self._get_head_record()
        LTHR = items[2]  # coherent/incoherent flag
        elast['S'] = {}

        if LTHR == 1:
            elast['type'] = 'coherent'
            params, sdata = self._get_tab1_record()
            temperature = params[0]
            LT = params[2]
            elast['S'][temperature] = sdata

            for t in range(LT):
                params, sdata = self._get_list_record()
                temperature = params[0]
                LT = params[2]
                elast['S'][temperature] = sdata

        elif LTHR == 2:
            elast['type'] = 'incoherent'
            params, wt = self._get_tab1_record()
            elast['bound_xs'] = params[0]
            elast['debye_waller'] = wt

    def _read_thermal_inelastic(self):
        self._print_info(7, 4)
        inel = self.thermal_inelastic

        # Get head record
        items = self._get_head_record()
        inel['temperature_used'] = 'actual' if items[3] == 0 else '0.0253 eV'  # Temperature flag
        inel['symmetric'] = (items[4] == 0)  # Symmetry flag
        header, B = self._get_list_record()
        inel['ln(S)'] = (header[2] == 1)
        inel['num_non_principal'] = header[5]
        inel['B'] = B
        if B[0] != 0.0:
            tab2 = self._get_tab2_record()
            n_beta = tab2.NBT[0]
            for i_beta in range(n_beta):
                #Read record for first temperature (always present)
                params, sab0 = self._get_tab1_record()
                n_temps = params[2] + 1

                # Create arrays on first pass through -- note that alphas and
                # temperatures only need to be stored on first beta
                if i_beta == 0:
                    alpha_values = sab0.x
                    n_alpha = alpha_values.shape[0]
                    beta_values = np.zeros(n_beta)
                    temp_values = np.zeros(n_temps)
                    temp_values[0] = params[0]
                    sab_values = np.zeros((n_alpha, n_beta, n_temps), order='F')

                # Store beta and S(a,b,0) for first beta
                beta_values[i_beta] = params[1]
                sab_values[:, i_beta, 0] = sab0.y

                for i_temp in range(1, n_temps):
                    # Read records for all the other temperatures
                    params, sabt = self._get_list_record()
                    if i_beta == 0:
                        temp_values[i_temp] = params[0]
                    sab_values[:, i_beta, i_temp] = sabt

            # Store arrays in dictionary
            inel['scattering_law'] = sab_values
            inel['alpha'] = alpha_values
            inel['beta'] = beta_values
            inel['temperature'] = temp_values

        params, teff = self._get_tab1_record()
        inel['teff'] = teff

        #Read teff curves for non-principal atoms with SCT approximation indicated:
        for i in range(len(B)//6-1):
            if B[6*(i+1)]==0.0:
                params, teff = self._get_tab1_record()
                inel['teff_%i'%(i+1)] = teff#Store as teff_1, teff_2, etc.


    def _read_independent_yield(self):
        self._print_info(8, 454)
        iyield = self.fission['yield_independent']

        # Initialize energies and yield dictionary
        iyield['energies'] = []
        iyield['data'] = {}
        iyield['interp'] = []

        items = self._get_head_record()
        LE = items[2] - 1  # Determine energy-dependence

        for i in range(LE + 1):
            items, values = self._get_list_record()
            E = items[0]  # Incident particle energy
            iyield['energies'].append(E)
            NFP = items[5]  # Number of fission product nuclide states
            if i > 0:
                iyield['interp'].append(items[2]) # Interpolation scheme

            # Get data for each yield
            iyield['data'][E] = {}
            iyield['data'][E]['zafp'] = [int(i) for i in values[0::4]] # ZA for fission products
            iyield['data'][E]['fps'] = values[1::4] # State designator
            iyield['data'][E]['yi'] = list(zip(values[2::4],values[3::4])) # Independent yield

        # Skip SEND record
        self._fh.readline()

    def _read_cumulative_yield(self):
        self._print_info(8, 459)
        cyield = self.fission['yield_cumulative']

        # Initialize energies and yield dictionary
        cyield['energies'] = []
        cyield['data'] = {}
        cyield['interp'] = []

        items = self._get_head_record()
        LE = items[2] - 1  # Determine energy-dependence

        for i in range(LE + 1):
            items, values = self._get_list_record()
            E = items[0]  # Incident particle energy
            cyield['energies'].append(E)
            NFP = items[5]  # Number of fission product nuclide states
            if i > 0:
                cyield['interp'].append(items[2]) # Interpolation scheme

            # Get data for each yield
            cyield['data'][E] = {}
            cyield['data'][E]['zafp'] = [int(i) for i in values[0::4]] # ZA for fission products
            cyield['data'][E]['fps'] = values[1::4] # State designator
            cyield['data'][E]['yc'] = list(zip(values[2::4],values[3::4])) # Cumulative yield

        # Skip SEND record
        self._fh.readline()

    def _read_decay(self):
        self._print_info(8, 457)
        decay = self.decay

        # Get head record
        items = self._get_head_record()
        decay['ZA'] = items[0]  # ZA identifier
        decay['awr'] = items[1]  # AWR
        decay['state']= items[2]  # State of the original nuclide
        decay['isomeric_state'] = items[3]  # Isomeric state for the original nuclide
        decay['stable'] = (items[4] == 1)  # Nucleus stability flag

        # Determine if radioactive/stable
        if not decay['stable']:
            NSP = items[5]  # Number of radiation types

            # Half-life and decay energies
            items, values = self._get_list_record()
            decay['half_life'] = (items[0], items[1])
            decay['NC'] = items[4]//2
            decay['energies'] = list(zip(values[0::2], values[1::2]))

            # Decay mode information
            items, values = self._get_list_record()
            decay['spin'] = items[0]  # Spin of the nuclide
            decay['parity'] = items[1]  # Parity of the nuclide
            NDK = items[5]  # Number of decay modes

            # Decay type (beta, gamma, etc.)
            decay['modes'] = []
            for i in range(NDK):
                mode = {}
                mode['type'] = _radiation_type(values[6*i])
                mode['isomeric_state'] = values[6*i + 1]
                mode['energy'] = tuple(values[6*i + 2:6*i + 4])
                mode['branching_ratio'] = tuple(values[6*i + 4:6*(i + 1)])
                decay['modes'].append(mode)

            discrete_type = {0.0: None, 1.0: 'allowed', 2.0: 'first-forbidden',
                             3.0: 'second-forbidden'}

            # Read spectra
            decay['spectra'] = {}
            for i in range(NSP):
                spectrum = {}

                items, values = self._get_list_record()
                # Decay radiation type
                spectrum['type'] = _radiation_type(items[1])
                # Continuous spectrum flag
                spectrum['continuous_flag'] = {0: 'discrete', 1: 'continuous',
                                               2: 'both'}[items[2]]
                spectrum['discrete_normalization'] = tuple(values[0:2])
                spectrum['energy_average'] = tuple(values[2:4])
                spectrum['continuous_normalization'] = tuple(values[4:6])

                NER = items[5]  # Number of tabulated discrete energies

                if not spectrum['continuous_flag'] == 'continuous':
                    # Information about discrete spectrum
                    spectrum['discrete'] = []
                    for j in range(NER):
                        items, values = self._get_list_record()
                        di = {}
                        di['energy'] = tuple(items[0:2])
                        di['from_mode'] = _radiation_type(values[0])
                        di['type'] = discrete_type[values[1]]
                        di['intensity'] = tuple(values[2:4])
                        if spectrum['type'] == 'ec/beta+':
                            di['positron_intensity'] = tuple(values[4:6])
                        elif spectrum['type'] == 'gamma':
                            di['internal_pair'] = tuple(values[4:6])
                            di['total_internal_conversion'] = tuple(values[6:8])
                            di['k_shell_conversion'] = tuple(values[8:10])
                            di['l_shell_conversion'] = tuple(values[10:12])
                        spectrum['discrete'].append(di)

                if not spectrum['continuous_flag'] == 'discrete':
                    # Read continuous spectrum
                    ci = {}
                    params, ci['probability'] = self._get_tab1_record()
                    ci['type'] = _radiation_type(params[0])

                    # Read covariance (Ek, Fk) table
                    LCOV = params[3]
                    if LCOV != 0:
                        items, values = self._get_list_record()
                        ci['covariance_lb'] = items[3]
                        ci['covariance'] = list(zip(values[0::2], values[1::2]))

                    spectrum['continuous'] = ci

                # Add spectrum to dictionary
                decay['spectra'][spectrum['type']] = spectrum

        else:
            items, values = self._get_list_record()
            items, values = self._get_list_record()
            decay['spin'] = items[0]
            decay['parity'] = items[1]

        # Skip SEND record
        self._fh.readline()

    def _read_radioactive_nuclide(self, MT):
        self._print_info(8, MT)

        # Get Reaction instance
        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(8)

        # Get head record
        items = self._get_head_record()
        NS = items[4]
        complete_chain = (items[5] == 0)

        if complete_chain:
            # If complete chain is specified here rather than in MF=8, get all
            # levels of radionuclide and corresponding decay modes
            rx.radionuclide_production = []
            for i in range(NS):
                items, values = self._get_list_record()
                radionuclide = {}
                radionuclide['za'] = items[0]
                radionuclide['excitation_energy'] = items[1]
                radionuclide['mf_multiplicity'] = items[2]
                radionuclide['level'] = items[3]
                radionuclide['modes'] = []
                for j in range(items[4]//6):
                    mode = {}
                    mode['half_life'] = values[6*j]
                    mode['type'] = _radiation_type(values[6*j + 1])
                    mode['za'] = values[6*j + 2]
                    mode['branching_ratio'] = values[6*j + 3]
                    mode['endpoint_energy'] = values[6*j + 4]
                    mode['chain_terminator'] = values[6*j + 5]
                    radionuclide['modes'].append(mode)
                rx.radionuclide_production.append(radionuclide)
        else:
            items = self._get_cont_record()
            rx.radionuclide_production = radionuclide = {}
            radionuclide['za'] = items[0]
            radionuclide['excitation_energy'] = items[1]
            radionuclide['mf_multiplicity'] = items[2]
            radionuclide['level'] = items[3]

    def _read_multiplicity(self, MT):
        self._print_info(9, MT)

        # Get Reaction instance
        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(9)

        # Get head record
        items = self._get_head_record()
        NS = items[4]  # Number of final states

        for i in range(NS):
            params, state = self._get_tab1_record()
            QM = params[0] # Mass difference Q value (eV)
            QI = params[1] # Reaction Q value (eV)
            IZAP = params[2] # 1000Z + A
            LFS = params[3] # Level number of the nuclide
            rx.multiplicities[LFS] = {'QM': QM, 'QI': QI, 'ZA': IZAP,
                                       'values': state}

    def _read_production_xs(self, MT):
        self._print_info(10, MT)

        # Get Reaction instance
        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(10)

        # Get head record
        items = self._get_head_record()
        NS = items[4]  # Number of final states

        for i in range(NS):
            params, state = self._get_tab1_record()
            QM = params[0]  # Mass difference Q value (eV)
            QI = params[1]  # Reaction Q value (eV)
            IZAP = params[2]  # 1000Z + A
            LFS = params[3]  # Level number of the nuclide
            rx.production[LFS] = {'QM': QM, 'QI': QI, 'ZA': IZAP,
                                   'values': state}

    def _read_photon_production_yield(self, MT):
        self._print_info(12, MT)

        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(12)
        rx.photon_production['yield'] = ppyield = {}

        # Determine option
        items = self._get_head_record()
        option = items[2]

        if option == 1:
            # Multiplicities given
            ppyield['type'] = 'multiplicity'
            n_discrete_photon = items[4]
            if n_discrete_photon > 1:
                items, ppyield['total'] = self._get_tab1_record()
            ppyield['discrete'] = []
            for k in range(n_discrete_photon):
                y = {}
                items, y['yield'] = self._get_tab1_record()
                y['energy_photon'] = items[0]
                y['energy_level'] = items[1]
                y['lp'] = items[2]
                y['law'] = items[3]
                ppyield['discrete'].append(y)

        elif option == 2:
            # Transition probability arrays given
            ppyield['type'] = 'transition'
            ppyield['transition'] = transition = {}

            # Determine whether simple (LG=1) or complex (LG=2) transitions
            lg = items[3]

            # Get transition data
            items, values = self._get_list_record()
            transition['energy_start'] = items[0]
            transition['energies'] = np.array(values[::lg + 1])
            transition['direct_probability'] = np.array(values[1::lg + 1])
            if lg == 2:
                # Complex case
                transition['conditional_probability'] = np.array(
                    values[2::lg + 1])

    def _read_photon_production_xs(self, MT):
        self._print_info(13, MT)

        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(13)
        rx.photon_production['xs'] = ppxs = {}

        # Determine option
        items = self._get_head_record()
        n_discrete_photon = items[4]
        if n_discrete_photon > 1:
            items, ppxs['total'] = self._get_tab1_record()
        ppxs['discrete'] = []
        for k in range(n_discrete_photon):
            xs = {}
            items, xs['xs'] = self._get_tab1_record()
            xs['energy_photon'] = items[0]
            xs['energy_level'] = items[1]
            xs['lp'] = items[2]
            xs['law'] = items[3]
            ppxs['discrete'].append(xs)

    def _read_photon_angular_distribution(self, MT):
        self._print_info(14, MT)

        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(14)
        rx.photon_production['angular_distribution'] = ppad = {}

        # Determine format for angular distributions
        items = self._get_head_record()
        ppad['isotropic'] = (items[2] == 1)

        if not ppad['isotropic']:
            ltt = items[3]
            n_discrete_photon = items[4]
            n_isotropic = items[5]
            ppad['discrete'] = []
            for i in range(n_isotropic):
                adist = AngularDistribution()
                adist.type = 'isotropic'

                items = self._get_cont_record()
                adist.energy_photon = items[0]
                adist.energy_level = items[1]
                ppad['discrete'].append(adist)

            if ltt == 1:
                # Legendre polynomial coefficients
                for i in range(n_isotropic, n_discrete_photon):
                    adist = AngularDistribution()
                    adist.type = 'legendre'

                    adist.tab2 = self._get_tab2_record()
                    adist.energy_photon = adist.tab2.params[0]
                    adist.energy_level = adist.tab2.params[1]
                    n_energy = adist.tab2.params[5]

                    adist.energy = np.zeros(n_energy)
                    adist.probability = []
                    for i in range(n_energy):
                        items, al = self._get_list_record()
                        adist.energy[i] = items[1]
                        coefficients = np.asarray([1.0] + al)
                        for i in range(len(coefficients)):
                            coefficients[i] *= (2.*i + 1.)/2.
                        adist.probability.append(Legendre(coefficients))
                    ppad['discrete'].append(adist)

            elif ltt == 2:
                # Tabulated probability distribution
                for i in range(n_isotropic, n_discrete_photon):
                    adist.type = 'tabulated'

                    adist.tab2 = self._get_tab2_record()
                    adist.energy_photon = adist.tab2.params[0]
                    adist.energy_level = adist.tab2.params[1]
                    n_energy = adist.tab2.params[5]

                    adist.energy = np.zeros(n_energy)
                    adist.probability = []
                    for i in range(n_energy):
                        params, f = self._get_tab1_record()
                        adist.energy[i] = params[1]
                        adist.probability.append(f)
                    ppad['discrete'].append(adist)

    def _read_photon_energy_distribution(self, MT):
        self._print_info(15, MT)

        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(15)
        rx.photon_production['energy_distribution'] = pped = []

        # Read HEAD record
        items = self._get_head_record()
        nc = items[4]

        for i in range(nc):
            edist = EnergyDistribution()

            # Read TAB1 record for p(E)
            params, applicability = self._get_tab1_record()
            lf = params[3]
            if lf == 1:
                # Arbitrary tabulated function -- only format currently
                # available in ENDF-102 for photon continuum
                tab2 = self._get_tab2_record()
                n_energies = tab2.params[5]

                energy = np.zeros(n_energies)
                pdf = []
                for j in range(n_energies):
                    params, func = self._get_tab1_record()
                    energy[j] = params[1]
                    pdf.append(func)
                edist = ArbitraryTabulated(energy, pdf)

            edist.applicability = applicability
            pped.append(edist)

    def _read_photon_interaction(self, MT):
        self._print_info(23, MT)

        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(23)

        # Skip HEAD record
        self._get_head_record()

        # Read cross section
        params, rx.xs = self._get_tab1_record()
        if MT >= 534 and MT <= 599:
            rx.subshell_binding_energy = params[0]
        if MT >= 534 and MT <= 572:
            rx.fluorescence_yield = params[1]

        # Skip SEND record
        self._fh.readline()

    def _read_electron_products(self, MT):
        self._print_info(26, MT)

        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(26)

        # Read HEAD record
        items = self._get_head_record()
        n_products = items[4]

        for i in range(n_products):
            product = {}
            rx.products.append(product)

            # Read TAB1 record for product yield
            params, product['yield'] = self._get_tab1_record()
            product['za'] = params[0]
            product['law'] = params[3]

            if product['law'] == 1:
                # Continuum energy-angle distribution
                tab2 = self._get_tab2_record()
                product['lang'] = tab2.params[2]
                product['lep'] = tab2.params[3]
                ne = tab2.params[5]
                product['energy'] = np.zeros(ne)
                product['n_discrete_energies'] = np.zeros(ne)
                product['energy_out'] = []
                product['b'] = []
                for i in range(ne):
                    items, values = self._get_list_record()
                    product['energy'][i] = items[1]
                    product['n_discrete_energies'][i] = items[2]
                    n_angle = items[3]
                    n_energy_out = items[5]
                    values = np.array(values)
                    values.shape = (n_energy_out, n_angle + 2)
                    product['energy_out'].append(values[:,0])
                    product['b'].append(values[:,1:])

            elif product['law'] == 2:
                # Discrete two-body scattering
                tab2 = self._get_tab2_record()
                ne = tab2.params[5]
                product['energy'] = np.zeros(ne)
                product['lang'] = np.zeros(ne, dtype=int)
                product['Al'] = []
                for i in range(ne):
                    items, values = self._get_list_record()
                    product['energy'][i] = items[1]
                    product['lang'][i] = items[2]
                    product['Al'].append(np.asarray(values))

            elif product['law'] == 8:
                # Energy transfer for excitation
                params, product['energy_transfer'] = self._get_tab1_record()

    def _read_scattering_functions(self, MT):
        self._print_info(27, MT)

        # Skip HEAD record
        self._get_head_record()

        # Get scattering function
        params, func = self._get_tab1_record()

        # Store in appropriate place
        if MT in (502, 504):
            rx = self.reactions[MT]
            rx.scattering_factor = func
        elif MT == 505:
            rx = self.reactions[502]
            rx.anomalous_scattering_imaginary = func
        elif MT == 506:
            rx = self.reactions[502]
            rx.anomalous_scattering_real = func

        # Skip SEND record
        self._fh.readline()

    def _read_atomic_relaxation(self):
        self._print_info(28, 533)

        # Read HEAD record
        params = self._get_head_record()
        n_subshells = params[4]

        # Read list of data
        subshells = {1: 'K', 2: 'L1', 3: 'L2', 4: 'L3', 5: 'M1',
                     6: 'M2', 7: 'M3', 8: 'M4', 9: 'M5', 10: 'N1',
                     11: 'N2', 12: 'N3', 13: 'N4', 14: 'N5', 15: 'N6',
                     16: 'N7', 17: 'O1', 18: 'O2', 19: 'O3', 20: 'O4',
                     21: 'O5', 22: 'O6', 23: 'O7', 24: 'O8', 25: 'O9',
                     26: 'P1', 27: 'P2', 28: 'P3', 29: 'P4', 30: 'P5',
                     31: 'P6', 32: 'P7', 33: 'P8', 34: 'P9', 35: 'P10',
                     36: 'P11', 37: 'Q1', 38: 'Q2', 39: 'Q3', 0: None}
        for i in range(n_subshells):
             params, list_items = self._get_list_record()
             subi = subshells[int(params[0])]
             n_transitions = int(params[5])
             ebi = list_items[0]
             eln = list_items[1]
             data = {'binding_energy': ebi, 'number_electrons': eln, 'transitions': []}
             for j in range(n_transitions):
                 subj = subshells[int(list_items[6*(j+1)])]
                 subk = subshells[int(list_items[6*(j+1) + 1])]
                 etr = list_items[6*(j+1) + 2]
                 ftr = list_items[6*(j+1) + 3]
                 data['transitions'].append((subj, subk, etr, ftr))
             self.atomic_relaxation[subi] = data

        # Skip SEND record
        self._fh.readline()

    def _get_text_record(self, line=None):
        if not line:
            line = self._fh.readline()
        if self._veryverbose:
            print('Get TEXT record')
        return line[0:66]

    def _get_cont_record(self, line=None, skipC=False):
        if self._veryverbose:
            print('Get CONT record')
        if not line:
            line = self._fh.readline()
        if skipC:
            C1 = None
            C2 = None
        else:
            C1 = endftod(line[:11])
            C2 = endftod(line[11:22])
        L1 = int(line[22:33])
        L2 = int(line[33:44])
        N1 = int(line[44:55])
        N2 = int(line[55:66])
        return [C1, C2, L1, L2, N1, N2]

    def _get_head_record(self, line=None):
        if not line:
            line = self._fh.readline()
        if self._veryverbose:
            print('Get HEAD record')
        ZA = int(endftod(line[:11]))
        AWR = endftod(line[11:22])
        L1 = int(line[22:33])
        L2 = int(line[33:44])
        N1 = int(line[44:55])
        N2 = int(line[55:66])
        return [ZA, AWR, L1, L2, N1, N2]

    def _get_list_record(self, onlyList=False):
        # determine how many items are in list
        if self._veryverbose:
            print('Get LIST record')
        items = self._get_cont_record()
        NPL = items[4]

        # read items
        itemsList = []
        m = 0
        for i in range((NPL-1)//6 + 1):
            line = self._fh.readline()
            toRead = min(6, NPL-m)
            for j in range(toRead):
                val = endftod(line[0:11])
                itemsList.append(val)
                line = line[11:]
            m = m + toRead
        if onlyList:
            return itemsList
        else:
            return (items, itemsList)

    def _get_tab1_record(self):
        if self._veryverbose:
            print('Get TAB1 record')
        return Tab1.from_file(self._fh)

    def _get_tab2_record(self):
        if self._veryverbose:
            print('Get TAB2 record')
        r = ENDFTab2Record()
        r.read(self._fh)
        return r


    def _print_info(self, MF, MT):
        if self._verbose:
            print('Reading MF={0}, MT={1} {2}'.format(MF, MT, label(MT)))

    def __repr__(self):
        try:
            name = self.target['zsymam'].replace(' ', '')
            library = self.info['library']
        except:
            name = 'Undetermined'
            library = 'None'
        return '<Evaluation: {0}, {1}>'.format(name, library)


class Tab1(object):
    """A one-dimensional tabulated function.

    Parameters
    ----------
    x : array_like
        Independent variable
    y : array_like
        Dependent variable
    nbt : array_like of int
        Breakpoints for interpolation regions
    interp : array_like of int
        Interpolation scheme identification number, e.g., 3 means y is linear in
        ln(x).

    """

    def __init__(self, x, y, nbt, interp):
        if len(nbt) == 0 and len(interp) == 0:
            self.n_regions = 1
            self.nbt = np.array([len(x)])
            self.interp = np.array([2])  # linear-linear by default
        else:
            # NR=0 implies linear-linear interpolation by default
            self.n_regions = len(nbt)
            self.nbt = np.asarray(nbt, dtype=int)
            self.interp = np.asarray(interp, dtype=int)

        self.n_pairs = len(x)
        self.x = np.asarray(x)  # Abscissa values
        self.y = np.asarray(y)  # Ordinate values

    @classmethod
    def from_file(cls, fh):
        """Create Tab1 object using open file handle.

        Parameters
        ----------
        fh : file
            File open and positioned to the start of a TAB1 record.

        """

        # Determine how many interpolation regions and total points there are
        line = fh.readline()
        C1 = endftod(line[:11])
        C2 = endftod(line[11:22])
        L1 = int(line[22:33])
        L2 = int(line[33:44])
        n_regions = int(line[44:55])
        n_pairs = int(line[55:66])
        params = [C1, C2, L1, L2]

        # Read the interpolation region data, namely NBT and INT
        nbt = np.zeros(n_regions)
        interp = np.zeros(n_regions)
        m = 0
        for i in range((n_regions - 1)//3 + 1):
            line = fh.readline()
            toRead = min(3, n_regions - m)
            for j in range(toRead):
                nbt[m] = int(line[0:11])
                interp[m] = int(line[11:22])
                line = line[22:]
                m += 1

        # Read tabulated pairs x(n) and y(n)
        x = np.zeros(n_pairs)
        y = np.zeros(n_pairs)
        m = 0
        for i in range((n_pairs - 1)//3 + 1):
            line = fh.readline()
            toRead = min(3, n_pairs - m)
            for j in range(toRead):
                x[m] = endftod(line[:11])
                y[m] = endftod(line[11:22])
                line = line[22:]
                m += 1

        return params, cls(x, y, nbt, interp)

    def __call__(self, x):
        # Check if input is array or scalar
        if isinstance(x, Iterable):
            iterable = True
            x = np.array(x)
        else:
            iterable = False
            x = np.array([x], dtype=float)

        # Create output array
        y = np.zeros_like(x)

        # Get indices for interpolation
        idx = np.searchsorted(self.x, x, side='right') - 1

        # Find lowest valid index
        i_low = np.searchsorted(idx, 0)

        for k in range(len(self.nbt)):
            # Determine which x values are within this interpolation range
            i_high = np.searchsorted(idx, self.nbt[k] - 1)

            # Get x values and bounding (x,y) pairs
            xk = x[i_low:i_high]
            xi = self.x[idx[i_low:i_high]]
            xi1 = self.x[idx[i_low:i_high] + 1]
            yi = self.y[idx[i_low:i_high]]
            yi1 = self.y[idx[i_low:i_high] + 1]

            if self.interp[k] == 1:
                # Histogram
                y[i_low:i_high] = yi

            elif self.interp[k] == 2:
                # Linear-linear
                y[i_low:i_high] = yi + (xk - xi)/(xi1 - xi)*(yi1 - yi)

            elif self.interp[k] == 3:
                # Linear-log
                y[i_low:i_high] = yi + np.log(xk/xi)/np.log(xi1/xi)*(yi1 - yi)

            elif self.interp[k] == 4:
                # Log-linear
                y[i_low:i_high] = yi*np.exp((xk - xi)/(xi1 - xi)*np.log(yi1/yi))

            elif self.interp[k] == 5:
                # Log-log
                y[i_low:i_high] = yi*np.exp(np.log(xk/xi)/np.log(xi1/xi)*np.log(yi1/yi))

            i_low = i_high

        # In some cases, the first/last point of x may be less than the first
        # value of self.x due only to precision, so we check if they're close
        # and set them equal if so. Otherwise, the interpolated value might be
        # out of range (and thus zero)
        if np.isclose(x[0], self.x[0], 1e-8):
            y[0] = self.y[0]
        if np.isclose(x[-1], self.x[-1], 1e-8):
            y[-1] = self.y[-1]

        return y if iterable else y[0]

    def integral(self):
        """Integral of the tabulated function over its tabulated range.

        Returns
        -------
        ndarray
            Array of same length as the tabulated data that represents partial
            integrals from the bottom of the range to each tabulated point.

        """

        # Create output array
        partial_sum = np.zeros(len(self.x) - 1)

        i_low = 0
        for k in range(len(self.nbt)):
            # Determine which x values are within this interpolation range
            i_high = self.nbt[k] - 1

            # Get x values and bounding (x,y) pairs
            x0 = self.x[i_low:i_high]
            x1 = self.x[i_low + 1:i_high + 1]
            y0 = self.y[i_low:i_high]
            y1 = self.y[i_low + 1:i_high + 1]

            if self.interp[k] == 1:
                # Histogram
                partial_sum[i_low:i_high] = y0*(x1 - x0)

            elif self.interp[k] == 2:
                # Linear-linear
                m = (y1 - y0)/(x1 - x0)
                partial_sum[i_low:i_high] = (y0 - m*x0)*(x1 - x0) + \
                                            m*(x1**2 - x0**2)/2

            elif self.interp[k] == 3:
                # Linear-log
                logx = np.log(x1/x0)
                m = (y1 - y0)/logx
                partial_sum[i_low:i_high] = y0 + m*(x1*(logx - 1) + x0)

            elif self.interp[k] == 4:
                # Log-linear
                m = np.log(y1/y0)/(x1 - x0)
                partial_sum[i_low:i_high] = y0/m*(np.exp(m*(x1 - x0)) - 1)

            elif self.interp[k] == 5:
                # Log-log
                m = np.log(y1/y0)/np.log(x1/x0)
                partial_sum[i_low:i_high] = y0/((m + 1)*x0**m)*(
                    x1**(m + 1) - x0**(m + 1))

            i_low = i_high

        return np.concatenate(([0.], np.cumsum(partial_sum)))


class EnergyDistribution(object):
    """Energy distribution from ENDF MF=5.

    Attributes
    ----------
    lf : int
        Flag specifying energy distribution.

    """
    def __init__(self):
        self.lf = None


class ArbitraryTabulated(EnergyDistribution):
    r"""Arbitrary tabulated function given in ENDF MF=5, LF=1 represented as

    .. math::
         f(E \rightarrow E') = g(E \rightarrow E')

    Attributes
    ----------
    energy : ndarray
        Array of incident neutron energies
    pdf : list of Tab1
        Tabulated outgoing energy distribution probability density functions

    """

    def __init__(self, energy, pdf):
        self.lf = 1
        self.energy = energy
        self.pdf = pdf


class GeneralEvaporation(EnergyDistribution):
    r"""General evaporation spectrum given in ENDF MF=5, LF=5 represented as

    .. math::
        f(E \rightarrow E') = g(E'/\theta(E))

    Attributes
    ----------
    theta : Tab1
        Tabulated function of incident neutron energy :math:`E`
    g : Tab1
        Tabulated function of :math:`x = E'/\theta(E)`
    u : float
        Constant introduced to define the proper upper limit for the final
        particle energy such that :math:`0 \le E' \le E - U`

    """

    def __init__(self, theta, g, u):
        self.lf = 5
        self.theta = theta
        self.g = g
        self.u = u


class Maxwellian(EnergyDistribution):
    r"""Simple Maxwellian fission spectrum given in ENDF MF=5, LF=7 represented
    as

    .. math::
        f(E \rightarrow E') = \frac{\sqrt{E'}}{I} e^{-E'/\theta(E)}

    Attributes
    ----------
    theta : Tab1
        Tabulated function of incident neutron energy
    u : float
        Constant introduced to define the proper upper limit for the final
        particle energy such that :math:`0 \le E' \le E - U`

    """

    def __init__(self, theta, u):
        self.lf = 7
        self.theta = theta
        self.u = u


class Evaporation(EnergyDistribution):
    r"""Evaporation spectrum given in ENDF MF=5, LF=9 represented as

    .. math::
        f(E \rightarrow E') = \frac{E'}{I} e^{-E'/\theta(E)}

    Attributes
    ----------
    theta : Tab1
        Tabulated function of incident neutron energy
    u : float
        Constant introduced to define the proper upper limit for the final
        particle energy such that :math:`0 \le E' \le E - U`

    """

    def __init__(self, theta, u):
        self.lf = 9
        self.theta = theta
        self.u = u


class Watt(EnergyDistribution):
    r"""Energy-dependent Watt spectrum given in ENDF MF=5, LF=11 represented as

    .. math::
        f(E \rightarrow E') = \frac{e^{-E'/a}}{I} \sinh \left ( \sqrt{bE'}
        \right )

    Attributes
    ----------
    a, b : Tab1
        Energy-dependent parameters tabulated as function of incident neutron
        energy
    u : float
        Constant introduced to define the proper upper limit for the final
        particle energy such that :math:`0 \le E' \le E - U`

    """

    def __init__(self, a, b, u):
        self.lf = 11
        self.a = a
        self.b = b
        self.u = u


class MadlandNix(EnergyDistribution):
    r"""Energy-dependent fission neutron spectrum (Madland and Nix) given in
    ENDF MF=5, LF=12 represented as

    .. math::
        f(E \rightarrow E') = \frac{1}{2} [ g(E', E_F(L)) + g(E', E_F(H))]

    where

    .. math::
        g(E',E_F) = \frac{1}{3\sqrt{E_F T_M}} \left [ u_2^{3/2} E_1 (u_2) -
        u_1^{3/2} E_1 (u_1) + \gamma \left ( \frac{3}{2}, u_2 \right ) - \gamma
        \left ( \frac{3}{2}, u_1 \right ) \right ] \\ u_1 = \left ( \sqrt{E'} -
        \sqrt{E_F} \right )^2 / T_M \\ u_2 = \left ( \sqrt{E'} + \sqrt{E_F}
        \right )^2 / T_M.

    Attributes
    ----------
    efl, efh : float
        Constants which represent the average kinetic energy per nucleon of the
        fission fragment (efl = light, efh = heavy)
    tm : Tab1
        Parameter tabulated as a function of incident neutron energy

    """

    def __init__(self, efl, efh, tm):
        self.lf = 12
        self.efl = efl
        self.efh = efh
        self.tm = tm


class ENDFTab2Record(object):
    def __init__(self):
        self.NBT = []
        self.INT = []

    def read(self, fh):
        # Determine how many interpolation regions and total points there are
        line = fh.readline()
        C1 = endftod(line[:11])
        C2 = endftod(line[11:22])
        L1 = int(line[22:33])
        L2 = int(line[33:44])
        NR = int(line[44:55])
        NZ = int(line[55:66])
        self.params = [C1, C2, L1, L2, NR, NZ]

        # Read the interpolation region data, namely NBT and INT
        m = 0
        for i in range((NR-1)//3 + 1):
            line = fh.readline()
            toRead = min(3,NR-m)
            for j in range(toRead):
                NBT = int(line[0:11])
                INT = int(line[11:22])
                self.NBT.append(NBT)
                self.INT.append(INT)
                line = line[22:]
            m = m + toRead


class AngularDistribution(object):
    pass


class Reaction(object):
    """Data for a single reaction including its cross section and secondary
    angle/energy distribution.

    Parameters
    ----------
    mt : int
        The MT number from the ENDF file.

    Attributes
    ----------
    angular_distribution : AngularDistribution
        Angular distribution represented as a tabulated function or as moments
        of a Legendre polynomial expansion from MF=4.
    complex_breakup_flag : int
        Complex breakup flag.
    xs : Tab1
        Tabulated cross section as a function of incident energy from MF=3.
    energy_distribution : list of EnergyDistribution
        List of partial energy distributions for the reaction from MF=5.
    files : list of int
        List of files (MF) that have been read for this reaction
    MT : int
        The MT number from the ENDF file.
    multiplicities : dict
        Multiplicities for production of radioactive nuclides as given in MF=9.
    product_distribution : list of dict
        Secondary energy or correlated energy-angle distribution from MF=6.
    production : dict
        Cross sections for production of radioactive nuclides as given in MF=10.
    Q_mass_difference : float
        Mass difference Q value in eV
    Q_reaction : float
        Reaction Q value in eV
    radionuclide_production : list of dict
        List of radioactive nuclides produced as given in MF=8.
    reference_frame : {'laboratory', 'center-of-mass', 'light-heavy'}
        Indicates what reference frame is used for outgoing energies and
        angles. Only relevant for product energy-angle distributions read in
        MF=6.
    subshell_binding_energy : float
        Subshell binding energy in eV as given in MF=23. Only relevant for
        photoelectric subshell ionization reactions in a photo-atomic
        sublibrary.
    fluorescence_yield : float
        Fluorescence yield (eV/photoionization) as given in MF=23. Only relevant
        for photoelectric subshell ionization reactions in a photo-atomic
        sublibrary.
    products : list of dict
        Secondary photon and electron distributions for electro-atomic
        reactions as given in MF=26.
    scattering_factor : Tab1
        Coherent or incoherent form factor as given in MF=27.
    anomalous_scattering_imaginary : Tab1
        Imaginary component of the anomalous scattering factor as given in
        MF=27.
    anomalous_scattering_real : Tab1
        Real component of the anomalous scattering factor as given in MF=27.

    """

    def __init__(self, mt):
        self.mt = mt
        self.files = []
        self.photon_production = {}
        self.xs = None
        self.Q_mass_difference = None
        self.Q_reaction = None
        self.complex_breakup_flag = None
        self.angular_distribution = None
        self.energy_distribution = []
        self.product_distribution = []
        self.reference_frame = None
        self.radionuclide_production = None
        self.multiplicities = {}
        self.production = {}
        self.subshell_binding_energy = None
        self.fluorescence_yield = None
        self.products = []
        self.scattering_factor = None
        self.anomalous_scattering_imaginary = None
        self.anomalous_scattering_real = None

    def __repr__(self):
        return '<ENDF Reaction: MT={0}, {1}>'.format(self.mt, label(self.mt))


class ResonanceRange(object):
    """Resolved resonance formalism data as given in MF=2, MT=151.

    Parameters
    ----------
    emin : float
        Minimum energy of the resolved resonance range in eV
    emax : float
        Maximum energy of the resolved resonance range in eV
    nro : int
        Flag designating energy-dependence of scattering radius (NRO). A value
        of 0 indicates it is energy-independent, a value of 1 indicates it is
        energy-dependent.
    naps : int
        Flag controlling use of the channel and scattering radius (NAPS)

    Attributes
    ----------
    energy_min : float
        Minimum energy of the resolved resonance range in eV
    energy_max : float
        Maximum energy of the resolved resonance range in eV
    nro : int
        Flag designating energy-dependence of scattering radius (NRO). A value
        of 0 indicates it is energy-independent, a value of 1 indicates it is
        energy-dependent.
    naps : int
        Flag controlling use of the channel and scattering radius (NAPS)

    """

    def __init__(self, emin, emax, nro, naps):
        self.energy_min = emin
        self.energy_max = emax
        self.nro = nro
        self.naps = naps
        self._prepared = False

    def read(self, ev):
        raise NotImplementedError


class MultiLevelBreitWigner(ResonanceRange):
    """Multi-level Breit-Wigner resolved resonance formalism data. This is
    identified by LRF=2 in the ENDF-6 format.

    Parameters
    ----------
    emin : float
        Minimum energy of the resolved resonance range in eV
    emax : float
        Maximum energy of the resolved resonance range in eV
    nro : int
        Flag designating energy-dependence of scattering radius (NRO). A value
        of 0 indicates it is energy-independent, a value of 1 indicates it is
        energy-dependent.
    naps : int
        Flag controlling use of the channel and scattering radius (NAPS)

    Attributes
    ----------
    spin : float
        Spin of the target nucleus
    scattering_radius_ind : float
        Scattering radius in units of 10^-12 cm
    resonances : list of Resonance
        List of resolved resonances
    l_values : ndarray of int
        Neutron orbital angular momentum values
    q : ndarray of float
        Q-value to be added to incident particle's center-of-mass energy to
        determine the channel energy for use in the penetrability factor. Given
        as a function of the l-value.
    competitive : ndarray of bool
        Whether a competitive width is given for each l-value.

    """

    def __init__(self, emin, emax, nro, naps):
        super(MultiLevelBreitWigner, self).__init__(emin, emax, nro, naps)
        self.spin = None
        self.scattering_radius_ind = None
        self.resonances = []
        self.l_values = None
        self.q = None
        self.competitive = None

    def read(self, ev):
        """Read multi-level Breit-Wigner resonance data from an ENDF evaluation.

        Parameters
        ----------
        ev : Evaluation
            ENDF evaluation positioned at the start of a resonance range
            subsection

        """

        # Read energy-dependent scattering radius if present
        if self.nro != 0:
            params, self.scattering_radius = ev._get_tab1_record()

        # Other scatter radius parameters
        items = ev._get_cont_record()
        self.spin = items[0]
        self.scattering_radius_ind = items[1]
        NLS = items[4]  # Number of l-values

        self.l_values = np.zeros(NLS, int)
        self.q = np.zeros(NLS)
        self.competitive = np.zeros(NLS, bool)

        # Read resonance widths, J values, etc
        for l in range(NLS):
            items, values = ev._get_list_record()
            self.q[l] = items[1]
            self.l_values[l] = items[2]
            self.competitive[l] = items[3]
            energy = values[0::6]
            spin = values[1::6]
            GT = values[2::6]
            GN = values[3::6]
            GG = values[4::6]
            GF = values[5::6]

            resonances = []
            for i, E in enumerate(energy):
                resonance = Resonance()
                resonance.energy = energy[i]
                resonance.spin = spin[i]
                resonance.width_total = GT[i]
                resonance.width_neutron = GN[i]
                resonance.width_gamma = GG[i]
                resonance.width_fissionA = GF[i]
                resonances.append(resonance)
            self.resonances.append(resonances)



class SingleLevelBreitWigner(MultiLevelBreitWigner):
    """Single-level Breit-Wigner resolved resonance formalism data. This is
    identified by LRF=1 in the ENDF-6 format.

    Parameters
    ----------
    emin : float
        Minimum energy of the resolved resonance range in eV
    emax : float
        Maximum energy of the resolved resonance range in eV
    nro : int
        Flag designating energy-dependence of scattering radius (NRO). A value
        of 0 indicates it is energy-independent, a value of 1 indicates it is
        energy-dependent.
    naps : int
        Flag controlling use of the channel and scattering radius (NAPS)

    """

    def __init__(self, emin, emax, nro, naps):
        super(SingleLevelBreitWigner, self).__init__(emin, emax, nro, naps)


class ReichMoore(ResonanceRange):
    """Reich-Moore resolved resonance formalism data. This is identified by LRF=3 in
    the ENDF-6 format.

    Parameters
    ----------
    emin : float
        Minimum energy of the resolved resonance range in eV
    emax : float
        Maximum energy of the resolved resonance range in eV
    nro : int
        Flag designating energy-dependence of scattering radius (NRO). A value
        of 0 indicates it is energy-independent, a value of 1 indicates it is
        energy-dependent.
    naps : int
        Flag controlling use of the channel and scattering radius (NAPS)

    Attributes
    ----------
    resonances : dict
        Dictionary whose keys are (l, J) pairs and values are instances are
        lists of Resonance instances.
    spin : float
        Spin of the target nucleus
    scattering_radius_ind : float
        Energy-independent scattering radius in units of 10^-12 cm
    scattering_radius : Tab1
        Scattering radius in units of 10^-12 cm as a function of energy
    LAD : int
        Indicate whether parameters can be used to compute angular distributions
    NLSC : int
        Number of l-values which must be used to converge the calculation
    apl : ndarray of float
        l-dependent scattering radius
    l_values : ndarray of int
        Neutron orbital angular momentum values

    """

    def __init__(self, emin, emax, nro, naps):
        super(ReichMoore, self).__init__(emin, emax, nro, naps)
        self.resonances = {}
        self.spin = None
        self.scattering_radius_ind = None
        self.scattering_radius = None
        self.LAD = None
        self.NLSC = None
        self.apl = None
        self.l_values = None

    def read(self, ev):
        """Read Reich-Moore resonance data from an ENDF evaluation.

        Parameters
        ----------
        ev : Evaluation
            ENDF evaluation positioned at the start of a resonance range
            subsection

        """
        # Read energy-dependent scattering radius if present
        if self.nro != 0:
            params, self.scattering_radius = ev._get_tab1_record()

        # Other scatter radius parameters
        items = ev._get_cont_record()
        self.spin = items[0]
        self.scattering_radius_ind = items[1]
        self.LAD = items[3]  # Flag for angular distribution
        NLS = items[4]  # Number of l-values
        self.NLSC = items[5]  # Number of l-values for convergence

        self.apl = np.zeros(NLS)
        self.l_values = np.zeros(NLS, int)

        # Read resonance widths, J values, etc
        for i in range(NLS):
            items, values = ev._get_list_record()
            self.apl[i] = items[1]
            self.l_values[i] = l = items[2]
            energy = values[0::6]
            spin = values[1::6]
            GN = values[2::6]
            GG = values[3::6]
            GFA = values[4::6]
            GFB = values[5::6]
            for i, E in enumerate(energy):
                resonance = Resonance()
                resonance.energy = energy[i]
                resonance.spin = J = spin[i]
                resonance.width_neutron = GN[i]
                resonance.width_gamma = GG[i]
                resonance.width_fissionA = GFA[i]
                resonance.width_fissionB = GFB[i]

                if not (l, abs(J)) in self.resonances:
                    self.resonances[l, abs(J)] = []
                self.resonances[l, abs(J)].append(resonance)


class AdlerAdler(ResonanceRange):
    """Adler-Adler resolved resonance formalism data. This is identified by LRF=4 in
    the ENDF-6 format.

    Parameters
    ----------
    emin : float
        Minimum energy of the resolved resonance range in eV
    emax : float
        Maximum energy of the resolved resonance range in eV
    nro : int
        Flag designating energy-dependence of scattering radius (NRO). A value
        of 0 indicates it is energy-independent, a value of 1 indicates it is
        energy-dependent.
    naps : int
        Flag controlling use of the channel and scattering radius (NAPS)

    Attributes
    ----------
    spin : float
        Spin of the target nucleus
    scattering_radius_ind : float
        Scattering radius in units of 10^-12 cm
    LI : float
        Flag indicating the kind of parameters
    AT, BT : ndarray
        Background constants for the total cross section
    AC, BC : ndarray
        Background constants for the radiative cross section
    AF, BF : ndarray
        Background constants for the fission cross section
    resonances : list of AdlerResonance
        List of resonances defined by Adler-Adler parameters

    """

    def __init__(self, emin, emax, nro, naps):
        super(AdlerAdler, self).__init__(emin, emax, nro, naps)

    def read(self, ev):
        """Read Adler-Adler resonance data from an ENDF evaluation.

        Parameters
        ----------
        ev : Evaluation
            ENDF evaluation positioned at the start of a resonance range
            subsection

        """
        # Read energy-dependent scattering radius if present
        if self.nro != 0:
            params, self.scattering_radius = ev._get_tab1_record()

        # Other scatter radius parameters
        items = ev._get_cont_record()
        self.spin = items[0]
        self.scattering_radius_ind = items[1]
        NLS = items[4]  # Number of l-values

        # Get AT, BT, AF, BF, AC, BC constants
        items, values = ev._get_list_record()
        self.LI = items[2]
        NX = items[5]
        self.AT = np.asarray(values[:4])
        self.BT = np.asarray(values[4:6])
        if NX == 2:
            self.AC = np.asarray(values[6:10])
            self.BC = np.asarray(values[10:12])
        elif NX == 3:
            self.AF = np.asarray(values[6:10])
            self.BF = np.asarray(values[10:12])
            self.AC = np.asarray(values[12:16])
            self.BC = np.asarray(values[16:18])

        self.resonances = []

        for ls in range(NLS):
            items = ev._get_cont_record()
            l_value = items[2]
            NJS = items[4]
            for j in range(NJS):
                items, values = ev._get_list_record()
                AJ = items[0]
                NLJ = items[5]
                for res in range(NLJ):
                    resonance = AdlerResonance()
                    resonance.L, resonance.J = l_value, AJ
                    resonance.DET = values[12*res]
                    resonance.DWT = values[12*res + 1]
                    resonance.DRT = values[12*res + 2]
                    resonance.DIT = values[12*res + 3]
                    resonance.DEF_ = values[12*res + 4]
                    resonance.DWF = values[12*res + 5]
                    resonance.GRF = values[12*res + 6]
                    resonance.GIF = values[12*res + 7]
                    resonance.DEC = values[12*res + 8]
                    resonance.DWC = values[12*res + 9]
                    resonance.GRC = values[12*res + 10]
                    resonance.GIC = values[12*res + 11]
                    self.resonances.append(resonance)


class RMatrixLimited(ResonanceRange):
    """R-Matrix Limited resolved resonance formalism data. This is identified by
    LRF=7 in the ENDF-6 format.

    Parameters
    ----------
    emin : float
        Minimum energy of the resolved resonance range in eV
    emax : float
        Maximum energy of the resolved resonance range in eV
    nro : int
        Flag designating energy-dependence of scattering radius (NRO). A value
        of 0 indicates it is energy-independent, a value of 1 indicates it is
        energy-dependent.
    naps : int
        Flag controlling use of the channel and scattering radius (NAPS)

    Attributes
    ----------
    IFG : int
        Flag indicating whether channel widths in eV or reduced-width amplitudes
        in eV^1/2 are given
    KRM : int
        Flag to specify which formulae for the R-matrix are to be used
    particle_pairs : list of dict
        List of particle pairs. Each particle pair is represented by a
        dictionary that contains the mass, atomic number, spin, and parity of
        each particle as well as other characteristics.
    spin_groups : list of dict
        List of spin groups. Each spin group is characterized by channels,
        resonance energies, and resonance widths.

    """

    def __init__(self, emin, emax, nro, naps):
        super(RMatrixLimited, self).__init__(emin, emax, nro, naps)
        self.IFG = None
        self.KRM = None
        self.particle_pairs = []
        self.spin_groups = []

    def read(self, ev):
        """Read R-Matrix limited resonance data from an ENDF evaluation.

        Parameters
        ----------
        ev : Evaluation
            ENDF evaluation positioned at the start of a resonance range
            subsection

        """
        items = ev._get_cont_record()
        self.IFG = items[2]  # reduced width amplitude?
        self.KRM = items[3]  # Specify which formulae are used
        n_spin_groups = items[4]  # Number of Jpi values (NJS)
        KRL = items[5]  # Flag for non-relativistic kinematics

        items, values = ev._get_list_record()
        n_pairs = items[5]//2  # Number of particle pairs (NPP)
        for i in range(n_pairs):
            pp = {'mass_a': values[12*i],
                  'mass_b': values[12*i + 1],
                  'z_a': values[12*i + 2],
                  'z_b': values[12*i + 3],
                  'spin_a': values[12*i + 4],
                  'spin_b': values[12*i + 5],
                  'q': values[12*i + 6],
                  'pnt': values[12*i + 7],
                  'shift': values[12*i + 8],
                  'MT': values[12*i + 9],
                  'parity_a': values[12*i + 10],
                  'parity_b': values[12*i + 11]}
            self.particle_pairs.append(pp)

        # loop over spin groups
        for i in range(n_spin_groups):
            sg = {'channels': []}
            self.spin_groups.append(sg)

            items, values = ev._get_list_record()
            J = items[0]
            parity = items[1]
            kbk = items[2]
            kps = items[3]
            n_channels = items[5]
            for j in range(n_channels):
                channel = {}
                channel['particle_pair'] = values[6*j]
                channel['l'] = values[6*j + 1]
                channel['spin'] = values[6*j + 2]
                channel['boundary'] = values[6*j + 3]
                channel['effective_radius'] = values[6*j + 4]
                channel['true_radius'] = values[6*j + 5]
                sg['channels'].append(channel)

            items, values = ev._get_list_record()
            n_resonances = items[3]
            sg['resonance_energies'] = np.zeros(n_resonances)
            sg['resonance_widths'] = np.zeros((n_channels, n_resonances))
            m = n_channels//6 + 1
            for j in range(n_resonances):
                sg['resonance_energies'][j] = values[m*j]
                for k in range(n_channels):
                    sg['resonance_widths'][k,j] = values[m*j + k + 1]

            # Optional extension (Background R-Matrix)
            if kbk > 0:
                items, values = ev._get_list_record()
                lbk = items[4]
                if lbk == 1:
                    params, rbr = ev._get_tab1_record()
                    params, rbi = ev._get_tab1_record()

            # Optional extension (Tabulated phase shifts)
            if kps > 0:
                items, values = ev._get_list_record()
                lps = items[4]
                if lps == 1:
                    params, psr = ev._get_tab1_record()
                    params, psi = ev._get_tab1_record()


class Unresolved(ResonanceRange):
    """Unresolved resonance parameters as identified by LRU=2 in MF=2.

    Parameters
    ----------
    emin : float
        Minimum energy of the unresolved resonance range in eV
    emax : float
        Maximum energy of the unresolved resonance range in eV
    nro : int
        Flag designating energy-dependence of scattering radius (NRO). A value
        of 0 indicates it is energy-independent, a value of 1 indicates it is
        energy-dependent.
    naps : int
        Flag controlling use of the channel and scattering radius (NAPS)

    Attributes
    ----------
    scattering_radius : float
        Scattering radius in units of 10^-12 cm
    LSSF : int
        Flag governing interpretation of file 3 cross sections
    spin : float
        Spin of the target nucleus
    l_values : ndarray of int
        Neutron orbital angular momentum values
    parameters : dict of dict
        Dictionary whose keys are l-values and whose values are dictionaries
        containing unresolved resonance parameters

    """

    def __init__(self, emin, emax, nro, naps):
        super(Unresolved, self).__init__(emin, emax, nro, naps)

    def read(self, ev):
        """Read unresolved resonance data from an ENDF evaluation.

        Parameters
        ----------
        ev : Evaluation
            ENDF evaluation positioned at the start of a resonance range
            subsection

        """
        # Read energy-dependent scattering radius if present
        if self.nro != 0:
            params, self.scattering_radius = ev._get_tab1_record()

        # Get SPI, AP, and LSSF
        if not (self.fission_widths and self.LRF == 1):
            items = ev._get_cont_record()
            self.spin = items[0]
            if self.nro == 0:
                self.scattering_radius = items[1]
            self.LSSF = items[2]

        if not self.fission_widths and self.LRF == 1:
            # Case A -- fission widths not given, all parameters are
            # energy-independent
            NLS = items[4]
            self.l_values = np.zeros(NLS)
            self.parameters = {}
            for ls in range(NLS):
                items, values = ev._get_list_record()
                l = items[2]
                NJS = items[5]
                self.l_values[ls] = l
                params = {}
                self.parameters[l] = params
                params['d'] = np.asarray(values[0::6])
                params['j'] = np.asarray(values[1::6])
                params['amun'] = np.asarray(values[2::6])
                params['gn0'] = np.asarray(values[3::6])
                params['gg'] = np.asarray(values[4::6])
                # params['gf'] = np.zeros(NJS)

        elif self.fission_widths and self.LRF == 1:
            # Case B -- fission widths given, only fission widths are
            # energy-dependent
            items, self.energies = ev._get_list_record()
            self.spin = items[0]
            if self.nro == 0:
                self.scatter_radius = items[1]
            self.LSSF = items[2]
            NE, NLS = items[4:6]
            self.l_values = np.zeros(NLS, int)
            self.parameters = {}
            for ls in range(NLS):
                items = ev._get_cont_record()
                l = items[2]
                NJS = items[4]
                self.l_values[ls] = l
                params = {}
                self.parameters[l] = params
                params['d'] = np.zeros(NJS)
                params['j'] = np.zeros(NJS)
                params['amun'] = np.zeros(NJS)
                params['gn0'] = np.zeros(NJS)
                params['gg'] = np.zeros(NJS)
                params['gf'] = []
                for j in range(NJS):
                    items, values = ev._get_list_record()
                    muf = items[3]
                    params['d'][j] = values[0]
                    params['j'][j] = values[1]
                    params['amun'][j] = values[2]
                    params['gn0'][j] = values[3]
                    params['gg'][j] = values[4]
                    params['gf'].append(np.asarray(values[6:]))

        elif self.LRF == 2:
            # Case C -- all parameters are energy-dependent
            NLS = items[4]
            self.l_values = np.zeros(NLS)
            self.parameters = {}
            for ls in range(NLS):
                items = ev._get_cont_record()
                l = items[2]
                NJS = items[4]
                self.l_values[ls] = l
                params = {}
                self.parameters[l] = params
                params['j'] = np.zeros(NJS)
                params['amux'] = np.zeros(NJS)
                params['amun'] = np.zeros(NJS)
                params['amug'] = np.zeros(NJS)
                params['amuf'] = np.zeros(NJS)
                params['energies'] = []
                params['d'] = []
                params['gx'] = []
                params['gn0'] = []
                params['gg'] = []
                params['gf'] = []
                for j in range(NJS):
                    items, values = ev._get_list_record()
                    ne = items[5]
                    params['j'][j] = items[0]
                    params['amux'][j] = values[2]
                    params['amun'][j] = values[3]
                    params['amug'][j] = values[4]
                    params['amuf'][j] = values[5]
                    params['energies'].append(np.asarray(values[6::6]))
                    params['d'].append(np.asarray(values[7::6]))
                    params['gx'].append(np.asarray(values[8::6]))
                    params['gn0'].append(np.asarray(values[9::6]))
                    params['gg'].append(np.asarray(values[10::6]))
                    params['gf'].append(np.asarray(values[11::6]))


class ScatteringRadius(ResonanceRange):
    """Energy range with no resonances and only a scattering radius.

    Parameters
    ----------
    emin : float
        Minimum energy of the resolved resonance range in eV
    emax : float
        Maximum energy of the resolved resonance range in eV
    nro : int
        Flag designating energy-dependence of scattering radius (NRO). A value
        of 0 indicates it is energy-independent, a value of 1 indicates it is
        energy-dependent.
    naps : int
        Flag controlling use of the channel and scattering radius (NAPS)

    """
    def __init__(self, emin, emax, nro, naps):
        super(ScatteringRadius, self).__init__(emin, emax, nro, naps)


_formalisms = {1: SingleLevelBreitWigner, 2: MultiLevelBreitWigner,
               3: ReichMoore, 4: AdlerAdler, 7: RMatrixLimited}


cdef class Resonance(object):
    # Attributes from ENDF file
    cdef public double energy
    cdef public double spin
    cdef public double width_total
    cdef public double width_neutron
    cdef public double width_gamma
    cdef public double width_fissionA
    cdef public double width_fissionB


class AdlerResonance(object):
    def __init__(self):
        pass
