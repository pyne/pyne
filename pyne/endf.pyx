#!/usr/bin/env python

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

import re
import os
from libc.stdlib cimport malloc, free

cimport numpy as np
import numpy as np

np.import_array()

from pyne cimport cpp_nucname

from math import e

from pyne import nucname
import pyne.rxdata as rx
from pyne.rxname import label
from pyne.utils import fromendf_tok, endftod

include "include/cython_version.pxi"
IF CYTHON_VERSION_MAJOR == 0 and CYTHON_VERSION_MINOR >= 17:
    from libc.stdlib cimport atof, atoi
    from libc.string cimport strtok, strcpy, strncpy
ELSE:
    from pyne._includes.libc.stdlib cimport atof, atoi
    from pyne._includes.libc.string cimport strtok, strcpy, strncpy

libraries = {0: "ENDF/B", 1: "ENDF/A", 2: "JEFF", 3: "EFF",
             4: "ENDF/B High Energy", 5: "CENDL", 6: "JENDL",
             31: "INDL/V", 32: "INDL/A", 33: "FENDL", 34: "IRDF",
             35: "BROND", 36: "INGDB-90", 37: "FENDL/A", 41: "BROND"}
FILE1_R = re.compile(r'1451 *\d{1,5}$')
CONTENTS_R = re.compile(' +\d{1,2} +\d{1,3} +\d{1,10} +')
SPACE66_R = re.compile(' {66}')
NUMERICAL_DATA_R = re.compile('[\d\-+. ]{80}\n$')
SPACE66_R = re.compile(' {66}')

class Library(rx.RxLib):
    "A class for a file which contains multiple ENDF evaluations."
    def __init__(self, fh):
        self.mts = {}
        self.structure = {}
        self.mat_dict = {}
        self.more_files = True
        self.intdict = {1: self._histogram, 2: self._linlin, 3: self._linlog, 4:
                        self._loglin, 5: self._loglog, 6:self._chargedparticles,
                        11: self._histogram, 12: self._linlin, 13: self._linlog,
                        14: self._loglin, 15: self._loglog, 21: self._histogram,
                        22: self._linlin, 23: self._linlog, 24: self._loglin,
                        25: self._loglog}
        self.chars_til_now = 0
        self.offset = 0
        self.fh = fh
        while self.more_files:
            self._read_headers()


    def load(self):
        """Read the ENDF file into a NumPy array.

        Returns
        --------
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
        # Skip the first line and get the material ID.
        fh.seek(self.chars_til_now)
        len_headline = len(fh.readline())
        self.offset += 81 - len_headline
        line = fh.readline()
        mat_id = int(line[66:70].strip() or -1)
        # originally in a float version of ZZAAA.M, ie 94242.1
        nuc = cpp_nucname.id(<int> (endftod(line[:11])*10))
        # Make a new dict in self.structure to contain the material data.
        if nuc not in self.structure:
            self.structure.update(
                {nuc:{'styles': "", 'docs': [], 'particles': [], 'data': {},
                         'matflags': {}}})
            self.mat_dict.update({nuc:{'end_line':[],
                                          'mfs':{}}})
        # Parse header (all lines with 1451)
        mf = 1
        stop = (self.chars_til_now+self.offset)/81
        while FILE1_R.search(line):
            # parse contents section
            if CONTENTS_R.match(line):
                # When MF and MT change, add offset due to SEND/FEND records.
                old_mf = mf
                mf, mt = int(line[22:33]), int(line[33:44])
                mt_length = int(line[44:55])
                if old_mf == mf:
                    start = stop + 1
                else:
                    start = stop + 2
                stop = start + mt_length
                self.mat_dict[nuc]['mfs'][mf,mt] = (81*start-self.offset,
                                                    81*stop-self.offset)
                line = fh.readline()
            # parse comment
            elif SPACE66_R.match(line):
                self.structure[nuc]['docs'].append(line[0:66])
                line = fh.readline()
            elif NUMERICAL_DATA_R.match(line):
                line = fh.readline()
                continue
            else:
                self.structure[nuc]['docs'].append(line[0:66])
                line = fh.readline()
        # Find where the end of the material is and then jump to it.
        self.chars_til_now = (stop + 4)*81 - self.offset
        fh.seek(self.chars_til_now)
        nextline = fh.readline()
        self.more_files = (nextline != '' and nextline[68:70] != "-1")
        # Update materials dict
        if mat_id != -1:
            self.mat_dict[nuc]['end_line'] = (self.chars_til_now+self.offset)/81
            setattr(self, "mat{0}".format(nuc), self.structure[nuc])
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
                    dict(zip(keys,line))))

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
        headlines = (len(headkeys)-1)/6 + 1
        arraylines = (npl-1)/6 + 1
        if len(itemkeys) == 1:
            array_len = npl - (headlines-1) * 6
            items={itemkeys[0]: lines[headlines:].flat[:array_len]}
        else:
            array_width = ((len(itemkeys)-1)/6 + 1)*6
            items_transposed = np.transpose(
                lines[headlines:headlines+arraylines].reshape(-1,
                                                              array_width))
            items = dict(zip(itemkeys, items_transposed))
        if 0 in items:
            del items[0]

        total_lines = 1+arraylines
        return head, items, total_lines

    def _get_tab1(self, headkeys, xykeys,lines):
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
        meta_len = (nr*2-1)/6 + 1
        data_len = (np_*2-1)/6 + 1
        intmeta = dict(zip(('intpoints','intschemes'),
                           (lines[1:1+meta_len].flat[:nr*2:2],
                            lines[1:1+meta_len].flat[1:nr*2:2])))
        intdata = dict(zip(xykeys,
            (lines[1+meta_len:1+meta_len+data_len].flat[:np_*2:2],
             lines[1+meta_len:1+meta_len+data_len].flat[1:np_*2:2])))
        intdata.update(intmeta)
        total_lines = 1 + meta_len + data_len
        return head, intdata, total_lines

    def _histogram(self, Eint, xs):
        dEint = float(Eint[-1]-Eint[0])
        return np.nansum((Eint[1:]-Eint[:-1]) * xs[:-1]/dEint)

    def _linlin(self, Eint, xs):
        dEint = float(Eint[-1]-Eint[0])
        return np.nansum((Eint[1:]-Eint[:-1])* (xs[1:]+xs[:-1])/2./dEint)

    def _linlog(self, Eint, xs):
        dEint = float(Eint[-1]-Eint[0])
        x1 = Eint[:-1]
        x2 = Eint[1:]
        y1 = xs[:-1]
        y2 = xs[1:]
        A = (y1-y2)/(np.log(x1/x2))
        B = y1-A*np.log(x1)
        return np.nansum(A*(x2*np.log(x2) - x1*np.log(x1)-x2+x1) + B*(x2-x1))/dEint

    def _loglin(self, Eint, xs):
        dEint = float(Eint[-1]-Eint[0])
        x1 = Eint[:-1]
        x2 = Eint[1:]
        y1 = xs[:-1]
        y2 = xs[1:]
        A = (np.log(y1)-np.log(y2))/(x1-x2)
        B = np.log(y1) - A*x1
        return np.nansum((y2-y1)/A)/dEint

    def _loglog(self, Eint, xs):
        dEint = float(Eint[-1]-Eint[0])
        x1 = Eint[:-1]
        x2 = Eint[1:]
        y1 = xs[:-1]
        y2 = xs[1:]
        A = - np.log(y2/y1)/np.log(x1/x2)
        B = - (np.log(y1)*np.log(x2) - np.log(y2)*np.log(x1))/np.log(x1/x2)
        return np.nansum(e**B / (A+1) * (x2**(A+1) - x1**(A+1))/dEint)

    def _chargedparticles(self, Eint, xs, flags=None):
        q = flags['Q']
        if q > 0:
            T = 0
        else:
            T = q
        dEint = float(Eint[-1]-Eint[0])
        x1 = Eint[:-1]
        x2 = Eint[1:]
        y1 = xs[:-1]
        y2 = xs[1:]
        B = np.log(y2*x2/(x1*y1)) / (1/(x1-T)**0.5 - 1/(x2-T)**0.5)
        A = e**(B/(x1-T)**0.5)*y1*x1
        # FIXME
        raise NotImplementedError("see docs for more details.")

    def integrate_tab_range(self, intscheme, Eint, xs):
        """Integrates across one tabulation range.

        Parameters
        ----------
        intscheme : int or float
            The interpolation scheme used in this range.
        Eint : array
            The energies at which we have xs data.
        xs : array
            The xs data corresponding to Eint.

        Returns
        -------
        sigma_g : float
            The group xs.
        """
        with np.errstate(divide="ignore", invalid="ignore"):
           return self.intdict[intscheme](Eint, xs)

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
        """Read the resonance data from one material in the library and updates
        self.structure.

        Parameters
        -----------
        mat_id: int
            Material id .
        """
        lrp = self.structure[mat_id]['matflags']['LRP']
        if (lrp == -1 or mat_id in (-1,0)):
            # If the LRP flag for the material is -1, there's no resonance data.
            # Also if the mat id is invalid.
            pass
        else:
            # Load the resonance data.
            mf2 = self.get_rx(mat_id,2,151).reshape(-1, 6)

            self.structure[mat_id]['matflags'].update(
                self._get_head(['ZA','AWR',0,0,'NIS',0], mf2[0]))
            total_lines = 1
            for isotope_num in range(
                    int(self.structure[mat_id]['matflags']['NIS'])):
                total_lines += self._read_nis(mf2[total_lines:], lrp, mat_id)
        for isotope in self.structure[mat_id]['data'].values():
            isotope['resolved'].sort()
            isotope['unresolved'].sort()

    def _read_nis(self, isotope_data, lrp, mat_id):
        """Read resonance data for a specific isotope.

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
        isotope_flags = self._get_cont(['ZAI','ABN',0,'LFW','NER',0],
                                       isotope_data[0])
        nuc_i = int(isotope_flags['ZAI']*10)
        self.structure[mat_id]['data'].update(
            {nuc_i:{'resolved':[],
                       'unresolved':[],
                       'datadocs':[],
                       'xs':{},
                       'output':{'channel1':[],
                                 'channel2':[]},
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
        range_flags = self._get_cont(('EL','EH','LRU','LRF','NRO','NAPS'),
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
                lbk_list_keys = {2: ('R0','R1','R2','S0','S1',0),
                                 3: ('R0','SO','GA',0,0,0)}
                aj_data['ch{}'.format(ch)] = {'LBK': lbk}
                ch_data = aj_data['ch{}'.format(ch)]
                if lbk == 0:
                    total_lines += 2
                elif lbk == 1:
                    total_lines += 2
                    rbr, rbr_size = self._get_tab1(
                        (0,0,0,0,'NR','NP'), ('Eint','RBR'),
                        subsection[total_lines:])[1:3]
                    total_lines += rbr_size
                    ch_data['RBR'] = rbr
                    rbi, rbi_size = self._get_tab1(
                        (0,0,0,0,'NR','NP'), ('Eint','RBI'),
                        (subsection[total_lines:]))[1:3]
                    total_lines += rbi_size
                    ch_data['RBI'] = rbi
                else:
                    ch_data, total_lines = self._cont_and_update(
                        ch_data, ('ED','EU',0,0,'LBK',0), subsection,
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
                        (0,0,0,0,'NR','NP'), ('Eint','PSR'),
                        subsection[total_lines:])[1:3]
                    total_lines += psr_size
                    ch_data['PSR'] = psr
                    psi, psi_size = self._get_tab1(
                        (0,0,0,0,'NR','NP'), ('Eint','PSI'),
                        (subsection[total_lines:]))[1:3]
                    total_lines += psi_size
                    ch_data['PSI'] = psi
                    total_lines += psi_size
            return total_lines

        lrf = int(range_flags['LRF'])
        subsection_dict = rx.DoubleSpinDict({})
        headers = [None,
                   ('SPI','AP',0,0,'NLS',0),
                   ('SPI','AP',0,0,'NLS',0),
                   ('SPI','AP','LAD',0,'NLS','NLSC'),
                   ('SPI','AP',0,0,'NLS',0),
                   None,
                   None,
                   (0,0,'IFG','KRM','NJS','KRL')]
        if range_flags['NRO'] > 0:
            intdata, total_lines = self._get_tab1((0,0,0,0,'NR','NP'),
                                                  ('E','AP'),
                                                  subsection)[1:3]
            subsection_dict['int'] = intdata
        else:
            total_lines = 0
        range_flags, total_lines = self._cont_and_update(
                range_flags, headers[lrf], subsection, total_lines)

        lrf_L_keys = [None,
                      ('AWRI','QX','L','LRX','6*NRS','NRS'),
                      ('AWRI','QX','L','LRX','6*NRS','NRS'),
                      ('AWRI','APL','L',0,'6*NRS','NRS'),
                      (0,0,'L',0,'NJS',0)]
        lrf_J_keys = [None, None, None, None, ('AJ',0,0,0,'12*NLJ','NLJ')]
        lrf_itemkeys = [None,
                        ('ER','AJ','GT','GN','GG','GF'),
                        ('ER','AJ','GT','GN','GG','GF'),
                        ('ER','AJ','GN','GG','GFA','GFB'),
                        ('DET','DWT','GRT','GIT','DEF','DWF','GRF','GIF','DEC',
                         'DWC','GRC','GIC')]
        if lrf == 4:
            # Adler-Adler
            bg_flags, bg, bg_size = self._get_list(
                ('AWRI',0,'LI',0,'6*NX','NX'),
                ('A1','A2','A3','A4','B1','B2'),
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
                (0,0,'NPP',0,'12*NPP','2*NPP'),
                ('MA','MB','ZA','ZB','IA','IB','Q','PNT','SHF','MT','PA','PB'),
                subsection[total_lines:])[1:3]
            total_lines += pp_size
            range_flags.update(particle_pair_data)
            for aj_section in range(int(range_flags['NJS'])):
                # Read first LIST record, with channel descriptions
                aj_flags, ch_items, ch_size = self._get_list(
                    ('AJ','PJ','KBK','KPS','6*NCH','NCH'),
                    ('IPP','L','SCH','BND','APE','APT'),
                    subsection[total_lines:])
                total_lines += ch_size
                # Second LIST record, with resonance energies and widths.
                er_flags, er_data, er_size = self._get_list(
                    (0,0,0,'NRS','6*NX','NX'), ('ER',), subsection[total_lines:])
                total_lines += er_size
                nch = int(aj_flags['NCH'])
                er_array_width = (nch/6+1)*6
                er_data = er_data['ER'].reshape(-1,er_array_width).transpose()
                aj_data = {'ER': er_data[0], 'GAM': er_data[1:1+nch].transpose()}
                aj_data.update(ch_items)
                aj = aj_flags['AJ']
                # Additional records
                if aj_flags['KBK'] > 0:
                    lbk_list_keys = ((),(),#('ED','EU',0,0,'LBK',0),
                                     ('R0','R1','R2','S0','S1',0),
                                     ('R0','SO','GA',0,0,0))
                    total_lines = read_kbks(nch, subsection, aj_data, total_lines)
                if aj_flags['KPS'] > 0:
                    total_lines = read_kpss(nch, subsection, aj_data, total_lines)
                subsection_dict[aj] = aj_data

        el, eh = range_flags['EL'], range_flags['EH']
        subsection_data = (el,eh,subsection_dict,range_flags)
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
        head_cont = ('SPI','AP','LSSF',0,'NLS',0)
        has_head_cont = {(0,1): True, (1,1): False, (0,2): True, (1,2): True}
        L_keys = {(0,1): ('AWRI',0,'L',0,'6*NJS','NJS'),
                  (1,1): ('AWRI',0,'L',0,'NJS',0),
                  (0,2): ('AWRI',0,'L',0,'NJS',0),
                  (1,2): ('AWRI',0,'L',0,'NJS',0)}
        j_keys = {(0,1): None,
                  (1,1): (0,0,'L','MUF','NE+6',0,'D','AJ','AMUN','GN0','GG',
                          0),
                  (0,2): ('AJ',0,'INT',0,'6*NE+6','NE',0,0,'AMUX','AMUN',
                      'AMUG','AMUF'),
                  (1,2): ('AJ',0,'INT',0,'6*NE+6','NE',0,0,'AMUX','AMUN',
                      'AMUG','AMUF')}
        itemkeys = {(0,1): ('D','AJ','AMUN','GN0','GG',0),
                    (1,1): ('GF',),
                    (0,2): ('ES','D','GX','GN0','GG','GF'),
                    (1,2): ('ES','D','GX','GN0','GG','GF')}

        lfw, lrf = int(isotope_flags['LFW']), int(range_flags['LRF'])
        subsection_dict = rx.DoubleSpinDict({})
        if range_flags['NRO'] > 0:
            tabhead,intdata,total_lines=self._get_tab1((0,0,0,0,'NR','NP'),
                                                       ('E','AP'),
                                                       subsection)
            subsection_dict['int']= intdata
        else:
            total_lines = 0
        if has_head_cont[(lfw, lrf)]:
            range_flags, total_lines = self._cont_and_update(
                range_flags, head_cont, subsection, total_lines)
        if (lfw, lrf) == (1,1):
            # Case B in ENDF manual p.70
            head_flags, es_array, lines = self._get_list(
                ('SPI','AP','LSSF',0,'NE','NLS'),
                ('ES',),
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
        subsection_data = (el,eh,subsection_dict,range_flags)
        isotope_dict = self.structure[mat_id]['data'][nuc_i]
        isotope_dict['unresolved'].append(subsection_data)
        return total_lines

    def _read_ap_only(self, subsection, range_flags, isotope_flags, mat_id,
                      nuc_i):
        "Read in scattering radius when it is the only resonance data given."
        subsection_dict = {}
        if range_flags['NRO'] > 0:
            tabhead,intdata,total_lines=self._get_tab1((0,0,0,0,'NR','NP'),
                                                       ('E','AP'),
                                                       subsection)
            subsection_dict['int']= intdata
        else:
            total_lines = 0
        range_flags, total_lines = self._cont_and_update(
            range_flags, ('SPI','AP',0,0,'NLS',0), subsection, total_lines)
        return total_lines

    def _read_xs(self, mat_id, mt, nuc_i=None):
        """Read in cross-section data. Read resonances with Library._read_res
        first.

        Parameters
        -----------
        mat_id: int
            ZZAAAM of material.
        mt: int
            Reaction number to find cross-section data of.
        nuc_i: int
            Isotope to find; if None, defaults to mat_id.
        """
        if nuc_i == None:
            nuc_i = mat_id
        xsdata = self.get_rx(mat_id, 3, mt).reshape(-1,6)
        total_lines = 0
        head_flags = self._get_head(('ZA','AWR',0,0,0,0),
                                    xsdata[total_lines])
        total_lines += 1
        int_flags, int_data, int_size = self._get_tab1(
            ('QM','QI',0,'LM','NR','NP'),
            ('Eint','xs'),
            xsdata[total_lines:])
        int_flags.update(head_flags)
        isotope_dict = self.structure[mat_id]['data'][nuc_i]
        isotope_dict['xs'].update({mt: (int_data, int_flags)})
        total_lines += int_size

    def get_xs(self, nuc, mt, nuc_i=None):
        """Grab cross-section data.

        Parameters
        -----------
        nuc: int
            ZZAAAM of nuclide to read.
        mt: int
            ENDF reaction number to read.
        nuc_i: int
            ZZAAAM of isotope to read. Defaults to nuc.

        Returns
        --------
        tuple
            Returns a tuple with xs data in tuple[0] and flags in tuple[1].
        """
        if not nuc_i:
            nuc_i = nuc
        if nuc not in self.structure:
            self._read_res(nuc)
        if nuc_i not in self.structure[nuc]['data'] or \
           mt not in self.structure[nuc]['data'][nuc_i]['xs']:
            self._read_xs(nuc, mt, nuc_i)
        return self.structure[nuc]['data'][nuc_i]['xs'][mt]

    def get_rx(self, nuc, mf, mt, lines=0):
        """Grab the data from one reaction type.

        Parameters
        -----------
        nuc: int
            ZZAAAM form of material to read from.
        mf: int
            ENDF file number (MF).
        mt: int
            ENDF reaction number (MT).
        lines: int
            Number of lines to read from this reaction, starting from the top.
            Default value is 0, which reads in the entire reaction.

        Returns
        --------
        data: NumPy array
            Contains the reaction data in an Nx6 array.
        """
        if nuc in self.structure:
            return self._read_nucmfmt(nuc, mf, mt, lines)
        else:
            raise ValueError("Material {} does not exist.".format(nuc))

    def _read_nucmfmt(self, nuc, mf, mt, lines):
        """Load in the data from one reaction into self.structure.

        Parameters
        -----------
        nuc : int
            ZZAAAM of nuclide.
        mf : int
            ENDF file number (MF).
        mt : int
            ENDF reaction number (MT).

        Returns
        --------
        array, 1d, float64
            1d, float64 NumPy array containing the reaction data.
        """
        opened_here = False
        if isinstance(self.fh, basestring):
            fh = open(self.fh, 'r')
            opened_here = True
        else:
            fh = self.fh
        start, stop = self.mat_dict[nuc]['mfs'][mf,mt]
        fh.readline()
        fh.seek(start)
        if lines == 0:
            s = fh.read(stop-start)
        else:
            s = fh.read(lines*81)
        if opened_here:
            fh.close
        return fromendf_tok(s)

class Evaluation(object):
    """
    Evaluation is the main class for an ENDF evaluation which contains a number
    of Files.
    """

    def __init__(self, filename, verbose=True):
        self.fh = open(filename, 'r')
        self.files = []
        self.verbose = verbose
        self.veryverbose = False

        # First we need to read MT=1, MT=451 which has a description of the ENDF
        # file and a list of what data exists in the file
        self._read_header()

    def read(self, reactions=None):
        if not reactions:
            if self.verbose:
                print 'No reaction given. Read all'
            reactions = []
            for r in self.reactionList[1:]:
                reactions.append(r[0:2])
        if isinstance(reactions, tuple):
            reactions = [reactions]
        # Start looping over the requested reactions
        # entry since it is the MT=451 block that we already read
        for rMF, rMT in reactions:
            found = False
            for MF, MT, NC, MOD in self.reactionList[1:]:
                if MF == rMF and MT == rMT:
                    found = True
                    # File 1 data
                    if MF == 1:
                        # Now read File1 (most likely unnecessary, but...)
                        file1 = self.find_file(1)
                        if not file1:
                            file1 = ENDFFile1()
                            self.files.append(file1)
                        # Number of total neutrons per fission
                        if MT == 452:
                            self._read_total_nu()
                        # Number of delayed neutrons per fission
                        elif MT == 455:
                            self._read_delayed_nu()
                        # Number of prompt neutrons per fission
                        elif MT == 456:
                            self._read_prompt_nu()
                        # Components of energy release due to fission
                        elif MT == 458:
                            self._read_fission_energy()
                        elif MT == 460:
                            self._read_delayed_photon()
                    elif MF == 2:
                        # Now read File2
                        file2 = self.find_file(2)
                        if not file2:
                            file2 = ENDFFile2()
                            self.files.append(file2)
                        if MT == 151:
                            self._read_resonances()
                    elif MF == 3:
                        file3 = self.find_file(3)
                        if not file3:
                            file3 = ENDFFile3()
                            self.files.append(file3)
                        self._read_reaction_xs(MT)
                    elif MF == 7:
                        # Now read File7
                        file7 = self.find_file(7)
                        if not file7:
                            file7 = ENDFFile7()
                            self.files.append(file7)
                        if MT == 2:
                            self._read_thermal_elastic()
                        if MT == 4:
                            self._read_thermal_inelastic()
                    elif MF == 8:
                        # Read File 8 -- decay and fission yield data
                        file8 = self.find_file(8)
                        if not file8:
                            file8 = ENDFFile8()
                            self.files.append(file8)
                        if MT == 454:
                            self._read_independent_yield()
                        elif MT == 459:
                            self._read_cumulative_yield()
                        elif MT == 457:
                            self._read_decay()
                    elif MF == 9:
                        # Read File 9 -- multiplicities
                        file9 = self.find_file(9)
                        if not file9:
                            file9 = ENDFFile9()
                            self.files.append(file9)
                        self._read_multiplicity(MT)
                    elif MF == 10:
                        # Read File 10 -- cross sections for production of
                        # radioactive nuclides
                        file10 = self.find_file(10)
                        if not file10:
                            file10 = ENDFFile10()
                            self.files.append(file10)
                        self._read_production_xs(MT)

            if not found:
                if self.verbose:
                    print 'Reaction not found'
                raise NotFound('Reaction')

    def _read_header(self):
        self.print_info(1, 451)

        # Find reaction
        self.seek_mfmt(1, 451)

        # Now read File1
        file1 = ENDFFile1()
        self.files.append(file1)

        # Create MT for description
        data = ENDFReaction(451)
        file1.reactions.append(data)

        # First HEAD record
        items = self._get_head_record()
        data.ZA = items[0]
        data.AWR = items[1]
        data.LRP = items[2]
        data.LFI = items[3]
        data.NLIB = items[4]
        data.NMOD = items[5]

        # Control record 1
        items = self._get_cont_record()
        data.ELIS = items[0]
        data.STA = int(items[1])
        data.LIS = items[2]
        data.LISO = items[3]
        data.NFOR = items[5]

        # Control record 2
        items = self._get_cont_record()
        data.AWI = items[0]
        data.EMAX = items[1]
        data.LREL = items[2]
        data.NSUB = items[4]
        data.NVER = items[5]

        # Control record 3
        items = self._get_cont_record()
        data.TEMP = items[0]
        data.LDRV = items[2]
        data.NWD = items[4]
        data.NXC = items[5]

        # Text record 1
        items = self._get_text_record()
        text = items[0]
        data.ZSYMAM = text[0:11]
        data.ALAB = text[11:22]
        data.EDATE = text[22:32]
        data.AUTH = text[32:66]

        # Text record 2
        items = self._get_text_record()
        text = items[0]
        data.REF = text[1:22]
        data.DDATE = text[22:32]
        data.RDATE = text[33:43]
        data.ENDATE = text[55:63]

        # Text record 3
        items = self._get_text_record()
        data.HSUB = items[0]

        # Now read descriptive records
        data.description = []
        for i in range(data.NWD - 3):
            line = self.fh.readline()[:66]
            data.description.append(line)

        # File numbers, reaction designations, and number of records
        self.reactionList = []
        while True:
            line = self.fh.readline()
            if line[72:75] == '  0':
                break
            items = self._get_cont_record(line, skipC=True)
            MF = items[2]
            MT = items[3]
            NC = items[4]
            MOD = items[5]
            self.reactionList.append((MF,MT,NC,MOD))

    def _read_total_nu(self):
        self.print_info(1, 452)

        # Find file 1
        file1 = self.find_file(1)

        # Find reaction
        self.seek_mfmt(1, 452)

        # Create total nu reaction
        nuTotal = ENDFReaction(452)
        file1.reactions.append(nuTotal)

        # Determine representation of total nu data
        items = self._get_head_record()
        nuTotal.LNU = items[3]

        # Polynomial representation
        if nuTotal.LNU == 1:
            nuTotal.coeffs = self._get_list_record()
        # Tabulated representation
        elif nuTotal.LNU == 2:
            nuTotal.value = self._get_tab1_record()

        # Skip SEND record
        self.fh.readline()

    def _read_delayed_nu(self):
        self.print_info(1, 455)

        # Find file 1
        file1 = self.find_file(1)

        # Find reaction
        self.seek_mfmt(1, 455)

        # Create delayed nu reaction
        nuDelay = ENDFReaction(455)
        file1.reactions.append(nuDelay)

        # Determine representation of delayed nu data
        items = self._get_head_record()
        nuDelay.LDG = items[2]
        nuDelay.LNU = items[3]

        # Nu tabulated and delayed-group constants are energy-independent
        if nuDelay.LNU == 2 and nuDelay.LDG == 0:
            nuDelay.decayConst = self._get_list_record(onlyList=True)
            nuDelay.value = self._get_tab1_record()
            self.fh.readline()
        elif nuDelay.LNU == 2 and nuDelay.LDG == 1:
            raise NotImplementedError
        elif nuDelay.LNU == 1 and nuDelay.LDG == 0:
            raise NotImplementedError
        elif nuDelay.LNU == 1 and nuDelay.LDG == 1:
            raise NotImplementedError

    def _read_prompt_nu(self):
        self.print_info(1, 456)

        # Create delayed nu reaction
        nuPrompt = ENDFReaction(456)
        self.find_file(1).reactions.append(nuPrompt)

        # Find reaction
        self.seek_mfmt(1, 456)

        # Determine representation of delayed nu data
        items = self._get_head_record()
        nuPrompt.LNU = items[3]

        # Tabulated values of nu
        if nuPrompt.LNU == 2:
            nuPrompt.value = self._get_tab1_record()
        # Spontaneous fission
        elif nuPrompt.LNU == 1:
            nuPrompt.value = self._get_list_record(onlyList=True)

        # Skip SEND record
        self.fh.readline()

    def _read_fission_energy(self):
        self.print_info(1, 458)

        # Create fission energy release reaction
        eRelease = ENDFReaction(458)
        self.find_file(1).reactions.append(eRelease)

        # Find reaction
        self.seek_mfmt(1, 458)

        # Skip HEAD record
        self._get_head_record()

        # Read LIST record containing components of fission energy release (or
        # coefficients)
        items, values = self._get_list_record()
        eRelease.NPLY = items[3]
        if eRelease.NPLY == 0:
            eRelease.fissProducts = (values[0], values[1])
            eRelease.promptNeuts = (values[2], values[3])
            eRelease.delayNeuts = (values[4], values[5])
            eRelease.promptGammas = (values[6], values[7])
            eRelease.delayGammas = (values[8], values[9])
            eRelease.delayBetas = (values[10], values[11])
            eRelease.neutrinos = (values[12], values[13])
            eRelease.pseudoQ = (values[14], values[15])
            eRelease.total = (values[16], values[17])
        elif eRelease.NPLY > 0:
            eRelease.fissProducts = zip(values[0::18], values[1::18])
            eRelease.promptNeuts = zip(values[2::18], values[3::18])
            eRelease.delayNeuts = zip(values[4::18], values[5::18])
            eRelease.promptGammas = zip(values[6::18], values[7::18])
            eRelease.delayGammas = zip(values[8::18], values[9::18])
            eRelease.delayBetas = zip(values[10::18], values[11::18])
            eRelease.neutrinos = zip(values[12::18], values[13::18])
            eRelease.pseudoQ = zip(values[14::18], values[15::18])
            eRelease.total = zip(values[16::18], values[16::18])

        # Skip SEND record
        self.fh.readline()

    def _read_reaction_xs(self, MT):
        self.print_info(3, MT)

        # Find file3
        file3 = self.find_file(3)

        # Create MT for reaction cross section
        xs = ENDFReaction(MT)
        file3.reactions.append(xs)

        # Find reaction
        self.seek_mfmt(3, MT)

        # Read HEAD record with ZA and atomic weight ratio
        items = self._get_head_record()
        xs.ZA = items[0]
        xs.AWR = items[1]

        # Read TAB1 record with reaction cross section
        xs.sigma = self._get_tab1_record()
        xs.QM = xs.sigma.params[0] # Mass difference Q value
        xs.QI = xs.sigma.params[1] # Reaction Q value
        xs.LR = xs.sigma.params[3] # Complex breakup flag

        # Skip SEND record
        self.fh.readline()

    def _read_delayed_photon(self):
        self.print_info(1, 460)

        # Create delayed photon data reaction
        dp = ENDFReaction(460)
        self.find_file(1).reactions.append(dp)

        # Find reaction
        self.seek_mfmt(1, 460)

        # Determine whether discrete or continuous representation
        items = self._get_head_record()
        dp.LO = items[2]
        dp.NG = items[4]

        # Discrete representation
        if dp.LO == 1:
            # Initialize lists for energies of photons and time dependence of
            # photon multiplicity
            dp.energy = []
            dp.multiplicity = []
            for i in range(dp.NG):
                # Read TAB1 record with multiplicity as function of time
                mult = self._get_tab1_record()
                dp.multiplicity.append(mult)

                # Determine energy
                E = mult.params[0]
                dp.energy.append(E)

        # Continuous representation
        elif dp.LO == 2:
            # Determine decay constant and number of precursor families
            dp.decayConst = self._get_list_record(onlyList=True)
            dp.NNF = len(dp.decayConst)

    def _read_resonances(self):
        self.print_info(2, 151)

        # Find reaction
        self.seek_mfmt(2, 151)

        # Now read File2
        file2 = self.find_file(2)

        # Create MT for resonances
        res = ENDFReaction(151)
        file2.reactions.append(res)
        res.resonances = []

        # Determine whether discrete or continuous representation
        items = self._get_head_record()
        res.NIS = items[4] # Number of isotopes

        for iso in range(res.NIS):
            items = self._get_cont_record()
            res.ABN = items[1] # isotopic abundance
            res.LFW = items[3] # fission widths present?
            res.NER = items[4] # number of resonance energy ranges

            for erange in range(res.NER):
                items = self._get_cont_record()
                res.EL = items[0] # lower limit of energy range
                res.EH = items[1] # upper limit of energy range
                res.LRU = items[2] # flag for resolved (1)/unresolved (2)
                res.LRF = items[3] # resonance representation
                res.NRO = items[4] # flag for energy dependence of scattering radius
                res.NAPS = items[5] # flag controlling use of channel/scattering radius

                # Only scattering radius specified
                if res.LRU == 0 and res.NRO == 0:
                    items = self._get_cont_record()
                    res.SPI = items[0]
                    res.AP = items[1]
                    res.NLS = items[4]
                # resolved resonance region
                elif res.LRU == 1:
                    self._read_resolved(res)
                # unresolved resonance region
                elif res.LRU == 2:
                    self._read_unresolved(res)

    def _read_resolved(self, res):
        # Single- or Multi-level Breit Wigner
        if res.LRF == 1 or res.LRF == 2:
            # Read energy-dependent scattering radius if present
            if res.NRO > 0:
                res.AP = self._get_tab1_record()

            # Other scatter radius parameters
            items = self._get_cont_record()
            res.SPI = items[0] # Spin, I, of the target nucleus
            if res.NRO == 0:
                res.AP = items[1]
            res.NLS = items[4] # Number of l-values

            # Read resonance widths, J values, etc
            for l in range(res.NLS):
                headerItems, items = self._get_list_record()
                QX, L, LRX = headerItems[1:4]
                energy = items[0::6]
                spin = items[1::6]
                GT = items[2::6]
                GN = items[3::6]
                GG = items[4::6]
                GF = items[5::6]
                for i, E in enumerate(energy):
                    resonance = BreitWigner()
                    resonance.QX = QX
                    resonance.L = L
                    resonance.LRX = LRX
                    resonance.E = energy[i]
                    resonance.J = spin[i]
                    resonance.GT = GT[i]
                    resonance.GN = GN[i]
                    resonance.GG = GG[i]
                    resonance.GF = GF[i]
                    res.resonances.append(resonance)

        # Reich-Moore
        elif res.LRF == 3:
            # Read energy-dependent scattering radius if present
            if res.NRO > 0:
                res.AP = self._get_tab1_record()

            # Other scatter radius parameters
            items = self._get_cont_record()
            res.SPI = items[0] # Spin, I, of the target nucleus
            if res.NRO == 0:
                res.AP = items[1]
            res.LAD = items[3] # Flag for angular distribution
            res.NLS = items[4] # Number of l-values
            res.NLSC = items[5] # Number of l-values for convergence

            # Read resonance widths, J values, etc
            for l in range(res.NLS):
                headerItems, items = self._get_list_record()
                APL, L = headerItems[1:3]
                energy = items[0::6]
                spin = items[1::6]
                GN = items[2::6]
                GG = items[3::6]
                GFA = items[4::6]
                GFB = items[5::6]
                for i, E in enumerate(energy):
                    resonance = ReichMoore()
                    resonance.APL = APL
                    resonance.L = L
                    resonance.E = energy[i]
                    resonance.J = spin[i]
                    resonance.GN = GN[i]
                    resonance.GG = GG[i]
                    resonance.GFA = GFA[i]
                    resonance.GFB = GFB[i]
                    res.resonances.append(resonance)

    def _read_unresolved(self, res):
        pass

    def _read_thermal_elastic(self):
        # Find file7
        file7 = self.find_file(7)

        # Create MT for resonances
        elast = ENDFReaction(2)
        self.print_info(7, 2)

        # Seek File 7
        self.seek_mfmt(7, 2)

        # Get head record
        items = self._get_head_record()
        elast.ZA = items[0] # ZA identifier
        elast.AWR = items[1] # AWR
        elast.LTHR = items[2] # coherent/incoherent flag
        if elast.LTHR == 1:
            if self.verbose:
                print 'Coherent elastic'
                temp = []
                eint = []
                set = []
                temp0 = self._get_tab1_record()
                # Save temp
                temp.append(temp0.params[0])
                # Save Eint
                eint = temp0.x
                # Save S(E, T0)
                set.append(temp0.y)
                elast.LT = temp0.params[2]
                if self.veryverbose:
                    print 'Number of temperatures:', elast.LT+1
                for t in range(elast.LT):
                    heads, s = self._get_list_record()
                    # Save S(E,T)
                    set.append(s)
                    # Save T
                    temp.append(heads[0])
                elast.temp = np.array(temp)
                elast.set = np.array(set)
                elast.eint = np.array(eint)
        elif elast.LTHR == 2:
            if self.verbose:
                print 'Incoherent elastic'
                temp = []
                eint = []
                set = []
                record = self._get_tab1_record()
                # Save cross section
                elast.sb = record.params[0]
                # Save Tint
                elast.t = np.array(record.x)
                # Save W(T)
                elast.w = np.array(record.y)
        else:
            print 'Invalid value of LHTR'
        file7.reactions.append(elast)

    def _read_thermal_inelastic(self):
        # Find file7
        file7 = self.find_file(7)

        # Create MT for resonances
        inel = ENDFReaction(4)
        self.print_info(7, 4)

        # Seek File 7
        self.seek_mfmt(7, 4)

        # Get head record
        items = self._get_head_record()
        inel.ZA = items[0] # ZA identifier
        inel.AWR = items[1] # AWR
        inel.LAT = items[3] # Temperature flag
        inel.LASYM = items[4] # Symmetry flag
        HeaderItems, B = self._get_list_record()
        inel.LLN = HeaderItems[2]
        inel.NS = HeaderItems[5]
        inel.B = B
        if B[0] == 0.0:
            if self.verbose:
                print 'No principal atom'
        else:
            nbeta = self._get_tab2_record()
            sabt = []
            beta = []
            for be in range(nbeta.NBT[0]):
                #Read record for first temperature (always present)
                sabt_temp = []
                temp = []
                temp0 = self._get_tab1_record()
                # Save S(be, 0, :)
                sabt_temp.append(temp0.y)
                # Save alpha(:)
                alpha = temp0.x
                # Save beta(be)
                beta.append(temp0.params[1])
                # Save temperature
                temp.append(temp0.params[0])
                inel.LT = temp0.params[2]
                if self.veryverbose:
                    print 'Number of temperatures:', inel.LT+1
                for t in range(inel.LT):
                    #Read records for all the other temperatures
                    headsab, sa = self._get_list_record()
                    # Save S(be, t+1, :)
                    sabt_temp.append(sa)
                    # Save temperature
                    temp.append(headsab[0])
                sabt.append(sabt_temp)
            # Prepare arrays for output
            sabt = np.array(sabt)
            inel.sabt = []
            for i in range(np.shape(sabt)[1]):
                inel.sabt.append(np.transpose(sabt[:, i, :]))
            inel.sabt = np.array(inel.sabt)
            inel.alpha = np.array(alpha)
            inel.beta = np.array(beta)
            inel.temp = np.array(temp)
        teffrecord = self._get_tab1_record()
        inel.teff = teffrecord.y
        file7.reactions.append(inel)

    def _read_independent_yield(self):
        self.print_info(8, 454)

        # Find file8
        file8 = self.find_file(8)

        # Create MT for resonances
        iyield = ENDFReaction(454)
        file8.reactions.append(iyield)

        # Seek File 8
        self.seek_mfmt(8, 454)

        # Initialize energies and yield dictionary
        iyield.energies = []
        iyield.data = {}
        iyield.interp = []

        items = self._get_head_record()
        iyield.ZA = items[0]
        iyield.AWR = items[1]
        LE = items[2] - 1 # Determine energy-dependence

        for i in range(LE):
            items, itemList = self._get_list_record()
            E = items[0] # Incident particle energy
            iyield.energies.append(E)
            NFP = items[5] # Number of fission product nuclide states
            if i > 0:
                iyield.interp.append(items[2]) # Interpolation scheme

            # Get data for each yield
            iyield.data[E] = {}
            iyield.data[E]['zafp'] = [int(i) for i in itemList[0::4]] # ZA for fission products
            iyield.data[E]['fps'] = itemList[1::4] # State designator
            iyield.data[E]['yi'] = zip(itemList[2::4],itemList[3::4]) # Independent yield

        # Skip SEND record
        self.fh.readline()

    def _read_cumulative_yield(self):
        self.print_info(8, 459)

        # Find file8
        file8 = self.find_file(8)

        # Create MT for resonances
        cyield = ENDFReaction(459)
        file8.reactions.append(cyield)

        # Seek File 8
        self.seek_mfmt(8, 459)

        # Initialize energies and yield dictionary
        cyield.energies = []
        cyield.data = {}
        cyield.interp = []

        items = self._get_head_record()
        cyield.ZA = items[0]
        cyield.AWR = items[1]
        LE = items[2] - 1 # Determine energy-dependence

        for i in range(LE):
            items, itemList = self._get_list_record()
            E = items[0] # Incident particle energy
            cyield.energies.append(E)
            NFP = items[5] # Number of fission product nuclide states
            if i > 0:
                cyield.interp.append(items[2]) # Interpolation scheme

            # Get data for each yield
            cyield.data[E] = {}
            cyield.data[E]['zafp'] = [int(i) for i in itemList[0::4]] # ZA for fission products
            cyield.data[E]['fps'] = itemList[1::4] # State designator
            cyield.data[E]['yc'] = zip(itemList[2::4],itemList[3::4]) # Cumulative yield

        # Skip SEND record
        self.fh.readline()

    def _read_decay(self):
        decay_type = []
        self.print_info(8, 457)

        # Find file8
        file8 = self.find_file(8)

        # Create MT for resonances
        decay = ENDFReaction(457)
        file8.reactions.append(decay)

        # Find reaction
        self.seek_mfmt(8, 457)

        # Get head record
        items = self._get_head_record()
        decay.ZA = items[0] # ZA identifier
        decay.AWR = items[1] # AWR
        decay.LIS = items[2] # State of the original nuclide
        decay.LISO = items[3] # Isomeric state for the original nuclide
        decay.NST = items[4] # Nucleus stability flag

        # Determine if radioactive (0)/stable (1)
        if decay.NST == 0:
            decay.NSP = items[5] # Number of radiation types

            # Half-life and decay energies
            items, itemList = self._get_list_record()
            decay.half_life = (items[0], items[1])
            decay.NC = items[4]/2
            decay.energies = zip(itemList[0::2], itemList[1::2])

            # Decay mode information
            items, itemList = self._get_list_record()
            decay.SPI = items[0] # Spin of the nuclide
            decay.PAR = items[1] # Parity of the nuclide
            decay.NDK = items[5] # Number of decay modes
            # Decay type (beta, gamma, etc.)
            decay.RTYP = []
            for i in itemList[0::6]:
                if i % 1.0 == 0:
                    decay.RTYP.append(decay_type[int(i)])
                else:
                    # TODO: Handle multiple decay
                    raise NotImplementedError
            decay.RFS = itemList[1::6] # Isomeric state for daughter
            decay.Q = zip(itemList[2::6], itemList[3::6]) # Total decay energy
            decay.BR = zip(itemList[4::6], itemList[5::6]) # Branching ratios

            # Read spectra
            for i in range(decay.NSP):
                items, itemList = self._get_list_record()
                STYP = decay_type[items[1]] # Decay radiation type
                LCON = items[2] # Continuous spectrum flag
                NER = items[5] # Number of tabulated discrete energies

                if LCON != 1:
                    for j in range(NER):
                        items, itemList = self._get_list_record()

                if LCON != 0:
                    r = self._get_tab1_record()
                    LCOV = r.params[3]

                    if LCOV != 0:
                        items, itemList = self._get_list_record()

        elif decay.NST == 1:
            items, itemList = self._get_list_record()
            items, itemList = self._get_list_record()
            decay.SPI = items[0]
            decay.PAR = items[1]

        # Skip SEND record
        self.fh.readline()

    def _read_multiplicity(self, MT):
        self.print_info(9, MT)

        # Find file9
        file9 = self.find_file(9)

        # Create MT for resonances
        mp = ENDFReaction(MT)
        file9.reactions.append(mp)

        # Find reaction
        self.seek_mfmt(9, MT)

        # Get head record
        items = self._get_head_record()
        mp.ZA = items[0]
        mp.AWR = items[1] # Atomic weight ratio
        mp.LIS = items[2] # Level number of the target
        mp.NS = items[4] # Number of final states

        mp.multiplicities = []
        for i in range(mp.NS):
            state = self._get_tab1_record()
            state.QM = state.params[0] # Mass difference Q value (eV)
            state.QI = state.params[1] # Reaction Q value (eV)
            state.IZAP = state.params[2] # 1000Z + A
            state.LFS = state.params[3] # Level number of the nuclide
            state.NR = state.params[4] # Number of energy ranges
            state.NP = state.params[5] # Number of energy points
            mp.multiplicities.append(state)

    def _read_production_xs(self, MT):
        self.print_info(10, MT)

        # Find file10
        file10 = self.find_file(10)

        # Create MT for resonances
        rxn = ENDFReaction(MT)
        file10.reactions.append(rxn)

        # Find reaction
        self.seek_mfmt(10, MT)

        # Get head record
        items = self._get_head_record()
        rxn.ZA = items[0]
        rxn.AWR = items[1] # Atomic weight ratio
        rxn.LIS = items[2] # Level number of the target
        rxn.NS = items[4] # Number of final states

        rxn.xs = []
        for i in range(rxn.NS):
            state = self._get_tab1_record()
            state.QM = state.params[0] # Mass difference Q value (eV)
            state.QI = state.params[1] # Reaction Q value (eV)
            state.IZAP = state.params[2] # 1000Z + A
            state.LFS = state.params[3] # Level number of the nuclide
            state.NR = state.params[4] # Number of energy ranges
            state.NP = state.params[5] # Number of energy points
            rxn.xs.append(state)

    def _get_text_record(self, line=None):
        if not line:
            line = self.fh.readline()
        if self.veryverbose:
            print 'Get TEXT record'
        HL = line[0:66]
        MAT = int(line[66:70])
        MF = int(line[70:72])
        MT = int(line[72:75])
        NS = int(line[75:80])
        return [HL, MAT, MF, MT, NS]

    def _get_cont_record(self, line=None, skipC=False):
        if self.veryverbose:
            print 'Get CONT record'
        if not line:
            line = self.fh.readline()
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
        MAT = int(line[66:70])
        MF = int(line[70:72])
        MT = int(line[72:75])
        NS = int(line[75:80])
        return [C1, C2, L1, L2, N1, N2, MAT, MF, MT, NS]

    def _get_head_record(self, line=None):
        if not line:
            line = self.fh.readline()
        if self.veryverbose:
            print 'Get HEAD record'
        ZA = int(endftod(line[:11]))
        AWR = endftod(line[11:22])
        L1 = int(line[22:33])
        L2 = int(line[33:44])
        N1 = int(line[44:55])
        N2 = int(line[55:66])
        MAT = int(line[66:70])
        MF = int(line[70:72])
        MT = int(line[72:75])
        NS = int(line[75:80])
        return [ZA, AWR, L1, L2, N1, N2, MAT, MF, MT, NS]

    def _get_list_record(self, onlyList=False):
        # determine how many items are in list
        if self.veryverbose:
            print 'Get LIST record'
        items = self._get_cont_record()
        NPL = items[4]

        # read items
        itemsList = []
        m = 0
        for i in range((NPL-1)/6 + 1):
            line = self.fh.readline()
            toRead = min(6,NPL-m)
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
        if self.veryverbose:
            print 'Get TAB1 record'
        r = ENDFTab1Record()
        r.read(self.fh)
        return r

    def _get_tab2_record(self):
        if self.veryverbose:
            print 'Get TAB2 record'
        r = ENDFTab2Record()
        r.read(self.fh)
        return r

    def find_file(self, fileNumber):
        for f in self.files:
            if f.fileNumber == fileNumber:
                return f
        return None

    def find_mt(self, MT):
        for f in self.files:
            for r in f.reactions:
                if r.MT == MT:
                    return r

#    def look_for_files(self):
#        files = set()
#        self.fh.seek(0)
#        for line in self.fh:
#            try:
#                fileNum = int(line[70:72])
#                if fileNum == 0:
#                    continue
#            except ValueError:
#                raise
#            except IndexError:
#                raise
#            files.add(fileNum)
#        print(files)

#    def seek_file(self, fileNum):
#        self.fh.seek(0)
#        fileString = '{0:2}'.format(fileNum)
#        while True:
#            position = self.fh.tell()
#            line = self.fh.readline()
#            if line == '':
#                # Reached EOF
#                print('Could not find File {0}'.format(fileNum))
#            if line[70:72] == fileString:
#                self.fh.seek(position)
#                break

    def seek_mfmt(self, MF, MT):
        self.fh.seek(0)
        searchString = '{0:2}{1:3}'.format(MF, MT)
        while True:
            position = self.fh.tell()
            line = self.fh.readline()
            if line == '':
                # Reached EOF
                if self.verbose:
                    print('Could not find MF={0}, MT={1}'.format(MF, MT))
                raise NotFound('Reaction')
            if line[70:75] == searchString:
                self.fh.seek(position)
                break

    def print_info(self, MF, MT):
        if self.verbose:
            print("Reading MF={0}, MT={1} {2}".format(MF, MT, label(MT)))

    def __iter__(self):
        for f in self.files:
            yield f

    def __repr__(self):
        try:
            name = libraries[self.files[0].NLIB]
            nuclide = self.files[0].ZA
        except:
            name = "Undetermined"
            nuclide = "None"
        return "<{0} Evaluation: {1}>".format(name, nuclide)

class ENDFTab1Record(object):
    def __init__(self):
        self.NBT = []
        self.INT = []
        self.x = []
        self.y = []

    def read(self, fh):
        # Determine how many interpolation regions and total points there are
        line = fh.readline()
        C1 = endftod(line[:11])
        C2 = endftod(line[11:22])
        L1 = int(line[22:33])
        L2 = int(line[33:44])
        NR = int(line[44:55])
        NP = int(line[55:66])
        self.params = [C1, C2, L1, L2, NR, NP]

        # Read the interpolation region data, namely NBT and INT
        m = 0
        for i in range((NR-1)/3 + 1):
            line = fh.readline()
            toRead = min(3,NR-m)
            for j in range(toRead):
                NBT = int(line[0:11])
                INT = int(line[11:22])
                self.NBT.append(NBT)
                self.INT.append(INT)
                line = line[22:]
            m = m + toRead

        # Read tabulated pairs x(n) and y(n)
        m = 0
        for i in range((NP-1)/3 + 1):
            line = fh.readline()
            toRead = min(3,NP-m)
            for j in range(toRead):
                x = endftod(line[:11])
                y = endftod(line[11:22])
                self.x.append(x)
                self.y.append(y)
                line = line[22:]
            m = m + toRead

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
        for i in range((NR-1)/3 + 1):
            line = fh.readline()
            toRead = min(3,NR-m)
            for j in range(toRead):
                NBT = int(line[0:11])
                INT = int(line[11:22])
                self.NBT.append(NBT)
                self.INT.append(INT)
                line = line[22:]
            m = m + toRead

class ENDFRecord(object):
    def __init__(self, fh):
        if fh:
            line = fh.readline()
            self.read(line)

class ENDFTextRecord(ENDFRecord):
    """
    An ENDFTextRecord is used either as the first entry on an ENDF tape (TPID)
    or to give comments in File 1.
    """

    def __init__(self, fh):
        super(ENDFTextRecord, self).__init__(fh)

    def read(self, line):
        HL = line[0:66]
        MAT = int(line[66:70])
        MF = int(line[70:72])
        MT = int(line[72:75])
        NS = int(line[75:80])
        self.items = [HL, MAT, MF, MT, NS]

class ENDFContRecord(ENDFRecord):
    """
    An ENDFContRecord is a control record.
    """

    def __init__(self, fh):
        super(ENDFContRecord, self).__init__(fh)

    def read(self, line):
        C1 = endftod(line[:11])
        C2 = endftod(line[11:22])
        L1 = int(line[22:33])
        L2 = int(line[33:44])
        N1 = int(line[44:55])
        N2 = int(line[55:66])
        MAT = int(line[66:70])
        MF = int(line[70:72])
        MT = int(line[72:75])
        NS = int(line[75:80])
        self.items = [C1, C2, L1, L2, N1, N2, MAT, MF, MT, NS]

#class ENDFEndRecord(ENDFRecord):
#    """
#    An ENDFEndRecord is an END record.
#    """
#
#    def __init__(self, fh):
#        super(ENDFEndRecord, self).__init__(fh)
#
#    def read(self, line):
#        MF = int(line[70:72])
#        MT = int(line[72:75])
#        NS = int(line[75:80])
#        self.items = [MF, MT, NS]
#
#class ENDFSendRecord(ENDFRecord):
#    """
#    An ENDFSendRecord is a SEND record.
#    """
#
#    def __init__(self, fh):
#        super(ENDFEndRecord, self).__init__(fh)
#
#    def read(self, line):
#        super(ENDFSendRecord, self).read(self.line)
#        if items[2] == 99999:
#            print 'SEND'
#        else:
#            raise NotFound('SEND')
#
#class ENDFFendRecord(ENDFRecord):
#    """
#    An ENDFFendRecord is a MEND record.
#    """
#
#    def __init__(self, fh):
#        super(ENDFEndRecord, self).__init__(fh)
#
#    def read(self, line):
#        super(ENDFFendRecord, self).read(self.line)
#        if (items[1] == 0) and (items[2] == 0):
#            print 'FEND'
#        else:
#            raise NotFound('FEND')
#
#class ENDFMendRecord(ENDFRecord):
#    """
#    An ENDFMendRecord is a MEND record.
#    """
#
#    def __init__(self, fh):
#        super(ENDFMEndRecord, self).__init__(fh)
#
#    def read(self, line):
#        super(ENDFMendRecord, self).read(self.line)
#        if (items[0] == 0) and (items[1] == 0) and (items[2] == 0):
#            print 'MEND'
#        else:
#            raise NotFound('MEND')
#
#class ENDFTendRecord(ENDFRecord):
#    """
#    An ENDFMendRecord is a MEND record.
#    """
#
#    def __init__(self, fh):
#        super(ENDFTEndRecord, self).__init__(fh)
#
#    def read(self, line):
#        super(ENDFTendRecord, self).read(self.line)
#        if (items[0] == -1) and (items[1] == 0) and (items[2] == 0):
#            print 'TEND'
#        else:
#            raise NotFound('TEND')
#

class ENDFHeadRecord(ENDFRecord):
    """
    An ENDFHeadRecord is the first in a section and has the same form as a
    control record, except that C1 and C2 fields always contain ZA and AWR,
    respectively.
    """

    def __init__(self, fh):
        super(ENDFHeadRecord, self).__init__(fh)

    def read(self, char * line):
        za_str = line[:11]
        ZA = int(endftod(za_str))
        awr_str = line[11:22]
        AWR = endftod(awr_str)
        L1 = int(line[22:33])
        L2 = int(line[33:44])
        N1 = int(line[44:55])
        N2 = int(line[55:66])
        MAT = int(line[66:70])
        MF = int(line[70:72])
        MT = int(line[72:75])
        NS = int(line[75:80])
        self.items = [ZA, AWR, L1, L2, N1, N2, MAT, MF, MT, NS]

class ENDFFile(object):
    """Abstract class for an ENDF file within an ENDF evaluation."""

    def __init__(self):
        self.reactions = []

    def __repr__(self):
        try:
            return "<ENDF File {0.fileNumber}: {0.ZA}>".format(self)
        except:
            return "<ENDF File {0.fileNumber}>".format(self)

class ENDFFile1(ENDFFile):
    """
    File1 contains general information about an evaluated data set such as
    descriptive data, number of neutrons per fission, delayed neutron data,
    number of prompt neutrons per fission, and components of energy release due
    to fission.
    """

    def __init__(self):
        super(ENDFFile1,self).__init__()

        self.fileNumber = 1

class ENDFFile2(ENDFFile):
    """
    File2 contains resonance parameters including resolved resonance parameters,
    unresolved resonance parameters, and prescribed procedures for resolved and
    unresolved resonances.
    """

    def __init__(self):
        super(ENDFFile2,self).__init__()

        self.fileNumber = 2

class ENDFFile3(ENDFFile):
    """
    File3 contains neutron cross sections and prescribed procedures for
    inciddent neutrons.
    """

    def __init__(self):
        super(ENDFFile3, self).__init__()
        self.fileNumber = 3

class ENDFFile4(ENDFFile):
    """
    File4 contains angular distributions of secondary particles.
    """

    def __init__(self):
        self.fileNumber = 4

class ENDFFile5(ENDFFile):
    """
    File5 contains energy distributions of secondary particles.
    """

    def __init__(self):
        self.fileNumber = 5

class ENDFFile6(ENDFFile):
    """
    File6 contains product energy-angle distributions.
    """

    def __init__(self):
        self.fileNumber = 6

class ENDFFile7(ENDFFile):
    """
    File7 contains thermal neutron scattering law data.
    """

    def __init__(self):
        super(ENDFFile7,self).__init__()
        self.fileNumber = 7

class ENDFFile8(ENDFFile):
    """
    File8 contains radioactive decay data.
    """

    def __init__(self):
        super(ENDFFile8,self).__init__()
        self.fileNumber = 8

class ENDFFile9(ENDFFile):
    """
    File9 contains multiplicites for production of radioactive elements.
    """

    def __init__(self):
        super(ENDFFile9,self).__init__()
        self.fileNumber = 9

class ENDFFile10(ENDFFile):
    """
    File10 contains cross sections for production of radioactive nuclides.
    """

    def __init__(self):
        super(ENDFFile10,self).__init__()
        self.fileNumber = 10

class ENDFReaction(ENDFFile):
    """A single MT record on an ENDF file."""

    def __init__(self, MT):
        self.MT = MT

    def __repr__(self):
        return "<ENDF Reaction: MT={0}, {1}>".format(self.MT, label(self.MT))

class Resonance(object):
    def __init__(self):
        pass

class BreitWigner(Resonance):
    def __init__(self):
        pass

    def __repr__(self):
        return "<Breit-Wigner Resonance: l={0.L} J={0.J} E={0.E}>".format(self)

class ReichMoore(Resonance):
    def __init__(self):
        pass

    def __repr__(self):
        return "<Reich-Moore Resonance: l={0.L} J={0.J} E={0.E}>".format(self)

class AdlerAdler(Resonance):
    def __init__(self):
        pass

class RMatrixLimited(Resonance):
    def __init__(self):
        pass


MTname = {1: "(n,total) Neutron total",
          2: "(z,z0) Elastic scattering",
          3: "(z,nonelas) Nonelastic neutron",
          4: "(z,n) One neutron in exit channel",
          5: "(z,anything) Miscellaneous",
          10: "(z,contin) Total continuum reaction",
          11: "(z,2nd) Production of 2n and d",
          16: "(z,2n) Production of 2n",
          17: "(z,3n) Production of 3n",
          18: "(z,fiss) Particle-induced fission",
          19: "(z,f) First-chance fission",
          20: "(z,nf) Second chance fission",
          21: "(z,2nf) Third-chance fission",
          22: "(z,na) Production of n and alpha",
          23: "(z,n3a) Production of n and 3 alphas",
          24: "(z,2na) Production of 2n and alpha",
          25: "(z,3na) Production of 3n and alpha",
          27: "(n,abs) Absorption",
          28: "(z,np) Production of n and p",
          29: "(z,n2a) Production of n and 2 alphas",
          30: "(z,2n2a) Production of 2n and 2 alphas",
          32: "(z,nd) Production of n and d",
          33: "(z,nt) Production of n and t",
          34: "(z,n3He) Production of n and He-3",
          35: "(z,nd2a) Production of n, d, and alpha",
          36: "(z,nt2a) Production of n, t, and 2 alphas",
          37: "(z,4n) Production of 4n",
          38: "(z,3nf) Fourth-chance fission",
          41: "(z,2np) Production of 2n and p",
          42: "(z,3np) Production of 3n and p",
          44: "(z,n2p) Production of n and 2p",
          45: "(z,npa) Production of n, p, and alpha",
          50: "(z,n0) Production of n, ground state",
          51: "(z,n1) Production of n, 1st excited state",
          52: "(z,n2) Production of n, 2nd excited state",
          53: "(z,n3) Production of n, 3rd excited state",
          54: "(z,n4) Production of n, 4th excited state",
          55: "(z,n5) Production of n, 5th excited state",
          56: "(z,n6) Production of n, 6th excited state",
          57: "(z,n7) Production of n, 7th excited state",
          58: "(z,n8) Production of n, 8th excited state",
          59: "(z,n9) Production of n, 9th excited state",
          60: "(z,n10) Production of n, 10th excited state",
          61: "(z,n11) Production of n, 11th excited state",
          62: "(z,n12) Production of n, 12th excited state",
          63: "(z,n13) Production of n, 13th excited state",
          64: "(z,n14) Production of n, 14th excited state",
          65: "(z,n15) Production of n, 15th excited state",
          66: "(z,n16) Production of n, 16th excited state",
          67: "(z,n17) Production of n, 17th excited state",
          68: "(z,n18) Production of n, 18th excited state",
          69: "(z,n19) Production of n, 19th excited state",
          70: "(z,n20) Production of n, 20th excited state",
          71: "(z,n21) Production of n, 21st excited state",
          72: "(z,n22) Production of n, 22nd excited state",
          73: "(z,n23) Production of n, 23rd excited state",
          74: "(z,n24) Production of n, 24th excited state",
          75: "(z,n25) Production of n, 25th excited state",
          76: "(z,n26) Production of n, 26th excited state",
          77: "(z,n27) Production of n, 27th excited state",
          78: "(z,n28) Production of n, 28th excited state",
          79: "(z,n29) Production of n, 29th excited state",
          80: "(z,n30) Production of n, 30th excited state",
          81: "(z,n31) Production of n, 31st excited state",
          82: "(z,n32) Production of n, 32nd excited state",
          83: "(z,n33) Production of n, 33rd excited state",
          84: "(z,n34) Production of n, 34th excited state",
          85: "(z,n35) Production of n, 35th excited state",
          86: "(z,n36) Production of n, 36th excited state",
          87: "(z,n37) Production of n, 37th excited state",
          88: "(z,n38) Production of n, 38th excited state",
          89: "(z,n39) Production of n, 39th excited state",
          90: "(z,n40) Production of n, 40th excited state",
          91: "(z,nc) Production of n in continuum",
          101: "(n,disap) Neutron disappeareance",
          102: "(z,gamma) Radiative capture",
          103: "(z,p) Production of p",
          104: "(z,d) Production of d",
          105: "(z,t) Production of t",
          106: "(z,3He) Production of He-3",
          107: "(z,a) Production of alpha",
          108: "(z,2a) Production of 2 alphas",
          109: "(z,3a) Production of 3 alphas",
          111: "(z,2p) Production of 2p",
          112: "(z,pa) Production of p and alpha",
          113: "(z,t2a) Production of t and 2 alphas",
          114: "(z,d2a) Production of d and 2 alphas",
          115: "(z,pd) Production of p and d",
          116: "(z,pt) Production of p and t",
          117: "(z,da) Production of d and a",
          151: "Resonance Parameters",
          201: "(z,Xn) Total neutron production",
          202: "(z,Xgamma) Total gamma production",
          203: "(z,Xp) Total proton production",
          204: "(z,Xd) Total deuteron production",
          205: "(z,Xt) Total triton production",
          206: "(z,X3He) Total He-3 production",
          207: "(z,Xa) Total alpha production",
          208: "(z,Xpi+) Total pi+ meson production",
          209: "(z,Xpi0) Total pi0 meson production",
          210: "(z,Xpi-) Total pi- meson production",
          211: "(z,Xmu+) Total anti-muon production",
          212: "(z,Xmu-) Total muon production",
          213: "(z,Xk+) Total positive kaon production",
          214: "(z,Xk0long) Total long-lived neutral kaon production",
          215: "(z,Xk0short) Total short-lived neutral kaon production",
          216: "(z,Xk-) Total negative kaon production",
          217: "(z,Xp-) Total anti-proton production",
          218: "(z,Xn-) Total anti-neutron production",
          251: "Average cosine of scattering angle",
          252: "Average logarithmic energy decrement",
          253: "Average xi^2/(2*xi)",
          451: "Desciptive data",
          452: "Total neutrons per fission",
          454: "Independent fission product yield",
          455: "Delayed neutron data",
          456: "Prompt neutrons per fission",
          457: "Radioactive decay data",
          458: "Energy release due to fission",
          459: "Cumulative fission product yield",
          460: "Delayed photon data",
          500: "Total charged-particle stopping power",
          501: "Total photon interaction",
          502: "Photon coherent scattering",
          504: "Photon incoherent scattering",
          505: "Imaginary scattering factor",
          506: "Real scattering factor",
          515: "Pair production, electron field",
          516: "Total pair production",
          517: "Pair production, nuclear field",
          522: "Photoelectric absorption",
          523: "Photo-excitation cross section",
          526: "Electro-atomic scattering",
          527: "Electro-atomic bremsstrahlung",
          528: "Electro-atomic excitation cross section",
          533: "Atomic relaxation data",
          534: "K (1s1/2) subshell",
          535: "L1 (2s1/2) subshell",
          536: "L2 (2p1/2) subshell",
          537: "L3 (2p3/2) subshell",
          538: "M1 (3s1/2) subshell",
          539: "M2 (3p1/2) subshell",
          540: "M3 (3p3/2) subshell",
          541: "M4 (3d1/2) subshell",
          542: "M5 (3d1/2) subshell",
          543: "N1 (4s1/2) subshell",
          544: "N2 (4p1/2) subshell",
          545: "N3 (4p3/2) subshell",
          546: "N4 (4d3/2) subshell",
          547: "N5 (4d5/2) subshell",
          548: "N6 (4f5/2) subshell",
          549: "N7 (4f7/2) subshell",
          550: "O1 (5s1/2) subshell",
          551: "O2 (5p1/2) subshell",
          552: "O3 (5p3/2) subshell",
          553: "O4 (5d3/2) subshell",
          554: "O5 (5d5/2) subshell",
          555: "O6 (5f5/2) subshell",
          556: "O7 (5f7/2) subshell",
          557: "O8 (5g7/2) subshell",
          558: "O9 (5g9/2) subshell",
          559: "P1 (6s1/2) subshell",
          560: "P2 (6p1/2) subshell",
          561: "P3 (6p3/2) subshell",
          562: "P4 (6d3/2) subshell",
          563: "P5 (6d5/2) subshell",
          564: "P6 (6f5/2) subshell",
          565: "P7 (6f7/2) subshell",
          566: "P8 (6g7/2) subshell",
          567: "P9 (6g9/2) subshell",
          568: "P10 (6h9/2) subshell",
          569: "P11 (6h11/2) subshell",
          570: "Q1 (7s1/2) subshell",
          571: "Q2 (7p1/2) subshell",
          572: "Q3 (7p3/2) subshell"}


class NotFound(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
