#!/usr/bin/env python

"""Module for parsing and manipulating data from ENDL evaluations.

For the moment, classes and functions in this module only implement the
specifications relevant to the EEDL (Evaluated Electron Data Library) and the
EPDL (Evaluated Photon Data Library). The formats are described by the
following documents:

"ENDL type formats for the Livermore Evaluated Photon Data Library, EPDL"
https://www.oecd-nea.org/dbdata/data/manual-endf/nds_eval_epdl.pdf

"ENDL type formats for the Livermore Evaluated Electron Data Library, EEDL"
https://www.oecd-nea.org/dbdata/data/manual-endf/nds_eval_eedl.pdf

For more information, contact Davide Mancusi <davide.mancusi@cea.fr>.
"""
from __future__ import print_function, division, unicode_literals

from warnings import warn
from pyne.utils import QAWarning

import pyne.rxdata as rx

import numpy as np
import re

warn(__name__ + ' is not yet QA compliant.', QAWarning)

END_OF_TABLE_RE = re.compile(' {71}1')


class Library(rx.RxLib):
    """A class for a file which contains multiple ENDL tables."""
    def __init__(self, fh):
        self.mts = {}
        self.structure = {}
        self.mat_dict = {}
        self.more_files = True
        self.intdict = {
                0: self._linlin,
                2: self._linlin,
                3: self._loglin,
                4: self._linlog,
                5: self._loglog,
                }
        self.fh = fh
        # read headers for all tables
        while self.more_files:
            self._read_headers()

#    def load(self):
#        """load()
#        Read the ENDL file into a NumPy array.
#
#        Returns
#        -------
#        data : np.array, 1d, float64
#            Returns a 1d float64 NumPy array.
#        """
#        opened_here = False
#        if isinstance(self.fh, basestring):
#            fh = open(self.fh, 'r')
#            opened_here = True
#        else:
#            fh = self.fh
#        fh.readline()
#        data = fromendf_tok(fh.read())
#        fh.seek(0)
#        if opened_here:
#            fh.close()
#        return data

    def _read_headers(self):
        opened_here = False
        if isinstance(self.fh, basestring):
            fh = open(self.fh, 'r')
            opened_here = True
        else:
            fh = self.fh

        while fh:
            # store the start of the table
            start = fh.tell()

            # get header lines
            line1 = fh.readline()
            line2 = fh.readline()

            # parse the first header line
            nuc_zzzaaa = line[0:6].strip()
            yi = int(line[7:9].strip() or -1)
            yo = int(line[10:12].strip() or -1)
            aw = pyne.utils.endftod(line[13:24].strip() or -1.)
            date = line[25:31].strip()
            iflag = int(line[31] or -1)

            # parse the second header line
            rdesc = int(line[0:2].strip() or -1)
            rprop = int(line[2:5].strip() or -1)
            rmod = int(line[5:8].strip() or -1)
            x1 = pyne.utils.endftod(line[21:32].strip() or -1.)

            # convert to Pyne native formats
            nuc = cpp_nucname.zzzaaa_to_id(nuc_zzzaaa)

            # make a new dict in self.structure to contain the table data.
            if nuc not in self.structure:
                self.structure.update(
                        { nuc: {
                            'particle_in': [],
                            'particle_out': [],
                            'data': {},
                            'start': start,
                            }})

            # skip to the end of the table
            line = ''
            while fh and not END_OF_TABLE_RE.match(line):
                line = fh.readline()



class NotFound(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
