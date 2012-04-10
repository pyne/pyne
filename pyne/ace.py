#!/usr/bin/env python

"""This module is for reading ACE-format cross sections. ACE stands for "A Compact
ENDF" format and originated from work on MCNP_. It is used in a number of other
Monte Carlo particle transport codes.

ACE-format cross sections are typically generated from ENDF_ files through a
cross section processing program like NJOY_. The ENDF data consists of tabulated
thermal data, ENDF/B resonance parameters, distribution parameters in the
unresolved resonance region, and tabulated data in the fast region. After the
ENDF data has been reconstructed and Doppler-broadened, the ACER module
generates ACE-format cross sections.

.. _MCNP: http://mcnp-green.lanl.gov/

.. _NJOY: http://t2.lanl.gov/codes/codes.html

.. _ENDF: http://www.nndc.bnl.gov/endf

.. moduleauthor:: Paul Romano <romano7@gmail.com>
"""

import struct

import numpy as np
from numpy import zeros, copy, meshgrid, interp, linspace, pi, arccos, concatenate
from bisect import bisect_right

from pyne import nucname

class Library(object):
    """A Library objects represents an ACE-formatted file which may contain
    multiple tables with data.

    Parameters
    ----------
    filename : str
        Path of the ACE library file to load.

    """

    def __init__(self, filename):
        # Determine whether file is ASCII or binary
        try:
            self.f = open(filename, 'r')
            # Grab 10 lines of the library
            s = ''.join([self.f.readline() for i in range(10)])

            # Try to decode it with ascii
            sd = s.decode('ascii')

            # No exception so proceed with ASCII
            self.f.seek(0)
            self.binary = False
        except UnicodeDecodeError:
            self.f.close()
            self.f = open(filename, 'rb')
            self.binary = True

        # Set verbosity
        self.verbose = False
        self.tables = {}

    def read(self, table_names=None):
        """Read through and parse the ACE-format library.

        Parameters
        ----------
        table_names : None, str, or iterable, optional
            Tables from the file to read in.  If None, reads in all of the 
            tables. If str, reads in only the single table of a matching name.
        """
        if isinstance(table_names, basestring):
            table_names = [table_names]

        if table_names is not None:
            table_names = set(table_names)

        if self.binary:
            self._read_binary(table_names)
        else:
            self._read_ascii(table_names)

    def _read_binary(self, table_names, recl_length=4096, entries=512):

        while True:
            start_position = self.f.tell()

            # Check for end-of-file
            if self.f.read(1) == '':
                return
            self.f.seek(start_position)

            # Read name, atomic weight ratio, temperature, date, comment, and
            # material
            name, awr, temp, date, comment, mat = \
                struct.unpack('=10sdd10s70s10s', self.f.read(116))
            name = name.strip()

            # Read ZAID/awr combinations
            data = struct.unpack('=' + 16*'id', self.f.read(192))

            # Read NXS
            NXS = list(struct.unpack('=16i', self.f.read(64)))

            # Determine length of XSS and number of records
            length = NXS[0]
            n_records = (length + entries - 1)/entries

            # verify that we are suppossed to read this table in
            if (table_names is not None) and (name not in table_names):
                self.f.seek(start_position + recl_length*(n_records + 1))
                continue

            # ensure we have a valid table type
            if 0 == len(name) or name[-1] not in table_types:
                # TODO: Make this a proper exception.
                print("Unsupported table: " + name)
                self.f.seek(start_position + recl_length*(n_records + 1))
                continue

            # get the table
            table = table_types[name[-1]](name, awr, temp)

            temp_in_K = round(temp * 1e6 / 8.617342e-5)
            if self.verbose:
                print("Loading nuclide {0} at {1} K ({2})".format(
                        nucname.serpent(name.partition('.')[0]), temp_in_K, name))
            self.tables[name] = table

            # Set NXS and read JXS
            table.NXS = NXS
            table.JXS = list(struct.unpack('=32i', self.f.read(128)))

            # Read XSS
            self.f.seek(start_position + recl_length)
            table.XSS = list(struct.unpack('={0}d'.format(length),
                                           self.f.read(length*8)))

            # Insert empty object at beginning of NXS, JXS, and XSS
            # arrays so that the indexing will be the same as
            # Fortran. This makes it easier to follow the ACE format
            # specification.
            table.NXS.insert(0, None)
            table.JXS.insert(0, None)
            table.XSS.insert(0, None)

            # Read all data blocks
            table._read_all()

            # Advance to next record
            self.f.seek(start_position + recl_length*(n_records + 1))

    def _read_ascii(self, table_names):

        lines = self.f.readlines()
        
        while True:
            # Check end condition
            if 0 == len(lines):
                break

            # Read name of table, atomic weight ratio, and temperature. If first
            # line is empty, we are at end of file
            words = lines[0].split()
            name = words[0]
            awr = float(words[1])
            temp = float(words[2])

            nxs = [int(i) for i in ''.join(lines[6:8]).split()]
            n_lines = (nxs[0] + 3)/4

            # verify that we are suppossed to read this table in
            if (table_names is not None) and (name not in table_names):
                lines = lines[12+n_lines:]
                continue

            # ensure we have a valid table type
            if 0 == len(name) or name[-1] not in table_types:
                # TODO: Make this a proper exception.
                print("Unsupported table: " + name)
                lines = lines[12+n_lines:]
                continue

            # get the table
            table = table_types[name[-1]](name, awr, temp)

            temp_in_K = round(temp * 1e6 / 8.617342e-5)
            if self.verbose:
                print("Loading nuclide {0} at {1} K ({2})".format(
                        nucname.serpent(name.partition('.')[0]), temp_in_K, name))
            self.tables[name] = table

            # Read comment
            table.comment = lines[1].strip()

            # Read NXS and JXS arrays
            table.NXS = nxs
            table.JXS = [int(i) for i in ''.join(lines[8:12]).split()]

            # Read XSS array
            #table.XSS = [float(i) for i in ''.join(lines[12:12+n_lines]).split()]
            datastr = ''.join(lines[12:12+n_lines])
            xss = np.fromstring(datastr, sep=' ')
            table.XSS = list(xss)

            # Insert empty object at beginning of NXS, JXS, and XSS
            # arrays so that the indexing will be the same as
            # Fortran. This makes it easier to follow the ACE format
            # specification.
            table.NXS.insert(0, None)
            table.JXS.insert(0, None)
            table.XSS.insert(0, None)

            # Read all data blocks
            table._read_all()
            lines = lines[12+n_lines:]

    def find_table(self, name):
        """Returns a cross-section table with a given name.

        Parameters
        ----------
        name : str
            Name of the cross-section table, e.g. 92235.70c

        """
        return self.tables.get(name, None)

    def __del__(self):
        self.f.close()


class AceTable(object):
    """Abstract superclass of all other classes for cross section tables."""

    def __init__(self, name, awr, temp):
        self.name = name
        self.awr = awr
        self.temp = temp

    def _read_all(self):
        raise NotImplementedError
        
        
class NeutronTable(AceTable):
    """A NeutronTable object contains continuous-energy neutron interaction data
    read from an ACE-formatted Type I table. These objects are not normally
    instantiated by the user but rather created when reading data using a
    Library object and stored within the ``tables`` attribute of a Library object.

    Parameters
    ----------
    name : str
        ZAID identifier of the table, e.g. 92235.70c.
    awr : float
        Atomic weight ratio of the target nuclide.
    temp : float
        Temperature of the target nuclide in eV.
    
    """

    def __init__(self, name, awr, temp):
        super(NeutronTable, self).__init__(name, awr, temp)

    def __repr__(self):
        if hasattr(self, 'name'):
            return "<ACE Continuous-E Neutron Table: {0}>".format(self.name)
        else:
            return "<ACE Continuous-E Neutron Table>"

    def _read_all(self):
        self._read_esz()
        self._read_nu()
        self._read_mtr()
        self._read_lqr()
        self._read_tyr()
        self._read_lsig()
        self._read_sig()
        self._read_land()
        self._read_and()
        self._read_ldlw()
        self._read_dlw()
        self._read_gpd()
        self._read_mtrp()
        self._read_lsigp()
        self._read_sigp()
        self._read_landp()
        self._read_andp()
        # Read LDLWP block
        # Read DLWP block
        # Read YP block
        # self._read_yp()
        self._read_fis()
        self._read_unr()

    def _read_esz(self):
        """Read ESZ block -- this block contains the energy grid, total
        xs, absorption xs, elastic scattering xs, and heating numbers.
        """
        
        NE = self.NXS[3]
        self.index = self.JXS[1]

        self.energy   = self._get_float(NE)
        self.sigma_t  = self._get_float(NE)
        self.sigma_a  = self._get_float(NE)
        sigma_el = self._get_float(NE)
        self.heating  = self._get_float(NE)

        # Create elastic scattering reaction
        MT = 2
        rxn = Reaction(MT, self)
        rxn.Q = 0.0
        rxn.IE = 1
        rxn.TY = 1
        rxn.sigma = sigma_el

        # Add elastic scattering to list of reactions
        self.reactions = []
        self.reactions.append(rxn)

        # Create photon-producting reaction list
        self.photonReactions = []

    def _read_nu(self):
        """Read the NU block -- this contains information on the prompt
        and delayed neutron precursor yields, decay constants, etc
        """

        JXS2 = self.JXS[2]

        # No NU block
        if JXS2 == 0:
            return

        # Either prompt nu or total nu is given
        if self.XSS[JXS2] > 0:
            KNU = JXS2
            LNU = int(self.XSS[KNU])

            # Polynomial function form of nu
            if LNU == 1:
                self.nu_t_type = "polynomial"
                NC = int(self.XSS[KNU+1])
                coeffs = self.XSS[KNU+2 : KNU+2+NC]
                
            # Tabular data form of nu
            elif LNU == 2:
                self.nu_t_type = "tabular"
                NR = int(self.XSS[KNU+1])
                if NR > 0:
                    interp_NBT = self.XSS[KNU+2    : KNU+2+NR  ]
                    interp_INT = self.XSS[KNU+2+NR : KNU+2+2*NR]
                NE = int(self.XSS[KNU+2+2*NR])
                self.nu_t_energy = self.XSS[KNU+3+2*NR    : KNU+3+2*NR+NE  ]
                self.nu_t_value  = self.XSS[KNU+3+2*NR+NE : KNU+3+2*NR+2*NE]

        # Both prompt nu and total nu
        elif self.XSS[JXS2] < 0:
            KNU = JXS2 + 1
            LNU = int(self.XSS[KNU])

            # Polynomial function form of nu
            if LNU == 1:
                self.nu_p_type = "polynomial"
                NC = int(self.XSS[KNU+1])
                coeffs = self.XSS[KNU+2 : KNU+2+NC]
                
            # Tabular data form of nu
            elif LNU == 2:
                self.nu_p_type = "tabular"
                NR = int(self.XSS[KNU+1])
                if NR > 0:
                    interp_NBT = self.XSS[KNU+2    : KNU+2+NR  ]
                    interp_INT = self.XSS[KNU+2+NR : KNU+2+2*NR]
                NE = int(self.XSS[KNU+2+2*NR])
                self.nu_p_energy = self.XSS[KNU+3+2*NR    : KNU+3+2*NR+NE  ]
                self.nu_p_value  = self.XSS[KNU+3+2*NR+NE : KNU+3+2*NR+2*NE]
                
            KNU = JXS2 + int(abs(self.XSS[JXS2])) + 1
            LNU = int(self.XSS[KNU])

            # Polynomial function form of nu
            if LNU == 1:
                self.nu_t_type = "polynomial"
                NC = int(self.XSS[KNU+1])
                coeffs = self.XSS[KNU+2 : KNU+2+NC]
                
            # Tabular data form of nu
            elif LNU == 2:
                self.nu_t_type = "tabular"
                NR = int(self.XSS[KNU+1])
                if NR > 0:
                    interp_NBT = self.XSS[KNU+2    : KNU+2+NR  ]
                    interp_INT = self.XSS[KNU+2+NR : KNU+2+2*NR]
                NE = int(self.XSS[KNU+2+2*NR])
                self.nu_t_energy = self.XSS[KNU+3+2*NR    : KNU+3+2*NR+NE  ]
                self.nu_t_value  = self.XSS[KNU+3+2*NR+NE : KNU+3+2*NR+2*NE]
    
        # Check for delayed nu data
        if self.JXS[24] > 0:
            KNU = self.JXS[24]
            NR = int(self.XSS[KNU+1])
            if NR > 0:
                interp_NBT = self.XSS[KNU+2    : KNU+2+NR  ]
                interp_INT = self.XSS[KNU+2+NR : KNU+2+2*NR]
            NE = int(self.XSS[KNU+2+2*NR])
            self.nu_d_energy = self.XSS[KNU+3+2*NR    : KNU+3+2*NR+NE  ]
            self.nu_d_value  = self.XSS[KNU+3+2*NR+NE : KNU+3+2*NR+2*NE]

            # Delayed neutron precursor distribution
            self.nu_d_precursor_const = {}
            self.nu_d_precursor_energy = {}
            self.nu_d_precursor_prob = {}
            i = self.JXS[25]
            n_group = self.NXS[8]
            for group in range(n_group):
                self.nu_d_precursor_const[group] = self.XSS[i]
                NR = int(self.XSS[i+1])
                if NR > 0:
                    interp_NBT = self.XSS[i+2    : i+2+NR]
                    interp_INT = self.XSS[i+2+NR : i+2+2*NR]
                NE = int(self.XSS[i+2+2*NR])
                self.nu_d_precursor_energy[group] = self.XSS[i+3+2*NR    : i+3+2*NR+NE  ]
                self.nu_d_precursor_prob[group]   = self.XSS[i+3+2*NR+NE : i+3+2*NR+2*NE]
                i = i+3+2*NR+2*NE

            # Energy distribution for delayed fission neutrons
            LED = self.JXS[26]
            LOCC = {}
            for group in range(n_group):
                LOCC[group] = self.XSS[LED + group]

    def _read_mtr(self):
        """Get the list of reaction MTs for this cross-section table. The
        MT values are somewhat arbitrary.
        """

        LMT = self.JXS[3]
        NMT = self.NXS[4]
        for i in range(NMT):
            MT = int(self.XSS[LMT+i])
            self.reactions.append(Reaction(MT,self))
            
    def _read_lqr(self):
        """Find Q-values for each reaction MT
        """

        JXS4 = self.JXS[4]
        for i, rxn in enumerate(self.reactions[1:]):
            rxn.Q = self.XSS[JXS4+i]

    def _read_tyr(self):
        """Find the neutron release for each reaction MT. A neutron
        release of 19 indicates fission. A neutron release greater
        than 100 signifies reactions other than fission taht have
        energy-dependent neutron multiplicities
        """

        JXS5 = self.JXS[5]
        for i, rxn in enumerate(self.reactions[1:]):
            rxn.TY = int(self.XSS[JXS5+i])

    def _read_lsig(self):
        """Determine location of cross sections for each reaction MT
        """

        LXS = self.JXS[6]
        for i, rxn in enumerate(self.reactions[1:]):
            rxn.LOCA = int(self.XSS[LXS+i])

    def _read_sig(self):
        """Read cross-sections for each reaction MT
        """

        JXS7 = self.JXS[7]
        for rxn in self.reactions[1:]:
            rxn.IE = int(self.XSS[JXS7+rxn.LOCA-1])
            NE = int(self.XSS[JXS7+rxn.LOCA])
            rxn.sigma = self.XSS[JXS7+rxn.LOCA+1 : JXS7+rxn.LOCA+1+NE]

    def _read_land(self):
        """Find locations for angular distributions
        """

        JXS8 = self.JXS[8]

        # Number of reactions is less than total since we only need
        # angular distribution for reactions with secondary
        # neutrons. Thus, MT > 100 are not included.
        NMT = self.NXS[5]

        # Need NMT + 1 since elastic scattering is included
        for i, rxn in enumerate(self.reactions[:NMT+1]):
            rxn.LOCB = int(self.XSS[JXS8+i])

    def _read_and(self):
        """Find the angular distribution for each reaction MT
        """

        JXS9 = self.JXS[9]
        NMT = self.NXS[5]

        # Angular distribution for all MT with secondary neutrons
        # including elastic scattering
        for rxn in self.reactions[:NMT+1]:

            # Check if angular distribution data exist 
            if rxn.LOCB == -1:
                # Angular distribution data are specified through LAWi
                # = 44 in the DLW block
                continue
            elif rxn.LOCB == 0:
                # No angular distribution data are given for this
                # reaction, isotropic scattering is asssumed (in CM if
                # TY < 0 and in LAB if TY > 0)
                continue

            self.index = JXS9 + rxn.LOCB - 1
            NE = self._get_int()
            rxn.ang_energy_in = self._get_float(NE)
            LC = self._get_int(NE)
            rxn.ang_location = LC

            rxn.ang_cos = {}
            rxn.ang_pdf = {}
            rxn.ang_cdf = {}
            for i, location in enumerate(LC):
                if location == 0:
                    # Isotropic angular distribution
                    continue
                elif location > 0:
                    # Equiprobable 32 bin distribution
                    # print([rxn,'equiprobable'])
                    rxn.ang_cos[i] = self._get_float(33)
                elif location < 0:
                    # Tabular angular distribution
                    JJ = self._get_int()
                    NP = self._get_int()
                    rxn.ang_cos[i] = self._get_float(NP)
                    rxn.ang_pdf[i] = self._get_float(NP)
                    rxn.ang_cdf[i] = self._get_float(NP)
        
    def _read_ldlw(self):
        """Find locations for energy distribution data for each reaction
        """

        LED = self.JXS[10]

        # Number of reactions is less than total since we only need
        # energy distribution for reactions with secondary
        # neutrons. Thus, MT > 100 are not included. Elastic
        # scattering is also not included.
        NMT = self.NXS[5]

        for i, rxn in enumerate(self.reactions[1:NMT+1]):
            rxn.LOCC = int(self.XSS[LED+i])

    def _read_dlw(self):
        """Determine the energy distribution for secondary neutrons for
        each reaction MT
        """

        LDIS = self.JXS[11]
        NMT = self.NXS[5]

        LOCC = [rxn.LOCC for rxn in self.reactions[1:NMT+1]]

        for irxn, rxn in enumerate(self.reactions[1:NMT+1]):
            self.index = LDIS + rxn.LOCC - 1
            LNW = self._get_int()
            LAW = self._get_int()
            # print([rxn,LAW])
            IDAT = self._get_int()
            NR = self._get_int()
            if NR > 0:
                interp_NBT = self._get_int(NR)
                interp_INT = self._get_int(NR)

            # Determine tabular energy points and probability of law
            # validity
            NE = self._get_int()
            rxn.e_dist_energy = self._get_float(NE)
            rxn.e_dist_pvalid = self._get_float(NE)

            rxn.e_dist_law = LAW
            self.index = LDIS + IDAT - 1

            # Tabular equiprobable energy bins (ENDF Law 1)
            if LAW == 1:
                NR = self._get_int()
                if NR > 0:
                    rxn.e_dist_NBT = self._get_int(NR)
                    rxn.e_dist_INT = self._get_int(NR)

                # Number of outgoing energies in each E_out table
                NE = self._get_int()
                rxn.e_dist_energy_in = self._get_float(NE)

                # Read E_out tables
                NET = self._get_int()
                self.e_dist_energy_out1 = self._get_float(NET)
                self.e_dist_energy_out2 = self._get_float(NET)
                self.e_dist_energy_outNE = self._get_float(NET)

            # Discrete photon energy
            elif LAW == 2:
                self.e_dist_LP = self._get_int()
                self.e_dist_EG = self._get_float()

            # Level scattering (ENDF Law 3)
            elif LAW == 3:
                rxn.e_dist_data = self._get_float(2)

            # Continuous tabular distribution (ENDF Law 1)
            elif LAW == 4:
                NR = self._get_int()
                if NR > 0:
                    rxn.e_dist_NBT = self._get_int(NR)
                    rxn.e_dist_INT = self._get_int(NR)

                # Number of outgoing energies in each E_out table
                NE = self._get_int()
                rxn.e_dist_energy_in = self._get_float(NE)
                L = self._get_float(NE)

                rxn.e_dist_intt = []        # Interpolation scheme (1=hist, 2=lin-lin)
                rxn.e_dist_energy_out = []  # Outgoing E grid for each incoming E
                rxn.e_dist_pdf = []         # Probability dist for " " "
                rxn.e_dist_cdf = []         # Cumulative dist for " " "
                for i in range(NE):
                    INTTp = self._get_int()
                    if INTTp > 10:
                        INTT = INTTp % 10
                        ND = (INTTp - INTT)/10
                    else:
                        INTT = INTTp
                        ND = 0
                    rxn.e_dist_intt.append(INTT)
                    if ND > 0:
                        print [rxn, ND, INTT]
                    

                    NP = self._get_int()
                    rxn.e_dist_energy_out.append(self._get_float(NP))
                    rxn.e_dist_pdf.append(self._get_float(NP))
                    rxn.e_dist_cdf.append(self._get_float(NP))
                 
            # General evaporation spectrum (ENDF-5 File 5 LF=5)
            elif LAW == 5:
                NR = self._get_int()
                if NR > 0:
                    rxn.e_dist_NBT = self._get_int(NR)
                    rxn.e_dist_INT = self._get_int(NR)
                
                NE = self._get_int()
                rxn.e_dist_energy_in = self._get_float(NE)
                rxn.e_dist_T = self._get_float(NE)
                NET = self._get_int()
                rxn.e_dist_X = self._get_float(NET)

            # Simple Maxwell fission spectrum (ENDF-6 File 5 LF=7) 
            elif LAW == 7:
                NR = self._get_int()
                if NR > 0:
                    rxn.e_dist_NBT = self._get_int(NR)
                    rxn.e_dist_INT = self._get_int(NR)

                NE = self._get_int()
                rxn.e_dist_energy_in = self._get_float(NE)
                rxn.e_dist_T = self._get_float(NE)
                rxn.e_dist_U = self._get_float()
                
            # Evaporation spectrum (ENDF-6 File 5 LF=9)
            elif LAW == 9:
                NR = self._get_int()
                if NR > 0:
                    rxn.e_dist_NBT = self._get_int(NR)
                    rxn.e_dist_INT = self._get_int(NR)

                NE = self._get_int()
                rxn.e_dist_energy_in = self._get_float(NE)
                rxn.e_dist_T = self._get_float(NE)
                rxn.e_dist_U = self._get_float()

            # Energy dependent Watt spectrum (ENDF-6 File 5 LF=11)
            elif LAW == 11:
                # Interpolation scheme between a's    
                NR = self._get_int()
                if NR > 0:
                    rxn.e_dist_NBTa = self._get_int(NR)
                    rxn.e_dist_INTa = self._get_int(NR)

                # Incident energy table and tabulated a's
                NE = self._get_int()
                rxn.e_dist_energya_in = self._get_float(NE)
                rxn.e_dist_a = self._get_float(NE)

                # Interpolation scheme between b's
                NR = self._get_int()
                if NR > 0:
                    rxn.e_dist_NBTb = self._get_int(NR)
                    rxn.e_dist_INTb = self._get_int(NR)

                # Incident energy table and tabulated b's
                NE = self._get_int()
                rxn.e_dist_energyb_in = self._get_float(NE)
                rxn.e_dist_b = self._get_float(NE)

                rxn.e_dist_U = self._get_float()

            # Tabular linear functions (UK Law 2)
            elif LAW == 22:
                # Interpolation scheme (not used in MCNP)
                NR = self._get_int()
                if NR > 0:
                    rxn.e_dist_NBT = self._get_int(NR)
                    rxn.e_dist_INT = self._get_int(NR)

                # Number of incident energies
                NE = self._get_int()
                rxn.e_dist_energy_in = self._get_float(NE)
                LOCE = self._get_int(NE)

                # Read linear functions
                rxn.e_dist_P = []
                rxn.e_dist_T = []
                rxn.e_dist_C = []
                for i in range(NE):
                    NF = self._get_int()
                    rxn.e_dist_P.append(self._get_float(NF))
                    rxn.e_dist_T.append(self._get_float(NF))
                    rxn.e_dist_C.append(self._get_float(NF))

            # From UK Law 6
            elif LAW == 24:
                # Interpolation scheme (not used in MCNP)
                NR = self._get_int()
                if NR > 0:
                    rxn.e_dist_NBT = self._get_int(NR)
                    rxn.e_dist_INT = self._get_int(NR)

                # Number of incident energies
                NE = self._get_int()
                rxn.e_dist_energy_in = self._get_float(NE)
                
                # Outgoing energy tables
                NET = self._get_int()
                rxn.e_dist_T = []
                for i in range(NE):
                    rxn.e_dist_T.append(self._get_float(NET))

            # Kalbach-87 Formalism (ENDF File 6 Law 1, LANG=2)
            elif LAW == 44:
                # Interpolation scheme
                NR = self._get_int()
                if NR > 0:
                    rxn.e_dist_NBT = self._get_int(NR)
                    rxn.e_dist_INT = self._get_int(NR)

                # Number of outgoing energies in each E_out table
                NE = self._get_int()
                rxn.e_dist_energy_in = self._get_float(NE)
                L = self._get_int(NE)

                rxn.e_dist_intt = []        # Interpolation scheme (1=hist, 2=lin-lin)
                rxn.e_dist_energy_out = []  # Outgoing E grid for each incoming E
                rxn.e_dist_pdf = []         # Probability dist for " " "
                rxn.e_dist_cdf = []         # Cumulative dist for " " "
                rxn.e_dist_frac = []        # Precompound fraction for " " "
                rxn.e_dist_ang = []         # Angular distribution slope for " " "
                for i in range(NE):
                    INTTp = self._get_int()
                    if INTTp > 10:
                        INTT = INTTp % 10
                        ND = (INTTp - INTT)/10
                    else:
                        INTT = INTTp
                    rxn.e_dist_intt.append(INTT)

                    NP = self._get_int()
                    rxn.e_dist_energy_out.append(self._get_float(NP))
                    rxn.e_dist_pdf.append(self._get_float(NP))
                    rxn.e_dist_cdf.append(self._get_float(NP))
                    rxn.e_dist_frac.append(self._get_float(NP))
                    rxn.e_dist_ang.append(self._get_float(NP))

            # Like 44, but tabular distribution instead of Kalbach-87
            elif LAW == 61:
                # Interpolation scheme
                NR = self._get_int()
                if NR > 0:
                    rxn.e_dist_NBT = self._get_int(NR)
                    rxn.e_dist_INT = self._get_int(NR)

                # Number of outgoing energies in each E_out table
                NE = self._get_int()
                rxn.e_dist_energy_in = self._get_float(NE)
                L = self._get_int(NE)

                rxn.e_dist_intt = []        # Interpolation scheme (1=hist, 2=lin-lin)
                rxn.e_dist_energy_out = []  # Outgoing E grid for each incoming E
                rxn.e_dist_pdf = []         # Probability dist for " " "
                rxn.e_dist_cdf = []         # Cumulative dist for " " "

                rxn.a_dist_intt = []
                rxn.a_dist_mu_out = [] # Cosine scattering angular grid
                rxn.a_dist_pdf = []    # Probability dist function
                rxn.a_dist_cdf = []
                for i in range(NE):
                    INTTp = self._get_int()
                    if INTTp > 10:
                        INTT = INTTp % 10
                        ND = (INTTp - INTT)/10
                    else:
                        INTT = INTTp
                    rxn.e_dist_intt.append(INTT)

                    # Secondary energy distribution
                    NP = self._get_int()
                    rxn.e_dist_energy_out.append(self._get_float(NP))
                    rxn.e_dist_pdf.append(self._get_float(NP))
                    rxn.e_dist_cdf.append(self._get_float(NP))

                    # Secondary angular distribution
                    LC = self._get_int(NP)
                    rxn.a_dist_intt.append([])
                    rxn.a_dist_mu_out.append([])
                    rxn.a_dist_pdf.append([])
                    rxn.a_dist_cdf.append([])
                    for j in range(NP):
                        rxn.a_dist_intt[-1].append(self._get_int())
                        NP = self._get_int()
                        rxn.a_dist_mu_out[-1].append(self._get_float(NP))
                        rxn.a_dist_pdf[-1].append(self._get_float(NP))
                        rxn.a_dist_cdf[-1].append(self._get_float(NP))

            # N-body phase space distribution (ENDF File 6 Law 6)
            elif LAW == 66:
                rxn.e_dist_nbodies = self._get_int()
                rxn.e_dist_massratio = self._get_float()

            # Laboratory angle-energy law (ENDF File 6 Law 7)
            elif LAW == 67:
                # Interpolation scheme
                NR = self._get_int()
                if NR > 0:
                    rxn.e_dist_NBT = self._get_int(NR)
                    rxn.e_dist_INT = self._get_int(NR)

                # Number of outgoing energies in each E_out table
                NE = self._get_int()
                rxn.e_dist_energy_in = self._get_float(NE)
                L = self._get_int(NE)

            if irxn+1 < len(LOCC):
                if self.index < LDIS + LOCC[irxn+1] - 1:
                    LNW = self._get_int()
                    LAW = self._get_int()
                    # print([rxn,rxn.e_dist_law,LAW])
                    
            # TODO: Read rest of data

    def _read_gpd(self):
        """Read total photon production cross section and secondary photon
        energies based on older 30x20 matrix formulation.

        """

        JXS12 = self.JXS[12]

        if JXS12 != 0:
            # Determine number of energies
            NE = self.NXS[3]

            # Read total photon production cross section
            self.index = JXS12
            self.sigma_photon = self._get_float(NE)

            # The following energies are the discrete incident neutron energies
            # for which the equiprobable secondary photon outgoing energies are
            # given
            self.e_in_photon_equi = [1.39e-10, 1.52e-7, 4.14e-7, 1.13e-6, 3.06e-6,
                                     8.32e-6,  2.26e-5, 6.14e-5, 1.67e-4, 4.54e-4,
                                     1.235e-3, 3.35e-3, 9.23e-3, 2.48e-2, 6.76e-2,
                                     0.184,    0.303,   0.500,   0.823,   1.353,
                                     1.738,    2.232,   2.865,   3.68,    6.07,
                                     7.79,     10.,     12.,     13.5,    15.]

            # Read equiprobable outgoing photon energies
            self.e_out_photon_equi = []
            for i in range(30):
                # Equiprobable outgoing photon energies for incident neutron
                # energy i
                self.e_out_photon_equi.append(self._get_float(20))

    def _read_mtrp(self):
        """
        Get the list of reaction MTs for photon-producing reactions for this
        cross-section table. The MT values are somewhat arbitrary.
        """

        LMT = self.JXS[13]
        NMT = self.NXS[6]
        for i in range(NMT):
            MT = int(self.XSS[LMT+i])
            self.photonReactions.append(Reaction(MT,self))

    def _read_lsigp(self):
        """Determine location of cross sections for each photon-producing reaction
        MT.
        """

        LXS = self.JXS[14]
        for i, rxn in enumerate(self.photonReactions):
            rxn.LOCA = int(self.XSS[LXS+i])

    def _read_sigp(self):
        """Read cross-sections for each photon-producing reaction MT.
        """

        JXS15 = self.JXS[15]
        for rxn in self.photonReactions:
            self.index = JXS15 + rxn.LOCA - 1
            MFTYPE = self._get_int()

            # Yield data taken from ENDF File 12 or 6
            if MFTYPE == 12 or MFTYPE == 16:
                MTMULT = self._get_int()

                # ENDF interpolation parameters
                NR = self._get_int()
                NBT = self._get_int(NR)
                INT = self._get_int(NR)

                # Energy-dependent yield
                NE = self._get_int()
                rxn.e_yield = self._get_float(NE)
                rxn.photon_yield = self._get_float(NE)

            # Cross-section data from ENDF File 13
            elif MFTYPE == 13:
                # Energy grid index at which data starts
                rxn.IE = self._get_int()

                # Cross sections
                NE = self._get_int()
                self.sigma = self._get_float(NE)
            else:
                raise

    def _read_landp(self):
        """Determine location of angular distribution for each photon-producing
        reaction MT.
        """

        JXS16 = self.JXS[16]
        for i, rxn in enumerate(self.photonReactions):
            rxn.LOCB = int(self.XSS[JXS16+i])

    def _read_andp(self):
        """Find the angular distribution for each photon-producing reaction
        MT."""

        JXS17 = self.JXS[17]
        for i, rxn in enumerate(self.photonReactions):
            if rxn.LOCB == 0:
                # No angular distribution data are given for this reaction,
                # isotropic scattering is asssumed in LAB
                continue

            self.index = JXS17 + rxn.LOCB - 1

            # Number of energies and incoming energy grid
            NE = self._get_int()
            self.a_dist_energy_in = self._get_float(NE)

            # Location of tables associated with each outgoing angle
            # distribution
            LC = self._get_int(NE)

            # 32 equiprobable cosine bins for each incoming energy
            self.a_dist_mu_out = {}
            for j, location in enumerate(LC):
                if location == 0:
                    continue
                self.index = JXS17 + location - 1
                self.a_dist_mu_out[j] = self._get_float(33)

    def _read_yp(self):
        """Read list of reactions required as photon production yield
        multipliers.
        """

        if self.NXS[6] != 0:
            self.index = self.JXS[20]
            NYP = self._get_int()
            self.MT_for_photon_yield = self._get_int(NYP)

    def _read_fis(self):
        """Read total fission cross-section data if present. Generally,
        this table is not provided since it is redundant.
        """

        # Check if fission block is present
        self.index = self.JXS[21]
        if self.index == 0:
            return

        # Read fission cross sections
        self.IE_fission = self._get_int()  # Energy grid index
        NE = self._get_int()
        self.sigma_f = self._get_float(NE)

    def _read_unr(self):
        """Read the unresolved resonance range probability tables if present
        """

        # Check if URR probability tables are present
        self.index = self.JXS[23]
        if self.index == 0:
            return

        N = self._get_int()   # Number of incident energies
        M = self._get_int()   # Length of probability table
        INT = self._get_int() # Interpolation parameter (2=lin-lin, 5=log-log)
        ILF = self._get_int() # Inelastic competition flag
        IOA = self._get_int() # Other absorption flag
        IFF = self._get_int() # Factors flag

        self.urr_energy = self._get_float(N) # Incident energies

        # Set up URR probability table
        self.urr_table = []
        for I in range(N):
            self.urr_table.append([])
            for J in range(6):
                self.urr_table[-1].append([None for K in range(M)])
        
        # Read values for URR probability tables
        for J in range(6):
            for K in range(M):
                for I in range(N):
                    self.urr_table[I][J][K] = self._get_float()

    def _get_float(self, n_values = 1):
        if n_values > 1:
            ind = self.index
            values = self.XSS[ind:ind+n_values]
            self.index += n_values
            return values
        else:
            value = self.XSS[self.index]
            self.index += 1
            return value
            
    def _get_int(self, n_values = 1):
        if n_values > 1:
            ind = self.index
            values = [int(i) for i in self.XSS[ind:ind+n_values]]
            self.index += n_values
            return values
        else:
            value = int(self.XSS[self.index])
            self.index += 1
            return value

    def find_reaction(self, MT):
        for r in self.reactions:
            if r.MT == MT:
                return r

    def __iter__(self):
        for r in self.reactions:
            yield r


class SabTable(AceTable):

    def __init__(self, name, awr, temp):
        super(SabTable, self).__init__(name, awr, temp)

    def _read_all(self):
        self._read_itie()
        self._read_itce()
        self._read_itxe()
        self._read_itca()

    def __repr__(self):
        if hasattr(self, 'name'):
            return "<ACE Thermal S(a,b) Table: {0}>".format(self.name)
        else:
            return "<ACE Thermal S(a,b) Table>"

    def plot_inelastic(self, index):
        E_out = [1e6 * E for E in self.inelastic_e_out[index]]
        mu_out = self.inelastic_mu_out[index]

        x = range(self.NXS[3]+1)
        y = E_out
        X,Y = meshgrid(x,y)
        Z = self._z_function(X,Y,E_out,mu_out)
        pcolor(X,Y,Z)

        # Plot options
        colorbar()
        xlim([0,self.NXS[3]])
        ylim([min(E_out),max(E_out)])
        xlabel('Equally likely cosine bin')
        ylabel('Outgoing Energy (eV)')

    def _z_function(self, x, y, E_out, mu_out):
        z = []
        for i in range(len(x)):
            x_ = x[i]
            y_ = y[i]
            z.append([mu_out[E_out.index(k)][j] for j,k in zip(x_,y_)])
        return z

    def _read_itie(self):
        """
        Read energy-dependent inelastic scattering cross sections
        """

        self.index = self.JXS[1]

        NE = self._get_int()
        self.inelastic_e_in = self._get_float(NE)
        self.inelastic_sigma = self._get_float(NE)

    def _read_itce(self):
        """
        Read energy-dependent elastic scattering cross sections
        """

        # Determine if ITCE block exists
        self.index = self.JXS[4]
        if self.index == 0:
            return

        # Read values
        NE = self._get_int()
        self.elastic_e_in = self._get_float(NE)
        self.elastic_P = self._get_float(NE)
        if self.NXS[5] == 4:
            self.elastic_type = 'sigma=P'
        else:
            self.elastic_type = 'sigma=P/E'

    def _read_itxe(self):
        """
        Read coupled energy/angle distributions for inelastic scattering
        """
        
        # Determine number of energies and angles
        NE_in = len(self.inelastic_e_in)
        NE_out = self.NXS[4]
        NMU = self.NXS[3]

        # Set index for reading values
        self.index = self.JXS[3]
        
        self.inelastic_e_out = []
        self.inelastic_mu_out = []
        # Loop over incoming energies 
        for i in range(NE_in):
            self.inelastic_e_out.append([])
            self.inelastic_mu_out.append([])

            # Loop over outgoing energies 
            for j in range(NE_out):
                # Read outgoing energy
                self.inelastic_e_out[-1].append(self._get_float())

                # Read discrete cosines for scattering from E_in to E_out
                self.inelastic_mu_out[-1].append(self._get_float(NMU+1))

    def _read_itca(self):
        """Read angular distributions for elastic scattering.
        """

        NMU = self.NXS[6]
        if self.JXS[4] == 0 or NMU == -1:
            return

        self.index = self.JXS[6]

        NE = len(self.elastic_e_in)
        self.elastic_mu_out = []
        for i in range(NE):
            self.elastic_mu_out.append(self._get_float(NMU))

    def _get_float(self, n_values = 1):
        if n_values > 1:
            values = self.XSS[self.index:self.index+n_values]
            self.index += n_values
            return values
        else:
            value = self.XSS[self.index]
            self.index += 1
            return value
            
    def _get_int(self, n_values = 1):
        if n_values > 1:
            values = [int(i) for i in 
                      self.XSS[self.index:self.index+n_values]]
            self.index += n_values
            return values
        else:
            value = int(self.XSS[self.index])
            self.index += 1
            return value

            
class Reaction(object):

    def __init__(self, MT, table=None):
        self.table = table # Reference to containing table
        self.MT = MT       # MT value
        self.Q = None      # Q-value
        self.TY = None     # Neutron release
        self.LOCA = None
        self.LOCB = None
        self.LOCC = None
        self.IE = None     # Energy grid index
        self.sigma = []    # Cross section values

    def broaden(self, T_high):
        pass        

    def threshold(self):
        """Return energy threshold for this reaction"""

        table = self.table
        return table.energy[self.IE]

    def __repr__(self):
        try:
            return "<ACE Reaction: MT={0} {1}>".format(
                self.MT, reaction_names[self.MT])
        except KeyError:
            return "<ACE Reaction: Unknown MT={0}>".format(self.MT)


class DosimetryTable(AceTable):

    def __init__(self, name, awr, temp):
        super(DosimetryTable, self).__init__(name, awr, temp)

    def __repr__(self):
        if hasattr(self, 'name'):
            return "<ACE Dosimetry Table: {0}>".format(self.name)
        else:
            return "<ACE Dosimetry Table>"
        

class NeutronDiscreteTable(AceTable):

    def __init__(self, name, awr, temp):
        super(NeutronDiscreteTable, self).__init__(name, awr, temp)

    def __repr__(self):
        if hasattr(self, 'name'):
            return "<ACE Discrete-E Neutron Table: {0}>".format(self.name)
        else:
            return "<ACE Discrete-E Neutron Table>"
        

class NeutronMGTable(AceTable):

    def __init__(self, name, awr, temp):
        super(NeutronMGTable, self).__init__(name, awr, temp)

    def __repr__(self):
        if hasattr(self, 'name'):
            return "<ACE Multigroup Neutron Table: {0}>".format(self.name)
        else:
            return "<ACE Multigroup Neutron Table>"
        

class PhotoatomicTable(AceTable):

    def __init__(self, name, awr, temp):
        super(PhotoatomicTable, self).__init__(name, awr, temp)

    def __repr__(self):
        if hasattr(self, 'name'):
            return "<ACE Continuous-E Photoatomic Table: {0}>".format(self.name)
        else:
            return "<ACE Continuous-E Photoatomic Table>"
        

class PhotoatomicMGTable(AceTable):

    def __init__(self, name, awr, temp):
        super(PhotoatomicMGTable, self).__init__(name, awr, temp)

    def __repr__(self):
        if hasattr(self, 'name'):
            return "<ACE Multigroup Photoatomic Table: {0}>".format(self.name)
        else:
            return "<ACE Multigroup Photoatomic Table>"
        

class ElectronTable(AceTable):

    def __init__(self, name, awr, temp):
        super(ElectronTable, self).__init__(name, awr, temp)

    def __repr__(self):
        if hasattr(self, 'name'):
            return "<ACE Electron Table: {0}>".format(self.name)
        else:
            return "<ACE Electron Table>"
        

class PhotonuclearTable(AceTable):

    def __init__(self, name, awr, temp):
        super(PhotonuclearTable, self).__init__(name, awr, temp)

    def __repr__(self):
        if hasattr(self, 'name'):
            return "<ACE Photonuclear Table: {0}>".format(self.name)
        else:
            return "<ACE Photonuclear Table>"

table_types = {
    "c": NeutronTable,
    "t": SabTable,
    "y": DosimetryTable,
    "d": NeutronDiscreteTable,
    "p": PhotoatomicTable,
    "m": NeutronMGTable,
    "g": PhotoatomicMGTable,
    "e": ElectronTable,
    "u": PhotonuclearTable}

reaction_names = {
    # TODO: This should be provided as part of the ENDF module functionality
    1: '(n,total)',
    2: '(n,elastic)',
    3: '(n,nonelastic)',
    4: '(n,inelastic)',
    5: '(misc)',
    10: '(n,continuum)',
    11: '(n,2n d)',
    16: '(n,2n)',
    17: '(n,3n)',
    18: '(n,fission)',
    19: '(n,f)',
    20: '(n,nf)',
    21: '(n,2nf)',
    22: '(n,na)',
    23: '(n,n3a)',
    24: '(n,2na)',
    25: '(n,3na)',
    28: '(n,np)',
    29: '(n,n2a)',
    30: '(n,2n2a)',
    32: '(n,nd)',
    33: '(n,nt)',
    34: '(n,n He-3)',
    35: '(n,nd3a)',
    36: '(n,nt2a)',
    37: '(n,4n)',
    38: '(n,3nf)',
    41: '(n,2np)',
    42: '(n,3np)',
    44: '(n,2np)',
    45: '(n,npa)',
    91: '(n,nc)',
    102: '(n,gamma)',
    103: '(n,p)',
    104: '(n,d)',
    105: '(n,t)',
    106: '(n,3He)',
    107: '(n,a)',
    108: '(n,2a)',
    109: '(n,3a)',
    111: '(n,2p)',
    112: '(n,pa)',
    113: '(n,t2a)',
    114: '(n,d2a)',
    115: '(n,pd)',
    116: '(n,pt)',
    117: '(n,da)',
    201: '(n,Xn)',
    202: '(n,Xgamma)',
    203: '(n,Xp)',
    204: '(n,Xd)',
    205: '(n,Xt)',
    206: '(n,X3He)',
    207: '(n,Xa)',
    444: '(damage)',
    649: '(n,pc)',
    699: '(n,dc)',
    749: '(n,tc)',
    799: '(n,3Hec)',
    849: '(n,ac)',
    }
"""Dictionary of MT reaction labels"""
reaction_names.update({mt: '(n,n{0})'.format(mt - 50) for mt in range(50, 91)})
reaction_names.update({mt: '(n,p{0})'.format(mt - 600) for mt in range(600, 649)})
reaction_names.update({mt: '(n,d{0})'.format(mt - 650) for mt in range(650, 699)})
reaction_names.update({mt: '(n,t{0})'.format(mt - 700) for mt in range(700, 649)})
reaction_names.update({mt: '(n,3He{0})'.format(mt - 750) for mt in range(750, 799)})
reaction_names.update({mt: '(n,a{0})'.format(mt - 800) for mt in range(700, 649)})


if __name__ == '__main__':
    # Might be nice to check environment variable DATAPATH to search for xsdir
    # and list files that could be read?
    pass
