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

.. moduleauthor:: Paul Romano <romano7@gmail.com>, Anthony Scopatz
"""

import struct
from warnings import warn
from collections import OrderedDict

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

    :attributes:
      **binary** : bool
        Identifies Whether the library is in binary format or not

      **tables** : dict
        Dictionary whose keys are the names of the ACE tables and whose values
        are the instances of subclasses of AceTable (e.g. NeutronTable)

      **verbose** : bool
        Determines whether output is printed to the stdout when reading a
        Library

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
            table.NXS.insert(0, 0)
            table.NXS = np.array(table.NXS, dtype=int)

            table.JXS.insert(0, 0)
            table.JXS = np.array(table.JXS, dtype=int)

            table.XSS.insert(0, 0.0)
            table.XSS = np.array(table.XSS, dtype=float)

            # Read all data blocks
            table._read_all()

            # Advance to next record
            self.f.seek(start_position + recl_length*(n_records + 1))

    def _read_ascii(self, table_names):

        lines = self.f.readlines()
        
        while 0 != len(lines):
            # Read name of table, atomic weight ratio, and temperature. If first
            # line is empty, we are at end of file
            words = lines[0].split()
            name = words[0]
            awr = float(words[1])
            temp = float(words[2])

            datastr = '0 ' + ' '.join(lines[6:8])
            nxs = np.fromstring(datastr, sep=' ', dtype=int)
            n_lines = (nxs[1] + 3)/4

            # verify that we are suppossed to read this table in
            if (table_names is not None) and (name not in table_names):
                lines = lines[12+n_lines:]
                continue

            # ensure we have a valid table type
            if 0 == len(name) or name[-1] not in table_types:
                warn("Unsupported table: " + name, RuntimeWarning)
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

            # Add NXS, JXS, and XSS arrays to table
            # Insert empty object at beginning of NXS, JXS, and XSS
            # arrays so that the indexing will be the same as
            # Fortran. This makes it easier to follow the ACE format
            # specification.
            table.NXS = nxs

            datastr = '0 ' + ' '.join(lines[8:12])
            table.JXS = np.fromstring(datastr, sep=' ', dtype=int)

            datastr = '0.0 ' + ' '.join(lines[12:12+n_lines])
            table.XSS = np.fromstring(datastr, sep=' ', dtype=float)

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
    Library object and stored within the ``tables`` attribute of a Library
    object.

    Parameters
    ----------
    name : str
        ZAID identifier of the table, e.g. '92235.70c'.
    awr : float
        Atomic weight ratio of the target nuclide.
    temp : float
        Temperature of the target nuclide in eV.
    
    :Attributes:
      **awr** : float
        Atomic weight ratio of the target nuclide.

      **energy** : list of floats
        The energy values (MeV) at which reaction cross-sections are tabulated.

      **name** : str
        ZAID identifier of the table, e.g. 92235.70c.

      **nu_p_energy** : list of floats
        Energies in MeV at which the number of prompt neutrons emitted per
        fission is tabulated.

      **nu_p_type** : str
        Indicates how number of prompt neutrons emitted per fission is
        stored. Can be either "polynomial" or "tabular".

      **nu_p_value** : list of floats
        The number of prompt neutrons emitted per fission, if data is stored in
        "tabular" form, or the polynomial coefficients for the "polynomial"
        form.

      **nu_t_energy** : list of floats
        Energies in MeV at which the total number of neutrons emitted per
        fission is tabulated.

      **nu_t_type** : str
        Indicates how total number of neutrons emitted per fission is
        stored. Can be either "polynomial" or "tabular".

      **nu_t_value** : list of floats
        The total number of neutrons emitted per fission, if data is stored in
        "tabular" form, or the polynomial coefficients for the "polynomial"
        form.

      **reactions** : list of Reactions
        A list of Reaction instances containing the cross sections, secondary
        angle and energy distributions, and other associated data for each
        reaction for this nuclide.

      **sigma_a** : list of floats
        The microscopic absorption cross section for each value on the energy
        grid.

      **sigma_t** : list of floats
        The microscopic total cross section for each value on the energy grid.

      **temp** : float
        Temperature of the target nuclide in eV.

    """

    def __init__(self, name, awr, temp):
        super(NeutronTable, self).__init__(name, awr, temp)
        self.reactions = OrderedDict()
        self.photon_reactions = OrderedDict()

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
        self._read_yp()
        self._read_fis()
        self._read_unr()

    def _read_esz(self):
        """Read ESZ block -- this block contains the energy grid, total
        xs, absorption xs, elastic scattering xs, and heating numbers.
        """
        
        NE = self.NXS[3]
        self.index = self.JXS[1]

        arr = self.XSS[self.index:self.index+NE*5]
        arr.shape = (5, NE)
        self.energy, self.sigma_t, self.sigma_a, sigma_el, self.heating = arr
        self.index += NE*5

        # Create elastic scattering reaction
        MT = 2
        rxn = Reaction(MT, self)
        rxn.Q = 0.0
        rxn.IE = 1
        rxn.TY = 1
        rxn.sigma = sigma_el
        self.reactions[MT] = rxn

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

            # FIXME The following code never will save LOCC on the object!
            # Energy distribution for delayed fission neutrons
            #LED = self.JXS[26]
            #LOCC = {}
            #for group in range(n_group):
            #    LOCC[group] = self.XSS[LED + group]

    def _read_mtr(self):
        """Get the list of reaction MTs for this cross-section table. The
        MT values are somewhat arbitrary.
        """
        LMT = self.JXS[3]
        NMT = self.NXS[4]
        mts = np.asarray(self.XSS[LMT:LMT+NMT], dtype=int)
        rxs = [(mt, Reaction(mt, self)) for mt in mts]
        self.reactions.update(rxs)
            
    def _read_lqr(self):
        """Find Q-values for each reaction MT
        """
        JXS4 = self.JXS[4]
        for i, rxn in enumerate(self.reactions.values()[1:]):
            rxn.Q = self.XSS[JXS4+i]

    def _read_tyr(self):
        """Find the neutron release for each reaction MT. A neutron
        release of 19 indicates fission. A neutron release greater
        than 100 signifies reactions other than fission taht have
        energy-dependent neutron multiplicities
        """
        NMT = self.NXS[4]
        JXS5 = self.JXS[5]
        tys = np.asarray(self.XSS[JXS5:JXS5+NMT], dtype=int)
        for ty, rxn in zip(tys, self.reactions.values()[1:]):
            rxn.TY = ty

    def _read_lsig(self):
        """Determine location of cross sections for each reaction MT
        """
        NMT = self.NXS[4]
        LXS = self.JXS[6]
        loca = np.asarray(self.XSS[LXS:LXS+NMT], dtype=int)
        for loc, rxn in zip(loca, self.reactions.values()[1:]):
            rxn.LOCA = loc

    def _read_sig(self):
        """Read cross-sections for each reaction MT
        """
        JXS7 = self.JXS[7]
        for rxn in self.reactions.values()[1:]:
            rxn.IE = int(self.XSS[JXS7+rxn.LOCA-1])
            NE = int(self.XSS[JXS7+rxn.LOCA])
            rxn.sigma = self.XSS[JXS7+rxn.LOCA+1:JXS7+rxn.LOCA+1+NE]

    def _read_land(self):
        """Find locations for angular distributions
        """
        JXS8 = self.JXS[8]

        # Number of reactions is less than total since we only need
        # angular distribution for reactions with secondary
        # neutrons. Thus, MT > 100 are not included.
        NMT = self.NXS[5]

        # Need NMT + 1 since elastic scattering is included
        locb = np.asarray(self.XSS[JXS8:JXS8+NMT+1], dtype=int)
        for loc, rxn in zip(locb, self.reactions.values()[:NMT+1]):
            rxn.LOCB = loc

    def _read_and(self):
        """Find the angular distribution for each reaction MT
        """
        JXS9 = self.JXS[9]
        NMT = self.NXS[5]

        # Angular distribution for all MT with secondary neutrons
        # including elastic scattering
        for rxn in self.reactions.values()[:NMT+1]:
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

            ind = JXS9 + rxn.LOCB - 1

            NE = int(self.XSS[ind])
            rxn.ang_energy_in = self.XSS[ind+1:ind+1+NE]
            LC = np.asarray(self.XSS[ind+1+NE:ind+1+2*NE], dtype=int)
            rxn.ang_location = LC
            ind = ind+1+2*NE

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
                    rxn.ang_cos[i] = self.XSS[ind:ind+33]
                    ind += 33
                elif location < 0:
                    # Tabular angular distribution
                    JJ = int(self.XSS[ind])
                    NP = int(self.XSS[ind+1])
                    ind += 2
                    ang_dat = self.XSS[ind:ind+3*NP]
                    ang_dat.shape = (3, NP)
                    rxn.ang_cos[i], rxn.ang_pdf[i], rxn.ang_cdf[i] = ang_dat
                    ind += 3 * NP

            self.index = ind
        
    def _read_ldlw(self):
        """Find locations for energy distribution data for each reaction
        """
        LED = self.JXS[10]

        # Number of reactions is less than total since we only need
        # energy distribution for reactions with secondary
        # neutrons. Thus, MT > 100 are not included. Elastic
        # scattering is also not included.
        NMT = self.NXS[5]
        locc = np.asarray(self.XSS[LED:LED+NMT], dtype=int)
        for loc, rxn in zip(locc, self.reactions.values()[1:NMT+1]):
            rxn.LOCC = loc

    def _read_dlw(self):
        """Determine the energy distribution for secondary neutrons for
        each reaction MT
        """
        LDIS = self.JXS[11]
        NMT = self.NXS[5]

        rxs = self.reactions.values()[1:NMT+1]
        for irxn, rxn in enumerate(rxs):
            ind = LDIS + rxn.LOCC - 1
            LNW = int(self.XSS[ind])
            LAW = int(self.XSS[ind+1])
            IDAT = int(self.XSS[ind+2])
            NR = int(self.XSS[ind+3])
            ind += 4
            if NR > 0:
                dat = np.asarray(self.XSS[ind:ind+2*NR], dtype=int)
                dat.shape = (2, NR)
                interp_NBT, interp_INT = dat
                ind += 2 * NR

            # Determine tabular energy points and probability of law
            # validity
            NE = int(self.XSS[ind])
            dat = self.XSS[ind+1:ind+1+2*NE]
            dat.shape = (2, NE)
            rxn.e_dist_energy, rxn.e_dist_pvalid = dat

            rxn.e_dist_law = LAW
            ind = LDIS + IDAT - 1
            self.index = ind

            if LAW == 1:
                # Tabular equiprobable energy bins (ENDF Law 1)
                NR = int(self.XSS[ind])
                ind += 1
                if NR > 0:
                    dat = np.asarray(self.XSS[ind:ind+2*NR], dtype=int)
                    dat.shape = (2, NR)
                    rxn.e_dist_NBT, rxn.e_dist_INT = dat
                    ind += 2 * NR                    

                # Number of outgoing energies in each E_out table
                NE = int(self.XSS[ind])
                rxn.e_dist_energy_in = self.XSS[ind+1:ind+1+NE]
                ind += 1 + NE

                # Read E_out tables
                NET = int(self.XSS[ind])
                dat = self.XSS[ind+1:ind+1+3*NET]
                dat.shape = (3, NET)
                self.e_dist_energy_out1, self.e_dist_energy_out2, \
                                         self.e_dist_energy_outNE = dat
                ind += 1 + 3 * NET
            elif LAW == 2:
                # Discrete photon energy
                self.e_dist_LP = int(self.XSS[ind])
                self.e_dist_EG = self.XSS[ind+1]
                ind += 2
            elif LAW == 3:
                # Level scattering (ENDF Law 3)
                rxn.e_dist_data = self.XSS[ind:ind+2]
                ind += 2
            elif LAW == 4:
                # Continuous tabular distribution (ENDF Law 1)
                NR = int(self.XSS[ind])
                ind += 1
                if NR > 0:
                    dat = np.asarray(self.XSS[ind:ind+2*NR], dtype=int)
                    dat.shape = (2, NR)
                    rxn.e_dist_NBT, rxn.e_dist_INT = dat
                    ind += 2 * NR                    

                # Number of outgoing energies in each E_out table
                NE = int(self.XSS[ind])
                rxn.e_dist_energy_in = self.XSS[ind+1:ind+1+NE]
                L = self.XSS[ind+1+NE:ind+1+2*NE]
                ind += 1 + 2*NE

                nps = []
                rxn.e_dist_intt = []        # Interpolation scheme (1=hist, 2=lin-lin)
                rxn.e_dist_energy_out = []  # Outgoing E grid for each incoming E
                rxn.e_dist_pdf = []         # Probability dist for " " "
                rxn.e_dist_cdf = []         # Cumulative dist for " " "
                for i in range(NE):
                    INTTp = int(self.XSS[ind])
                    if INTTp > 10:
                        INTT = INTTp % 10
                        ND = (INTTp - INTT)/10
                    else:
                        INTT = INTTp
                        ND = 0
                    rxn.e_dist_intt.append(INTT)
                    #if ND > 0:
                    #    print [rxn, ND, INTT]

                    NP = int(self.XSS[ind+1])
                    nps.append(NP)
                    dat = self.XSS[ind+2:ind+2+3*NP]
                    dat.shape = (3, NP)
                    rxn.e_dist_energy_out.append(dat[0])
                    rxn.e_dist_pdf.append(dat[1])
                    rxn.e_dist_cdf.append(dat[2])
                    ind += 2 + 3*NP

                # convert to arrays if possible
                rxn.e_dist_intt = np.array(rxn.e_dist_intt)
                nps = np.array(nps)
                if all((nps[1:] - nps[:-1]) == 0):
                    rxn.e_dist_energy_out = np.array(rxn.e_dist_energy_out)
                    rxn.e_dist_pdf = np.array(rxn.e_dist_pdf)
                    rxn.e_dist_cdf = np.array(rxn.e_dist_cdf)
            elif LAW == 5:
                # General evaporation spectrum (ENDF-5 File 5 LF=5)
                NR = int(self.XSS[ind])
                ind += 1
                if NR > 0:
                    dat = np.asarray(self.XSS[ind:ind+2*NR], dtype=int)
                    dat.shape = (2, NR)
                    rxn.e_dist_NBT, rxn.e_dist_INT = dat
                    ind += 2 * NR                    
                
                NE = int(self.XSS[ind])
                rxn.e_dist_energy_in = self.XSS[ind+1:ind+1+NE]
                rxn.e_dist_T = self.XSS[ind+1+NE:ind+1+2*NE]
                ind += 1+ 2*NE

                NET = int(self.XSS[ind])
                rxn.e_dist_X = self.XSS[ind+1:ind+1+NET]
                ind += 1 + NET
            elif LAW == 7:
                # Simple Maxwell fission spectrum (ENDF-6 File 5 LF=7) 
                NR = int(self.XSS[ind])
                ind += 1
                if NR > 0:
                    dat = np.asarray(self.XSS[ind:ind+2*NR], dtype=int)
                    dat.shape = (2, NR)
                    rxn.e_dist_NBT, rxn.e_dist_INT = dat
                    ind += 2 * NR                    

                NE = int(self.XSS[ind])
                rxn.e_dist_energy_in = self.XSS[ind+1:ind+1+NE]
                rxn.e_dist_T = self.XSS[ind+1+NE:ind+1+2*NE]
                rxn.e_dist_U = self.XSS[ind+1+2*NE]
                ind += 2 + 2*NE
            elif LAW == 9:
                # Evaporation spectrum (ENDF-6 File 5 LF=9)
                NR = int(self.XSS[ind])
                ind += 1
                if NR > 0:
                    dat = np.asarray(self.XSS[ind:ind+2*NR], dtype=int)
                    dat.shape = (2, NR)
                    rxn.e_dist_NBT, rxn.e_dist_INT = dat
                    ind += 2 * NR                    

                NE = int(self.XSS[ind])
                rxn.e_dist_energy_in = self.XSS[ind+1:ind+1+NE]
                rxn.e_dist_T = self.XSS[ind+1+NE:ind+1+2*NE]
                rxn.e_dist_U = self.XSS[ind+1+2*NE]
                ind += 2 + 2*NE
            elif LAW == 11:
                # Energy dependent Watt spectrum (ENDF-6 File 5 LF=11)
                # Interpolation scheme between a's    
                NR = int(self.XSS[ind])
                ind += 1
                if NR > 0:
                    dat = np.asarray(self.XSS[ind:ind+2*NR], dtype=int)
                    dat.shape = (2, NR)
                    rxn.e_dist_NBTa, rxn.e_dist_INTa = dat
                    ind += 2 * NR                    

                # Incident energy table and tabulated a's
                NE = int(self.XSS[ind])
                rxn.e_dist_energya_in = self.XSS[ind+1:ind+1+NE]
                rxn.e_dist_a = self.XSS[ind+1+NE:ind+1+2*NE]
                ind += 1 + 2*NE

                # Interpolation scheme between b's
                NR = int(self.XSS[ind])
                ind += 1
                if NR > 0:
                    dat = np.asarray(self.XSS[ind:ind+2*NR], dtype=int)
                    dat.shape = (2, NR)
                    rxn.e_dist_NBTb, rxn.e_dist_INTb = dat
                    ind += 2 * NR                    

                # Incident energy table and tabulated b's
                NE = int(self.XSS[ind])
                rxn.e_dist_energyb_in = self.XSS[ind+1:ind+1+NE]
                rxn.e_dist_b = self.XSS[ind+1+NE:ind+1+2*NE]

                rxn.e_dist_U = self.XSS[ind+1+2*NE]
                ind += 2 + 2*NE
            elif LAW == 22:
                # Tabular linear functions (UK Law 2)
                # Interpolation scheme (not used in MCNP)
                NR = int(self.XSS[ind])
                ind += 1
                if NR > 0:
                    dat = np.asarray(self.XSS[ind:ind+2*NR], dtype=int)
                    dat.shape = (2, NR)
                    rxn.e_dist_NBT, rxn.e_dist_INT = dat
                    ind += 2 * NR                    

                # Number of incident energies
                NE = int(self.XSS[ind])
                rxn.e_dist_energy_in = self.XSS[ind+1:ind+1+NE]
                LOCE = np.asarray(self.XSS[ind+1+NE:ind+1+2*NE], dtype=int)
                ind += 1 + 2*NE

                # Read linear functions
                nfs = []
                rxn.e_dist_P = []
                rxn.e_dist_T = []
                rxn.e_dist_C = []
                for i in range(NE):
                    NF = int(self.XSS[ind])
                    nfs.append(NF)
                    dat = self.XSS[ind+1:ind+1+3*NF]
                    dat.shape = (3, NF)
                    rxn.e_dist_P.append(dat[0])
                    rxn.e_dist_T.append(dat[1])
                    rxn.e_dist_C.append(dat[2])
                    ind += 1 + 3*NF

                # convert to arrays if possible
                nfs = np.array(nfs)
                if all((nfs[1:] - nfs[:-1]) == 0):
                    rxn.e_dist_P = np.array(rxn.e_dist_P)
                    rxn.e_dist_T = np.array(rxn.e_dist_T)
                    rxn.e_dist_C = np.array(rxn.e_dist_C)
            elif LAW == 24:
                # From UK Law 6
                # Interpolation scheme (not used in MCNP)
                NR = int(self.XSS[ind])
                ind += 1
                if NR > 0:
                    dat = np.asarray(self.XSS[ind:ind+2*NR], dtype=int)
                    dat.shape = (2, NR)
                    rxn.e_dist_NBT, rxn.e_dist_INT = dat
                    ind += 2 * NR                    

                # Number of incident energies
                NE = int(self.XSS[ind])
                rxn.e_dist_energy_in = self.XSS[ind+1:ind+1+NE]
                ind += 1 + NE
                
                # Outgoing energy tables
                NET = int(self.XSS[ind])
                rxn.e_dist_T = self.XSS[ind+1:ind+1+NE*NET]
                rxn.e_dist_T.shape = (NE, NET)
                ind += 1 + NE*NET
            elif LAW == 44:
                # Kalbach-87 Formalism (ENDF File 6 Law 1, LANG=2)
                # Interpolation scheme
                NR = int(self.XSS[ind])
                ind += 1
                if NR > 0:
                    dat = np.asarray(self.XSS[ind:ind+2*NR], dtype=int)
                    dat.shape = (2, NR)
                    rxn.e_dist_NBT, rxn.e_dist_INT = dat
                    ind += 2 * NR                    

                # Number of outgoing energies in each E_out table
                NE = int(self.XSS[ind])
                rxn.e_dist_energy_in = self.XSS[ind+1:ind+1+NE]
                L = np.asarray(self.XSS[ind+1+NE:ind+1+2*NE], dtype=int)
                ind += 1 + 2*NE

                nps = []
                rxn.e_dist_intt = []        # Interpolation scheme (1=hist, 2=lin-lin)
                rxn.e_dist_energy_out = []  # Outgoing E grid for each incoming E
                rxn.e_dist_pdf = []         # Probability dist for " " "
                rxn.e_dist_cdf = []         # Cumulative dist for " " "
                rxn.e_dist_frac = []        # Precompound fraction for " " "
                rxn.e_dist_ang = []         # Angular distribution slope for " " "
                for i in range(NE):
                    INTTp = int(self.XSS[ind])
                    if INTTp > 10:
                        INTT = INTTp % 10
                        ND = (INTTp - INTT)/10
                    else:
                        INTT = INTTp
                    rxn.e_dist_intt.append(INTT)

                    NP = int(self.XSS[ind+1])
                    nps.append(NP)
                    ind += 2

                    dat = self.XSS[ind:ind+5*NP]
                    dat.shape = (5, NP)
                    rxn.e_dist_energy_out.append(dat[0])
                    rxn.e_dist_pdf.append(dat[1])
                    rxn.e_dist_cdf.append(dat[2])
                    rxn.e_dist_frac.append(dat[3])
                    rxn.e_dist_ang.append(dat[4])
                    ind += 5*NP

                # convert to arrays if possible
                rxn.e_dist_intt = np.array(rxn.e_dist_intt)
                nps = np.array(nps)
                if all((nps[1:] - nps[:-1]) == 0):
                    rxn.e_dist_energy_out = np.array(rxn.e_dist_energy_out)
                    rxn.e_dist_pdf = np.array(rxn.e_dist_pdf)
                    rxn.e_dist_cdf = np.array(rxn.e_dist_cdf)
            elif LAW == 61:
                # Like 44, but tabular distribution instead of Kalbach-87
                # Interpolation scheme
                NR = int(self.XSS[ind])
                ind += 1
                if NR > 0:
                    dat = np.asarray(self.XSS[ind:ind+2*NR], dtype=int)
                    dat.shape = (2, NR)
                    rxn.e_dist_NBT, rxn.e_dist_INT = dat
                    ind += 2 * NR                    

                # Number of outgoing energies in each E_out table
                NE = int(self.XSS[ind])
                rxn.e_dist_energy_in = self.XSS[ind+1:ind+1+NE]
                L = np.asarray(self.XSS[ind+1+NE:ind+1+2*NE], dtype=int)
                ind += 1 + 2*NE

                npes = []
                rxn.e_dist_intt = []        # Interpolation scheme (1=hist, 2=lin-lin)
                rxn.e_dist_energy_out = []  # Outgoing E grid for each incoming E
                rxn.e_dist_pdf = []         # Probability dist for " " "
                rxn.e_dist_cdf = []         # Cumulative dist for " " "

                npas = []
                rxn.a_dist_intt = []
                rxn.a_dist_mu_out = [] # Cosine scattering angular grid
                rxn.a_dist_pdf = []    # Probability dist function
                rxn.a_dist_cdf = []
                for i in range(NE):
                    INTTp = int(self.XSS[ind])
                    if INTTp > 10:
                        INTT = INTTp % 10
                        ND = (INTTp - INTT)/10
                    else:
                        INTT = INTTp
                    rxn.e_dist_intt.append(INTT)

                    # Secondary energy distribution
                    NPE = int(self.XSS[ind+1])
                    npes.append(NPE)
                    dat = self.XSS[ind+2:ind+2+4*NPE]
                    dat.shape = (4, NPE)
                    rxn.e_dist_energy_out.append(dat[0])
                    rxn.e_dist_pdf.append(dat[1])
                    rxn.e_dist_cdf.append(dat[2])
                    LC = np.asarray(dat[3], dtype=int)
                    ind += 2 + 4*NPE

                    # Secondary angular distribution
                    rxn.a_dist_intt.append([])
                    rxn.a_dist_mu_out.append([])
                    rxn.a_dist_pdf.append([])
                    rxn.a_dist_cdf.append([])
                    for j in range(NPE):
                        rxn.a_dist_intt[-1].append(int(self.XSS[ind]))
                        NPA = int(self.XSS[ind+1])
                        npas.append(NPA)
                        dat = self.XSS[ind+2:ind+2+3*NPA]
                        dat.shape = (3, NPA)
                        rxn.a_dist_mu_out[-1].append(dat[0])
                        rxn.a_dist_pdf[-1].append(dat[1])
                        rxn.a_dist_cdf[-1].append(dat[2])
                        ind += 2 + 3*NPA

                # convert to arrays if possible
                rxn.e_dist_intt = np.array(rxn.e_dist_intt)
                npes = np.array(npes)
                npas = np.array(npas)
                if all((npes[1:] - npes[:-1]) == 0):
                    rxn.e_dist_energy_out = np.array(rxn.e_dist_energy_out)
                    rxn.e_dist_pdf = np.array(rxn.e_dist_pdf)
                    rxn.e_dist_cdf = np.array(rxn.e_dist_cdf)

                    rxn.a_dist_intt = np.array(rxn.a_dist_intt)
                    if all((npas[1:] - npas[:-1]) == 0):
                        rxn.a_dist_mu_out = np.array(rxn.a_dist_mu_out)
                        rxn.a_dist_pdf = np.array(rxn.a_dist_pdf)
                        rxn.a_dist_cdf = np.array(rxn.a_dist_cdf)
            elif LAW == 66:
                # N-body phase space distribution (ENDF File 6 Law 6)
                rxn.e_dist_nbodies = int(self.XSS[ind])
                rxn.e_dist_massratio = self.XSS[ind+1]
                ind += 2
            elif LAW == 67:
                # Laboratory angle-energy law (ENDF File 6 Law 7)
                # Interpolation scheme
                NR = int(self.XSS[ind])
                ind += 1
                if NR > 0:
                    dat = np.asarray(self.XSS[ind:ind+2*NR], dtype=int)
                    dat.shape = (2, NR)
                    rxn.e_dist_NBT, rxn.e_dist_INT = dat
                    ind += 2 * NR                    

                # Number of outgoing energies in each E_out table
                NE = int(self.XSS[ind])
                rxn.e_dist_energy_in = self.XSS[ind+1:ind+1+NE]
                L = np.asarray(self.XSS[ind+1+NE:ind+1+2*NE], dtype=int)
                ind += 1 + 2*NE


            # Bump up index for next loop
            if irxn+1 < NMT:
                if ind < LDIS + rxs[irxn+1].LOCC - 1:
                    LNW = int(self.XSS[ind])
                    LAW = int(self.XSS[ind+1])
                    ind += 2
                    
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
            ind = JXS12
            self.sigma_photon = self.XSS[ind:ind+NE]
            ind += NE

            # The following energies are the discrete incident neutron energies
            # for which the equiprobable secondary photon outgoing energies are
            # given
            self.e_in_photon_equi = np.array(
                                    [1.39e-10, 1.52e-7, 4.14e-7, 1.13e-6, 3.06e-6,
                                     8.32e-6,  2.26e-5, 6.14e-5, 1.67e-4, 4.54e-4,
                                     1.235e-3, 3.35e-3, 9.23e-3, 2.48e-2, 6.76e-2,
                                     0.184,    0.303,   0.500,   0.823,   1.353,
                                     1.738,    2.232,   2.865,   3.68,    6.07,
                                     7.79,     10.,     12.,     13.5,    15.])

            # Read equiprobable outgoing photon energies
            # Equiprobable outgoing photon energies for incident neutron
            # energy i
            self.e_out_photon_equi = self.XSS[ind:ind+600]
            self.e_out_photon_equi.shape = (30, 20)

    def _read_mtrp(self):
        """Get the list of reaction MTs for photon-producing reactions for this
        cross-section table. The MT values are somewhat arbitrary.
        """
        LMT = self.JXS[13]
        NMT = self.NXS[6]
        mts = np.asarray(self.XSS[LMT:LMT+NMT], dtype=int)
        rxs = [(mt, Reaction(mt, self)) for mt in mts]
        self.photon_reactions.update(rxs)

    def _read_lsigp(self):
        """Determine location of cross sections for each photon-producing reaction
        MT.
        """
        LXS = self.JXS[14]
        NMT = self.NXS[6]
        loca = np.asarray(self.XSS[LXS:LXS+NMT], dtype=int)
        for loc, rxn in zip(loca, self.photon_reactions.values()):
            rxn.LOCA = loc

    def _read_sigp(self):
        """Read cross-sections for each photon-producing reaction MT.
        """
        JXS15 = self.JXS[15]
        for rxn in self.photon_reactions.values():
            ind = JXS15 + rxn.LOCA - 1
            MFTYPE = int(self.XSS[ind])
            ind += 1

            if MFTYPE == 12 or MFTYPE == 16:
                # Yield data taken from ENDF File 12 or 6
                MTMULT = int(self.XSS[ind])
                ind += 1
    
                # ENDF interpolation parameters
                NR = int(self.XSS[ind])
                dat = np.asarray(self.XSS[ind+1:ind+1+2*NR], dtype=int)
                dat.shape = (2, NR)
                NBT, INT = dat
                ind += 1 + 2*NR

                # Energy-dependent yield
                NE = int(self.XSS[ind])
                dat = self.XSS[ind+1:ind+1+2*NE]
                dat.shape = (2, NE)
                rxn.e_yield, rxn.photon_yield = dat
                ind += 1 + 2*NE
            elif MFTYPE == 13:
                # Cross-section data from ENDF File 13
                # Energy grid index at which data starts
                rxn.IE = int(self.XSS[ind])

                # Cross sections
                NE = int(self.XSS[ind+1])
                self.sigma = self.XSS[ind+2:ind+2+NE]
                ind += 2 + NE
            else:
                raise ValueError("MFTYPE must be 12, 13, 16. Got {}".format(MFTYPE))

    def _read_landp(self):
        """Determine location of angular distribution for each photon-producing
        reaction MT.
        """
        JXS16 = self.JXS[16]
        NMT = self.NXS[6]
        locb = np.asarray(self.XSS[JXS16:JXS16+NMT], dtype=int)
        for loc, rxn in zip(locb, self.photon_reactions.values()):
            rxn.LOCB = loc

    def _read_andp(self):
        """Find the angular distribution for each photon-producing reaction
        MT."""
        JXS17 = self.JXS[17]
        for i, rxn in enumerate(self.photon_reactions.values()):
            if rxn.LOCB == 0:
                # No angular distribution data are given for this reaction,
                # isotropic scattering is asssumed in LAB
                continue

            ind = JXS17 + rxn.LOCB - 1

            # Number of energies and incoming energy grid
            NE = int(self.XSS[ind])
            self.a_dist_energy_in = self.XSS[ind+1:ind+1+NE]
            ind += 1 + NE

            # Location of tables associated with each outgoing angle
            # distribution
            LC = np.asarray(self.XSS[ind:ind+NE], dtype=int)

            # 32 equiprobable cosine bins for each incoming energy
            a_dist_mu_out = {}
            for j, location in enumerate(LC):
                if location == 0:
                    continue
                ind = JXS17 + location - 1
                a_dist_mu_out[j] = self.XSS[ind:ind+33]
            self.a_dist_mu_out = a_dist_mu_out

    def _read_yp(self):
        """Read list of reactions required as photon production yield
        multipliers.
        """
        if self.NXS[6] != 0:
            ind = self.JXS[20]
            NYP = int(self.XSS[ind])
            if NYP > 0:
                dat = np.asarray(self.XSS[ind+1:ind+1+NYP], dtype=int)
                self.MT_for_photon_yield = dat

    def _read_fis(self):
        """Read total fission cross-section data if present. Generally,
        this table is not provided since it is redundant.
        """
        # Check if fission block is present
        ind = self.JXS[21]
        if ind == 0:
            return

        # Read fission cross sections
        self.IE_fission = int(self.XSS[ind])  # Energy grid index
        NE = int(self.XSS[ind+1])
        self.sigma_f = self.XSS[ind+2:ind+2+NE]

    def _read_unr(self):
        """Read the unresolved resonance range probability tables if present.
        """
        # Check if URR probability tables are present
        ind = self.JXS[23]
        if ind == 0:
            return

        N = int(self.XSS[ind])     # Number of incident energies
        M = int(self.XSS[ind+1])   # Length of probability table
        INT = int(self.XSS[ind+2]) # Interpolation parameter (2=lin-lin, 5=log-log)
        ILF = int(self.XSS[ind+3]) # Inelastic competition flag
        IOA = int(self.XSS[ind+4]) # Other absorption flag
        IFF = int(self.XSS[ind+5]) # Factors flag
        ind += 6

        self.urr_energy = self.XSS[ind:ind+N] # Incident energies
        ind += N

        # Set up URR probability table
        urr_table = self.XSS[ind:ind+N*6*M]
        urr_table.shape = (N, 6, M)
        urr_table.strides = (8, 8*M*N, 8*N)  # Needed to get in N, 6, M order
        self.urr_table = urr_table

    def _get_float(self, n_values=1):
        if n_values > 1:
            ind = self.index
            values = self.XSS[ind:ind+n_values]
            self.index += n_values
            return values
        else:
            value = self.XSS[self.index]
            self.index += 1
            return value
            
    def _get_int(self, n_values=1):
        if n_values > 1:
            ind = self.index
            #values = [int(i) for i in self.XSS[ind:ind+n_values]]
            values = np.asarray(self.XSS[ind:ind+n_values], dtype=int)
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
    """A SabTable object contains thermal scattering data as represented by
    an S(alpha, beta) table.

    :Attributes:
      **awr** : float
        Atomic weight ratio of the target nuclide.

      **elastic_e_in** : list of floats
        Incoming energies in MeV for which the elastic cross section is
        tabulated.

      **elastic_P** : list of floats
        Elastic scattering cross section for data derived in the incoherent
        approximation, or Bragg edge parameters for data derived in the coherent
        approximation.

      **elastic_type** : str
        Describes the behavior of the elastic cross section, i.e. whether it was
        derived in the incoherent or coherent approximation.

      **inelastic_e_in** : list of floats
        Incoming energies in MeV for which the inelastic cross section is
        tabulated.

      **inelastic_sigma** : list of floats
        Inelastic scattering cross section in barns at each energy.

      **name** : str
        ZAID identifier of the table, e.g. 92235.70c.

      **temp** : float
        Temperature of the target nuclide in eV.

    Parameters
    ----------
    name : str
        ZAID identifier of the table, e.g. lwtr.10t.
    awr : float
        Atomic weight ratio of the target nuclide.
    temp : float
        Temperature of the target nuclide in eV.

    """
    

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
    """A Reaction object represents a single reaction channel for a nuclide with
    an associated cross section and, if present, a secondary angle and energy
    distribution. These objects are stored within the ``reactions`` attribute on
    subclasses of AceTable, e.g. NeutronTable.

    Parameters
    ----------
    MT : int
        The ENDF MT number for this reaction. On occasion, MCNP uses MT numbers
        that don't correspond exactly to the ENDF specification.
    table : AceTable
        The ACE table which contains this reaction. This is useful if data on
        the parent nuclide is needed (for instance, the energy grid at which
        cross sections are tabulated)

    :Attributes:
      **ang_energy_in** : list of floats
        Incoming energies in MeV at which angular distributions are tabulated.

      **ang_energy_cos** : list of floats
        Scattering cosines corresponding to each point of the angular distribution
        functions.

      **ang_energy_pdf** : list of floats
        Probability distribution function for angular distribution.

      **ang_energy_cdf** : list of floats
        Cumulative distribution function for angular distribution.

      **e_dist_energy** : list of floats
        Incoming energies in MeV at which energy distributions are tabulated.

      **e_dist_law** : int
        ACE law used for secondary energy distribution.

      **IE** : int
        The index on the energy grid corresponding to the threshold of this
        reaction.

      **MT** : int
        The ENDF MT number for this reaction. On occasion, MCNP uses MT numbers
        that don't correspond exactly to the ENDF specification.

      **Q** : float
        The Q-value of this reaction in MeV.

      **sigma** : list of floats
        Microscopic cross section for this reaction at each point on the energy
        grid above the threshold value.

      **TY** : int
        An integer whose absolute value is the number of neutrons emitted in
        this reaction. If negative, it indicates that scattering should be
        performed in the center-of-mass system. If positive, scattering should
        be preformed in the laboratory system.

    """

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
        return self.table.energy[self.IE]

    def __repr__(self):
        name = reaction_names.get(self.MT, None)
        if name is not None:
            rep = "<ACE Reaction: MT={0} {1}>".format(self.MT, name)
        else:
            rep = "<ACE Reaction: Unknown MT={0}>".format(self.MT)
        return rep


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
