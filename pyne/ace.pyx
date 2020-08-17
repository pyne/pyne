"""This module is for reading ACE-format cross sections. ACE stands for "A
Compact ENDF" format and originated from work on MCNP_. It is used in a number
of other Monte Carlo particle transport codes.

ACE-format cross sections are typically generated from ENDF_ files through a
cross section processing program like NJOY_. The ENDF data consists of tabulated
thermal data, ENDF/B resonance parameters, distribution parameters in the
unresolved resonance region, and tabulated data in the fast region. After the
ENDF data has been reconstructed and Doppler-broadened, the ACER module
generates ACE-format cross sections.

.. _MCNP: https://laws.lanl.gov/vhosts/mcnp.lanl.gov/
.. _NJOY: http://t2.lanl.gov/nis/codes.shtml
.. _ENDF: http://www.nndc.bnl.gov/endf

.. moduleauthor:: Paul Romano <paul.k.romano@gmail.com>, Anthony Scopatz <scopatz@gmail.com>
"""

from __future__ import division, unicode_literals
import io
import struct
from warnings import warn
from pyne.utils import QA_warn

try:
    from collections.abc import OrderedDict
except ImportError:
    from collections import OrderedDict

cimport numpy as np
import numpy as np
from bisect import bisect_right

from pyne cimport nucname
from pyne import nucname
from pyne.rxname import label

# fromstring func should depend on numpy verison
from pyne._utils import fromstring_split, fromstring_token
cdef bint NP_LE_V15 = int(np.__version__.split('.')[1]) <= 5 and np.__version__.startswith('1')

QA_warn(__name__)

def ascii_to_binary(ascii_file, binary_file):
    """Convert an ACE file in ASCII format (type 1) to binary format (type 2).

    Parameters
    ----------
    ascii_file : str
        Filename of ASCII ACE file
    binary_file : str
        Filename of binary ACE file to be written

    """

    # Open ASCII file
    ascii = open(ascii_file, 'r')

    # Set default record length
    record_length = 4096

    # Read data from ASCII file
    lines = ascii.readlines()
    ascii.close()

    # Open binary file
    binary = open(binary_file, 'wb')

    idx = 0

    while idx < len(lines):
        #check if it's a > 2.0.0 version header
        if lines[idx].split()[0][1] == '.':
            if lines[idx + 1].split()[3] == '3':
                idx = idx + 3
            else:
                raise NotImplementedError('Only backwards compatible ACE'
                                          'headers currently supported')
        # Read/write header block
        hz = lines[idx][:10].encode('UTF-8')
        aw0 = float(lines[idx][10:22])
        tz = float(lines[idx][22:34])
        hd = lines[idx][35:45].encode('UTF-8')
        hk = lines[idx + 1][:70].encode('UTF-8')
        hm = lines[idx + 1][70:80].encode('UTF-8')
        binary.write(struct.pack(str('=10sdd10s70s10s'), hz, aw0, tz, hd, hk, hm))

        # Read/write IZ/AW pairs
        data = ' '.join(lines[idx + 2:idx + 6]).split()
        iz = list(map(int, data[::2]))
        aw = list(map(float, data[1::2]))
        izaw = [item for sublist in zip(iz, aw) for item in sublist]
        binary.write(struct.pack(str('=' + 16*'id'), *izaw))

        # Read/write NXS and JXS arrays. Null bytes are added at the end so
        # that XSS will start at the second record
        nxs = list(map(int, ' '.join(lines[idx + 6:idx + 8]).split()))
        jxs = list(map(int, ' '.join(lines[idx + 8:idx + 12]).split()))
        binary.write(struct.pack(str('=16i32i{0}x'.format(record_length - 500)),
                                 *(nxs + jxs)))

        # Read/write XSS array. Null bytes are added to form a complete record
        # at the end of the file
        n_lines = (nxs[0] + 3)//4
        xss = list(map(float, ' '.join(lines[
            idx + 12:idx + 12 + n_lines]).split()))
        extra_bytes = record_length - ((len(xss)*8 - 1) % record_length + 1)
        binary.write(struct.pack(str('={0}d{1}x'.format(nxs[0], extra_bytes)),
                                 *xss))

        # Advance to next table in file
        idx += 12 + n_lines

    # Close binary file
    binary.close()


class Library(object):
    """
    A Library objects represents an ACE-formatted file which may contain
    multiple tables with data.

    Parameters
    ----------
    filename : str
        Path of the ACE library file to load.

    Attributes
    ----------
    binary : bool
        Identifies Whether the library is in binary format or not
    tables : dict
        Dictionary whose keys are the names of the ACE tables and whose
        values are the instances of subclasses of AceTable (e.g. NeutronTable)
    verbose : bool
        Determines whether output is printed to the stdout when reading a
        Library
    """

    def __init__(self, filename):
        # Determine whether file is ASCII or binary
        self.f = None
        try:
            self.f = io.open(filename, 'rb')
            # Grab 10 lines of the library
            sb = b''.join([self.f.readline() for i in range(10)])

            # Try to decode it with ascii
            sd = sb.decode('ascii')

            # No exception so proceed with ASCII - reopen in non-binary
            self.f.close()
            self.f = io.open(filename, 'r')
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
        """read(table_names=None)

        Read through and parse the ACE-format library.

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
            if len(self.f.read(1)) == 0:
                return
            self.f.seek(start_position)

            # Read name, atomic mass ratio, temperature, date, comment, and
            # material
            name, awr, temp, date, comment, mat = \
                struct.unpack(str('=10sdd10s70s10s'), self.f.read(116))
            name = name.strip()

            # Read ZAID/awr combinations
            data = struct.unpack(str('=' + 16*'id'), self.f.read(192))

            # Read NXS
            nxs = list(struct.unpack(str('=16i'), self.f.read(64)))

            # Determine length of XSS and number of records
            length = nxs[0]
            n_records = (length + entries - 1)//entries

            # name is bytes, make it a string
            name = name.decode()
            # verify that we are supposed to read this table in
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

            if self.verbose:
                temp_in_K = round(temp * 1e6 / 8.617342e-5)
                print("Loading nuclide {0} at {1} K".format(name, temp_in_K))
            self.tables[name] = table

            # Read JXS
            table.jxs = list(struct.unpack(str('=32i'), self.f.read(128)))

            # Read XSS
            self.f.seek(start_position + recl_length)
            table.xss = list(struct.unpack(str('={0}d'.format(length)),
                                           self.f.read(length*8)))

            # Insert empty object at beginning of NXS, JXS, and XSS
            # arrays so that the indexing will be the same as
            # Fortran. This makes it easier to follow the ACE format
            # specification.
            table.nxs = nxs
            table.nxs.insert(0, 0)
            table.nxs = np.array(table.nxs, dtype=int)

            table.jxs.insert(0, 0)
            table.jxs = np.array(table.jxs, dtype=int)

            table.xss.insert(0, 0.0)
            table.xss = np.array(table.xss, dtype=float)

            # Read all data blocks
            table._read_all()

            # Advance to next record
            self.f.seek(start_position + recl_length*(n_records + 1))

    def _read_ascii(self, table_names):
        cdef list lines, rawdata

        f = self.f
        tables_seen = set()

        cdef int i
        lines = [f.readline() for i in range(13)]

        while (0 != len(lines)) and (lines[0] != ''):
            # Read name of table, atomic mass ratio, and temperature. If first
            # line is empty, we are at end of file

            # check if it's a 2.0 style header
            if lines[0].split()[0][1] == '.':
                words = lines[0].split()
                version = words[0]
                name = words[1]
                if len(words) == 3:
                    source = words[2]
                words = lines[1].split()
                awr = float(words[0])
                temp = float(words[1])
                commentlines = int(words[3])
                for i in range(commentlines):
                    lines.pop(0)
                    lines.append(f.readline())
            else:
                words = lines[0].split()
                name = words[0]
                awr = float(words[1])
                temp = float(words[2])

            datastr = '0 ' + ' '.join(lines[6:8])
            nxs = fromstring_split(datastr, dtype=int)

            n_lines = (nxs[1] + 3)//4
            n_bytes = len(lines[-1]) * (n_lines - 2) + 1

            # Ensure that we have more tables to read in
            if (table_names is not None) and (table_names < tables_seen):
                break
            tables_seen.add(name)

            # verify that we are suppossed to read this table in
            if (table_names is not None) and (name not in table_names):
                cur = f.tell()
                f.seek(cur + n_bytes)
                f.readline()
                lines = [f.readline() for i in range(13)]
                continue

            # ensure we have a valid table type
            if 0 == len(name) or name[-1] not in table_types:
                warn("Unsupported table: " + name, RuntimeWarning)
                cur = f.tell()
                f.seek(cur + n_bytes)
                f.readline()
                lines = [f.readline() for i in range(13)]
                continue

            # read and and fix over-shoot
            lines += f.readlines(n_bytes)
            if 12+n_lines < len(lines):
                goback = sum([len(line) for line in lines[12+n_lines:]])
                lines = lines[:12+n_lines]
                cur = f.tell()
                f.seek(cur - goback)

            # get the table
            table = table_types[name[-1]](name, awr, temp)

            if self.verbose:
                temp_in_K = round(temp * 1e6 / 8.617342e-5)
                print("Loading nuclide {0} at {1} K".format(name, temp_in_K))
            self.tables[name] = table

            # Read comment
            table.comment = lines[1].strip()

            # Add NXS, JXS, and XSS arrays to table
            # Insert empty object at beginning of NXS, JXS, and XSS
            # arrays so that the indexing will be the same as
            # Fortran. This makes it easier to follow the ACE format
            # specification.
            table.nxs = nxs

            datastr = '0 ' + ' '.join(lines[8:12])
            table.jxs = fromstring_split(datastr, dtype=int)

            datastr = '0.0 ' + ''.join(lines[12:12+n_lines])
            if NP_LE_V15:
                #table.xss = np.fromstring(datastr, sep=" ")
                table.xss = fromstring_split(datastr, dtype=float)
            else:
                table.xss = fromstring_token(datastr, inplace=True, maxsize=4*n_lines+1)

            # Read all data blocks
            table._read_all()
            lines = [f.readline() for i in range(13)]

        f.seek(0)

    def find_table(self, name):
        """find_table(name)

        Returns a cross-section table with a given name.

        Parameters
        ----------
        name : str
            Name of the cross-section table, e.g. 92235.70c

        """
        return self.tables.get(name, None)

    def __del__(self):
        if self.f is not None:
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
        Atomic mass ratio of the target nuclide.
    temp : float
        Temperature of the target nuclide in eV.

    Attributes
    ----------
    awr : float
        Atomic mass ratio of the target nuclide.

    energy : list of floats
        The energy values (MeV) at which reaction cross-sections are tabulated.

    name : str
        ZAID identifier of the table, e.g. 92235.70c.

    nu_p_energy : list of floats
        Energies in MeV at which the number of prompt neutrons emitted per
        fission is tabulated.

    nu_p_type : str
        Indicates how number of prompt neutrons emitted per fission is
        stored. Can be either "polynomial" or "tabular".

    nu_p_value : list of floats
        The number of prompt neutrons emitted per fission, if data is stored in
        "tabular" form, or the polynomial coefficients for the "polynomial"
        form.

    nu_t_energy : list of floats
        Energies in MeV at which the total number of neutrons emitted per
        fission is tabulated.

    nu_t_type : str
        Indicates how total number of neutrons emitted per fission is
        stored. Can be either "polynomial" or "tabular".

    nu_t_value : list of floats
        The total number of neutrons emitted per fission, if data is stored in
        "tabular" form, or the polynomial coefficients for the "polynomial"
        form.

    reactions : list of Reactions
        A list of Reaction instances containing the cross sections, secondary
        angle and energy distributions, and other associated data for each
        reaction for this nuclide.

    sigma_a : list of floats
        The microscopic absorption cross section for each value on the energy
        grid.

    sigma_t : list of floats
        The microscopic total cross section for each value on the energy grid.

    temp : float
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
        self._read_cross_sections()
        self._read_nu()
        self._read_angular_distributions()
        self._read_energy_distributions()
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

    def _read_cross_sections(self):
        """Reads and parses the ESZ, MTR, LQR, TRY, LSIG, and SIG blocks. These
        blocks contain the energy grid, all reaction cross sections, the total
        cross section, average heating numbers, and a list of reactions with
        their Q-values and multiplicites.
        """

        cdef int n_energies, n_reactions, loc

        # Determine number of energies on nuclide grid and number of reactions
        # excluding elastic scattering
        n_energies = self.nxs[3]
        n_reactions = self.nxs[4]

        # Read energy grid and total, absorption, elastic scattering, and
        # heating cross sections -- note that this appear separate from the rest
        # of the reaction cross sections
        arr = self.xss[self.jxs[1]:self.jxs[1] + 5*n_energies]
        arr.shape = (5, n_energies)
        self.energy, self.sigma_t, self.sigma_a, sigma_el, self.heating = arr

        # Create elastic scattering reaction
        elastic_scatter = Reaction(2, self)
        elastic_scatter.Q = 0.0
        elastic_scatter.IE = 0
        elastic_scatter.multiplicity = 1
        elastic_scatter.sigma = sigma_el
        self.reactions[2] = elastic_scatter

        # Create all other reactions with MT values
        mts = np.asarray(self.xss[self.jxs[3]:self.jxs[3] + n_reactions], dtype=int)
        qvalues = np.asarray(self.xss[self.jxs[4]:self.jxs[4] +
                                      n_reactions], dtype=float)
        tys = np.asarray(self.xss[self.jxs[5]:self.jxs[5] + n_reactions], dtype=int)

        # Create all reactions other than elastic scatter
        reactions = [(mt, Reaction(mt, self)) for mt in mts]
        self.reactions.update(reactions)

        # Loop over all reactions other than elastic scattering
        for i, reaction in enumerate(list(self.reactions.values())[1:]):
            # Copy Q values and multiplicities and determine if scattering
            # should be treated in the center-of-mass or lab system
            reaction.Q = qvalues[i]
            reaction.multiplicity = abs(tys[i])
            reaction.center_of_mass = (tys[i] < 0)

            # Get locator for cross-section data
            loc = int(self.xss[self.jxs[6] + i])

            # Determine starting index on energy grid
            reaction.IE = int(self.xss[self.jxs[7] + loc - 1]) - 1

            # Determine number of energies in reaction
            n_energies = int(self.xss[self.jxs[7] + loc])

            # Read reaction cross section
            reaction.sigma = self.xss[self.jxs[7] + loc + 1:
                                          self.jxs[7] + loc + 1 + n_energies]

    def _read_nu(self):
        """Read the NU block -- this contains information on the prompt
        and delayed neutron precursor yields, decay constants, etc
        """
        cdef int ind, i, jxs2, KNU, LNU, NR, NE, NC

        jxs2 = self.jxs[2]

        # No NU block
        if jxs2 == 0:
            return

        # Either prompt nu or total nu is given
        if self.xss[jxs2] > 0:
            KNU = jxs2
            LNU = int(self.xss[KNU])

            # Polynomial function form of nu
            if LNU == 1:
                self.nu_t_type = "polynomial"
                NC = int(self.xss[KNU+1])
                coeffs = self.xss[KNU+2 : KNU+2+NC]

            # Tabular data form of nu
            elif LNU == 2:
                self.nu_t_type = "tabular"
                NR = int(self.xss[KNU+1])
                if NR > 0:
                    interp_NBT = self.xss[KNU+2    : KNU+2+NR  ]
                    interp_INT = self.xss[KNU+2+NR : KNU+2+2*NR]
                NE = int(self.xss[KNU+2+2*NR])
                self.nu_t_energy = self.xss[KNU+3+2*NR    : KNU+3+2*NR+NE  ]
                self.nu_t_value  = self.xss[KNU+3+2*NR+NE : KNU+3+2*NR+2*NE]
        # Both prompt nu and total nu
        elif self.xss[jxs2] < 0:
            KNU = jxs2 + 1
            LNU = int(self.xss[KNU])

            # Polynomial function form of nu
            if LNU == 1:
                self.nu_p_type = "polynomial"
                NC = int(self.xss[KNU+1])
                coeffs = self.xss[KNU+2 : KNU+2+NC]

            # Tabular data form of nu
            elif LNU == 2:
                self.nu_p_type = "tabular"
                NR = int(self.xss[KNU+1])
                if NR > 0:
                    interp_NBT = self.xss[KNU+2    : KNU+2+NR  ]
                    interp_INT = self.xss[KNU+2+NR : KNU+2+2*NR]
                NE = int(self.xss[KNU+2+2*NR])
                self.nu_p_energy = self.xss[KNU+3+2*NR    : KNU+3+2*NR+NE  ]
                self.nu_p_value  = self.xss[KNU+3+2*NR+NE : KNU+3+2*NR+2*NE]

            KNU = jxs2 + int(abs(self.xss[jxs2])) + 1
            LNU = int(self.xss[KNU])

            # Polynomial function form of nu
            if LNU == 1:
                self.nu_t_type = "polynomial"
                NC = int(self.xss[KNU+1])
                coeffs = self.xss[KNU+2 : KNU+2+NC]

            # Tabular data form of nu
            elif LNU == 2:
                self.nu_t_type = "tabular"
                NR = int(self.xss[KNU+1])
                if NR > 0:
                    interp_NBT = self.xss[KNU+2    : KNU+2+NR  ]
                    interp_INT = self.xss[KNU+2+NR : KNU+2+2*NR]
                NE = int(self.xss[KNU+2+2*NR])
                self.nu_t_energy = self.xss[KNU+3+2*NR    : KNU+3+2*NR+NE  ]
                self.nu_t_value  = self.xss[KNU+3+2*NR+NE : KNU+3+2*NR+2*NE]

        # Check for delayed nu data
        if self.jxs[24] > 0:
            KNU = self.jxs[24]
            NR = int(self.xss[KNU+1])
            if NR > 0:
                interp_NBT = self.xss[KNU+2    : KNU+2+NR  ]
                interp_INT = self.xss[KNU+2+NR : KNU+2+2*NR]
            NE = int(self.xss[KNU+2+2*NR])
            self.nu_d_energy = self.xss[KNU+3+2*NR    : KNU+3+2*NR+NE  ]
            self.nu_d_value  = self.xss[KNU+3+2*NR+NE : KNU+3+2*NR+2*NE]

            # Delayed neutron precursor distribution
            self.nu_d_precursor_const = {}
            self.nu_d_precursor_energy = {}
            self.nu_d_precursor_prob = {}
            i = self.jxs[25]
            n_group = self.nxs[8]
            for group in range(n_group):
                self.nu_d_precursor_const[group] = self.xss[i]
                NR = int(self.xss[i+1])
                if NR > 0:
                    interp_NBT = self.xss[i+2    : i+2+NR]
                    interp_INT = self.xss[i+2+NR : i+2+2*NR]
                NE = int(self.xss[i+2+2*NR])
                self.nu_d_precursor_energy[group] = self.xss[i+3+2*NR    : i+3+2*NR+NE  ]
                self.nu_d_precursor_prob[group]   = self.xss[i+3+2*NR+NE : i+3+2*NR+2*NE]
                i = i+3+2*NR+2*NE

            # Energy distribution for delayed fission neutrons
            LED = self.jxs[26]
            self.nu_d_energy_dist = []
            for group in range(n_group):
                location_start = self.xss[LED + group]
                energy_dist = self._get_energy_distribution(
                    location_start, delayed_n=True)
                self.nu_d_energy_dist.append(energy_dist)


    def _read_angular_distributions(self):
        """Find the angular distribution for each reaction MT
        """
        cdef int ind, i, j, n_reactions, n_energies, n_bins
        cdef dict ang_cos, ang_pdf, ang_cdf

        # Number of reactions with secondary neutrons (including elastic
        # scattering)
        n_reactions = self.nxs[5] + 1

        # Angular distribution for all reactions with secondary neutrons
        for i, reaction in enumerate(list(self.reactions.values())[:n_reactions]):
            loc = int(self.xss[self.jxs[8] + i])

            # Check if angular distribution data exist
            if loc == -1:
                # Angular distribution data are specified through LAWi
                # = 44 in the DLW block
                continue
            elif loc == 0:
                # No angular distribution data are given for this
                # reaction, isotropic scattering is asssumed (in CM if
                # TY < 0 and in LAB if TY > 0)
                continue

            ind = self.jxs[9] + loc

            # Number of energies at which angular distributions are tabulated
            n_energies = int(self.xss[ind - 1])

            # Incoming energy grid
            reaction.ang_energy_in = self.xss[ind:ind + n_energies]
            ind += n_energies

            # Read locations for angular distributions
            locations = np.asarray(self.xss[ind:ind + n_energies], dtype=int)
            ind += n_energies

            ang_cos = {}
            ang_pdf = {}
            ang_cdf = {}
            for j, location in enumerate(locations):
                if location > 0:
                    # Equiprobable 32 bin distribution
                    # print([reaction,'equiprobable'])
                    ang_cos[j] = self.xss[ind:ind + 33]
                    ind += 33
                elif location < 0:
                    # Tabular angular distribution
                    JJ = int(self.xss[ind])
                    n_bins = int(self.xss[ind + 1])
                    ind += 2
                    ang_dat = self.xss[ind:ind + 3*n_bins]
                    ang_dat.shape = (3, n_bins)
                    ang_cos[j], ang_pdf[j], ang_cdf[j] = ang_dat
                    ind += 3 * n_bins
                else:
                    # Isotropic angular distribution
                    ang_cos = np.array([-1., 0., 1.])
                    ang_pdf = np.array([0.5, 0.5, 0.5])
                    ang_cdf = np.array([0., 0.5, 1.])

            reaction.ang_cos = ang_cos
            reaction.ang_pdf = ang_pdf
            reaction.ang_cdf = ang_cdf

    def _read_energy_distributions(self):
        """Determine the energy distribution for secondary neutrons for
        each reaction MT
        """
        cdef int i

        # Number of reactions with secondary neutrons other than elastic
        # scattering. For elastic scattering, the outgoing energy can be
        # determined from kinematics.
        n_reactions = self.nxs[5]

        for i, reaction in enumerate(list(self.reactions.values())[1:n_reactions + 1]):
            # Determine locator for ith energy distribution
            location_start = int(self.xss[self.jxs[10] + i])

            # Read energy distribution data
            reaction.energy_dist = self._get_energy_distribution(location_start)

    def _get_energy_distribution(self, location_start, delayed_n=False):
        """Returns an EnergyDistribution object from data read in starting at
        location_start.
        """

        cdef int ind, i, n_reactions, NE, n_regions, location_next_law, law, location_data, NPE, NPA

        # Create EnergyDistribution object
        edist = EnergyDistribution()

        # Determine location of energy distribution
        if delayed_n:
            location_dist = self.jxs[27]
        else:
            location_dist = self.jxs[11]

        # Set starting index for energy distribution
        ind = location_dist + location_start - 1

        location_next_law = int(self.xss[ind])
        law = int(self.xss[ind+1])
        location_data = int(self.xss[ind+2])

        # Number of interpolation regions for law applicability regime
        n_regions = int(self.xss[ind+3])
        ind += 4
        if n_regions > 0:
            dat = np.asarray(self.xss[ind:ind + 2*n_regions], dtype=int)
            dat.shape = (2, n_regions)
            interp_NBT, interp_INT = dat
            ind += 2 * n_regions

        # Determine tabular energy points and probability of law
        # validity
        NE = int(self.xss[ind])
        dat = self.xss[ind+1:ind+1+2*NE]
        dat.shape = (2, NE)
        edist.energy, edist.pvalid = dat

        edist.law = law
        ind = location_dist + location_data - 1

        if law == 1:
            # Tabular equiprobable energy bins (ENDF Law 1)
            n_regions = int(self.xss[ind])
            ind += 1
            if n_regions > 0:
                dat = np.asarray(self.xss[ind:ind+2*n_regions], dtype=int)
                dat.shape = (2, n_regions)
                edist.NBT, edist.INT = dat
                ind += 2 * n_regions

            # Number of outgoing energies in each E_out table
            NE = int(self.xss[ind])
            edist.energy_in = self.xss[ind+1:ind+1+NE]
            ind += 1 + NE

            # Read E_out tables
            NET = int(self.xss[ind])
            dat = self.xss[ind+1:ind+1+3*NET]
            dat.shape = (3, NET)
            self.e_dist_energy_out1, self.e_dist_energy_out2, \
                                     self.e_dist_energy_outNE = dat
            ind += 1 + 3 * NET
        elif law == 2:
            # Discrete photon energy
            self.e_dist_LP = int(self.xss[ind])
            self.e_dist_EG = self.xss[ind+1]
            ind += 2
        elif law == 3:
            # Level scattering (ENDF Law 3)
            edist.data = self.xss[ind:ind+2]
            ind += 2
        elif law == 4:
            # Continuous tabular distribution (ENDF Law 1)
            n_regions = int(self.xss[ind])
            ind += 1
            if n_regions > 0:
                dat = np.asarray(self.xss[ind:ind+2*n_regions], dtype=int)
                dat.shape = (2, n_regions)
                edist.NBT, edist.INT = dat
                ind += 2 * n_regions

            # Number of outgoing energies in each E_out table
            NE = int(self.xss[ind])
            edist.energy_in = self.xss[ind+1:ind+1+NE]
            L = self.xss[ind+1+NE:ind+1+2*NE]
            ind += 1 + 2*NE

            nps = []
            edist.intt = []        # Interpolation scheme (1=hist, 2=lin-lin)
            edist.energy_out = []  # Outgoing E grid for each incoming E
            edist.pdf = []         # Probability dist for " " "
            edist.cdf = []         # Cumulative dist for " " "
            for i in range(NE):
                INTTp = int(self.xss[ind])
                if INTTp > 10:
                    INTT = INTTp % 10
                    ND = (INTTp - INTT)/10
                else:
                    INTT = INTTp
                    ND = 0
                edist.intt.append(INTT)
                #if ND > 0:
                #    print [reaction, ND, INTT]

                NP = int(self.xss[ind+1])
                nps.append(NP)
                dat = self.xss[ind+2:ind+2+3*NP]
                dat.shape = (3, NP)
                edist.energy_out.append(dat[0])
                edist.pdf.append(dat[1])
                edist.cdf.append(dat[2])
                ind += 2 + 3*NP

            # convert to arrays if possible
            edist.intt = np.array(edist.intt)
            nps = np.array(nps)
            if all((nps[1:] - nps[:-1]) == 0):
                edist.energy_out = np.array(edist.energy_out)
                edist.pdf = np.array(edist.pdf)
                edist.cdf = np.array(edist.cdf)
        elif law == 5:
            # General evaporation spectrum (ENDF-5 File 5 LF=5)
            n_regions = int(self.xss[ind])
            ind += 1
            if n_regions > 0:
                dat = np.asarray(self.xss[ind:ind+2*n_regions], dtype=int)
                dat.shape = (2, n_regions)
                edist.NBT, edist.INT = dat
                ind += 2 * n_regions

            NE = int(self.xss[ind])
            edist.energy_in = self.xss[ind+1:ind+1+NE]
            edist.T = self.xss[ind+1+NE:ind+1+2*NE]
            ind += 1+ 2*NE

            NET = int(self.xss[ind])
            edist.X = self.xss[ind+1:ind+1+NET]
            ind += 1 + NET
        elif law == 7:
            # Simple Maxwell fission spectrum (ENDF-6 File 5 LF=7)
            n_regions = int(self.xss[ind])
            ind += 1
            if n_regions > 0:
                dat = np.asarray(self.xss[ind:ind+2*n_regions], dtype=int)
                dat.shape = (2, n_regions)
                edist.NBT, edist.INT = dat
                ind += 2 * n_regions

            NE = int(self.xss[ind])
            edist.energy_in = self.xss[ind+1:ind+1+NE]
            edist.T = self.xss[ind+1+NE:ind+1+2*NE]
            edist.U = self.xss[ind+1+2*NE]
            ind += 2 + 2*NE
        elif law == 9:
            # Evaporation spectrum (ENDF-6 File 5 LF=9)
            n_regions = int(self.xss[ind])
            ind += 1
            if n_regions > 0:
                dat = np.asarray(self.xss[ind:ind+2*n_regions], dtype=int)
                dat.shape = (2, n_regions)
                edist.NBT, edist.INT = dat
                ind += 2 * n_regions

            NE = int(self.xss[ind])
            edist.energy_in = self.xss[ind+1:ind+1+NE]
            edist.T = self.xss[ind+1+NE:ind+1+2*NE]
            edist.U = self.xss[ind+1+2*NE]
            ind += 2 + 2*NE
        elif law == 11:
            # Energy dependent Watt spectrum (ENDF-6 File 5 LF=11)
            # Interpolation scheme between a's
            n_regions = int(self.xss[ind])
            ind += 1
            if n_regions > 0:
                dat = np.asarray(self.xss[ind:ind+2*n_regions], dtype=int)
                dat.shape = (2, n_regions)
                edist.NBTa, edist.INTa = dat
                ind += 2 * n_regions

            # Incident energy table and tabulated a's
            NE = int(self.xss[ind])
            edist.energya_in = self.xss[ind+1:ind+1+NE]
            edist.a = self.xss[ind+1+NE:ind+1+2*NE]
            ind += 1 + 2*NE

            # Interpolation scheme between b's
            n_regions = int(self.xss[ind])
            ind += 1
            if n_regions > 0:
                dat = np.asarray(self.xss[ind:ind+2*n_regions], dtype=int)
                dat.shape = (2, n_regions)
                edist.NBTb, edist.INTb = dat
                ind += 2 * n_regions

            # Incident energy table and tabulated b's
            NE = int(self.xss[ind])
            edist.energyb_in = self.xss[ind+1:ind+1+NE]
            edist.b = self.xss[ind+1+NE:ind+1+2*NE]

            edist.U = self.xss[ind+1+2*NE]
            ind += 2 + 2*NE
        elif law == 22:
            # Tabular linear functions (UK Law 2)
            # Interpolation scheme (not used in MCNP)
            n_regions = int(self.xss[ind])
            ind += 1
            if n_regions > 0:
                dat = np.asarray(self.xss[ind:ind+2*n_regions], dtype=int)
                dat.shape = (2, n_regions)
                edist.NBT, edist.INT = dat
                ind += 2 * n_regions

            # Number of incident energies
            NE = int(self.xss[ind])
            edist.energy_in = self.xss[ind+1:ind+1+NE]
            LOCE = np.asarray(self.xss[ind+1+NE:ind+1+2*NE], dtype=int)
            ind += 1 + 2*NE

            # Read linear functions
            nfs = []
            edist.P = []
            edist.T = []
            edist.C = []
            for i in range(NE):
                NF = int(self.xss[ind])
                nfs.append(NF)
                dat = self.xss[ind+1:ind+1+3*NF]
                dat.shape = (3, NF)
                edist.P.append(dat[0])
                edist.T.append(dat[1])
                edist.C.append(dat[2])
                ind += 1 + 3*NF

            # convert to arrays if possible
            nfs = np.array(nfs)
            if all((nfs[1:] - nfs[:-1]) == 0):
                edist.P = np.array(edist.P)
                edist.T = np.array(edist.T)
                edist.C = np.array(edist.C)
        elif law == 24:
            # From UK Law 6
            # Interpolation scheme (not used in MCNP)
            n_regions = int(self.xss[ind])
            ind += 1
            if n_regions > 0:
                dat = np.asarray(self.xss[ind:ind+2*n_regions], dtype=int)
                dat.shape = (2, n_regions)
                edist.NBT, edist.INT = dat
                ind += 2 * n_regions

            # Number of incident energies
            NE = int(self.xss[ind])
            edist.energy_in = self.xss[ind+1:ind+1+NE]
            ind += 1 + NE

            # Outgoing energy tables
            NET = int(self.xss[ind])
            edist.T = self.xss[ind+1:ind+1+NE*NET]
            edist.T.shape = (NE, NET)
            ind += 1 + NE*NET
        elif law == 44:
            # Kalbach-87 Formalism (ENDF File 6 Law 1, LANG=2)
            # Interpolation scheme
            n_regions = int(self.xss[ind])
            ind += 1
            if n_regions > 0:
                dat = np.asarray(self.xss[ind:ind+2*n_regions], dtype=int)
                dat.shape = (2, n_regions)
                edist.NBT, edist.INT = dat
                ind += 2 * n_regions

            # Number of outgoing energies in each E_out table
            NE = int(self.xss[ind])
            edist.energy_in = self.xss[ind+1:ind+1+NE]
            L = np.asarray(self.xss[ind+1+NE:ind+1+2*NE], dtype=int)
            ind += 1 + 2*NE

            nps = []
            edist.intt = []        # Interpolation scheme (1=hist, 2=lin-lin)
            edist.energy_out = []  # Outgoing E grid for each incoming E
            edist.pdf = []         # Probability dist for " " "
            edist.cdf = []         # Cumulative dist for " " "
            edist.frac = []        # Precompound fraction for " " "
            edist.ang = []         # Angular distribution slope for " " "
            for i in range(NE):
                INTTp = int(self.xss[ind])
                if INTTp > 10:
                    INTT = INTTp % 10
                    ND = (INTTp - INTT)/10
                else:
                    INTT = INTTp
                edist.intt.append(INTT)

                NP = int(self.xss[ind+1])
                nps.append(NP)
                ind += 2

                dat = self.xss[ind:ind+5*NP]
                dat.shape = (5, NP)
                edist.energy_out.append(dat[0])
                edist.pdf.append(dat[1])
                edist.cdf.append(dat[2])
                edist.frac.append(dat[3])
                edist.ang.append(dat[4])
                ind += 5*NP

            # convert to arrays if possible
            edist.intt = np.array(edist.intt)
            nps = np.array(nps)
            if all((nps[1:] - nps[:-1]) == 0):
                edist.energy_out = np.array(edist.energy_out)
                edist.pdf = np.array(edist.pdf)
                edist.cdf = np.array(edist.cdf)
        elif law == 61:
            # Like 44, but tabular distribution instead of Kalbach-87
            # Interpolation scheme
            n_regions = int(self.xss[ind])
            ind += 1
            if n_regions > 0:
                dat = np.asarray(self.xss[ind:ind+2*n_regions], dtype=int)
                dat.shape = (2, n_regions)
                edist.NBT, edist.INT = dat
                ind += 2 * n_regions

            # Number of outgoing energies in each E_out table
            NE = int(self.xss[ind])
            edist.energy_in = self.xss[ind+1:ind+1+NE]
            L = np.asarray(self.xss[ind+1+NE:ind+1+2*NE], dtype=int)
            ind += 1 + 2*NE

            npes = []
            edist.intt = []        # Interpolation scheme (1=hist, 2=lin-lin)
            edist.energy_out = []  # Outgoing E grid for each incoming E
            edist.pdf = []         # Probability dist for " " "
            edist.cdf = []         # Cumulative dist for " " "

            npas = []
            edist.a_dist_intt = []
            edist.a_dist_mu_out = [] # Cosine scattering angular grid
            edist.a_dist_pdf = []    # Probability dist function
            edist.a_dist_cdf = []
            for i in range(NE):
                INTTp = int(self.xss[ind])
                if INTTp > 10:
                    INTT = INTTp % 10
                    ND = (INTTp - INTT)/10
                else:
                    INTT = INTTp
                edist.intt.append(INTT)

                # Secondary energy distribution
                NPE = int(self.xss[ind+1])
                npes.append(NPE)
                dat = self.xss[ind+2:ind+2+4*NPE]
                dat.shape = (4, NPE)
                edist.energy_out.append(dat[0])
                edist.pdf.append(dat[1])
                edist.cdf.append(dat[2])
                LC = np.asarray(dat[3], dtype=int)
                ind += 2 + 4*NPE

                # Secondary angular distribution
                edist.a_dist_intt.append([])
                edist.a_dist_mu_out.append([])
                edist.a_dist_pdf.append([])
                edist.a_dist_cdf.append([])
                for j in range(NPE):
                    edist.a_dist_intt[-1].append(int(self.xss[ind]))
                    NPA = int(self.xss[ind+1])
                    npas.append(NPA)
                    dat = self.xss[ind+2:ind+2+3*NPA]
                    dat.shape = (3, NPA)
                    edist.a_dist_mu_out[-1].append(dat[0])
                    edist.a_dist_pdf[-1].append(dat[1])
                    edist.a_dist_cdf[-1].append(dat[2])
                    ind += 2 + 3*NPA

            # convert to arrays if possible
            edist.intt = np.array(edist.intt)
            npes = np.array(npes)
            npas = np.array(npas)
            if all((npes[1:] - npes[:-1]) == 0):
                edist.energy_out = np.array(edist.energy_out)
                edist.pdf = np.array(edist.pdf)
                edist.cdf = np.array(edist.cdf)

                edist.a_dist_intt = np.array(edist.a_dist_intt)
                if all((npas[1:] - npas[:-1]) == 0):
                    edist.a_dist_mu_out = np.array(edist.a_dist_mu_out)
                    edist.a_dist_pdf = np.array(edist.a_dist_pdf)
                    edist.a_dist_cdf = np.array(edist.a_dist_cdf)
        elif law == 66:
            # N-body phase space distribution (ENDF File 6 Law 6)
            edist.nbodies = int(self.xss[ind])
            edist.massratio = self.xss[ind+1]
            ind += 2
        elif law == 67:
            # Laboratory angle-energy law (ENDF File 6 Law 7)
            # Interpolation scheme
            n_regions = int(self.xss[ind])
            ind += 1
            if n_regions > 0:
                dat = np.asarray(self.xss[ind:ind+2*n_regions], dtype=int)
                dat.shape = (2, n_regions)
                edist.NBT, edist.INT = dat
                ind += 2 * n_regions

            # Number of outgoing energies in each E_out table
            NE = int(self.xss[ind])
            edist.energy_in = self.xss[ind+1:ind+1+NE]
            L = np.asarray(self.xss[ind+1+NE:ind+1+2*NE], dtype=int)
            ind += 1 + 2*NE

        # TODO: Read rest of data

        return edist


    def _read_gpd(self):
        """Read total photon production cross section.
        """
        cdef int ind, jxs12, NE

        jxs12 = self.jxs[12]
        if jxs12 != 0:
            # Determine number of energies
            NE = self.nxs[3]

            # Read total photon production cross section
            ind = jxs12
            self.sigma_photon = self.xss[ind:ind+NE]

            # The MCNP manual also specifies that this block contains secondary
            # photon energies based on a 30x20 matrix formulation. However, the
            # ENDF/B-VII.0 libraries distributed with MCNP as well as other
            # libraries do not contain this 30x20 matrix.

            # # The following energies are the discrete incident neutron energies
            # # for which the equiprobable secondary photon outgoing energies are
            # # given
            # self.e_in_photon_equi = np.array(
            #                         [1.39e-10, 1.52e-7, 4.14e-7, 1.13e-6, 3.06e-6,
            #                          8.32e-6,  2.26e-5, 6.14e-5, 1.67e-4, 4.54e-4,
            #                          1.235e-3, 3.35e-3, 9.23e-3, 2.48e-2, 6.76e-2,
            #                          0.184,    0.303,   0.500,   0.823,   1.353,
            #                          1.738,    2.232,   2.865,   3.68,    6.07,
            #                          7.79,     10.,     12.,     13.5,    15.])

            # # Read equiprobable outgoing photon energies
            # # Equiprobable outgoing photon energies for incident neutron
            # # energy i
            # e_out_photon_equi = self.xss[ind:ind+600]
            # if len(e_out_photon_equi) == 600:
            #     self.e_out_photon_equi = e_out_photon_equi
            #     self.e_out_photon_equi.shape = (30, 20)

    def _read_mtrp(self):
        """Get the list of reaction MTs for photon-producing reactions for this
        cross-section table. The MT values are somewhat arbitrary.
        """
        LMT = self.jxs[13]
        NMT = self.nxs[6]
        mts = np.asarray(self.xss[LMT:LMT+NMT], dtype=int)
        rxs = [(mt, Reaction(mt, self)) for mt in mts]
        self.photon_reactions.update(rxs)

    def _read_lsigp(self):
        """Determine location of cross sections for each photon-producing reaction
        MT.
        """
        LXS = self.jxs[14]
        NMT = self.nxs[6]
        loca = np.asarray(self.xss[LXS:LXS+NMT], dtype=int)
        for loc, rxn in zip(loca, self.photon_reactions.values()):
            rxn.LOCA = loc

    def _read_sigp(self):
        """Read cross-sections for each photon-producing reaction MT.
        """
        cdef int ind, jxs15, MFTYPE, n_regions, NE

        jxs15 = self.jxs[15]
        for rxn in self.photon_reactions.values():
            ind = jxs15 + rxn.LOCA - 1
            MFTYPE = int(self.xss[ind])
            ind += 1

            if MFTYPE == 12 or MFTYPE == 16:
                # Yield data taken from ENDF File 12 or 6
                MTMULT = int(self.xss[ind])
                ind += 1

                # ENDF interpolation parameters
                n_regions = int(self.xss[ind])
                dat = np.asarray(self.xss[ind+1:ind+1+2*n_regions], dtype=int)
                dat.shape = (2, n_regions)
                NBT, INT = dat
                ind += 1 + 2*n_regions

                # Energy-dependent yield
                NE = int(self.xss[ind])
                dat = self.xss[ind+1:ind+1+2*NE]
                dat.shape = (2, NE)
                rxn.e_yield, rxn.photon_yield = dat
                ind += 1 + 2*NE
            elif MFTYPE == 13:
                # Cross-section data from ENDF File 13
                # Energy grid index at which data starts
                rxn.IE = int(self.xss[ind]) - 1

                # Cross sections
                NE = int(self.xss[ind+1])
                self.sigma = self.xss[ind+2:ind+2+NE]
                ind += 2 + NE
            else:
                raise ValueError("MFTYPE must be 12, 13, 16. Got {}".format(MFTYPE))

    def _read_landp(self):
        """Determine location of angular distribution for each photon-producing
        reaction MT.
        """
        jxs16 = self.jxs[16]
        NMT = self.nxs[6]
        locb = np.asarray(self.xss[jxs16:jxs16+NMT], dtype=int)
        for loc, rxn in zip(locb, self.photon_reactions.values()):
            rxn.LOCB = loc

    def _read_andp(self):
        """Find the angular distribution for each photon-producing reaction
        MT."""
        cdef int ind, i, j, jxs17, NE

        jxs17 = self.jxs[17]
        for i, rxn in enumerate(self.photon_reactions.values()):
            if rxn.LOCB == 0:
                # No angular distribution data are given for this reaction,
                # isotropic scattering is asssumed in LAB
                continue

            ind = jxs17 + rxn.LOCB - 1

            # Number of energies and incoming energy grid
            NE = int(self.xss[ind])
            self.a_dist_energy_in = self.xss[ind+1:ind+1+NE]
            ind += 1 + NE

            # Location of tables associated with each outgoing angle
            # distribution
            LC = np.asarray(self.xss[ind:ind+NE], dtype=int)

            # 32 equiprobable cosine bins for each incoming energy
            a_dist_mu_out = {}
            for j, location in enumerate(LC):
                if location == 0:
                    continue
                ind = jxs17 + location - 1
                a_dist_mu_out[j] = self.xss[ind:ind+33]
            self.a_dist_mu_out = a_dist_mu_out

    def _read_yp(self):
        """Read list of reactions required as photon production yield
        multipliers.
        """
        if self.nxs[6] != 0:
            ind = self.jxs[20]
            NYP = int(self.xss[ind])
            if NYP > 0:
                dat = np.asarray(self.xss[ind+1:ind+1+NYP], dtype=int)
                self.MT_for_photon_yield = dat

    def _read_fis(self):
        """Read total fission cross-section data if present. Generally,
        this table is not provided since it is redundant.
        """
        # Check if fission block is present
        ind = self.jxs[21]
        if ind == 0:
            return

        # Read fission cross sections
        self.IE_fission = int(self.xss[ind]) - 1  # Energy grid index
        NE = int(self.xss[ind+1])
        self.sigma_f = self.xss[ind+2:ind+2+NE]

    def _read_unr(self):
        """Read the unresolved resonance range probability tables if present.
        """
        cdef int ind, N, M, INT, ILF, IOA, IFF

        # Check if URR probability tables are present
        ind = self.jxs[23]
        if ind == 0:
            return

        N = int(self.xss[ind])     # Number of incident energies
        M = int(self.xss[ind+1])   # Length of probability table
        INT = int(self.xss[ind+2]) # Interpolation parameter (2=lin-lin, 5=log-log)
        ILF = int(self.xss[ind+3]) # Inelastic competition flag
        IOA = int(self.xss[ind+4]) # Other absorption flag
        IFF = int(self.xss[ind+5]) # Factors flag
        ind += 6

        self.urr_energy = self.xss[ind:ind+N] # Incident energies
        ind += N

        # Set up URR probability table
        urr_table = self.xss[ind:ind+N*6*M]
        urr_table.shape = (N, 6, M)
        self.urr_table = urr_table

    def find_reaction(self, mt):
        return self.reactions.get(mt, None)

    def __iter__(self):
        # Generators not supported in Cython
        #for r in self.reactions.values():
        #    yield r
        return iter(self.reactions.values())

class EnergyDistribution(object):
    def __init__(self):
        pass

class SabTable(AceTable):
    """A SabTable object contains thermal scattering data as represented by
    an S(alpha, beta) table.

    Parameters
    ----------
    name : str
        ZAID identifier of the table, e.g. lwtr.10t.
    awr : float
        Atomic mass ratio of the target nuclide.
    temp : float
        Temperature of the target nuclide in eV.

    Attributes
    ----------
    awr : float
        Atomic mass ratio of the target nuclide.

    elastic_e_in : list of floats
        Incoming energies in MeV for which the elastic cross section is
        tabulated.

    elastic_P : list of floats
        Elastic scattering cross section for data derived in the incoherent
        approximation, or Bragg edge parameters for data derived in the coherent
        approximation.

    elastic_type : str
        Describes the behavior of the elastic cross section, i.e. whether it was
        derived in the incoherent or coherent approximation.

    inelastic_e_in : list of floats
        Incoming energies in MeV for which the inelastic cross section is
        tabulated.

    inelastic_sigma : list of floats
        Inelastic scattering cross section in barns at each energy.

    name : str
        ZAID identifier of the table, e.g. 92235.70c.

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
        """Read energy-dependent inelastic scattering cross sections.
        """
        ind = self.jxs[1]
        NE = int(self.xss[ind])
        self.inelastic_e_in = self.xss[ind+1:ind+1+NE]
        self.inelastic_sigma = self.xss[ind+1+NE:ind+1+2*NE]

    def _read_itce(self):
        """Read energy-dependent elastic scattering cross sections.
        """
        # Determine if ITCE block exists
        ind = self.jxs[4]
        if ind == 0:
            return

        # Read values
        NE = int(self.xss[ind])
        self.elastic_e_in = self.xss[ind+1:ind+1+NE]
        self.elastic_P = self.xss[ind+1+NE:ind+1+2*NE]

        if self.nxs[5] == 4:
            self.elastic_type = 'sigma=P'
        else:
            self.elastic_type = 'sigma=P/E'

    def _read_itxe(self):
        """Read coupled energy/angle distributions for inelastic scattering.
        """
        # Determine number of energies and angles
        NE_in = len(self.inelastic_e_in)
        NE_out = self.nxs[4]
        NMU = self.nxs[3]
        ind = self.jxs[3]

        self.inelastic_e_out = self.xss[ind:ind+NE_in*NE_out*(NMU+2):NMU+2]
        self.inelastic_e_out.shape = (NE_in, NE_out)

        self.inelastic_mu_out = self.xss[ind:ind+NE_in*NE_out*(NMU+2)]
        self.inelastic_mu_out.shape = (NE_in, NE_out, NMU+2)
        self.inelastic_mu_out = self.inelastic_mu_out[:,:,1:]

    def _read_itca(self):
        """Read angular distributions for elastic scattering.
        """
        NMU = self.nxs[6]
        if self.jxs[4] == 0 or NMU == -1:
            return
        ind = self.jxs[6]

        NE = len(self.elastic_e_in)
        self.elastic_mu_out = self.xss[ind:ind+NE*NMU]
        self.elastic_mu_out.shape = (NE, NMU)


class Reaction(object):
    """Reaction(MT, table=None)

    A Reaction object represents a single reaction channel for a nuclide with
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

    Attributes
    ----------
    ang_energy_in : list of floats
        Incoming energies in MeV at which angular distributions are tabulated.

    ang_energy_cos : list of floats
        Scattering cosines corresponding to each point of the angular distribution
        functions.

    ang_energy_pdf : list of floats
        Probability distribution function for angular distribution.

    ang_energy_cdf : list of floats
        Cumulative distribution function for angular distribution.

    e_dist_energy : list of floats
        Incoming energies in MeV at which energy distributions are tabulated.

    e_dist_law : int
        ACE law used for secondary energy distribution.

    IE : int
        The index on the energy grid corresponding to the threshold of this
        reaction.

    MT : int
        The ENDF MT number for this reaction. On occasion, MCNP uses MT numbers
        that don't correspond exactly to the ENDF specification.

    Q : float
        The Q-value of this reaction in MeV.

    sigma : list of floats
        Microscopic cross section for this reaction at each point on the energy
        grid above the threshold value.

    TY : int
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
        self.IE = 0        # Energy grid index
        self.sigma = []    # Cross section values

    def broaden(self, T_high):
        pass

    def threshold(self):
        """threshold()

        Return energy threshold for this reaction.
        """
        return self.table.energy[self.IE]

    def __repr__(self):
        name = label(self.MT)
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

if __name__ == '__main__':
    # Might be nice to check environment variable DATAPATH to search for xsdir
    # and list files that could be read?
    pass
