#!/usr/bin/env python

"""
The CCCC module contains a number of classes for reading various cross section,
flux, geometry, and data files with specifications given by the Committee for
Computer Code Coordination. The following types of files can be read using
classes from this module: ISOTXS, DLAYXS, BRKOXS, RTFLUX, ATFLUX, RZFLUX, MATXS,
and SPECTR.

The ISOTXS reader was originally derived from Professor James Holloway's
open-source C++ classes from the University of Michigan and later expanded by
Nick Touran for work on his PhD thesis. DLAYXS was later added by Paul Romano.
RTFLUX was done by Elliott Biondo.

A description of several CCCC formats are available online for ISOTXS_, MATXS_,
RTFLUX_, and RZFLUX_. Other format specifications can be found in Los Alamos
Report LA-5324-MS_.

.. _ISOTXS: http://t2.lanl.gov/codes/transx-hyper/isotxs.html

.. _MATXS: http://t2.lanl.gov/codes/transx-hyper/matxs.html

.. _RTFLUX: http://t2.lanl.gov/codes/transx-hyper/rtflux.html

.. _RZFLUX: http://t2.lanl.gov/codes/transx-hyper/rzflux.html

.. _LA-5324-MS: http://www.osti.gov/bridge/servlets/purl/5369298-uIcX6p/

"""

from __future__ import division
from warnings import warn
import numpy as np

from pyne.utils import QAWarning
from pyne.binaryreader import _BinaryReader, _FortranRecord

warn(__name__ + " is not yet QA compliant.", QAWarning)


class Isotxs(_BinaryReader):
    """An Isotxs object represents a binary ISOTXS file written according to the
    CCCC specifications.

    :Attributes:
      **chi** : list of floats
        Fission yields by group.

      **emax** : list of floats
        Maximum energy bound for each group

      **emin** : float
        Minimum energy bound of set

      **fc** : dict
        Dictionary with file-control information

      **fileVersion** : int
        Version of the ISOTXS file.

      **label** : str
        File identification string

      **nuclides** : list of _Nuclides
        List of individual nuclides in the ISOTXS file.

      **vel** : float
        Mean neutron velocity in each group.

    Parameters
    ----------
    filename : str
        Path of the ISOTXS file to load.

    """

    def __init__(self, filename):
        super(Isotxs, self).__init__(filename)

        # Initialize attributes
        self.fc = {}       # file control info
        self.nuclides = []  # : List of nuclides in ISOTXS file.

    def read(self):
        """Read through and parse the ISOTXS file."""

        self._read_file_ID()
        self._read_file_control()
        self._read_file_data()

        # Read file-wide chi-distribution matrix if present. Note that if
        # file-wide chi is given as a vector, it will be read during
        # the read_file_data method.
        if self.fc['ichidst'] > 1:
            self._read_chi_data()

        # Read nuclide data
        for nucName in self.nucNames:
            # Create nuclide object
            nuc = _Nuclide(nucName)

            # Read nuclide name and global data
            self._read_nuclide_data(nuc)

            # Read nuclide cross sections
            self._read_nuclide_xs(nuc)

            # Read nuclide chi data if present
            if nuc.libParams['chiFlag'] > 1:
                self._read_nuclide_chi(nuc)

            # Read nuclide scattering matrix
            for block in range(self.fc['nscmax']):
                for subBlock in range(self.fc['nsblok']):
                    if nuc.libParams['ords'][block] > 0:
                        self._read_nuclide_scatter(nuc, block, subBlock)

            # Add nuclide to dictionary
            self.nuclides.append(nuc)

    def _read_file_ID(self):
        """Reads the file identification block. This block is always present in
        the ISOTXS format and contains a label and file version number.
        """

        # Get first record from file
        fileID = self.get_fortran_record()

        # Read data from file identification record
        self.label = fileID.get_string(24)[0]
        self.fileVersion = fileID.get_int()[0]

    def _read_file_control(self):
        """Reads the file control block. This block is always present and gives
        many parameters for the file including number of energy groups, number
        of isotopes, etc.
        """

        # Get file control record
        fc = self.get_fortran_record()

        # Read data from file control record
        self.fc['ngroup'] = fc.get_int()[0]  # Number of energy groups in file
        self.fc['niso'] = fc.get_int()[0]  # Number of isotopes in file
        # Maximum number of upscatter groups
        self.fc['maxup'] = fc.get_int()[0]
        # Maximum number of downscatter groups
        self.fc['maxdown'] = fc.get_int()[0]
        self.fc['maxord'] = fc.get_int()[0]  # Maximum scattering order
        self.fc['ichidst'] = fc.get_int()[0]  # File-wide fission spectrum flag
        # Max blocks of scatter data (seems to be actual number)
        self.fc['nscmax'] = fc.get_int()[0]
        self.fc['nsblok'] = fc.get_int()[0]  # Number of subblocks

    def _read_file_data(self):
        """Reads the file data block. This block is always present and contains
        isotope names, global chi distribution, energy group structure, and
        locations of each nuclide record.
        """

        # Get file data record
        fileData = self.get_fortran_record()

        # Skip identification label of file
        fileData.get_string(12*8)

        # Read nuclide label for each nuclide
        self.nucNames = fileData.get_string(8, self.fc['niso'])
        self.nucNames = [name.strip() for name in self.nucNames]

        # Read file-wide chi distribution vector
        if self.fc['ichidst'] == 1:
            self.chi = fileData.get_float(self.fc['ngroup'])

        #: Mean neutron velocity in each group
        self.vel = fileData.get_float(self.fc['ngroup'])

        # Read maximum energy bound of each group
        self.emax = fileData.get_float(self.fc['ngroup'])

        # Read minimum energy bound of set
        self.emin = fileData.get_float()[0]

        # Read number of records to be skipped to read data for a given nuclide
        self.locs = fileData.get_int(self.fc['niso'])

    def _read_chi_data(self):
        """Reads file-wide chi-distribution matrix. In most cases, chi will be
        given as a vector, not a matrix, and thus in such cases this routine is
        not needed.
        """

        raise NotImplementedError

    def _read_nuclide_data(self, nuc):
        """Read the following individual nuclide XS record. Load data into nuc.
        This record contains non-mg data like atomic mass, temperature, and some
        flags.
        """

        # Get nuclide data record
        r = self.get_fortran_record()

        # Read nuclide data
        nuc.libParams['nuclide'] = r.get_string(
            8)[0].strip()  # absolute nuclide label
        nuc.libParams['libName'] = r.get_string(
            8)[0]  # library name (ENDFV, etc. )
        nuc.libParams['isoIdent'] = r.get_string(8)[0]
        nuc.libParams['amass'] = r.get_float()[0]  # gram atomic mass
        # thermal energy yield/fission
        nuc.libParams['efiss'] = r.get_float()[0]
        # thermal energy yield/capture
        nuc.libParams['ecapt'] = r.get_float()[0]
        nuc.libParams['temp'] = r.get_float()[0]  # nuclide temperature (K)
        # potential scattering (b/atom)
        nuc.libParams['sigPot'] = r.get_float()[0]
        # density of nuclide (atom/b-cm)
        nuc.libParams['adens'] = r.get_float()[0]
        nuc.libParams['classif'] = r.get_int()[0]  # nuclide classification
        nuc.libParams['chiFlag'] = r.get_int()[0]  # fission spectrum flag
        nuc.libParams['fisFlag'] = r.get_int()[0]  # (n,f) cross section flag
        nuc.libParams['nalph'] = r.get_int()[0]  # (n,alpha) cross section flag
        nuc.libParams['np'] = r.get_int()[0]  # (n,p) cross section flag
        nuc.libParams['n2n'] = r.get_int()[0]  # (n,2n) cross section flag
        nuc.libParams['nd'] = r.get_int()[0]  # (n,d) cross section flag
        nuc.libParams['nt'] = r.get_int()[0]  # (n,t) cross section flag
        nuc.libParams['ltot'] = r.get_int()[0]  # number of moments of total xs
        # number of moments of transport xs
        nuc.libParams['ltrn'] = r.get_int()[0]
        # number of coord directions for transport xs
        nuc.libParams['strpd'] = r.get_int()[0]

        # Read scattering matrix type identifications for each scatter
        # block. Could be total, inelastic, elastic, n2n
        nuc.libParams['scatFlag'] = r.get_int(self.fc['nscmax'])

        # Read number of scattering orders in each scatter block.
        nuc.libParams['ords'] = r.get_int(self.fc['nscmax'])

        # Read number of groups that scatter into group j, including
        # self-scatter, in scatter block n.
        nuc.libParams['jband'] = {}
        for n in range(self.fc['nscmax']):
            for j in range(self.fc['ngroup']):
                nuc.libParams['jband'][j, n] = r.get_int()[0]

        # Read position of in-group scattering cross section for group j,
        # scattering block n, counted from first word of group j data
        nuc.libParams['jj'] = {}
        for n in range(self.fc['nscmax']):
            for j in range(self.fc['ngroup']):
                nuc.libParams['jj'][j, n] = r.get_int()[0]

    def _read_nuclide_xs(self, nuc):
        """Reads principal microscopic multigroup cross-section data for a
        single nuclide.
        """

        # Get cross section record
        r = self.get_fortran_record()

        # PL-weighted transport cross section in group g for Legendre order l
        for l in range(nuc.libParams['ltrn']):
            for g in range(self.fc['ngroup']):
                nuc.micros['transport', g, l] = r.get_float()[0]

        # PL-weighted total cross section in group g for Legendre order l
        for l in range(nuc.libParams['ltot']):
            for g in range(self.fc['ngroup']):
                nuc.micros['total', g, l] = r.get_float()[0]

        # Microscopic (n,gamma) cross section in group g
        for g in range(self.fc['ngroup']):
            nuc.micros['n,g', g] = r.get_float()[0]

        # Read fission data if present
        if nuc.libParams['fisFlag'] > 0:

            # Microscopic (n,fission) cross section in group g
            for g in range(self.fc['ngroup']):
                nuc.micros['fis', g] = r.get_float()[0]

            # Total number of neutrons/fission in group g
            for g in range(self.fc['ngroup']):
                nuc.micros['nu', g] = r.get_float()[0]

        # Read fission spectrum vector if present
        if nuc.libParams['chiFlag'] == 1:
            # Nuclide chi in group g
            for g in range(self.fc['ngroup']):
                nuc.micros['chi', g] = r.get_float()[0]
        else:
            if nuc.libParams['fisFlag'] > 0:
                # Make sure file-wide chi exists
                assert self.fc['ichidst'] == 1, "Fissile nuclide %s in library but no individual or global chi!" % nuc

                # Set the chi to the file-wide chi distribution if this nuclide
                # has a fission cross section
                for g in range(self.fc['ngroup']):
                    nuc.micros['chi', g] = self.chi[g]

        # Read some other important cross sections, if they exist
        for xstype in ['nalph', 'np', 'n2n', 'nd', 'nt']:
            if nuc.libParams[xstype]:
                for g in range(self.fc['ngroup']):
                    nuc.micros[xstype, g] = r.get_float()[0]

        # Read coordinate direction transport cross section (for various
        # coordinate directions)
        if nuc.libParams['strpd'] > 0:
            for i in range(nuc.libParams['strpd']):
                for g in range(self.fc['ngroup']):
                    nuc.micros['strpd', g, i] = r.get_float()[0]

    def _read_nuclide_chi(self, nuc):
        """Reads nuclide-level fission spectrum matrix. In most cases, chi will
        be given as a vector, not a matrix, and thus in such cases this routine
        is not needed.
        """

        raise NotImplementedError

    def _read_nuclide_scatter(self, nuc, block, subBlock):
        """Read nuclide scattering matrix.

        In some versions of the specification, the written description of the
        scattering matrix is wrong! The person who was typing that version had
        shifted their right hand one key to the right on the keyboard resulting
        in gibberish. The CCCC-IV pdf has the correct specification.
        """

        # Get record
        r = self.get_fortran_record()

        # Copy values for number of groups and number of subblocks
        ng = self.fc['ngroup']
        nsblok = self.fc['nsblok']

        # Make sure blocks and subblocks are indexed starting from 1
        m = subBlock + 1
        n = block + 1

        # Determine number of scattering orders in this block
        lordn = nuc.libParams['ords'][block]

        # This is basically how many scattering cross sections there are for
        # this scatter type for this nuclide
        jl = (m - 1)*((ng - 1)//nsblok + 1) + 1
        jup = m*((ng - 1)//nsblok + 1)
        ju = min(ng, jup)

        # Figure out kmax for this sub-block.
        kmax = 0
        for j in range(jl, ju+1):
            g = j - 1  # convert to groups starting at 0
            kmax += nuc.libParams['jband'][g, block]
            # scattering from group j

        for order in range(lordn):
            # for k in range(kmax):
            for j in range(jl, ju+1):
                # There are JBAND values for scattering into group j listed in
                # order of the "from" group as from j+jup to j, from j+jup-1 to
                # j, ...,from j to j, from j-1 to j, j-2 to j, ... , j-down to j
                # anything listed to the left of j represents
                # upscatter. anything to the right is downscatter. n,2n on
                # MC**2-2 ISOTXS scatter matrix are reaction based and need to
                # be multiplied by 2 to get the correct neutron balance.
                g = j-1
                assert g >= 0, "loading negative group in ISOTXS."
                jup = nuc.libParams['jj'][g, block] - 1
                jdown = nuc.libParams['jband'][g, block] - \
                    nuc.libParams['jj'][g, block]
                fromgroups = list(range(j-jdown, j+jup+1))
                fromgroups.reverse()
                for k in fromgroups:
                    fromg = k-1
                    nuc.micros['scat', block, g, fromg, order] = r.get_float()[
                        0]

    def find_nuclide(self, name):
        """Returns a nuclide with a given name.

        Parameters
        ----------
        name : str
            Path of the ISOTXS file to load.

        Returns
        -------
        nuc : Nuclide
            Object containing microscopic cross sections and other data.

        """

        for nuc in self:
            if nuc.name == name:
                return nuc
        return None

    def __iter__(self):
        for nuc in self.nuclides:
            yield nuc

    def __repr__(self):
        return "<ISOTXS File: {0}>".format(self.f.name)


class Dlayxs(_BinaryReader):
    """A Dlayxs object represents the data stored in a CCCC-format DLAYXS
    file. This file contains delayed neutron precursor yields, emission spectra,
    and decay constants reduced to multigroup form. Typically, the data in a
    DLAYXS file would be related to cross-section files in ISOTXS and GRUPXS.

    :Attributes:
      **isotopes** : list of strs
        Names of the isotopes in the DLAYXS file.

      **isotopeFamily** : dict
        Dictionary whose keys are the isotope names and whose values are

      **decay** : dict
        Dictionary whose keys are names of nuclides and whose values are decay
        constants for each delayed neutron family.

      **spectrum** : dict

      **nGroups** : int
        Number of energy groups

      **nIsotopes** : int
        Number of isotopes

      **nFamilies** : int
        Number of delayed neutron families

      **nu** : dict


    Parameters
    ----------
    filename : str
        Path of the DLAYXS file to load.

    """

    def __init__(self, filename):
        super(Dlayxs, self).__init__(filename)

        self.isotopeFamily = {}
        self.decay = {}
        self.spectrum = {}
        self.nu = {}

    def read(self):
        """Read through and parse data in the DLAYXS file."""

        self._read_file_ID()
        self._read_file_control()
        (decay, spectrum) = self._read_spectra()
        self._read_yield()

        for isotope in self.isotopes:
            self.decay[isotope] = {}
            self.spectrum[isotope] = {}
            for gDelay in [1, 2, 3, 4, 5, 6]:
                family = self.isotopeFamily[isotope][gDelay-1]
                self.decay[isotope][gDelay] = decay[family]
                self.spectrum[isotope][gDelay] = spectrum[family]

    def _read_file_ID(self):
        """Read file ID block"""

        id = self.get_fortran_record()
        self.label = id.get_string(24)[0]
        fileID = id.get_int()[0]

    def _read_file_control(self):
        """Read file control block."""

        fileControl = self.get_fortran_record()
        self.nGroups = fileControl.get_int()[0]
        self.nIsotopes = fileControl.get_int()[0]
        self.nFamilies = fileControl.get_int()[0]

    def _read_spectra(self):
        """Read the decay constants and delayed neutron spectra"""

        fileData = self.get_fortran_record()
        self.isotopes = fileData.get_string(8, self.nIsotopes)

        # Read decay constants for each family. We will follow the convention
        # of the CCCC files that the families are indexed starting from 1.
        decay = {}
        for family in range(1, self.nFamilies+1):
            decay[family] = fileData.get_float()[0]

        # Read the delayed neutron spectra for each family
        spectra = {}
        for family in range(1, self.nFamilies+1):
            spectra[family] = fileData.get_float(self.nGroups)

        # This reads the maximum E for each energy group in eV as well as the
        # minimum energy bound of the set in eV.

        self.energySpectra = fileData.get_float(self.nGroups)
        self.minEnergy = fileData.get_float()[0]

        # Determine the number of families to which fission each isotope
        # contributes to delayed neutron precursors and the number of records
        # to be skipped to read data for each isotope

        ## nFamilies = fileData.get_int(self.nIsotopes)
        ## nSkip     = fileData.get_int(self.nIsotopes)

        return decay, spectra

    def _read_yield(self):
        """Read the delayed neutron precursor yields"""

        for isotope in self.isotopes:
            yieldData = self.get_fortran_record()
            self.nu[isotope] = {}
            for gDelay in [1, 2, 3, 4, 5, 6]:
                self.nu[isotope][gDelay] = yieldData.get_float(self.nGroups)
            self.isotopeFamily[isotope] = yieldData.get_int(6)


class Brkoxs(_BinaryReader):
    """A Brkoxs object represents data stored in a BRKOXS file from the CCCC
    format specification. This file is given in conjunction with an ISOTXS (or
    GRUPXS) file when the Bondarenko self-shielding method is to be used.

    Parameters
    ----------
    filename : str
        Path of the BRKOXS file to read.

    """

    def __init__(self, filename):
        super(Brkoxs, self).__init__(filename)


class Rtflux(object):
    """An Rtflux object represents data stored in a RTFLUX file from the CCCC
    format specification. This file contains regular (i.e. not adjoint) total
    fluxes. Attribute names mirror those described in the CCCC specification,
    found here:

    http://t2.lanl.gov/nis/codes/transx-hyper/rtflux.html

    Attributes:
    -----------
    hname: str
        Name of file ("rtflux" or "atflux")
    huse: str
        User identification string
    ivers: int
        File version
    ndim: int
        Number of dimenstions
    ngroup: int
        Number of energy groups
    ninti: int
        Number of fine mesh intervals in the first dimension
    nintj: int
        Number of fine mesh intervals in the second dimension
    nintk: int
        Number of fine mesh intervals in the third dimension
    iter: int
        Outer interation number
    effk: float
        Effective multiplication (k)
    nblok: int
        Number of Fortran data blocks
    flux: ndarray
        Fluxes in the form flux(i, j) where i is interval and j is energy group
    adjoint: bool
        Specify if fluxes are adjoint (e.g. for an atflux file)
    """

    def __init__(self, filename):
        """
        Parameters
        ----------
        filename : str
            Path to the RTFLUX file to be read.
        """

        b = _BinaryReader(filename)
        fr = b.get_fortran_record()

        # read file identification
        self.hname = fr.get_string(8)[0].strip()
        self.huse = fr.get_string(8)[0].strip()
        self.ivers = fr.get_string(8)[0].strip()
        mult = fr.get_int(1)

        if self.hname == "rtflux":
            self.adjoint = False
        elif self.hname == "atflux":
            self.adjoint = True

        # read specifcations
        fr = b.get_fortran_record()
        self.ndim, self.ngroup, self.ninti, self.nintj, self.nintk, self.niter \
            = fr.get_int(6)
        self.effk = fr.get_float(1)[0]
        if not self.adjoint:
            self.power = fr.get_float(1)[0]
        else:
            fr.get_float(1)
        self.nblok = fr.get_int(1)[0]

        # read fluxes
        flux = []

        # This is the 1D binary spec, specified by CCCC.
        # It does not work the the PyNE binary reader, but using the 3D format
        # does work, as tested.
        #
        # if self.ndim == 1:
        #    for m in range(1, self.nblok + 1):
        #        fr = b.get_fortran_record()
        #        print fr.num_bytes
        #        jl = (m - 1)*((self.ngroup - 1)/self.nblok + 1) + 1
        #        jup = m*((self.ngroup -1)/self.nblok + 1)
        #        ju = min(self.ngroup, jup)
        #        flux += fr.get_double(int(self.ninti*(ju-jl+1)))

        # 3D binary spec
        for l in range(1, self.ngroup + 1):
            for k in range(1, self.nintk + 1):
                for m in range(1, self.nblok + 1):
                    fr = b.get_fortran_record()
                    jl = (m - 1)*((self.nintj - 1)/self.nblok + 1) + 1
                    jup = m*((self.nintj - 1)/self.nblok + 1)
                    ju = min(self.nintj, jup)
                    flux += fr.get_double(int(self.ninti*(ju-jl+1)))

        flux2 = []
        num_intervals = self.ninti*self.nintj*self.nintk
        for i in range(self.ngroup):
            if not self.adjoint:
                flux2.insert(0, flux[i*num_intervals:(i+1)*num_intervals])
            else:
                flux2.append(flux[i*num_intervals:(i+1)*num_intervals])

        flux2 = np.array(flux2)
        flux2 = flux2.transpose()

        self.flux = flux2
        b.close()

    def to_mesh(self, m, tag_name):
        """This member function tags supplied PyNE Mesh object with the fluxes
        contained in the rtflux file.

        Parameters
        ----------
        m: PyNE Mesh
            A PyNE Mesh object with same x, y, z intervals used to generate
            the rtflux file.
        tag_name: str
             The tag name to use to tag the fluxes onto the mesh.
        """

        from pyne.mesh import HAVE_PYMOAB
        if HAVE_PYMOAB:
            from pyne.mesh import Mesh, NativeMeshTag
        else:
            warn("The PyMOAB optional dependency could not be imported. "
                 "All aspects of the partisn module are not imported.",
                 QAWarning)

        if not m.structured:
            raise ValueError("Only structured mesh is supported.")

        mesh_dims = [len(x) - 1 for x in m.structured_coords]
        if mesh_dims != [self.ninti, self.nintj, self.nintk]:
            raise ValueError("Supplied mesh does not comform to rtflux bounds")

        temp = m.structured_ordering
        m.structured_ordering = 'zyx'
        m.tag = NativeMeshTag(self.ngroup, float, name=tag_name)
        m.tag[:] = self.flux
        m.structured_ordering = temp


class Atflux(Rtflux):
    """An Atflux object represents data stored in a ATFLUX file from the CCCC
    format specification. This file contains adjoint total fluxes. Note that
    this is the same format as RTFLUX. See Rtflux class for a complete list of
    atrributes. The RTFLUX/ATFLUX binary specification is found here:

    http://t2.lanl.gov/nis/codes/transx-hyper/rtflux.html
    """

    def __init__(self, filename):
        """
        Parameters
        ----------
        filename : str
        Path to the ATFLUX file to be read.
        """
        super(Atflux, self).__init__(filename)


class Rzflux(_BinaryReader):
    """A Rzflux object represents data stored in a RZFLUX file from the CCCC
    format specification. This file contains volumetric averages of fluxes by
    broad energy groups for different geometric zones.

    Parameters
    ----------
    filename : str
        Path to the RZFLUX file to be read.

    """

    def __init__(self, filename):
        super(Rzflux, self).__init__(filename)


class Matxs(_BinaryReader):
    """A Matxs object represents data stored in a MATXS file. This file contains
    generalized cross-sections.

    Parameters
    ----------
    filename : str
        Path to the MATXS file to be read.

    """

    def __init__(self, filename):
        super(Matxs, self).__init__(filename)


class Spectr(_BinaryReader):
    """Reads ultra-fine group spectrum file from MC**2"""

    def __init__(self, filename):
        super(SPECTR, self).__init__(filename)
        self.fc = {}
        self.read1D()
        self.flux = self.read2D()

    def read1D(self):
        t1 = self.get_fortran_record()
        self.fc['eig'] = t1.get_float()[0]
        self.fc['buck'] = t1.get_float()[0]
        self.fc['emax'] = t1.get_float()[0]
        self.fc['deltau'] = t1.get_float()[0]
        self.fc['ngrp'] = t1.get_int()[0]
        self.fc['mgcsd'] = t1.get_int()[0]
        self.fc['ncsd'] = t1.get_int()[0]

    def read2D(self):
        t2 = self.get_fortran_record()
        flux = []
        for g in range(self.fc['ngrp']):
            flux.append(t2.get_float()[0])
        return flux


class _Nuclide(object):
    """Contains data about a single nuclide in an ISOTXS file. Originally,
    Touran had his own Nuclide class so this one is provided to supply the basic
    capabilities needed.
    """

    def __init__(self, name):
        self.name = name
        self.libParams = {}
        self.micros = {}

    def __repr__(self):
        return "<Nuclide: {0}>".format(self.name)


if __name__ == '__main__':
    lib = Isotxs('ISOTXS')
