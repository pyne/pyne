#!/usr/bin/env python

"""
The CCCC module contains a number of classes for reading various cross section,
flux, geometry, and data files with specifications given by the Committee for
Computer Code Coordination. The following types of files can be read using
classes from this module: ISOTXS, DLAYXS, BRKOXS, RTFLUX, ATFLUX, RZFLUX, MATXS,
and SPECTR.

The ISOTXS reader was derived from Prof. James Holloway's open-source C++
classes from the University of Michigan and later expanded by Nick Touran for
work on his PhD thesis. DLAYXS was later added by Paul Romano.

A description of the various CCCC formats are available online at:
    http://t2.lanl.gov/codes/transx-hyper/isotxs.html
    http://t2.lanl.gov/codes/transx-hyper/matxs.html
    http://t2.lanl.gov/codes/transx-hyper/rtflux.html
    http://t2.lanl.gov/codes/transx-hyper/rzflux.html

Other format specifications can be found in Los Alamos report LA-5324-MS,
"Report of the Subcommittee on Standard Interface Files." at:
    http://www.osti.gov/bridge/servlets/purl/5369298-uIcX6p/

"""

import struct
import math

from binaryreader import _BinaryReader, _FortranRecord

class Isotxs(_BinaryReader):
    """An Isotxs object represents a binary ISOTXS file written according to the
    CCCC specifications.

    Parameters
    ----------
    filename : str
        Path of the ISOTXS file to load.
    
    """
    
    def __init__(self, filename):
        super(Isotxs, self).__init__(filename)

        # Initialize attributes
        self.fc = {}       # file control info
        self.nuclides = [] #: List of nuclides in ISOTXS file.

    def read(self):
        """Read through and parse the ISOTXS file."""

        self.read_file_ID()
        self.read_file_control()
        self.read_file_data()

        # Read file-wide chi-distribution matrix if present. Note that if
        # file-wide chi is given as a vector, it will be read during
        # the read_file_data method.
        if self.fc['ichidst']>1:
            self.read_chi_data()

        # Read nuclide data
        for nucName in self.nucNames:
            # Create nuclide object
            nuc = _Nuclide(nucName)
            
            # Read nuclide name and global data
            self.read_nuclide_data(nuc)

            # Read nuclide cross sections
            self.read_nuclide_xs(nuc)   

            # Read nuclide chi data if present
            if nuc.libParams['chiFlag']>1:
                self.read_nuclide_chi(nuc)
                
            # Read nuclide scattering matrix
            for block in range(self.fc['nscmax']):
                for subBlock in range(self.fc['nsblok']):
                    if nuc.libParams['ords'][block] > 0:
                        self.read_nuclide_scatter(nuc, block, subBlock)

            # Add nuclide to dictionary
            self.nuclides.append(nuc)
        
    def read_file_ID(self):
        """Reads the file identification block. This block is always present in
        the ISOTXS format and contains a label and file version number.
        """

        # Get first record from file
        fileID = self.get_fortran_record()

        # Read data from file identification record
        self.label = fileID.get_string(24)
        self.fileVersion = fileID.get_int()
        
    def read_file_control(self):
        """Reads the file control block. This block is always present and gives
        many parameters for the file including number of energy groups, number
        of isotopes, etc.
        """
        
        # Get file control record
        fc = self.get_fortran_record()

        # Read data from file control record
        self.fc['ngroup'] = fc.get_int() # Number of energy groups in file
        self.fc['niso'] = fc.get_int() # Number of isotopes in file
        self.fc['maxup'] = fc.get_int() # Maximum number of upscatter groups
        self.fc['maxdown'] = fc.get_int() # Maximum number of downscatter groups
        self.fc['maxord'] = fc.get_int() # Maximum scattering order
        self.fc['ichidst'] = fc.get_int() # File-wide fission spectrum flag
        self.fc['nscmax'] = fc.get_int() # Max blocks of scatter data (seems to be actual number)
        self.fc['nsblok'] = fc.get_int() # Number of subblocks
        
    def read_file_data(self):
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
        if self.fc['ichidst']==1:
            self.chi = fileData.get_float(self.fc['ngroup'])
        
        #: Mean neutron velocity in each group
        self.vel = fileData.get_float(self.fc['ngroup'])

        # Read maximum energy bound of each group
        self.emax = fileData.get_float(self.fc['ngroup'])
        
        # Read minimum energy bound of set
        self.emin = fileData.get_float()
        
        # Read number of records to be skipped to read data for a given nuclide
        self.locs = fileData.get_int(self.fc['niso'])
            
    def read_chi_data(self):
        """Reads file-wide chi-distribution matrix. In most cases, chi will be
        given as a vector, not a matrix, and thus in such cases this routine is
        not needed.
        """

        raise NotImplementedError
    
    def read_nuclide_data(self, nuc):
        """Read the following individual nuclide XS record. Load data into nuc.
        This record contains non-mg data like atomic mass, temperature, and some
        flags.
        """

        # Get nuclide data record
        r = self.get_fortran_record()

        # Read nuclide data
        nuc.libParams['nuclide'] = r.get_string(8).strip() # absolute nuclide label
        nuc.libParams['libName'] = r.get_string(8) # library name (ENDFV, etc. )
        nuc.libParams['isoIdent'] = r.get_string(8)
        nuc.libParams['amass'] = r.get_float() # gram atomic weight
        nuc.libParams['efiss'] = r.get_float() # thermal energy yield/fission
        nuc.libParams['ecapt'] = r.get_float() # thermal energy yield/capture
        nuc.libParams['temp'] = r.get_float() # nuclide temperature (K)
        nuc.libParams['sigPot'] = r.get_float() # potential scattering (b/atom)
        nuc.libParams['adens'] = r.get_float() # density of nuclide (atom/b-cm)
        nuc.libParams['classif'] = r.get_int() # nuclide classification
        nuc.libParams['chiFlag'] = r.get_int() # fission spectrum flag
        nuc.libParams['fisFlag'] = r.get_int() # (n,f) cross section flag
        nuc.libParams['nalph'] = r.get_int() # (n,alpha) cross section flag
        nuc.libParams['np'] = r.get_int() # (n,p) cross section flag
        nuc.libParams['n2n'] = r.get_int() # (n,2n) cross section flag
        nuc.libParams['nd'] = r.get_int() # (n,d) cross section flag
        nuc.libParams['nt'] = r.get_int() # (n,t) cross section flag
        nuc.libParams['ltot'] = r.get_int() # number of moments of total xs
        nuc.libParams['ltrn'] = r.get_int() # number of moments of transport xs
        nuc.libParams['strpd'] = r.get_int() # number of coord directions for transport xs
        
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
                nuc.libParams['jband'][j,n] = r.get_int()
                
        # Read position of in-group scattering cross section for group j,
        # scattering block n, counted from first word of group j data
        nuc.libParams['jj'] = {}
        for n in range(self.fc['nscmax']):
            for j in range(self.fc['ngroup']):
                nuc.libParams['jj'][j,n] = r.get_int()
        
    def read_nuclide_xs(self, nuc):
        """Reads principal microscopic multigroup cross-section data for a
        single nuclide.
        """

        # Get cross section record
        r = self.get_fortran_record()
        
        # PL-weighted transport cross section in group g for Legendre order l
        for l in range(nuc.libParams['ltrn']):
            for g in range(self.fc['ngroup']):
                nuc.micros['transport',g,l] = r.get_float()
        
        # PL-weighted total cross section in group g for Legendre order l
        for l in range(nuc.libParams['ltot']):
            for g in range(self.fc['ngroup']):
                nuc.micros['total',g,l] = r.get_float()
        
        # Microscopic (n,gamma) cross section in group g
        for g in range(self.fc['ngroup']):
            nuc.micros['n,g',g] = r.get_float()
    
        # Read fission data if present
        if nuc.libParams['fisFlag'] > 0:
            
            # Microscopic (n,fission) cross section in group g
            for g in range(self.fc['ngroup']):
                nuc.micros['fis',g] = r.get_float()
        
            # Total number of neutrons/fission in group g
            for g in range(self.fc['ngroup']):
                nuc.micros['nu',g] = r.get_float()
        
        # Read fission spectrum vector if present
        if nuc.libParams['chiFlag'] == 1:
            # Nuclide chi in group g
            for g in range(self.fc['ngroup']):
                nuc.micros['chi',g]=r.get_float()
        else:
            if nuc.libParams['fisFlag'] > 0:
                # Make sure file-wide chi exists
                assert self.fc['ichidst'] == 1, "Fissile nuclide %s in library but no individual or global chi!" % nuc
                
                # Set the chi to the file-wide chi distribution if this nuclide
                # has a fission cross section
                for g in range(self.fc['ngroup']):
                    nuc.micros['chi',g]=self.chi[g]
        
        # Read some other important cross sections, if they exist
        for xstype in ['nalph','np','n2n','nd','nt']:
            if nuc.libParams[xstype]:
                for g in range(self.fc['ngroup']):
                    nuc.micros[xstype,g]=r.get_float()
        
        # Read coordinate direction transport cross section (for various
        # coordinate directions)
        if nuc.libParams['strpd'] > 0:
            for i in range(nuc.libParams['strpd']):
                for g in range(self.fc['ngroup']):
                    nuc.micros['strpd',g,i] = r.get_float()
        
    def read_nuclide_chi(self, nuc):
        """Reads nuclide-level fission spectrum matrix. In most cases, chi will
        be given as a vector, not a matrix, and thus in such cases this routine
        is not needed.
        """

        raise NotImplementedError
    
    def read_nuclide_scatter(self, nuc, block, subBlock):
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
        jl = (m - 1)*((ng - 1)/nsblok + 1) + 1
        jup = m*((ng - 1)/nsblok + 1)
        ju = min(ng, jup)

        # Figure out kmax for this sub-block. 
        kmax = 0
        for j in range(jl, ju+1):
            g = j - 1  # convert to groups starting at 0
            kmax += nuc.libParams['jband'][g,block] 
            # scattering from group j
        
        for order in range(lordn):
            #for k in range(kmax):
            for j in range(jl, ju+1):
                # There are JBAND values for scattering into group j listed in
                # order of the "from" group as from j+jup to j, from j+jup-1 to
                # j, ...,from j to j, from j-1 to j, j-2 to j, ... , j-down to j
                # anything listed to the left of j represents
                # upscatter. anything to the right is downscatter. n,2n on
                # MC**2-2 ISOTXS scatter matrix are reaction based and need to
                # be multiplied by 2 to get the correct neutron balance.
                g = j-1 
                assert g>=0, "loading negative group in ISOTXS."
                jup   = nuc.libParams['jj'][g,block] - 1
                jdown = nuc.libParams['jband'][g,block] - nuc.libParams['jj'][g,block]
                fromGroups = range(j-jdown,j+jup+1)
                fromGroups.reverse()
                for k in fromGroups:
                    fromG = k-1
                    nuc.micros['scat',block,g,fromG,order] = r.get_float()

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
        
        self.read_file_ID()
        self.read_file_control()
        (decay, spectrum) = self.read_spectra()
        self.read_yield()
        
        for isotope in self.isotopes:
            self.decay[isotope]    = {}
            self.spectrum[isotope] = {}
            for gDelay in [1,2,3,4,5,6]:
                family = self.isotopeFamily[isotope][gDelay-1]
                self.decay[isotope][gDelay]    = decay[family]
                self.spectrum[isotope][gDelay] = spectrum[family]
        
    def read_file_ID(self):
        """Read file ID block"""
        
        id = self.get_fortran_record()
        self.label = id.get_string(24)
        fileID = id.get_int()
        
    def read_file_control(self):
        """Read file control block."""
        
        fileControl = self.get_fortran_record()
        self.nGroups = fileControl.get_int()
        self.nIsotopes = fileControl.get_int()
        self.nFamilies = fileControl.get_int()
        
    def read_spectra(self):
        """Read the decay constants and delayed neutron spectra"""
        
        fileData = self.get_fortran_record()
        self.isotopes = fileData.get_string(8, self.nIsotopes)
        
        # Read decay constants for each family. We will follow the convention
        # of the CCCC files that the families are indexed starting from 1.
        decay = {}
        for family in range(1, self.nFamilies+1):
            decay[family] = fileData.get_float()
           
        # Read the delayed neutron spectra for each family
        spectra = {}
        for family in range(1, self.nFamilies+1):
            spectra[family] = fileData.get_float(self.nGroups)
            
        # This reads the maximum E for each energy group in eV as well as the
        # minimum energy bound of the set in eV. 
        
        self.energySpectra = fileData.get_float(self.nGroups)
        self.minEnergy = fileData.get_float()
        
        # Determine the number of families to which fission each isotope
        # contributes to delayed neutron precursors and the number of records
        # to be skipped to read data for each isotope
        
        ## nFamilies = fileData.get_int(self.nIsotopes)
        ## nSkip     = fileData.get_int(self.nIsotopes)
        
        return decay, spectra
    
    def read_yield(self):
        """Read the delayed neutron precursor yields"""
        
        for isotope in self.isotopes:
            yieldData = self.get_fortran_record()
            self.nu[isotope] = {}
            for gDelay in [1,2,3,4,5,6]:
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

class Rtflux(_BinaryReader):
    """A Rtflux object represents data stored in a RTFLUX file from the CCCC
    format specification. This file contains regular total fluxes.

    Parameters
    ----------
    filename : str
        Path to the RTFLUX file to be read.

    """

    def __init__(self, filename):
        super(Rtflux, self).__init__(filename)


class Atflux(_BinaryReader):
    """A Atflux object represents data stored in a ATFLUX file from the CCCC
    format specification. This file contains adjoint total fluxes.

    Parameters
    ----------
    filename : str
        Path to the ATFLUX file to be read.

    """

    def __init__(self, filename):
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
        self.fc['eig'] = t1.get_float()
        self.fc['buck'] = t1.get_float()
        self.fc['emax'] = t1.get_float()
        self.fc['deltau'] = t1.get_float()
        self.fc['ngrp'] = t1.get_int()
        self.fc['mgcsd'] = t1.get_int()
        self.fc['ncsd'] = t1.get_int()
        
    def read2D(self):
        t2 = self.get_fortran_record()
        flux = []
        for g in range(self.fc['ngrp']):
            flux.append(t2.get_float())
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


if __name__=='__main__':
    lib = Isotxs('ISOTXS')
