#!/usr/bin/env python

"""
ISOTXS Reader Based on the CCCC standard. This code was derived from
Prof. James Holloway's open-source C++ classes from the University of Michigan.

Initial work on this module by Nick Touran (for his Thesis mostly). DLAYXS added
by Paul Romano.

Some description of ISOTXS available at
http://t2.lanl.gov/codes/transx-hyper/isotxs.html.

TODO: Find out Nick's thoughts on CrossSection class.
"""

import struct
import math

class CCCCRecord(object):
    """A single formatted record from a CCCC file."""

    def __init__(self,data):
        self.data=data
        self.pos=0
        self.intSize   = struct.calcsize('i')
        self.floatSize = struct.calcsize('f')
    
    def getInt(self):
        """Returns the next integer in the record."""

        (i,) = struct.unpack('i',self.data[self.pos:self.pos+self.intSize])
        self.pos += self.intSize
        return i
                             
    def getFloat(self):
        """Returns the next float in the record."""

        (f,) =  struct.unpack('f',self.data[self.pos:self.pos+self.floatSize])
        self.pos += self.floatSize
        return f
    
    def getDouble(self):
        """Returns the next double in the record."""

        (d,) =  struct.unpack('d',self.data[self.pos:self.pos+self.floatSize*2])
        self.pos += self.floatSize*2
        return d
    
    def getString(self,length):
        """Returns the next string of a specified length in the record."""

        relevantData = self.data[self.pos:self.pos+length]
        #print len(relevantData)
        (s,) = struct.unpack('%ds'%length,relevantData)
        self.pos += length
        return s
    
    def getList(self, type, length, strLength = 0):
        """Returns a list of integers, floats, or strings."""

        if type == "int":
            results = [self.getInt() for i in range(length)]
        elif type == "float":
            results = [self.getFloat() for i in range(length)]
        elif type == "double":
            results = [self.getDouble() for i in range(length)]
        elif type == "string":
            results = [self.getString(strLength).strip() for i in range(length)]
        else:
            print("Do not recognize type: {0}".format(type))
            return None
        return results

    def __repr__(self):
        return "<CCCC Record>"


class CCCCReader(object):
    """Reads a binary file according to CCCC standards. This was created
    following Prof. James Paul Holloway's (hagar@umich.edu) alpha release of
    ccccutils written in C++ from 2001.
    """
    
    def __init__(self,fName='ISOTXS'):
        self.intSize   = struct.calcsize('i')
        self.floatSize = struct.calcsize('d')
        self.f = open(fName,'rb')
        
    def getInt(self):
        (i,) = struct.unpack('i',self.f.read(self.intSize))
        return i
                             
    def getFloat(self):
        (f,) =  struct.unpack('d',self.f.read(self.floatSize))
        return f

    def getRecord(self):
        """CCCC records start with an int and end with the same int. This int
        represents the number of bytes that the record is. That makes it easy to
        read.
        """

        numBytes = self.getInt()
        if numBytes % self.intSize:
            print('Error: numBytes %d is not a multiple of byte size: %d' % (numBytes,self.intSize))
            return

        # Read numBytes from the record
        rec = self.f.read(numBytes)
        
        # now read end of record
        numBytes2 = self.getInt()
        if numBytes2 != numBytes:
            print('Error: numBytes2 %d is not a equal to original byte count: %d' % (numBytes2,numBytes))
            return
        
        return CCCCRecord(rec)


class ISOTXS(CCCCReader):
    """Reads a binary ISOTXS file according to the CCCC specifications."""
    
    def __init__(self,fName,debug=False):
        super(ISOTXS,self).__init__(fName)

        # Initialize attributes
        self.fc = {}       # file control info
        self.nuclides = [] # actual nuclides

    def read(self):
        """Read through and parse the ISOTXS file."""

        self.readFileID()
        self.readFileControl()
        self.readFileData()

        # Read file-wide chi-distribution matrix if present. Note that if
        # file-wide chi is given as a vector, it will be read during
        # the readFileData method.
        if self.fc['ichidst']>1:
            self.readChiData()

        # Read nuclide data
        for nucName in self.nucNames:
            # Create nuclide object
            nuc = Nuclide(nucName)
            
            # Read nuclide name and global data
            self.readNuclideData(nuc)

            # Read nuclide cross sections
            self.readNuclideXS(nuc)   

            # Read nuclide chi data if present
            if nuc.libParams['chiFlag']>1:
                self.readNuclideChi(nuc)
                
            # Read nuclide scattering matrix
            for block in range(self.fc['nscmax']):
                for subBlock in range(self.fc['nsblok']):
                    if nuc.libParams['ords'][block]>0:
                        self.readNuclideScatter(nuc,block,subBlock)

            # Add nuclide to dictionary
            self.nuclides.append(nuc)
        
    def readFileID(self):
        """Reads the file identification block. This block is always present in
        the ISOTXS format and contains a label and file version number.
        """

        # Get first record from file
        fileID = self.getRecord()

        # Read data from file identification record
        self.label       = fileID.getString(24)
        self.fileVersion = fileID.getInt()
        
    def readFileControl(self):
        """Reads the file control block. This block is always present and gives
        many parameters for the file including number of energy groups, number
        of isotopes, etc.
        """
        
        # Get file control record
        fc = self.getRecord()

        # Read data from file control record
        self.fc['ngroup']  = fc.getInt() # Number of energy groups in file
        self.fc['niso']    = fc.getInt() # Number of isotopes in file
        self.fc['maxup']   = fc.getInt() # Maximum number of upscatter groups
        self.fc['maxdown'] = fc.getInt() # Maximum number of downscatter groups
        self.fc['maxord']  = fc.getInt() # Maximum scattering order
        self.fc['ichidst'] = fc.getInt() # File-wide fission spectrum flag
        self.fc['nscmax']  = fc.getInt() # Max blocks of scatter data (seems to be actual number)
        self.fc['nsblok']  = fc.getInt() # Number of subblocks
        
    def readFileData(self):
        """Reads the file data block. This block is always present and contains
        isotope names, global chi distribution, energy group structure, and
        locations of each nuclide record.
        """

        # Get file data record
        fileData = self.getRecord()

        # Skip identification label of file
        fileData.getString(12*8)
        
        # Read nuclide label for each nuclide
        self.nucNames = fileData.getList('string', self.fc['niso'], 8)
        self.nucNames = [name.strip() for name in self.nucNames]
            
        # Read file-wide chi distribution vector
        if self.fc['ichidst']==1:
            self.chi = fileData.getList('float', self.fc['ngroup'])
        
        # Read mean neutron velocity in each group
        self.vel = fileData.getList('float', self.fc['ngroup'])

        # Read maximum energy bound of each group
        self.emax = fileData.getList('float', self.fc['ngroup'])
        
        # Read minimum energy bound of set
        self.emin = fileData.getFloat()
        
        # Read number of records to be skipped to read data for a given nuclide
        self.locs = fileData.getList('int', self.fc['niso'])
            
    def readChiData(self):
        """Reads file-wide chi-distribution matrix. In most cases, chi will be
        given as a vector, not a matrix, and thus in such cases this routine is
        not needed.
        """

        raise NotImplementedError
    
    def readNuclideData(self, nuc):
        """Read the following individual nuclide XS record. Load data into nuc.
        This record contains non-mg data like atomic mass, temperature, and some
        flags.
        """

        # Get nuclide data record
        r = self.getRecord()

        # Read nuclide data
        nuc.libParams['nuclide']  = r.getString(8).strip() # absolute nuclide label
        nuc.libParams['libName']  = r.getString(8) # library name (ENDFV, etc. )
        nuc.libParams['isoIdent'] = r.getString(8)
        nuc.libParams['amass']    = r.getFloat()   # gram atomic weight
        nuc.libParams['efiss']    = r.getFloat()   # thermal energy yield/fission
        nuc.libParams['ecapt']    = r.getFloat()   # thermal energy yield/capture
        nuc.libParams['temp']     = r.getFloat()   # nuclide temperature (K)
        nuc.libParams['sigPot']   = r.getFloat()   # potential scattering (b/atom)
        nuc.libParams['adens']    = r.getFloat()   # density of nuclide (atom/b-cm)
        nuc.libParams['classif']  = r.getInt()     # nuclide classification
        nuc.libParams['chiFlag']  = r.getInt()     # fission spectrum flag
        nuc.libParams['fisFlag']  = r.getInt()     # (n,f) cross section flag
        nuc.libParams['nalph']    = r.getInt()     # (n,alpha) cross section flag
        nuc.libParams['np']       = r.getInt()     # (n,p) cross section flag
        nuc.libParams['n2n']      = r.getInt()     # (n,2n) cross section flag
        nuc.libParams['nd']       = r.getInt()     # (n,d) cross section flag
        nuc.libParams['nt']       = r.getInt()     # (n,t) cross section flag
        nuc.libParams['ltot']     = r.getInt()     # number of moments of total xs
        nuc.libParams['ltrn']     = r.getInt()     # number of moments of transport xs
        nuc.libParams['strpd']    = r.getInt()     # number of coord directions for transport xs
        
        # Read scattering matrix type identifications for each scatter
        # block. Could be total, inelastic, elastic, n2n
        nuc.libParams['scatFlag'] = r.getList('int', self.fc['nscmax'])

        # Read number of scattering orders in each scatter block.
        nuc.libParams['ords'] = r.getList('int', self.fc['nscmax'])

        # Read number of groups that scatter into group j, including
        # self-scatter, in scatter block n.
        nuc.libParams['jband'] = {}
        for n in range(self.fc['nscmax']):
            for j in range(self.fc['ngroup']):
                nuc.libParams['jband'][j,n] = r.getInt()
                
        # Read position of in-group scattering cross section for group j,
        # scattering block n, counted from first word of group j data
        nuc.libParams['jj'] = {}
        for n in range(self.fc['nscmax']):
            for j in range(self.fc['ngroup']):
                nuc.libParams['jj'][j,n] = r.getInt()
                
        
    def readNuclideXS(self, nuc):
        """Reads principal microscopic multigroup cross-section data for a
        single nuclide.
        """

        # Get cross section record
        r = self.getRecord()
        
        # PL-weighted transport cross section in group g for Legendre order l
        for l in range(nuc.libParams['ltrn']):
            for g in range(self.fc['ngroup']):
                nuc.micros['transport',g,l] = r.getFloat()
        
        # PL-weighted total cross section in group g for Legendre order l
        for l in range(nuc.libParams['ltot']):
            for g in range(self.fc['ngroup']):
                nuc.micros['total',g,l] = r.getFloat()
        
        # Microscopic (n,gamma) cross section in group g
        for g in range(self.fc['ngroup']):
            nuc.micros['n,g',g] = r.getFloat()
    
        # Read fission data if present
        if nuc.libParams['fisFlag'] > 0:
            
            # Microscopic (n,fission) cross section in group g
            for g in range(self.fc['ngroup']):
                nuc.micros['fis',g] = r.getFloat()
        
            # Total number of neutrons/fission in group g
            for g in range(self.fc['ngroup']):
                nuc.micros['nu',g] = r.getFloat()
        
        # Read fission spectrum vector if present
        if nuc.libParams['chiFlag'] == 1:
            # Nuclide chi in group g
            for g in range(self.fc['ngroup']):
                nuc.micros['chi',g]=r.getFloat()
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
                    nuc.micros[xstype,g]=r.getFloat()
        
        # Read coordinate direction transport cross section (for various
        # coordinate directions)
        if nuc.libParams['strpd'] > 0:
            for i in range(nuc.libParams['strpd']):
                for g in range(self.fc['ngroup']):
                    nuc.micros['strpd',g,i] = r.getFloat()
        
    def readNuclideChi(self, nuc):
        """Reads nuclide-level fission spectrum matrix. In most cases, chi will
        be given as a vector, not a matrix, and thus in such cases this routine
        is not needed.
        """

        raise NotImplementedError
    
    def readNuclideScatter(self, nuc, block, subBlock):
        """Read nuclide scattering matrix.
        TODO: Tidy this method up a bit.
        """

        r = self.getRecord()
        ng     = self.fc['ngroup']
        nsblok = self.fc['nsblok']
        m=subBlock+1 # fix starting at zero problem and use same indices as CCCC specification
        n=block+1
        kmax=0
        # be careful with starting indices at 0 here!!
        lordn =nuc.libParams['ords'][block]
        # figure out kmax for this sub-block. 
        # this is basically how many scattering cross sections there are for this scatter type for this nuclide 
        jl=(m-1)*((ng-1)/nsblok+1)+1
        jup= m*((ng-1)/nsblok+1)
        ju = min(ng,jup)
        for j in range(jl,ju+1):
            g = j-1  # convert to groups starting at 0
            kmax+=nuc.libParams['jband'][g,block]
            # scattering from group j
        
        #print('jl, jup, ju, kmax: %d %d %d %d' %(jl,jup,ju,kmax))
        
        for order in range(lordn):
            #for k in range(kmax):
            for j in range(jl,ju+1):       # so close. just gotta finish decifering that last paragraph in 7D record description. 
                # god, I love the debugger. 
                # hmm, seems the ISOTXS description online at LANL is messed up. The cccc-iv pdf got me going. 
                # XS are listed by group 
                # there are JBAND values for scattering into group j listed in order of the "from" group as
                # from j+jup to j, from j+jup-1 to j, ...,from j to j, from j-1 to j, j-2 to j, ... , j-down to j
                # anything listed to the left of j represents upscatter. anything to the right is downscatter. yay.
                # n,2n on MC**2-2 ISOTXS scatter matrix are reaction based and need to be multiplied by 2 to get
                # the correct neutron balance.  
                g = j-1 
                assert g>=0, "loading negative group in ISOTXS."
                jup   = nuc.libParams['jj'][g,block] - 1
                jdown = nuc.libParams['jband'][g,block] - nuc.libParams['jj'][g,block]
                fromGroups = range(j-jdown,j+jup+1)
                fromGroups.reverse()
                for k in fromGroups:
                    fromG = k-1
                    nuc.micros['scat',block,g,fromG,order] = r.getFloat()
                #print nuc.micros['scat'][k,l]

    def findNuclide(self, name):
        for nuc in self:
            if nuc.name == name:
                return nuc
        return None
    
    def collapse(self,emaxList,spectr):
        """Given a SPECTR file, this method further condesnses a broad group
        library to the energy structure in emaxList where each value is in
        MeV. To condense to a single group, pass the argument
        emaxList=[1e7]. This method is not complete yet.

        TODO: Find out more about this from Nick.
        """
        
        # u = ln(E/Emax)
        # E = Emax*e^u
        import math
        for hfg,flux in enumerate(spectr.flux):
            # spectr.flux is constant lethargy. 
            #what is the hyper fine energy?
            hfLethargy = hfg*spectr.fc['deltau'] # this hyper-fine group's lethargy
            hfEnergy = spectr.fc['emax']*math.exp(hfLethargy)
            # which broad group are we in?
            ebgMin=0
            for bg,ebgMax in enumerate(self.emax):
                if ebgMin<hfEnergy<ebgMax:
                    break
                ebgMin=ebgMax
            
            # which new broader group are we in?
            nbgMin=0
            for nbg,nbgMax in enumerate(emaxList):
                if nbgMin<=hfEnergy<nbgMax:
                    break
                nbgMin=nbgMax

    def __iter__(self):
        for nuc in self.nuclides:
            yield nuc
            
    def __repr__(self):
        return "<ISOTXS File: {0}>".format(self.f.name)

                
class DLAYXS(CCCCReader):
    """Reads a binary DLAYXS file according to CCCC specification."""
    
    def __init__(self, filename):
        super(DLAYXS,self).__init__(filename)
        
        self.isotopeFamily = {}
        self.decay = {}
        self.spectrum = {}
        self.nu = {}
    
    def read(self):
        
        self.readFileID()
        self.readFileControl()
        (decay, spectrum) = self.readSpectra()
        self.readYield()
        
        for isotope in self.isotopes:
            self.decay[isotope]    = {}
            self.spectrum[isotope] = {}
            for gDelay in [1,2,3,4,5,6]:
                family = self.isotopeFamily[isotope][gDelay-1]
                self.decay[isotope][gDelay]    = decay[family]
                self.spectrum[isotope][gDelay] = spectrum[family]
        
    def readFileID(self):
        """read file ID block"""
        
        id = self.getRecord()
        self.label = id.getString(24)
        fileID     = id.getInt()
        
    def readFileControl(self):
        
        fileControl    = self.getRecord()
        self.nGroups   = fileControl.getInt()
        self.nIsotopes = fileControl.getInt()
        self.nFamilies = fileControl.getInt()
        
    def readSpectra(self):
        """Read the decay constants and delayed neutron spectra"""
        
        fileData = self.getRecord()
        self.isotopes = fileData.getList('string', self.nIsotopes, 8)
        
        # Read decay constants for each family. We will follow the convention
        # of the CCCC files that the families are indexed starting from 1.
        decay = {}
        for family in range(1, self.nFamilies+1):
            decay[family] = fileData.getFloat()
           
        # Read the delayed neutron spectra for each family
        spectra = {}
        for family in range(1, self.nFamilies+1):
            spectra[family] = fileData.getList('float', self.nGroups)
            
        # This reads the maximum E for each energy group in eV as well as the
        # minimum energy bound of the set in eV. 
        
        self.energySpectra = fileData.getList('float', self.nGroups)
        self.minEnergy = fileData.getFloat()
        
        # Determine the number of families to which fission each isotope
        # contributes to delayed neutron precursors and the number of records
        # to be skipped to read data for each isotope
        
        ## nFamilies = fileData.getList('int', self.nIsotopes)
        ## nSkip     = fileData.getList('int', self.nIsotopes)
        
        return decay, spectra
    
    def readYield(self):
        """Read the delayed neutron precursor yields"""
        
        for isotope in self.isotopes:
            yieldData = self.getRecord()
            self.nu[isotope] = {}
            for gDelay in [1,2,3,4,5,6]:
                self.nu[isotope][gDelay] = yieldData.getList('float', self.nGroups)
            self.isotopeFamily[isotope] = yieldData.getList('int', 6)


class SPECTR(CCCCReader):
    """Reads ultra-fine group spectrum file from MC**2"""

    def __init__(self,fName,debug=False):
        CCCCReader.__init__(self,fName)
        self.fc={}
        self.read1D()
        self.flux=self.read2D()
        
    def read1D(self):
        t1 = self.getRecord()
        self.fc['eig']=t1.getFloat()
        self.fc['buck']=t1.getFloat()
        self.fc['emax']=t1.getFloat()
        self.fc['deltau']=t1.getFloat()
        self.fc['ngrp']=t1.getInt()
        self.fc['mgcsd']=t1.getInt()
        self.fc['ncsd']=t1.getInt()
        
    def read2D(self):
        t2 = self.getRecord()
        flux=[]
        for g in range(self.fc['ngrp']):
            flux.append(t2.getFloat())
        return flux


class Nuclide(object):
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
    lib = ISOTXS('ISOTXS')
