#!/usr/bin/env python

"""
Module for parsing and manipulating data from ENDF evaluations. Currently, it
only can read several MTs from File 1, but with time it will be expanded to
include the entire ENDF format.

All the classes and functions in this module are based on document ENDF-102
titled "Data Formats and Procedures for the Evaluated Nuclear Data File
ENDF-6". The latest version from June 2009 can be found at
http://www-nds.iaea.org/ndspub/documents/endf/endf102/endf102.pdf

For more information on this module, contact Paul Romano <romano7@gmail.com>
"""

import re
import numpy as np
import matplotlib.pyplot as plt

number = " (\d.\d+(?:\+|\-)\d)"

libraries = {0: "ENDF/B", 1: "ENDF/A", 2: "JEFF", 3: "EFF",
             4: "ENDF/B High Energy", 5: "CENDL", 6: "JENDL",
             31: "INDL/V", 32: "INDL/A", 33: "FENDL", 34: "IRDF",
             35: "BROND", 36: "INGDB-90", 37: "FENDL/A", 41: "BROND"}

class Evaluation(object):
    """
    Evaluation is the main class for an ENDF evaluation which contains a number
    of Files.
    """

    def __init__(self, filename):
        self.fh = open(filename, 'r')
        self.files = []
        self.verbose = True

    def read(self):
        # First we need to read MT=1, MT=451 which has a description of the ENDF
        # file and a list of what data exists in the file
        self.readHeader()

        # Now we can start looping over the list of data - we can skip the first
        # entry since it is the MT=451 block that we already read
        for MF, MT, NC, MOD in self.reactionList[1:]:
            # File 1 data
            if MF == 1:
                # Number of total neutrons per fission
                if MT == 452:
                    self.readTotalNu()
                # Number of delayed neutrons per fission
                elif MT == 455:
                    self.readDelayedNu()
                # Number of prompt neutrons per fission
                elif MT == 456:
                    self.readPromptNu()
                # Components of energy release due to fission
                elif MT == 458:
                    self.readFissionEnergy()
                elif MT == 460:
                    self.readDelayedPhoton()
            elif MF == 2:
                if MT == 151:
                    self.readResonances()
                    
    def readHeader(self):
        if self.verbose:
            print("Reading File 1...")
        self.printInfo(451)

        # Find File 1 in evaluation
        self.seekFile(1)

        # Now read File1
        file1 = ENDFFile1()
        self.files.append(file1)

        # Create MT for description
        data = ENDFReaction(451)
        file1.reactions.append(data)

        # First HEAD record
        items = self.getHeadRecord()
        data.ZA   = items[0]
        data.AWR  = items[1]
        data.LRP  = items[2]
        data.LFI  = items[3]
        data.NLIB = items[4]
        data.NMOD = items[5]

        # Control record 1
        items = self.getContRecord()
        data.ELIS = items[0]
        data.STA  = int(items[1])
        data.LIS  = items[2]
        data.LISO = items[3]
        data.NFOR = items[5]

        # Control record 2
        items = self.getContRecord()
        data.AWI  = items[0]
        data.EMAX = items[1]
        data.LREL = items[2]
        data.NSUB = items[4]
        data.NVER = items[5]

        # Control record 3
        items = self.getContRecord()
        data.TEMP = items[0]
        data.LDRV = items[2]
        data.NWD  = items[4]
        data.NXC  = items[5]

        # Text record 1
        items = self.getTextRecord()
        text = items[0]
        data.ZSYMAM = text[0:11]
        data.ALAB   = text[11:22]
        data.EDATE  = text[22:32]
        data.AUTH   = text[32:66]

        # Text record 2
        items = self.getTextRecord()
        text = items[0]
        data.REF    = text[1:22]
        data.DDATE  = text[22:32]
        data.RDATE  = text[33:43]
        data.ENDATE = text[55:63]

        # Text record 3
        items = self.getTextRecord()
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
            items = self.getContRecord(line, skipC=True)
            MF  = items[2]
            MT  = items[3]
            NC  = items[4]
            MOD = items[5]
            self.reactionList.append((MF,MT,NC,MOD))

    def readTotalNu(self):
        self.printInfo(452)

        # Find file 1
        file1 = self.findFile(1)

        # Create total nu reaction
        nuTotal = ENDFReaction(452)
        file1.reactions.append(nuTotal)

        # Determine representation of total nu data
        items = self.getHeadRecord()
        nuTotal.LNU = items[3]

        # Polynomial representation
        if nuTotal.LNU == 1:
            nuTotal.coeffs = self.getListRecord()
        # Tabulated representation
        elif nuTotal.LNU == 2:
            nuTotal.value = self.getTab1Record()

        # Skip SEND record
        self.fh.readline()

    def readDelayedNu(self):
        self.printInfo(455)

        # Find file 1
        file1 = self.findFile(1)

        # Create delayed nu reaction
        nuDelay = ENDFReaction(455)
        file1.reactions.append(nuDelay)

        # Determine representation of delayed nu data
        items = self.getHeadRecord()
        nuDelay.LDG = items[2]
        nuDelay.LNU = items[3]

        # Nu tabulated and delayed-group constants are energy-independent
        if nuDelay.LNU == 2 and nuDelay.LDG == 0:
            nuDelay.decayConst = self.getListRecord(onlyList=True)
            nuDelay.value = self.getTab1Record()
            self.fh.readline()
        elif nuDelay.LNU == 2 and nuDelay.LDG == 1:
            raise NotImplementedError
        elif nuDelay.LNU == 1 and nuDelay.LDG == 0:
            raise NotImplementedError
        elif nuDelay.LNU == 1 and nuDelay.LDG == 1:
            raise NotImplementedError

    def readPromptNu(self):
        self.printInfo(456)

        # Create delayed nu reaction
        nuPrompt = ENDFReaction(456)
        self.findFile(1).reactions.append(nuPrompt)

        # Determine representation of delayed nu data
        items = self.getHeadRecord()
        nuPrompt.LNU = items[3]

        # Tabulated values of nu
        if nuPrompt.LNU == 2:
            nuPrompt.value = self.getTab1Record()
        # Spontaneous fission
        elif nuPrompt.LNU == 1:
            nuPrompt.value = self.getListRecord(onlyList=True)

        # Skip SEND record
        self.fh.readline()

    def readFissionEnergy(self):
        self.printInfo(458)

        # Create fission energy release reaction
        eRelease = ENDFReaction(458)
        self.findFile(1).reactions.append(eRelease)

        # Skip HEAD record
        self.getHeadRecord()

        # Read LIST record containing components of fission energy release (or
        # coefficients)
        items, values = self.getListRecord()
        NPLY = items[3]
        if NPLY == 0:
            eRelease.fissProducts = (values[0], values[1])
            eRelease.promptNeuts  = (values[2], values[3])
            eRelease.delayNeuts   = (values[4], values[5])
            eRelease.promptGammas = (values[6], values[7])
            eRelease.delayGammas  = (values[8], values[9])
            eRelease.delayBetas   = (values[10], values[11])
            eRelease.neutrinos    = (values[12], values[13])
            eRelease.pseudoQ      = (values[14], values[15])
            eRelease.total        = (values[16], values[17])
        elif NPLY > 0:
            raise NotImplementedError

        # Skip SEND record
        self.fh.readline()

    def readDelayedPhoton(self):
        self.printInfo(460)

        # Create delayed photon data reaction
        dp = ENDFReaction(460)
        self.findFile(1).reactions.append(dp)

        # Determine whether discrete or continuous representation
        items = self.getHeadRecord()
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
                mult = self.getTab1Record()
                dp.multiplicity.append(mult)

                # Determine energy
                E = mult.params[0]
                dp.energy.append(E)

        # Continuous representation
        elif dp.LO == 2:
            # Determine decay constant and number of precursor families
            dp.decayConst = self.getListRecord(onlyList=True)
            dp.NNF = len(dp.decayConst)

    def readResonances(self):
        if self.verbose:
            print("Reading File 2...")
        self.printInfo(151)

        # Find File 2 in evaluation
        self.seekFile(2)

        # Now read File1
        file2 = ENDFFile2()
        self.files.append(file2)

        # Create MT for resonances
        res = ENDFReaction(151)
        file2.reactions.append(res)
        res.resonances = []

        # Determine whether discrete or continuous representation
        items = self.getHeadRecord()
        res.NIS = items[4] # Number of isotopes
        
        for iso in range(res.NIS):
            items = self.getContRecord()
            res.ABN = items[1] # isotopic abundance
            res.LFW = items[3] # fission widths present?
            res.NER = items[4] # number of resonance energy ranges

            for erange in range(res.NER):
                items = self.getContRecord()
                res.EL   = items[0] # lower limit of energy range
                res.EH   = items[1] # upper limit of energy range
                res.LRU  = items[2] # flag for resolved (1)/unresolved (2)
                res.LRF  = items[3] # resonance representation
                res.NRO  = items[4] # flag for energy dependence of scattering radius
                res.NAPS = items[5] # flag controlling use of channel/scattering radius

                # Only scattering radius specified
                if res.LRU == 0 and res.NRO == 0:
                    items = self.getContRecord()
                    res.SPI = items[0]
                    res.AP = items[1]
                    res.NLS = items[4]
                # Resolved resonance region
                elif res.LRU == 1:
                    self.readResolved(res)
                # Unresolved resonance region
                elif res.LRU == 2:
                    self.readUnresolved(res)

    def readResolved(self, res):
        # Single- or Multi-level Breit Wigner
        if res.LRF == 1 or res.LRF == 2:
            # Read energy-dependent scattering radius if present
            if res.NRO > 0:
                res.AP = self.getTab1Record()
        
            # Other scatter radius parameters
            items = self.getContRecord()
            res.SPI = items[0] # Spin, I, of the target nucleus
            if res.NRO == 0:
                res.AP = items[1]
            res.NLS = items[4] # Number of l-values

            # Read resonance widths, J values, etc
            for l in range(res.NLS):
                headerItems, items = self.getListRecord()
                QX, L, LRX = headerItems[1:4]
                energy = items[0::6]
                spin   = items[1::6]
                GT     = items[2::6]
                GN     = items[3::6]
                GG     = items[4::6]
                GF     = items[5::6]
                for i, E in enumerate(energy):
                    resonance = BreitWigner()
                    resonance.QX  = QX
                    resonance.L   = L
                    resonance.LRX = LRX
                    resonance.E   = energy[i]
                    resonance.J   = spin[i]
                    resonance.GT  = GT[i]
                    resonance.GN  = GN[i]
                    resonance.GG  = GG[i]
                    resonance.GF  = GF[i]
                    res.resonances.append(resonance)

        # Reich-Moore
        elif res.LRF == 3:
            # Read energy-dependent scattering radius if present
            if res.NRO > 0:
                res.AP = self.getTab1Record()
        
            # Other scatter radius parameters
            items = self.getContRecord()
            res.SPI  = items[0] # Spin, I, of the target nucleus
            if res.NRO == 0:
                res.AP = items[1]
            res.LAD  = items[3] # Flag for angular distribution
            res.NLS  = items[4] # Number of l-values
            res.NLSC = items[5] # Number of l-values for convergence
            
            # Read resonance widths, J values, etc
            for l in range(res.NLS):
                headerItems, items = self.getListRecord()
                APL, L = headerItems[1:3]
                energy = items[0::6]
                spin   = items[1::6]
                GN     = items[2::6]
                GG     = items[3::6]
                GFA    = items[4::6]
                GFB    = items[5::6]
                for i, E in enumerate(energy):
                    resonance = ReichMoore()
                    resonance.APL = APL
                    resonance.L   = L
                    resonance.E   = energy[i]
                    resonance.J   = spin[i]
                    resonance.GN  = GN[i]
                    resonance.GG  = GG[i]
                    resonance.GFA = GFA[i]
                    resonance.GFB = GFB[i]
                    res.resonances.append(resonance)

    def readUnresolved(self, res):
        pass

    def getTextRecord(self, line=None):
        if not line:
            line = self.fh.readline()
        HL  = line[0:66]
        MAT = int(line[66:70])
        MF  = int(line[70:72])
        MT  = int(line[72:75])
        NS  = int(line[75:80])
        return [HL, MAT, MF, MT, NS]

    def getContRecord(self, line=None, skipC=False):
        if not line:
            line = self.fh.readline()
        if skipC:
            C1 = None
            C2 = None
        else:
            C1 = convert(line[:11])
            C2 = convert(line[11:22])
        L1  = int(line[22:33])
        L2  = int(line[33:44])
        N1  = int(line[44:55])
        N2  = int(line[55:66])
        MAT = int(line[66:70])
        MF  = int(line[70:72])
        MT  = int(line[72:75])
        NS  = int(line[75:80])
        return [C1, C2, L1, L2, N1, N2, MAT, MF, MT, NS]

    def getHeadRecord(self, line=None):
        if not line:
            line = self.fh.readline()
        ZA  = int(convert(line[:11]))
        AWR = convert(line[11:22])
        L1  = int(line[22:33])
        L2  = int(line[33:44])
        N1  = int(line[44:55])
        N2  = int(line[55:66])
        MAT = int(line[66:70])
        MF  = int(line[70:72])
        MT  = int(line[72:75])
        NS  = int(line[75:80])
        return [ZA, AWR, L1, L2, N1, N2, MAT, MF, MT, NS]

    def getListRecord(self, onlyList=False):
        # determine how many items are in list
        items = self.getContRecord()
        NPL = items[4]
        
        # read items
        itemsList = []
        m = 0
        for i in range((NPL-1)/6 + 1):
            line = self.fh.readline()
            toRead = min(6,NPL-m)
            for j in range(toRead):
                val = convert(line[0:11])
                itemsList.append(val)
                line = line[11:]
            m = m + toRead
        if onlyList:
            return itemsList
        else:
            return (items, itemsList)

    def getTab1Record(self):
        r = ENDFTab1Record()
        r.read(self.fh)
        return r

    def findFile(self, fileNumber):
        for f in self.files:
            if f.fileNumber == fileNumber:
                return f

    def findMT(self, MT):
        for f in self.files:
            for r in f.reactions:
                if r.MT == MT:
                    return r

    def lookForFiles(self):
        files = set()
        self.fh.seek(0)
        for line in self.fh:
            try:
                fileNum = int(line[70:72])
                if fileNum == 0:
                    continue
            except ValueError:
                raise
            except IndexError:
                raise
            files.add(fileNum)
        print(files)

    def seekFile(self, fileNum):
        self.fh.seek(0)
        fileString = '{0:2}'.format(fileNum)
        while True:
            position = self.fh.tell()
            line = self.fh.readline()
            if line == '':
                # Reached EOF
                print('Could not find File {0}'.format(fileNum))
            if line[70:72] == fileString:
                self.fh.seek(position)
                break

    def printInfo(self, MT):
        if self.verbose:
            print("   MT={0} {1}".format(MT, MTname[MT]))

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
        C1 = convert(line[:11])
        C2 = convert(line[11:22])
        L1  = int(line[22:33])
        L2  = int(line[33:44])
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
                x = convert(line[:11])
                y = convert(line[11:22])
                self.x.append(x)
                self.y.append(y)
                line = line[22:]
            m = m + toRead

    def plot(self):
        plt.plot(self.x, self.y)

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
        HL  = line[0:66]
        MAT = int(line[66:70])
        MF  = int(line[70:72])
        MT  = int(line[72:75])
        NS  = int(line[75:80])
        self.items = [HL, MAT, MF, MT, NS]

class ENDFContRecord(ENDFRecord): 
    """
    An ENDFContRecord is a control record.
    """

    def __init__(self, fh):
        super(ENDFContRecord, self).__init__(fh)

    def read(self, line):
        C1  = convert(line[:11])
        C2  = convert(line[11:22])
        L1  = int(line[22:33])
        L2  = int(line[33:44])
        N1  = int(line[44:55])
        N2  = int(line[55:66])
        MAT = int(line[66:70])
        MF  = int(line[70:72])
        MT  = int(line[72:75])
        NS  = int(line[75:80])
        self.items = [C1, C2, L1, L2, N1, N2, MAT, MF, MT, NS]

class ENDFHeadRecord(ENDFRecord): 
    """
    An ENDFHeadRecord is the first in a section and has the same form as a
    control record, except that C1 and C2 fields always contain ZA and AWR,
    respectively.
    """

    def __init__(self, fh):
        super(ENDFHeadRecord, self).__init__(fh)

    def read(self, line):
        ZA  = int(convert(line[:11]))
        AWR = convert(line[11:22])
        L1  = int(line[22:33])
        L2  = int(line[33:44])
        N1  = int(line[44:55])
        N2  = int(line[55:66])
        MAT = int(line[66:70])
        MF  = int(line[70:72])
        MT  = int(line[72:75])
        NS  = int(line[75:80])
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
        self.fileNumber = 7

class ENDFFile8(ENDFFile):
    """
    File8 contains radioactive decay data.
    """

    def __init__(self):
        self.fileNumber = 8

class ENDFFile9(ENDFFile):
    """
    File9 contains multiplicites for production of radioactive elements.
    """

    def __init__(self):
        self.fileNumber = 9

class ENDFFile10(ENDFFile):
    """
    File10 contains cross sections for production of radioactive nuclides.
    """

    def __init__(self):
        self.fileNumber = 10

class ENDFReaction(ENDFFile):
    """A single MT record on an ENDF file."""

    def __init__(self, MT):
        self.MT = MT

    def __repr__(self):
        return "<ENDF Reaction: MT={0}, {1}>".format(self.MT, MTname[self.MT])

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


def convert(string):
    """
    This function converts a number listed on an ENDF tape into a float or int
    depending on whether an exponent is present.
    """

    if string[-2] in '+-':
        return float(string[:-2] + 'e' + string[-2:])
    else:
        return eval(string)

MTname = {451: "Desciptive Data",
          452: "Total Neutrons per Fission",
          455: "Delayed Neutron Data",
          456: "Prompt Neutrons per Fission",
          458: "Energy Release Due to Fission",
          460: "Delayed Photon Data",
          151: "Resonance Parameters",
          1: "(n,total) Neutron total",
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
          60: "(z,n10) Production of n, 10th excited state"}
