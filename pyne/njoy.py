#!/usr/bin/env python

"""
Purpose:
  Automatic generation of Njoy input data, including dragr data.
  Generation of DRAGLIB and ACELIB.

Copyright:
  This code was originally distirbuted under the LGPL license (below).
  However, we have been given written permission from the author Alain
  Herbert to redistribute it under PyNE's BSD license.
  
Original Copyright:
  Copyright (C) 2003 Ecole Polytechnique de Montreal
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version

.. moduleauthor:: A. Hebert, R. Karthikeyan
"""

import os
import time

class PyNjoyError(Exception):
    """Exception indicating an error in PyNjoy."""
    pass

class Njoy99(object):
    def __init__(self):
        self.iwt = 4
        self.legendre = 1
        self.legendregg = 6
        self.scatteringLaw = None
        self.eFiss = None
        self.branchingNG = None
        self.branchingN2N = None
        self.gstr = 0
        self.oldlib = None
        self.purr = None

    def pendf(self, eaf=0):
        """Generate a pointwise PENDF file from the ENDF raw data using the
        MODER, RECONR, BROADR, PURR (if dilutions present), and THERMR
        modules. This is the starting point for all other datatype generations
        including DRAGLIB, ACE, WIMSD, etc.

        Parameters
        ----------
        eaf : int
            If eaf is 1, simplified processing is performed to be compatible
            with the EAF nuclear library.

        """

        print " --- make pendf for " + self.hmat + " ---"
        myCwd = os.getcwd()
        myNjoy = myCwd + '/' + self.execDir + "/xnjoy<file_data"
        if not os.path.isfile(os.path.expandvars(self.evaluationFile)):
            raise PyNjoyError("evaluation file " + self.evaluationFile + 
                              " not found")
        if not os.path.isdir(self.evaluationName):
            os.mkdir(self.evaluationName)
        os.chdir(self.evaluationName)
        textDil=""
        if self.dilutions:
            nbDil = len(self.dilutions)
            if nbDil > 10:
                raise PyNjoyError("cannot have more than 10 dilutions")
            for dil in self.dilutions:
                textDil += " %E" % dil
        else:
            nbDil = 0
        nbTmp = len(self.temperatures)
        if nbTmp > 10:
            raise PyNjoyError("cannot have more than 10 temperatures")
        textTmp=""
        for tmp in self.temperatures:
            textTmp += " %E" % tmp

        matsab_inc = 221
        nbAtoms = 1
        elasOpt = 0
        if self.scatteringLaw:
            if self.scatteringMat == 1:
                # Hydrogen in H20
                nbAtoms = 2
                matsab_inc = 222
            elif self.scatteringMat == 7:
                # Hydrogen in ZrH
                nbAtoms = 2
                elasOpt = 1
                matsab_inc = 225
            elif self.scatteringMat == 11:
                # Hydrogen in D20
                nbAtoms = 2
                matsab_inc = 228
            elif self.scatteringMat == 26:
                # Beryllium metal
                elasOpt = 1
                matsab_inc = 231
            elif self.scatteringMat == 27:
                # Beryllium in BeO 
                nbAtoms = 2
                elasOpt = 1
                matsab_inc = 233
            elif self.scatteringMat == 31:
                # Graphite
                elasOpt = 1
                matsab_inc = 229
            elif self.scatteringMat == 37:
                # Hydrogen in CH2
                nbAtoms = 2
                matsab_inc = 223
            elif self.scatteringMat == 40:
                # Hydrogen in C6H6 (benzine)
                nbAtoms = 2
                matsab_inc = 227
            elif self.scatteringMat == 58:
                # Zirconium in ZrH
                nbAtoms = 2
                elasOpt = 1
                matsab_inc = 235
            unitLaw = -27
            matLaw = self.scatteringMat
            typeLaw = 4
            os.system("ln -s " + self.scatteringLaw + " tape26")
        else:
            unitLaw = 0
            matLaw = 0
            typeLaw = 1
            nbAtoms = 1
            elasOpt = 0
        htime = time.ctime(time.time())
        self.__dict__.update({"textDil": textDil, "nbDil": nbDil,
                              "textTmp": textTmp, "nbTmp": nbTmp,  
                              "unitLaw": unitLaw, "matLaw": matLaw,
                              "typeLaw": typeLaw, "nbAtoms": nbAtoms,
                              "elasOpt": elasOpt, "htime": htime,
                              "matsab_inc" : matsab_inc})

        if self.scatteringLaw:
            text_data = "moder\n20 -21\nmoder\n26 -27\n"
        else:
            text_data = "moder\n20 -21\n"
        text_data += (
            "reconr\n-21 -22\n'pendf tape from %(evaluationName)s'/\n"
            "%(mat)d 1/\n0.001  0.  0.005/\n"
            "'%(hmat)s from %(evaluationName)s at %(htime)s' /\n0/\n"
            "broadr\n-21 -22 -23\n%(mat)d %(nbTmp)d/\n0.001/\n"
            "%(textTmp)s/\n0/\n" % self.__dict__)
        if self.dilutions and self.purr:
            text_data +=  (
                "purr\n-21 -23 -24\n%(mat)d %(nbTmp)d %(nbDil)d 20 32/\n"
                "%(textTmp)s/\n%(textDil)s/\n0/\n" % self.__dict__)
        elif self.dilutions:
            text_data +=  (
                "unresr\n-21 -23 -24\n%(mat)d %(nbTmp)d %(nbDil)d 1/\n"
                "%(textTmp)s/\n%(textDil)s/\n0/\n" % self.__dict__)
        if self.dilutions:
            text_data +=  (
                "thermr\n0 -24 -35\n0 %(mat)d 16 %(nbTmp)d %(typeLaw)d 0 "
                "%(nbAtoms)d 221 0\n%(textTmp)s/\n0.001 4.0\nmoder\n-35 29\n"
                "stop\n" % self.__dict__)
        else:
            if self.scatteringLaw:
                text_data +=  (
                    "thermr\n%(unitLaw)d -23 -35\n%(matLaw)d %(mat)d 16 "
                    "%(nbTmp)d %(typeLaw)d %(elasOpt)d %(nbAtoms)d "
                    "%(matsab_inc)d 0/\n%(textTmp)s/\n0.001 4.0\n"
                    "moder\n-35 29\nstop\n" % self.__dict__)
            elif eaf == 0:
                text_data +=  (
                    "thermr\n0 -23 -35\n0 %(mat)d 16 %(nbTmp)d "
                    "%(typeLaw)d 0 %(nbAtoms)d 221 0\n%(textTmp)s/\n"
                    "0.001 4.0\nmoder\n-35 29\nstop\n" % self.__dict__)
            else:
                text_data +=  (
                    "moder\n-23 29\nstop\n" % self.__dict__)
        file_data = open("file_data",'w')
        file_data.write(text_data)
        file_data.close()
        os.system("ln -s " + self.evaluationFile + " tape20")
        os.system(myNjoy)
        os.system("mv tape29 pendf" + self.hmat)
        os.system("mv file_data file_data_pendf" + self.hmat)
        os.system("mv output out_pendf_" + self.hmat)
        os.system("chmod 644 out_pendf_" + self.hmat)
        for fileName in os.listdir(os.getcwd()):
            if fileName[:4] == 'tape':
                os.remove(fileName)
        os.chdir(myCwd)

    def gendf(self, eaf=0):
        """Generate a multigroup GENDF file using the MODER and GROUPR
        modules. This requires that the pendf() method has already been called.

        Parameters
        ----------
        eaf : int
            If eaf is 1, simplified processing is performed to be compatible
            with the EAF nuclear library.

        """

        print " --- make gendf for " + self.hmat + " ---"
        myCwd = os.getcwd()
        myNjoy = myCwd + '/' + self.execDir + "/xnjoy<file_data"
        if not os.path.isfile(os.path.expandvars(self.evaluationFile)):
            raise PyNjoyError("evaluation file " + self.evaluationFile +
                              " not found")
        os.chdir(self.evaluationName)
        matsab_inc = 221
        matsab_coh = 0
        if self.scatteringLaw:
            if self.scatteringMat == 1:
                matsab_inc = 222
            elif self.scatteringMat == 7:
                matsab_inc = 225        
                matsab_coh = 226        
            elif self.scatteringMat == 11:
                matsab_inc = 228
            elif self.scatteringMat == 26:
                matsab_inc = 231
                matsab_coh = 232
            elif self.scatteringMat == 27:
                matsab_inc = 233
                matsab_coh = 234
            elif self.scatteringMat == 31:
                matsab_inc = 229
                matsab_coh = 230
            elif self.scatteringMat == 37:
                matsab_inc = 223
                matsab_coh = 224
            elif self.scatteringMat == 40:
                matsab_inc = 227
            elif self.scatteringMat == 58:
                matsab_inc = 235             
                matsab_coh = 236             
        if self.iwt:
            newiwt = self.iwt
        else:
            newiwt = 4
        if self.dilutions:
            newiwt = -abs(newiwt)
            nbDil = len(self.dilutions)
            if nbDil > 10:
                raise PyNjoyError("cannot have more than 10 dilutions")
            textDil = ""
            for dil in self.dilutions:
                textDil += " %E" % dil
        else:
            newiwt = abs(newiwt)
            nbDil = 1
            textDil = "1.0e10"
        nbTmp = len(self.temperatures)
        if nbTmp > 10:
            raise PyNjoyError("cannot have more than 10 temperatures")
        textTmp = ""
        for tmp in self.temperatures:
            textTmp += " %E" % tmp
        htime = time.ctime(time.time())
        self.__dict__.update({"textDil": textDil,
                              "nbDil"  : nbDil,
                              "textTmp": textTmp,
                              "nbTmp"  : nbTmp,
                              "htime"  : htime,
                              "newiwt" : newiwt,
                              "matsab_inc" : matsab_inc,
                              "matsab_coh" : matsab_coh,
                              "autosup": self.autolib[1]})

        text_data = ("moder\n20 -21\nmoder\n29 -25\ngroupr"
                     "-21 -25 0 -26\n%(mat)d %(nstr)d %(gstr)d %(newiwt)d "
                     "%(legendre)d %(nbTmp)d %(nbDil)d 1\n'%(hmat)s from "
                     "%(evaluationName)s at %(htime)s' /\n%(textTmp)s/\n"
                     "%(textDil)s/\n" % self.__dict__)
        if self.dilutions:
            text_data += ("%(autosup)f %(potential)f 20000 / "
                          "Homog. Flux Calc.Param" % self.__dict__)
        if newiwt == 1 or newiwt == -1:
            text_data += self.wght
        elif newiwt == 4 or newiwt == -4:
            text_data += "0.2 0.0253 820.3e3 1.40e6 / iwt=4 parameters"
        for tmp in self.temperatures:
            if eaf == 0:
                if matsab_coh != 0:
                    text_data += ("3/\n3 %(matsab_inc)d /\n"
                                  "3 %(matsab_coh)d /\n" % self.__dict__)
                else:
                    text_data += "3/\n3 %(matsab_inc)d /\n" % self.__dict__
            else:
              text_data += "3/\n"
            if self.fission:
                text_data += "3 452 /\n"
            if self.fission == 2:
                text_data += "3 455 /\n5 455 /\n"
            if self.gstr != 0:
                text_data += "16 / photon interaction matrices\n"
            if eaf == 0:
                if matsab_coh != 0:
                    text_data += ("6 /\n6 %(matsab_inc)d /\n6 "
                                  "%(matsab_coh)d /\n0/\n" % self.__dict__)
                else:
                    text_data += "6 /\n6 %(matsab_inc)d /\n0/\n" % self.__dict__
            else:
                text_data += "6 /\n0/\n"
        text_data += "0/\nmoder\n-26 30\nstop"
        file_data = open("file_data",'w')
        file_data.write(text_data)
        file_data.close()
        os.system("ln -s " + self.evaluationFile + " tape20")
        os.system("ln -s pendf" + self.hmat + " tape29")
        os.system(myNjoy)
        os.system("mv file_data file_data_gendf" + self.hmat)
        os.system("mv tape30 gendf" + self.hmat)
        os.system("mv output out_gendf_" + self.hmat)
        os.system("chmod 644 out_gendf_" + self.hmat)
        for fileName in os.listdir(os.getcwd()):
            if fileName[:4] == 'tape':
                os.remove(fileName)
        os.chdir(myCwd)

    def gamma(self):
        """Generate photo-atomic (gamma) group ENDF file using the MODER,
        RECONR, and GAMINR modules.
        """

        print " --- make gamma gendf for " + self.hmatgg + " ---"
        myCwd = os.getcwd()
        myNjoy = myCwd + '/' + self.execDir + "/xnjoy<file_data"
        if not os.path.isfile(os.path.expandvars(self.evaluationFile)):
          raise PyNjoyError("evaluation file " + self.evaluationFile +
                            " not found")
        if not os.path.isdir(self.evaluationName):
            os.mkdir(self.evaluationName)
        os.chdir(self.evaluationName)
        htime = time.ctime(time.time())
        self.__dict__.update({"htime"  : htime})

        text_data = ("moder\n40 -41\nreconr\n-41 -42\n"
                     "'pendf tape from %(evaluationName)s'/\n"
                     "%(matgg)d 1/\n0.002 /\n'%(hmatgg)s (gamma) from "
                     "%(evaluationName)s at %(htime)s' /\n0/\ngaminr\n"
                     "-41 -42 0 -43\n%(matgg)d %(gstr)d 3 %(legendregg)d 1\n"
                     "'%(hmatgg)s (gamma) from %(evaluationName)s at "
                     "%(htime)s' /\n-1 /\n0 /\nmoder\n-43 44\nstop\n" 
                     % self.__dict__)
        file_data = open("file_data", 'w')
        file_data.write(text_data)
        file_data.close()
        os.system("ln -s " + self.evaluationFile + " tape40")
        os.system(myNjoy)
        os.system("mv tape44 gamma" + self.hmatgg)
        os.system("mv file_data file_data_gamma" + self.hmatgg)
        os.system("mv output out_gamma_" + self.hmatgg)
        os.system("chmod 644 out_gamma_" + self.hmatgg)
        for fileName in os.listdir(os.getcwd()):
          if fileName[:4] == 'tape':
              os.remove(fileName)
        os.chdir(myCwd)

    def draglib(self, fp=0):
        """Generate a DRAGLIB file using the MODER and DRAGR modules and
        add/update the new isotopic data in the DRAGLIB file.

        Parameters
        ----------
        fp : int
            If fp is 1, the scattering information are stored as diagonal
            matrices in the DRAGLIB.

        """

        myCwd = os.getcwd()
        myNjoy = myCwd + '/' + self.execDir + "/xnjoy<file_data"
        if not os.path.isfile(os.path.expandvars(self.evaluationFile)):
          raise PyNjoyError("evaluation file " + self.evaluationFile +
                            " not found")
        evaluationNameBase = os.path.basename(self.evaluationName)
        if self.oldlib:
            os.system("mv ../" + self.oldlib + " draglib" + evaluationNameBase)
        os.chdir(self.evaluationName)
        if os.path.isfile("drag"): os.remove("drag")
        if os.path.isfile("draglib" + evaluationNameBase):
            iold = 29
            name1 = "draglib" + evaluationNameBase + ".bis.gz"
            if not os.path.isfile(name1):
                os.system("cp draglib" + evaluationNameBase +
                          " draglib" + evaluationNameBase + ".bis")
                os.system("gzip draglib" + evaluationNameBase + ".bis")
            else:
                name2 = "draglib" + evaluationNameBase + ".bis2.gz"
                os.system("cp draglib" + evaluationNameBase +
                          " draglib" + evaluationNameBase + ".bis2")
                os.system("gzip draglib" + evaluationNameBase + ".bis2")
                if os.path.isfile(name2):
                    len1 = os.stat(name1).st_size
                    len2 = os.stat(name2).st_size
                    if len2 > len1:
                        os.remove(name1)
                        os.rename(name2, name1)
                    else:
                        os.remove(name2)
            os.system("mv draglib" + evaluationNameBase + " tape29")
            print " append data for " + self.hmat + " to existing draglib file"
        else:
            iold = 0
            print " create a new draglib file for " + self.hmat
        htime = time.ctime(time.time())

        if self.dilutions:
            self.__dict__.update({"htime": htime,
                                  "iold" : iold,
                                  "ss0"  : self.ss[0],
                                  "ss1"  : self.ss[1],
                                  "auto0": self.autolib[0],
                                  "auto1": self.autolib[1],
                                  "auto2": self.autolib[2],
                                  "fp"   : fp})
            text_data = ("moder\n20 -21\nmoder\n22 -23\nmoder\n24 -25\n"
                         "dragr\n-21 -23 -25 0 0 %(iold)d 30 %(fp)d/\n" 
                         % self.__dict__)
            if iold == 0:
                text_data += ("'draglib from %(evaluationName)s at "
                              "%(htime)s'/\n" % self.__dict__)
            text_data += ("%(mat)d %(hmat)s /\n'%(hmat)s from "
                          "%(evaluationName)s (%(mat)d) at %(htime)s' /\n"
                          "%(ss0)E %(ss1)E /\n%(auto0)E %(auto1)E %(auto2)E /\n"
                          "0/\nstop\n" % self.__dict__)
        else:
            self.__dict__.update({"htime": htime, "iold": iold, "fp": fp})
            text_data = ("moder\n20 -21\nmoder\n24 -25\ndragr\n-21 0 -25 0 0 "
                         "%(iold)d 30 %(fp)d/\n" % self.__dict__)
            if iold == 0:
                text_data += ("'draglib from %(evaluationName)s at %(htime)s'/\n"
                              % self.__dict__)
            text_data += ("%(mat)d %(hmat)s /\n'%(hmat)s from "
                          "%(evaluationName)s (%(mat)d) at %(htime)s' /\n"
                          "0.1 1.0E10 /\n0/\nstop\n" % self.__dict__)
        file_data = open("file_data",'w')
        file_data.write(text_data)
        file_data.close()
        os.system("ln -s " + self.evaluationFile + " tape20")
        if self.dilutions:
            os.system("ln -s pendf" + self.hmat + " tape22")
        os.system("ln -s gendf" + self.hmat + " tape24")
        os.system(myNjoy)
        os.system("mv file_data file_data_dendf" + self.hmat)
        stats = os.stat("tape30") # get the stats of the file
        size = stats[6] # extract the file size in bytes from the stats list
        if size > 8 :
            os.system("mv tape30 draglib" + evaluationNameBase)
        else:
            os.system("mv tape29 draglib" + evaluationNameBase)
            raise PyNjoyError("draglib file for " + self.hmat + " not created")
        file_in = open('output','r')
        file_out = open("out_draglib_" + self.hmat,'w')
        while 1:
            line=file_in.readline()
            if not line: break
            ind=line.find('???????????')
            if ind != -1:
                if not self.eFiss:
                    raise PyNjoyError("self.eFiss instance variable not set")
                line = line[:ind] + "%E" % self.eFiss + line[ind+11:]
                self.eFiss = None
            if self.branchingNG:
                ind=line.find(' ng ')
                if ind != -1:
                  jnd = line[ind+3:].find(' 0.000 ')
                  if jnd == -1:
                      raise PyNjoyError("unable to set the isomeric ng "
                                        "branching ratio")
                  line = line[:ind+jnd+4] + "%5.3f" % self.branchingNG \
                      + line[ind+jnd+9:]
                  self.branchingNG = None
            if self.branchingN2N:
                ind=line.find(' n2n ')
                if ind != -1:
                    jnd = line[ind+4:].find(' 0.000 ')
                    if jnd == -1:
                        raise PyNjoyError("unable to set the isomeric n2n "
                                          "branching ratio")
                    line = line[:ind+jnd+5] + "%5.3f"%self.branchingN2N \
                        + line[ind+jnd+10:]
                    self.branchingN2N = None
            file_out.writelines(line)
        file_out.close()
        file_in.close()
        os.remove('output')
        os.system("chmod 644 out_draglib_" + self.hmat)
        for fileName in os.listdir(os.getcwd()):
            if fileName[:4] == 'temp' or fileName[:4] == 'tape':
                os.remove(fileName)
        os.chdir(myCwd)

    def matxs(self):
        """Generate an ASCII MATXS file using the MODER and MATXSR modules."""

        print " --- make matxs for " + self.hmat + " ---"
        myCwd = os.getcwd()
        myNjoy = myCwd + '/' + self.execDir + "/xnjoy<file_data"
        os.chdir(self.evaluationName)
        listGro = [0, 239, 30, 27, 50, 68, 100, 35, 69, 187, 70, 620, 80, 100,
                   640, 174, 175, 172, 33, 1968, 315, 172, 175, 281, 349, 89]
        listGro2 = [0, 94, 12, 21, 22, 48, 24, 36, 38, 42]
        htime = time.ctime(time.time())
        self.__dict__.update({"nbGro": listGro[self.nstr - 1],
                              "htime": htime})
        if self.gstr == 0:
            text_data = ("moder\n24 -25\nmatxsr\n-25 0 28/\n1 '%(hmat)s from "
                         "%(evaluationName)s (%(mat)d) at %(htime)s'/\n"
                         "1 2 1 1\n'neutron library'/\n'n'\n%(nbGro)d\n"
                         "'nscat' 'ntherm'/\n1 1\n1 1\n%(hmat)s %(mat)d /\n"
                         "stop\n" % self.__dict__)
        else:
            self.__dict__.update({"nbGro2": listGro2[self.gstr - 1]})
            text_data = ("moder\n24 -25\nmoder\n26 -27\nmatxsr\n-25 -27 28/\n"
                         "1 '%(hmat)s coupled-set from %(evaluationName)s "
                         "(%(mat)d+%(matgg)d) at %(htime)s'/\n2 3 1 1\n"
                         "'neutron-gamma library'/\n'n' 'g'\n%(nbGro)d "
                         "%(nbGro2)d\n'nscat' 'ng' 'gscat' 'ntherm'/\n1 1 2 1\n"
                         "1 2 2 1\n%(hmat)s %(mat)d %(matgg)d/\nstop\n" 
                         % self.__dict__)
            os.system("ln -s gamma" + self.hmatgg + " tape26")
        file_data = open("file_data", 'w')
        file_data.write(text_data)
        file_data.close()
        os.system("ln -s gendf" + self.hmat + " tape24")
        os.system(myNjoy)
        os.system("mv file_data file_data_matxs" + self.hmat)
        os.system("mv tape28 matxs" + self.hmat)
        os.system("mv output out_matxs_" + self.hmat)
        os.system("chmod 644 out_matxs_" + self.hmat)
        for fileName in os.listdir(os.getcwd()):
            if fileName[:4] == 'tape':
                os.remove(fileName)
        os.chdir(myCwd)

    def makeFp(self, eaf=0):
        """Creates a PENDF, GENDF, and DRAGLIB file for a single fission
        product.

        Parameters
        ----------
        eaf : int
            If eaf is 1, simplified processing is performed to be compatible
            with the EAF nuclear library.

        """

        self.scatteringLaw = None
        self.fission = None
        self.dilutions = None
        keeplegendre = self.legendre
        self.legendre = 0
        self.pendf(eaf)
        self.gendf(eaf)
        self.draglib(fp=1)
        self.legendre = keeplegendre

    def burnup(self):
        """Process burnup data for the complete library. This requires a file
        whose name starts with chain, e.g. chaincandu, that contains information
        about the energy from all isotopes generated using single DRAGR
        runs. The 'chain' file is generated automatically.

        """

        myCwd = os.getcwd()
        myNjoy = myCwd + '/' + self.execDir + "/xnjoy<tempFile"
        evaluationNameBase =  os.path.basename(self.evaluationName)
        os.chdir(self.evaluationName)
        if os.path.isfile("drag"):
            os.remove("drag")
        if os.path.isfile("draglib" + evaluationNameBase):
            iold = 29
            name1 = "draglib" + evaluationNameBase + ".bis.gz"
            if not os.path.isfile(name1):
                os.system("cp draglib" + evaluationNameBase +
                          " draglib" + evaluationNameBase + ".bis")
                os.system("gzip draglib" + evaluationNameBase + ".bis")
            else:
                name2 = "draglib" + evaluationNameBase + ".bis2.gz"
                os.system("cp draglib" + evaluationNameBase +
                          " draglib" + evaluationNameBase + ".bis2")
                os.system("gzip draglib" + evaluationNameBase + ".bis2")
                if os.path.isfile(name2):
                    len1 = os.stat(name1).st_size
                    len2 = os.stat(name2).st_size
                    if len2 > len1:
                        os.remove(name1)
                        os.rename(name2, name1)
                    else:
                        os.remove(name2)
            os.system("mv draglib" + evaluationNameBase + " tape29")
        else:
            iold = 0
        htime = time.ctime(time.time())
        self.__dict__.update({"iold": iold, "htime": htime})
        text_data = "dragr\n0 0 0 23 24 %(iold)d 30/\n" % self.__dict__
        if iold == 0:
          text_data += ("'draglib from %(evaluationName)s at "
                        "%(htime)s'/\n" % self.__dict__)
        file_data = open("file_data", 'w')
        file_data.write(text_data)
        file_data.close()
        yield_file = os.path.expandvars(self.fissionFile)
        if os.path.isdir(yield_file):
            tape23 = open("tape23",'w')
            tape23.write("LIBRARY, DUMMY TAPE HEADER\n")
            for fileName in os.listdir(yield_file):
                file = open(yield_file + fileName, 'r')
                lines = file.readlines()
                n = len(lines)
                tape23.writelines(lines[1:n-1])
                file.close()
            tape23.writelines(lines[n-1:n])
            tape23.close()
        else:
            os.system("ln -s " + self.fissionFile + " tape23")
        os.system("ln -s " + self.decayFile + " tape24")
        chainFileName = 'chain' + evaluationNameBase
        listFiles = os.listdir(os.getcwd())
        if chainFileName not in listFiles:
            print 'Make the burnup chain file named',chainFileName
            data_dict = {}
            mat_dict = {}
            for fileName in listFiles:
                if fileName[:11] == 'out_draglib':
                    file = open(fileName,'r')
                    while True:
                        line = file.readline()
                        if not line:
                            break
                        if (line[:41] == 
                            ' isotopic specification line for material'):
                            mat = line[42:49]
                            line = file.readline()
                            line = file.readline()
                            info = ''
                            key = line[:8]
                            while line[:6] != ' -----':
                                info = info + line
                                line = file.readline()
                            data_dict[key] = info[:-1]
                            mat_dict[key] = 10*int(mat)
                            if key.find('_') != -1:
                                mat_dict[key] = mat_dict[key] - 1
            dictKeys = data_dict.keys()
            dictKeys.sort(lambda a, b: mat_dict[a] - mat_dict[b])
            chainFile = open(chainFileName, 'w')
            for key in dictKeys:
                line = data_dict[key]
                pos = 0
                while pos != -1:
                    pos = line.find('\n')
                    chainFile.write(line[:pos-1] + '\n')
                    line = line[pos+1:]
            chainFile.write("end /\n")
            chainFile.write("stop\n")
            chainFile.close()
        else:
            print 'Use existing burnup chain file named',chainFileName
        os.system("cat file_data " + chainFileName + " > tempFile")
        os.system(myNjoy)
        os.system("mv tape30 draglib" + evaluationNameBase)
        os.system("mv output out_draglib_burnup")
        os.system("chmod 644 out_draglib_burnup")
        os.system("mv tempFile file_data_burnup")
        for fileName in os.listdir(os.getcwd()):
            if fileName[:4] == 'temp' or fileName[:4] == 'tape':
                os.remove(fileName)
        os.remove("file_data")
        os.chdir(myCwd)

  def acer(self):
    myCwd = os.getcwd()
    myNjoy = myCwd + '/' + self.execDir + "/xnjoy<file_data"
    if not os.path.isfile(os.path.expandvars(self.evaluationFile)):
      raise PyNjoyError("evaluation file " + self.evaluationFile + " not found")
    textTmp=""
    for tmp in self.tempace:
      textTmp = textTmp + " " + "%E"%tmp
    os.chdir(self.evaluationName)
    self.__dict__.update({"textTmp": textTmp})
    text_data = """
    moder
    20 -21
    moder
    29 -25
    acer
    -21 -25 0 38 39
    1 0 1 %(suff)f/
    'pendf tape from %(evaluationName)s'/
    %(mat)d  %(textTmp)s /
    1 1/
    0.001/
    acer / Check ACE files
    0 38 0 40 41
    7 1 1 -1/
    /
    """%self.__dict__
    if self.scatteringLaw:
      matsab_inc = 221
      matsab_coh = 0
      nbAtoms = 1
      elasOpt = 0
      if self.scatteringMat == 1:    # H1_H20
        nbAtoms = 2
        matsab_inc = 222               # h2o
      elif self.scatteringMat == 7:  # H_ZRH
        nbAtoms = 2
        elasOpt = 1
        matsab_inc = 225               # zrhyd
        matsab_coh = 226               # zrhyd$
      elif self.scatteringMat == 11: # H2_D20
        nbAtoms = 2
        matsab_inc = 228               # d2o
      elif self.scatteringMat == 26: # Be
        elasOpt = 1
        matsab_inc = 231               # be
        matsab_coh = 232               # be$
      elif self.scatteringMat == 27: # BeO 
        nbAtoms = 2
        elasOpt = 1
        matsab_inc = 233               # beo
        matsab_inc = 234               # beo$
      elif self.scatteringMat == 31: # C_GRAPH
        elasOpt = 1
        matsab_inc = 229               # graph
        matsab_inc = 230               # graph$
      elif self.scatteringMat == 37: # H_CH2
        nbAtoms = 2
        matsab_inc = 223               # poly
      elif self.scatteringMat == 40: # H_C6H6 (benzine)
        nbAtoms = 2
        matsab_inc = 227               # benz
      elif self.scatteringMat == 58: # Zr_ZRH
        nbAtoms = 2
        elasOpt = 1
        matsab_inc = 235               # zrhyd
        matsab_inc = 236               # zrhyd$
      text_data = text_data + """
    acer
    -21 -25 0 48 49
    2 0 1 %(suff)f/
    'pendf tape from %(evaluationName)s'/
    %(mat)d  %(textTmp)s %(scatName)s/
    %(za)d 0 0 /
    %(matsab_inc)d 16 %(matsab_coh)d %(elasOpt)d %(nbAtoms)d 4.0 0 0 /
    acer / Check ACE files
    0 48 0 50 51
    7 1 1 -1/
    /
    stop
    """%self.__dict__
    else:
      text_data = text_data + """
    stop
    """
    file_data = open("file_data",'w')
    file_data.write(text_data)
    file_data.close()
    os.system("ln -s " + self.evaluationFile + " tape20")
    os.system("ln -s pendf" + self.hmat + " tape29")
    os.system(myNjoy)
    if os.path.isfile("tape48"):
     inp=open("tape48","r")
     outp=open("acecandu",'a')
     line = inp.readlines()
     outp.writelines(line)
     inp1=open("tape49","r")
     outp1=open("acexsdir",'a')
     line = inp1.readlines()
     outp1.writelines(line)
    if os.path.isfile("tape38"):
     inp=open("tape38","r")
     outp=open("acecandu",'a')
     line = inp.readlines()
     outp.writelines(line)
     inp1=open("tape39","r")
     outp1=open("acexsdir",'a')
     line = inp1.readlines()
     outp1.writelines(line)
     print " ace file for" + self.hmat + "created"
     if not os.path.isdir(self.dirName): os.mkdir(self.dirName)
     os.system("mv tape38 " + '   ' + self.dirName + "/"+self.hmat +'_'+str(int(self.tempace[0]))+".ace")
    else: 
     raise PyNjoyError("ace file for " + self.hmat + " not created")
#
    for fileName in os.listdir(os.getcwd()):
      if fileName[:4] == 'tape': os.remove(fileName)
    os.remove("file_data")
    os.chdir(myCwd)
