.. _usersguide_njoy:

===============
NJOY Automation
===============

.. currentmodule:: pyne.njoy

.. automodule:: pyne.njoy

The purpose of this module is to enable simple generation of cross section
libraries and avoiding the need to write NJOY input files by hand. A complete
processing of an evaluated datafile is done by invoking methods of the Njoy99
class. This class contains the various instance variables and methods required
to use NJOY in the simplest possible way.

For a complete specification for the classes in the ``njoy`` module, please
refer to the Library Reference entry for :ref:`pyne_njoy`.

-----------------
Njoy99 Attributes
-----------------

The generation of cross section libraries is controlled by a number of
attributes of the Njoy99 class. These attributes are passed as input parameters
to NJOY for setting options, temperatures, tolerances, etc.

self.evaluationName
  Path of directory where you want to store the pendf, gendf, draglib and acelib
  files. The path can be prefixed by /tmp/ to force the files to be created
  locally on the /tmp directory.

self.legendre
  Order of Legendre polynomials for neutrons (= 1 for linear anisotropy in LAB -
  default)

self.legendregg
  Order of Legendre polynomials for gamma particles (default: = 6)

self.nstr
  Option for a particular neutron group structure (= 22 for the XMAS 172–group
  structure)

self.gstr
  Option for a particular gamma group structure, for producing neutron–gamma
  coupled sets (equal to zero by default)

self.iwt
  Type of flux weighting in groupr (= 1 for user-defined spectra; = 3 for 1/E
  weighting; = 4 recommended/default)

self.wght
  User-defined weighting spectra to be used if self.iwt = 1. Given as
  """–delimited string.

self.autolib
  Three-component tuple containing the energy limits of the autolibs (must
  correspond to energy-group limits) and the elementary lethargy width of the
  autolibs

self.temperatures
  Value of temperatures at which the cross sections are generated

self.hmat
  Material name that is included in the DRAGLIB – User dependent

self.mat
  mat number of nuclide under consideration

self.hmatgg
  Photo-atomic element name that is included in the DRAGLIB – User dependent

self.matgg
  Photo-atomic mat number of element under consideration

self.za
  ZA number, mainly required for generation of S(α, β) cross section in acer

self.scatName
  Name of S(α, β) cross section identifier for inclusion in xsdir

self.suff
  The suffix to be attached to nuclide in ACELIB - User dependent

self.evaluationFile
  Path of the evaluated datafile.

self.scatteringLaw
  Path of file having thermal scattering data (default = None)

self.scatteringMat
  mat number in scattering data file

self.fission
  Choice for including delayed neutron fission data in groupr module

self.ss
  Two-component tuple containing energy limits in eV for the self-shielding
  domain

self.potential
  Value of the potential cross section used in the flux calculator of groupr

self.dilutions
  Tuple containing the dilution values that need to be considered for
  calculation of resonance integrals and probability tables

self.dirName
  Directory name to store data for independent verification of ACELIB.

self.tempace
  Temperature at which ACELIB needs to be generated.

self.eFiss
  Fission energy in MeV. Used in cases where this value is not available in the
  evaluation.

self.branchingNG
  Radiative capture isomeric branching ratio (default = None). If you use this
  value, don’t forget to reset it to None after the isotope is completed.

self.branchingN2N
  N2N isomeric branching ratio (default = None). If you use this value, don’t
  forget to reset it to None after the isotope is completed.

self.purr
  Set to 1 to use purr module. By default, use unresr.

self.oldlib
  Name of an existing DRAGLIB file that is modified by self.draglib().

------------------
DRAGLIB Generation
------------------

The modules that will be used in the generation of DRAGLIB are MODER, RECONR,
BROADR, PURR (if dilutions present), THERMR, GROUPR and DRAGR. Figure 3 gives
the flow chart for generation of DRAGLIB formatted library. It is important to
identify the nuclides that are needed to be included as part of library and
corresponding evaluated datafiles from respective data centers need to be
compiled in a particular directory. This will help in cross verification and
assessment at any stage of library generation.  In this section specific
examples of elements will be provided to understand the nuances of DRAGLIB
library generation. The examples will be such that all the reactor type
materials will be chosen. They are scattering material (heavy water), structural
material (Zirconium), fission product (Xe-135), Actinide (U-235). A special
example for burnup dependent data will also be provided. For isotopes that have
resonances and whose presence in fuel can alter the flux in the energy region
between 2.76792 eV and 1.66156e4 eV, cross sections are generated at specific
dilution values. The choice of dilution values is shown in Table-1 and shown in
Figure 4. The value for potential scattering cross section is obtained using the
Fortran code getmf2.f which is provided by IAEA. Using this code, and the
evaluated datafile for the particular element, one can obtain the value for
self.potential. In general we cannot provide more than ten temperature values
and ten dilution values. If one has to generate cross sections for more than ten
dilution values, it has to be split, as shown in example for U-235. After each
"self.dragr()" run, one will obtain a file "out draglib elementname". Energy
information is recovered from this ascii file by method self.burnup() and is
collected in a file named 'chain' + self.evaluationName (stored on directory
self.evaluationName). Sometimes, it is likely that fission energy is not
provided in the tape.  In that case one obtains "????????" in the "out draglib
elementname" file. In that case, one has to obtain the energy value from some
other source (typically, from another evaluation) and provide it using the
self.eFiss instance variable, even if one does not use the tape for generation
of multigroup data.

File 8 of ENDF evaluation contains half-lives, decay modes, decay energies, and
radiation spectra for most isotopes. Information concerning the decay of the
reaction products is also given in this file.  In addition, fission product
yield data (MT=454 and 459) for fissionable materials and spontaneous
radioactive decay data (MT=457) for the nucleus are included. 

File 8 information is processed by module DRAGR. A large number of fission
products are included in the evaluated file for each element capable of
undergoing fission. For example, in the fission product yield data file included
in ENDF/B-VI rel. 8, one can notice that there are information of 1232 fission
products for 0.0253 eV fission of 233 U, 1247 fission products for 0.0253 eV
fission of 235 U etc. But the evaluations are not available for all the
nuclides, as most of them have very short half-lives and in the reactor context,
can be considered insignificant. They are subsequently lumped by a procedure
that is built in DRAGR. If there are nuclides with long half lives, but are not
available as evaluated files, a warning is provided before lumping the
corresponding element. The DRAGR user has the complete control over the lumping
process. DRAGR currently has no capability to produce pseudo fission product,
i.e., custom library isotopes made from the combination of many minor ENDF
fission products. All the isotopes missing in the ’chain’
+ self.evaluationName file are tested against a lumping criterion and are lumped. The criterion for
lumping a depleting isotope is a half life less than thirty days and a fission yield less than 0.01%. If this
criterion is not met, this isotope is lumped and a warning message is issued.

Information on energies for various reaction types like (n, :math:`\gamma`), (n,
f), (n, 2n), (n, 3n), (n, 4n), (n, :math:`\alpha`), (n, p), (n, 2
:math:`\alpha`), (n, np), (n, d), (n, t) are recovered from earlier DRAGR
single-isotope calculations and used for inclusion in relevant depletion data in
DRAGLIB format. The fission energy (n, f) is obtained from MF1 MT458 and the
energy from delayed betas and gammas are subtracted from it. Information
regarding energies for other reactions are derived by DRAGR from MF3. The
corresponding MT numbers for the above mentioned reactions (other than (n, f))
are 102, 16, 17, 37, 107, 103, 108, 28, 104, 105 respectively.  The complete
information required to do the depletion calculations is provided in ten
specific records of the DRAGLIB file.

Heavy Water
-----------

Heavy water is used in CANDU reactors as moderator and coolant. So it is important to generate
consistent data for efficient analysis of CANDU lattices. The instance variables and methods that will
be used are

.. code-block:: python

    candu.hmat = "H2_D2O"
    candu.temperatures = (293.6, 323.6, 573.6)
    candu.mat = 128
    candu.evaluationFile = "/home/user/Tripoli4/JEF2/eval/jef2.neutron.H2.bcd"
    candu.scatteringLaw = "/home/user/Tripoli4/JEF2/eval/jef2.neutron.D_D2O.bcd.therm"
    candu.scatteringMat = 11
    candu.fission = None
    candu.dilutions = None
    candu.pendf()
    candu.gendf()
    candu.draglib()

Zirconium
---------

Zirconium is used in CANDU reactors as cladding material and also for pressure
tube and calandria tube. Zirconium has some resonances and as a result of this
it is important to generate the cross sections of Zirconium for certain dilution
values. The instance variables and methods that will be used are

.. code-block:: python

    candu.hmat = "Zr0"
    candu.temperatures = (293.6, 323.6, 573.6)
    candu.mat = 4000
    candu.evaluationFile = "/home/user/Njoy99/evaluations/database/jendl306.asc"
    candu.fission = None
    candu.ss = (2.76792, 1.66156e4)
    candu.potential = 6.5144
    candu.dilutions = (1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, 158.676849, \
        56.3173141, 19.9880447, 7.09412289, 2.51783395)
    candu.pendf()
    candu.gendf()
    candu.draglib()

Xenon-135
---------

Xenon is a very important fission product in nuclear reactors. It is very
important to estimate the number density of this nuclide as a function of
burnup. This will help in estimating reactivity variations due to change in
concentration of the nuclide. By using makeFp option, we avoid generating
scattering matrices in (NG x NG) format, where NG is number of groups. We
instead generate only the scattering matrices along the diagonal.

.. code-block:: python

    candu.hmat = "Xe135"
    candu.temperatures = (293.6, 323.6, 573.6)
    candu.scatteringLaw = None
    candu.legendre = 0
    candu.fission = None
    candu.dilutions = None
    candu.mat = 5458
    candu.evaluationFile = "/home/user/Njoy99/evaluations/database/jendl310.asc"
    candu.makeFp()

Uranium-235
-----------

U-235 is one of the most prevalently used fissile material for energy
production. In case of CANDU reactors, natural uranium is used as fuel, where
weight (%) of U-235 is 0.711.

.. code-block:: python

    candu.hmat = "U235"
    candu.temperatures = ( 293.6, 323.6, 573.6, )
    candu.mat = 9228
    candu.evaluationFile = "/home/user/Njoy99/evaluations/database/U-235"
    candu.fission = 2
    candu.ss = (2.76792, 1.22773e5)
    candu.potential = 11.6070
    candu.dilutions = (1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
        11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5)
    candu.pendf()
    candu.gendf()
    candu.draglib()
    candu.dilutions = (1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, 1259.67004, \
        750.448669, 447.079956, 266.347961, 158.676849)
    candu.pendf()
    candu.gendf()
    candu.draglib()

Process Burnup
--------------

It is important to identify the tapes provided by evaluators which contains
information on fission yields and decay chains. These files are provided as
mentioned below, along with the chain(self) file mentioned already.

.. code-block:: python

    candu.fissionFile = "/home/user/Njoy99/evaluations/database/TAPE.107"
    candu.decayFile = "/home/user/Njoy99/evaluations/database/TAPE.106"
    candu.burnup()

-------------------------
ACELIB Library Generation
-------------------------

The modules that will be used in the generation of ACELIB are MODER, RECONR,
BROADR, PURR (if dilutions present), THERMR, and ACER. Figure 3 gives the flow
chart for generation of ACE formatted library. In this section specific examples
of elements will be provided to understand the nuances of ACELIB library
generation. The examples will be along the same lines as that for DRAGLIB
library generation, i.e scattering material (heavy water), structural material
(Zirconium), fission product (Xe135), Actinide (U-235). The prsent script is
such that the ACELIBs are appended in a single file named "acecandu" and is
available in the same directory as the draglib file. The other important file
that is generated is the "acexsdir", which contains the information about the
nuclides for which the cross sections are generated and the temperature at which
the ACELIB is generated. A small code has been written-"append.f", which will
read the file acecandu and acexsdir and create the file "myxsdir". This file has
to be appended to existing xsdir file provided with MCNP5 data. It is important
to provide suffix values "self.suff" for different temperatures. This will be
automatically appended to the ZA value and written in main ACELIB and xsdir
file. Care should be taken not to repeat the ".suff" value already used in xsdir
file for other evaluations. For each temperature provide different
"self.dirName" so that all the required data for comparison with PENDF tape
generated using PREPRO code is made possible. The "self.comp()" does the task of
verifying the ACELIBs generated, and is described in the next section.  The
example provided here helps in generating ACELIB alone. In case a DRAGLIB file
also is to be generated, please refer to Figure 2.

Heavy Water
-----------

.. code-block:: python

    candu.hmat = "H2_D2O"
    candu.temperatures = (293.6, 323.6, 573.6)
    candu.mat = 128
    candu.za = 1002
    candu.scatName = "hwtr"
    candu.evaluationFile = "/home/user/Tripoli4/JEF2/eval/jef2.neutron.H2.bcd"
    candu.scatteringLaw = "/home/user/Tripoli4/JEF2/eval/jef2.neutron.D_D2O.bcd.therm"
    candu.scatteringMat = 11
    candu.fission = None
    candu.dilutions = None
    candu.pendf()
    candu.dirName = "D2O-1"
    candu.tempace = (293.6,)
    candu.suff = 0.20
    candu.acer()
    candu.comp()
    candu.dirName = "D2O-2"
    candu.tempace = (323.6,)
    candu.suff = 0.21
    candu.acer()
    candu.comp()
    candu.dirName = "D2O-3"
    candu.tempace = (573.6,)
    candu.suff = 0.22
    candu.acer()
    candu.comp()

Zirconium
---------

.. code-block:: python

    candu.hmat = "Zr0"
    candu.temperatures = (293.6, 323.6, 573.6)
    candu.mat = 4000
    candu.za = 40000
    candu.evaluationFile = "/home/user/Njoy99/evaluations/database/jendl306.asc"
    candu.fission = None
    candu.dilutions = (1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
        158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395)
    candu.pendf()
    candu.dirName = "Zr-1"
    candu.tempace = (293.6,)
    candu.suff = 0.20
    candu.acer()
    candu.comp()
    candu.dirName = "Zr-2"
    candu.tempace = (323.6,)
    candu.suff = 0.21
    candu.acer()
    candu.comp()
    candu.dirName = "Zr-3"
    candu.tempace = (573.6,)
    candu.suff = 0.22
    candu.acer()
    candu.comp()

Xenon-135
---------

.. code-block:: python

    candu.hmat = "Xe135"
    candu.temperatures = ( 293.6, 323.6, 573.6, )
    candu.scatteringLaw = None
    candu.legendre = 0
    candu.fission = None
    candu.dilutions = None
    candu.mat = 5458
    candu.za = 54135
    candu.evaluationFile = "/home/user/Njoy99/evaluations/database/jendl310.asc"
    candu.makeFp()
    candu.dirName = "Xe135-1"
    candu.tempace = (293.6,)
    candu.suff = 0.20
    candu.acer()
    candu.comp()
    candu.dirName = "Xe135-2"
    candu.tempace = (323.6,)
    candu.suff = 0.21
    candu.acer()
    candu.comp()
    candu.dirName = "Xe135-3"
    candu.tempace = (573.6,)
    candu.suff = 0.22
    candu.acer()
    candu.comp()

Uranium-235
-----------

.. code-block:: python

    candu.hmat = "U235"
    candu.temperatures = (293.6, 323.6, 573.6)
    candu.mat = 9228
    candu.za = 92235
    candu.scatteringLaw = None
    candu.legendre = 0
    candu.evaluationFile = "/home/user/Njoy99/evaluations/database/U-235"
    candu.fission = 2
    candu.dilutions = (1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
        1259.67004, 750.448669, 447.079956, 266.347961, 158.676849)
    candu.dirName = "U235-1"
    candu.tempace = (293.6,)
    candu.suff = 0.20
    candu.acer()
    candu.comp()
    candu.dirName = "U235-2"
    candu.tempace = (323.6,)
    candu.suff = 0.21
    candu.acer()
    candu.comp()
    candu.dirName = "U235-3"
    candu.tempace = (573.6,)
    candu.suff = 0.22
    candu.acer()
    candu.comp()
