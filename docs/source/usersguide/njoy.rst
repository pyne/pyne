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

*****************
Njoy99 Attributes
*****************

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

self.wght - User-defined weighting spectra to be used if self.iwt = 1. Given as """–delimited string.

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
