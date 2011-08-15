.. _usersguide_cccc:

============
CCCC Formats
============

.. currentmodule:: pyne.cccc

The CCCC module contains a number of classes for reading various cross section,
flux, geometry, and data files with specifications given by the Committee for
Computer Code Coordination. The following types of files can be read using
classes from this module: ISOTXS, DLAYXS, BRKOXS, RTFLUX, ATFLUX, RZFLUX, MATXS,
and SPECTR.

The ISOTXS reader was originally derived from Professor James Holloway's
open-source C++ classes from the University of Michigan and later expanded by
Nick Touran for work on his PhD thesis. DLAYXS was later added by Paul Romano.

A description of several CCCC formats are available online for ISOTXS_, MATXS_,
RTFLUX_, and RZFLUX_. Other format specifications can be found in Los Alamos
Report LA-5324-MS_.
    
.. _ISOTXS: http://t2.lanl.gov/codes/transx-hyper/isotxs.html
.. _MATXS: http://t2.lanl.gov/codes/transx-hyper/matxs.html
.. _RTFLUX: http://t2.lanl.gov/codes/transx-hyper/rtflux.html
.. _RZFLUX: http://t2.lanl.gov/codes/transx-hyper/rzflux.html
.. _LA-5324-MS: http://www.osti.gov/bridge/servlets/purl/5369298-uIcX6p/

For a complete specification for the classes in the ``cccc`` module, please
refer to the Library Reference entry for :ref:`pyne_cccc`.

***************************
Example Use of Isotxs Class
***************************

To load data from an ISOTXS file, one needs to simply initialize an instance of
the Isotxs class specifying the path to the ISOTXS file.

.. code-block:: ipython

   In [1]: from pyne import cccc

   In [2]: isoFile = cccc.Isotxs('ISOTXS')

   In [3]: isoFile.read()

After the file has been read, the data from the ISOTXS is now accessible through
the attributes of the ``isoFile`` object. Some of the attributes are as follows:

:chi: 
  Normalized fission spectrum, giving the fraction of neutrons emitted in each
  energy group.

:emax:
  Maximum energies for each of the energy groups in units of eV.

:emin:
  Minimum energy for the lowest (thermal) energy group, typically somewhere
  around 0.00001 eV.

:fc:
  Dictionary containing several important parameters. For example
  ``isofile.fc['niso']`` tells you the number of nuclides in the file and
  ``isoFile.fc['ngroup']`` tells you the number of energy groups.

:label:
  Identification label describing the data in the file.

:nucNames:
  A list of the labels for each nuclide.

:nuclides:
  A list of nuclides present in the file. Each nuclide is an instance of
  the ``Nuclide`` class which contains data on cross-sections.

:vel:
  Mean neutron velocity in each energy group.
