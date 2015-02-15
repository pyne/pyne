.. _theorymanual_decay:

===============================
Data
===============================

.. currentmodule:: pyne.data

The :py:mod:`pyne.data` module implements access to various different sets of 
fundamental nuclear data stored in the `nuc_data.h5` HDF5 database. These 
include the following distinct datasets:

* Fundamental constants common to nuclear engineering problems
* Atomic mass data (JAEA)
* Natural isotopic abundance data (IUPAC2009)
* Energy/fission (ORIGEN-S) http://web.ornl.gov/~webworks/cppr/y2001/rpt/97914.pdf
* Gamma Energy/fission (ORIGEN-S)
* Simple cross sections (KAERI)
* Dose factors for air, inhalation, soil, ingestion, uptake in fluids and lungs (Hanford report)
* Coherent/incoherent/total bound scattering lengths (NIST)
* Fission product yield from WIMSD and NDS (IAEA)
* X-ray conversion coefficients (BNL)
* ENSDF decay and level data (BNL)

Module allows the grabbing of q_values (energy per disintegration) for the
calculation of decay heat. This currently consists of the nuclide, it's
q_value, and the percent of q coming from gammas. This data is from
'ORIGEN-S DECAY DATA LIBRARY AND HALF-LIFE UNCERTAINTIES'
(http://web.ornl.gov/~webworks/cppr/y2001/rpt/97914.pdf)

Module allows the grabbing of dose rate factors for the calculation of radiotoxicity. There are four dose rates provided:
 1. external from air (mrem/h per Ci/m^3)
     Table includes: nuclide, air dose rate factor, ratio to inhalation dose (All EPA values)
 2. external from 15 cm of soil (mrem/h per Ci/m^2)
     Table includes: nuclide, GENII, EPA, DOE, GENII/EPA, DOE/EPA
 3. ingestion (mrem/pCi)
     Table includes: nuclide, f1 (fraction of the activity ingested that enters body fluids.), GENII, EPA, DOE, GENII/EPA, DOE/EPA
 4. inhalation (mrem/pCi)
     Table includes: nuclide, lung model*, GENII, EPA, DOE, GENII/EPA, DOE/EPA 

This data is from: 
[Exposure Scenarios and Unit Dose Factors for the Hanford
Immobilized Low-Activity Tank Waste Performance Assessment, ref.
HNF-SD-WM-TI-707 Rev. 1 December 1999] Appendix O of HNF-5636 [DATA PACKAGES
FOR THE HANFORD IMMOBILIZED LOW-ACTIVITY TANK WASTE PERFORMANCE ASSESSMENT:
2001 VERSION]

This module provides a way to grab and store raw data for neutron scattering 
lengths.  This data comes from Neutron News, Vol. 3, No. 3, 1992, pp. 29-37 via 
a NIST webpage (http://www.ncnr.nist.gov/resources/n-lengths/list.html).  Please
contact Alan Munter, <alan.munter@nist.gov> for more information.

**********************
Additional Information
**********************
For further information, please see:

* `The ENSDF manual <https://>`_, 
* `IAEA nuclear data section <https://>`_, and
* `the benchmark study <http://>`_.



Note that the benchmark study shows quite high agreement between this method
and ORIGEN v2.2.
