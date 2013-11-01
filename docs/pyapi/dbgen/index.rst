.. _pyne_dbgen:

============================================
Nuclear Data Generation -- :mod:`pyne.dbgen`
============================================
Pyne provides an easy-to-use, repeatable aggregation utility for nuclear data.  This command line
utility is called ``nuc_data_make`` builds and installs an HDF5 file named ``nuc_data.h5`` to the 
current PyNE install.  Nuclear data is gathered from a vareity of sources, including the web 
and the data files for other programs installed on your system (such as MCNP).

All of the code to produce the ``nuc_data.h5`` file is found in the ``dbgen`` sub-package.  This
package was designed to be modular.  Therfore ``nuc_data_make`` may be run such that only an available 
subset of ``nuc_data.h5`` is produced.  

As this is the library refence portion of the documention, the underlying functionality for each 
module is displayed rather than how to use the end product.  However, most modules here are 
divided conceptually into three parts, run in series:

#. Gather raw data and place it in a build directory.
#. Parse the raw data and put it in a form suitable for storage in ``nuc_data.h5``.
#. Write the parsed data to ``nuc_data.h5``.

**Database Generation Modules**

.. toctree::
    :maxdepth: 1

    nuc_data_make
    kaeri
    atomic_mass
    atomic_weight
    decay
    scattering_lengths
    simple_xs
    cinder
    eaf
    materials_library

++++++++++++++++
Attribution
++++++++++++++++
The PyNE development team would like to thank the following people and organizations
for allowing us to redistribute open nuclear data in binary form as part of 
``nuc_data_make``.

* **KAERI** and **Jonghwa Chang** for the data available at the excellent 
  http://atom.kaeri.re.kr website.
* **NIST** and **Alan Munter** for `neutron scattering length data`_.
* **The Atomic Mass Data Center**, **Georges Audi**, and **Wang Meng**  for 
  `atomic mass related data`_.
* **PNNL** and **Ronald J. McConn Jr** for providing us with the 
  `materials compendium data`_.
  

.. _neutron scattering length data: http://www.ncnr.nist.gov/resources/n-lengths/list.html
.. _atomic mass related data: http://amdc.in2p3.fr/
.. _materials compendium data: http://www.pnnl.gov/main/publications/external/technical_reports/PNNL-15870Rev1.pdf

++++++++++++++++++++++++++++
Prebuilt Nuclear Data
++++++++++++++++++++++++++++
For developers who wish to generate the open nuclear data file that is distributed along with PyNE, 
please run the following command::

    nuc_data_make --fetch-prebuilt False --make-open-only True -o prebuilt_nuc_data.h5

