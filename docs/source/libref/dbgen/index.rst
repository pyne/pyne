.. _pyne_dbgen:

============================================
Nuclear Data Generation -- :mod:`pyne.dbgen`
============================================
Pyne provides an easy-to-use, repeatable aggregation utility for nuclear data.  This command line
utility called ``nuc_data_make`` builds and installs an HDF5 file named ``nuc_data.h5`` to the 
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
    atomic_weight
    decay
    scattering_lengths
    simple_xs
    cinder
