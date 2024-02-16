.. _conda:

^^^^^^^^^^^^^^^^^^^^^^^^^^
Conda Package Manager
^^^^^^^^^^^^^^^^^^^^^^^^^^
The first step is to make sure that you have 
the conda package manager installed. 
You can download and install either anaconda or miniconda from 
`the downloads page <https://www.anaconda.com/distribution/#download-section>`_.
Make sure that you follow the full instructions and that the 
conda command line utility is on your PATH.  This is the default 
option when installing conda.

--------------------------
Binary Package (For Users)
--------------------------

To list all of the versions of PyNE available on `conda-forge
<https://conda-forge.github.io/>`_ channel, in your terminal window run::

    conda search -c conda-forge pyne

Binary distributions of the latest release for linux (64-bit) 
using the conda package manager can be installed by running the command::

    conda install -c conda-forge pyne

If you want to install PyNE with the correct package specification, try
``pkg_name=version=build_string``.

For example, if you want to install ``pyne version=0.7.1`` with build option ``moab_openmc``, you would enter::

    conda install -c conda-forge pyne=0.7.1=moab_openmc*

where version should be replaced with the version number to be installed.

Conda binaries do not have moab/pymoab/mesh support (yet).
