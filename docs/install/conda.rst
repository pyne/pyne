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

To list all of the versions of `pyne` available on `conda-forge
<https://conda-forge.github.io/>`_ channel, in your terminal window run:

    conda search -c conda-forge pyne

Binary distributions of the latest release for mac and linux (64-bit) 
using the conda package manager can be installed by running the command:

    conda install -c conda-forge pyne

If you want to install PyNE with the correct package specification, try
``pkg_name=version=build_string``.

For example, if you want to install `pyne` `version=0.7.0` with build option `moab_openmc_py36he21c9c4_1`, you would enter:

    conda install -c conda-forge pyne=0.7.0=moab_openmc_py36h094fa6c_0

where version should be replaced with the version number to be installed.

Conda binaries do not have moab/pymoab/mesh support (yet).
