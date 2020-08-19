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
Binary distributions of the latest release for mac and linux (64-bit) 
using the conda package manager can be installed by running the command::

    conda install -c conda-forge pyne

Conda binaries do not have moab/pymoab/mesh support (yet).
