.. _conda:

^^^^^^^^^^^^^^^^^^^^^^^^^^
Conda Package Manager
^^^^^^^^^^^^^^^^^^^^^^^^^^
The first step is to make sure that you have 
the conda package manager installed. 
You can download and install either anaconda or miniconda from 
`the Continuum downloads page <http://continuum.io/downloads>`_.
Make sure that you follow the full instructions and that the 
conda command line utility is on your PATH.  This is the default 
option when installing conda.

--------------------------
Binary Package (For Users)
--------------------------
Binary distributions of the latest release for mac and linux (64-bit) 
using the conda package manager can be installed by running the command::

    conda install -c conda-forge pyne

where VERSION should be replaced with the version number to be installed.

A windows 32-bit binary is also available on conda via the same command but
it is highly experimental and likely broken. Conda binaries do not have 
moab/pymoab/mesh support (yet).

----------------------------------
Create a Package (For Developer)
----------------------------------
In a new shell, run the following::

    conda install conda-build jinja2 nose setuptools pytables hdf5 scipy

on linux you may also need to run::

    conda install patchelf

Then dowload the latest conda-recipes `here 
<https://github.com/conda/conda-recipes/archive/master.zip>`_

cd to the conda-recipes directory and run::

    conda build pyne
    conda install $(conda build --output pyne)
    nuc_data_make

You should now be able to update the conda package that lives on binstar.
