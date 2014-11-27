.. _conda:

^^^^^^^^^^^^^^^^^^^^^^^^^^
Conda Package manager
^^^^^^^^^^^^^^^^^^^^^^^^^^

On mac and linux PyNE can be installed via the package manager conda. 
After installing anaconda or miniconda from 
`the Continuum downloads page <http://continuum.io/downloads>`_ 
add conda's binary directory to your bash profile by adding::

    export PATH=/path/to/anaconda/bin:$PATH

to your .bashrc or .bash_profile. Then in a new shell::

    conda install conda-build jinja2 nose setuptools pytables hdf5 scipy

on linux you may also need to run::

    conda install patchelf

Then dowload the latest conda-recipes `here 
<https://github.com/conda/conda-recipes/archive/master.zip>`_

cd to the conda-recipes directory and run::

    conda build pyne
    conda install $(conda build --output pyne)
    nuc_data_make

------
Binary
------
Binary distributions of the latest release (0.4) for mac and linux (64-bit) 
using the conda package manager can be installed by running the command::

    conda install -c https://conda.binstar.org/pyne pyne

A windows 32-bit binary is also available on conda via the same command but
it is highly experimental and likely broken. Conda binaries do not have 
moab/pytaps/mesh support (yet).

