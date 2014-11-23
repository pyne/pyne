.. _osx_source:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Source install - Mac OSX
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
These instructions are based on using the homebrew http://brew.sh/ package manager
Install command line tools from https://developer.apple.com/downloads/
you will need to create an account in order to download::

    ruby -e "$(curl -fsSL https://raw.github.com/mxcl/homebrew/go/install)"
    brew doctor
    brew tap homebrew/science
    brew install hdf5
    brew install cmake
    brew install python

Add::

    export PATH=/usr/local/bin:$PATH
    export PATH=/usr/local/share/python:$PATH

to ~/.bash_profile, then::

    source ~/.bash_profile
    sudo pip install numpy
    sudo chown -R $(whoami) /usr/local
    brew install gfortran
    pip install scipy
    pip install cython
    pip install numexpr
    pip install tables

download pyne-staging cd to that directory::

    cd Downloads/pyne-staging
    python setup.py install


Once those lines have been added, run the following command before running 
``nuc_data_make``::

    source ~/.bashrc
