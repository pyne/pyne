.. _osx_source:

-----------------------------
Source install - Mac OSX
-----------------------------

These instructions are based on using the homebrew http://brew.sh/ package manager
Install command line tools from https://developer.apple.com/downloads/
you will need to create an account in order to download::

    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
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

then follow the source install documentations :ref:`linux_source`
