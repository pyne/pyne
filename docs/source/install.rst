.. _install:

============
Installation
============
-------------
Dependencies
-------------
PyNE has the following dependencies:

   #. `NumPy <http://numpy.scipy.org/>`_
   #. `SciPy <http://www.scipy.org/>`_
   #. `Cython <http://cython.org/>`_
   #. `HDF5 <http://www.hdfgroup.org/HDF5/>`_
   #. `PyTables <http://www.pytables.org/>`_

------
Binary
------
A binary distribution of PyNE is hopefully coming soon.  Until then, please
install from source.


.. _install_source:

------
Source
------
Installing PyNE from source is a two-step process.  First, download and 
unzip the source (`zip`_, `tar`_).  Then run the following commands from 
the unzipped directory::

    cd pyne/
    python setup.py install --user
    nuc_data_make

The ``setup.py`` command compiles and installs the PyNE source code.
The ``nuc_data_make`` builds and installs a database of nuclear data.
Unfortunately, this must be done as a second step because most nuclear 
data is under some form of license restriction or export control which 
prevents the developers from distributing it with PyNE.  However, the 
``nuc_data_make`` program (which is installed by ``setup.py``) will
do its best to find relevant nuclear data elsewhere on your machine
or from public sources on the internet.  


.. _win_install:

********************
Windows Installation
********************
Depending on the current state of your system, installing on Windows may 
be more or less involved.  We recommend the following procedure.  This 
ensures that all dependencies are installed correctly and PyNE has been 
built and tested using this setup.

#. Install the Enthought Python Distribution (`EPD`_).
#. Determine your HDF5 version by running the following command::

    python -c "import tables; print tables.getHDF5Version()"

#. Download the `HDF5 Windows binaries`_ for your version.
   Navigate to something like ``http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-{h5ver}/bin/windows/``
   and select the appropriate 32- or 64-bit file.  Do not download the source-only files.
#. Unzip HDF5 to the C-drive (``C:\\hdf5-{h5ver}``).
#. Download and unzip the source (`zip`_). 
#. Move into the source directory and run the PyNE setup command with the ``--hdf5`` option::

    cd pyne\
    python setup.py install --user --hdf5=C:\\hdf5-{h5ver}
           
And voila, everything will have installed correctly.  Should this still fail, 
please report your problem to pyne-dev@googlegroups.com.

********************
Linux + EPD
********************
Assuming you are on some flavor of Linux and you primarily use Python 
through the Enthought Python Distribution (`EPD`_), you can install PyNE
to be based off of the EPD packages.

First, you'll need to know where on your system EPD is installed.
Call this variable ``EPD_DIR``; for example if you have installed it 
to your home directory then ``EPD_DIR=$HOME/epd``.  You'll then need
to add the following lines to your :file:`~/.bashrc` file *after* 
installing EPD but *prior to* installing PyNE:

.. code-block:: bash

    export PATH=$EPD_DIR/bin:$HOME/.local/bin:$PATH
    export CPATH=$EPD_DIR/include:$CPATH
    export LD_LIBRARY_PATH=$EPD_DIR/lib:$LD_LIBRARY_PATH

Or as in the example:

.. code-block:: bash

    export PATH=$HOME/epd/bin:$HOME/.local/bin:$PATH
    export CPATH=$HOME/epd/include:$CPATH
    export LD_LIBRARY_PATH=$HOME/epd/lib:$LD_LIBRARY_PATH

You may now proceed with the PyNE install :ref:`as above <install_source>`.

.. _zip: https://github.com/pyne/pyne/zipball/0.1-rc
.. _tar: https://github.com/pyne/pyne/tarball/0.1-rc

.. _EPD: http://www.enthought.com/products/epd.php
.. _HDF5 Windows binaries: http://www.hdfgroup.org/ftp/HDF5/releases/
