PyNE: The Nuclear Engineering Toolkit
=====================================
The pyne project aims to provide a common set of tools for nuclear 
science and engineering needs.

If you are interested in the package itself, or would like to help
and contribute, please let us know either on the mailing list 
(pyne-dev@googlegroups.com) or `github`_.

.. _github: https://github.com/pyne/pyne

.. install-start

.. _install:

============
Installation
============
-------------
Dependencies
-------------
PyNE has the following dependencies:

   #. `CMake <http://www.cmake.org/>`_ (>= 2.8.5)
   #. `NumPy <http://www.numpy.org/>`_
   #. `SciPy <http://www.scipy.org/>`_
   #. `Cython <http://cython.org/>`_
   #. `HDF5 <http://www.hdfgroup.org/HDF5/>`_
   #. `PyTables <http://www.pytables.org/>`_
   #. `Python 2.7 <http://www.python.org/>`_

Additionally, building the documentation requires the following:

   #. `Sphinx <http://sphinx-doc.org/>`_
   #. `SciSphinx <https://github.com/numfocus/scisphinx>`_
   #. `breathe <http://michaeljones.github.io/breathe/>`_ 

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
    scripts/nuc_data_make

The ``setup.py`` command compiles and installs the PyNE source code.
The ``nuc_data_make`` builds and installs a database of nuclear data.
Unfortunately, this must be done as a second step because most nuclear 
data is under some form of license restriction or export control which 
prevents the developers from distributing it with PyNE.  However, the 
``nuc_data_make`` program (which is installed by ``setup.py``) will
do its best to find relevant nuclear data elsewhere on your machine
or from public sources on the internet.  

On MacOSX, it may be necessary to add the pyne library path to the 
``DYLD_FALLBACK_LIBRARY_PATH`` environment variable *before* running 
``nuc_data_make``. To do this, add the following lines to your 
``~/.bashrc`` file where ``/path/to/pyne/lib`` is the absolute path to the 
directory containing libpyne.dylib :: 

    DYLD_FALLBACK_LIBRARY_PATH="${DYLD_FALLBACK_LIBRARY_PATH}:/path/to/pyne/lib"
    export DYLD_FALLBACK_LIBRARY_PATH

Once those lines have been added, run the following command before running 
``nuc_data_make`` ::

    source ~/.bashrc

.. install-end


============
Contributing
============
We highly encourage contributions to PyNE! If you would like to contribute, 
it is as easy as forking the repository on GitHub, making your changes, and 
issuing a pull request. If you have any questions about this process don't 
hesitate to ask the mailing list (pyne-dev@googlegroups.com).


.. _zip: https://github.com/pyne/pyne/zipball/0.3
.. _tar: https://github.com/pyne/pyne/tarball/0.3

