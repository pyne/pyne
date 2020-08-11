.. _osx_source:

^^^^^^^^^^^^^^^^^^^^^^^^
Source install - Mac OSX
^^^^^^^^^^^^^^^^^^^^^^^^

=============
Dependencies:
=============

Many, if not all of these dependencies can be installed using a package manager
(i.e. `MacPorts <https://www.macports.org/>`__, `Anaconda
<https://www.anaconda.com/>`__, `Homebrew <https://brew.sh/>`__, etc.)

----
PyNE
----
#. `Python <https://www.python.org/>`__
#. `GCC <https://gcc.gnu.org/>`__
#. `CMake <http://www.cmake.org/>`_ (>= 2.8.5)
#. `NumPy <http://www.numpy.org/>`_ (>= 1.8.0)
#. `SciPy <https://www.scipy.org/docs.html>`__
#. `Cython <https://cython.org>`__ (>= 0.19.1)
#. `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`__
#. `PyTables <https://www.pytables.org/>`__
#. `LAPACK <http://www.netlib.org/lapack>`__
#. `BLAS <http://www.netlib.org/blas/#_documentation>`__
#. `Numexpr <https://numexpr.readthedocs.io/en/latest/user_guide.html>`__

--------------------
Website Dependencies
--------------------
These are only necessary if you intend to contribute to the PyNE website.

#. `Sphinx <https://www.sphinx-doc.org/en/master/>`__
#. `Sphinxcontrib-bibtex <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/>`__
#. `PrettyTable <https://code.google.com/archive/p/prettytable/>`__
#. `Cloud Sphinx <https://foss.heptapod.net/doc-utils/cloud_sptheme>`__
#. `Jupyter <https://jupyter.readthedocs.io/en/latest/install.html>`__

---------------------
Optional Dependencies
---------------------
#. `MOAB <https://press3.mcs.anl.gov/sigma/moab-library>`__
#. `PyTAPS <https://pythonhosted.org/PyTAPS/index.html>`__


=============
Installation:
=============

Download and unzip
the source (`zip`_, `tar`_) or checkout a version from the PyNE repository
(`Github`_).  Then run the following commands from the directory above the
unzipped directory

Add::

    export PATH=/usr/local/bin:$PATH
    export PATH=/usr/local/share/python:$PATH

to ~/.bash_profile, then::

    source ~/.bash_profile
    sudo chown -R $(whoami)/usr/local

download pyne cd to that directory::

    cd Downloads/pyne
    python setup.py install


Once those lines have been added, run the following command before running
``scripts/nuc_data_make``::

    source ~/.bashrc

.. _zip: https://github.com/pyne/pyne/zipball/0.5.1
.. _tar: https://github.com/pyne/pyne/tarball/0.5.1
.. _GitHub: http://github.com/pyne/pyne
