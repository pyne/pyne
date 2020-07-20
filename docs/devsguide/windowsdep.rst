.. _windowsdep:

=============================
Windows Dependencies for PyNE
=============================

Something to note before you begin the installation process is the version of Python 
that you are using. PyNE is currently supported from Python 2 through Python 3.6 
(NOT 3.7, which is coming in future updates).


-------------
PyNE Library:
-------------

PyNE is what is known as a library in Python (as with many languages, 
libraries can be built and installed to add functionalities to the language), 
and it is designed to assist in Nuclear Engineering projects. 


Python
''''''

Python is a popular computer language for building applications, 
data science, and visualization. 

``To Install:``

Anaconda is a comprehensive Python package 
which provides a free installation for MacOS, Windows, and Linux. Follow the instructions 
for your device's platform `here <https://docs.anaconda.com/anaconda/install/>`_ if you would 
like to use Anaconda. Installing Python through Anaconda will result in many of the packages and 
dependencies for PyNE being satisfied.

If you want to download Python without the rest of Anaconda, use 
the links in the resources below. Keep in-mind that PyNE is currently 
not supported on Python 3.7 or later, but you can follow the instructions
for Anaconda users in the New Developer's Guide.

``Resources:``

#. `Python Download (Not Anaconda) <https://www.python.org/downloads/>`_
#. `Python Wiki <https://wiki.python.org/moin/>`_
#. `Beginner’s Guide <https://wiki.python.org/moin/BeginnersGuide>`_


GCC
'''

GCC stands for GNU Compiler Collection, and GNU is a Unix-Like operating system 
that allocates resources and communicates with hardware. GNU is often used with 
a kernel called Linux. The Compiler Collection is one element that allows for 
PyNE to work as a Python library, even though it is built with several languages.

``To Install:``

To download GCC, visit `Installing GCC <https://gcc.gnu.org/install/index.html>`_

``Resources:``

#. `Documentation <https://gcc.gnu.org/onlinedocs/gfortran/#toc-Compiler-Characteristics-1>`_
#. `GCC FAQ <https://gcc.gnu.org/wiki/FAQ>`_
#. `GCC Glossary <https://gcc.gnu.org/wiki/GCC_glossary>`_


CMake
'''''

CMake is designed to build, test, and compile packages programs across platforms.

``To Install:``

To download CMake, you can visit `Downloads <https://cmake.org/download/>`_.

``Resources:``

#. `Wiki <https://gitlab.kitware.com/cmake/community/-/wikis/home>`_
#. `Tutorial <https://cmake.org/cmake/help/latest/guide/tutorial/index.html>`_


Numpy
'''''

Numpy is a library in Python that adds support to matrices and arrays, 
as well as mathematical functions acting on those arrays.

``To Install:``

First, install PIP 

#. Download this file `get-pip.py <https://bootstrap.pypa.io/get-pip.py>`_

#. Launch your Command Prompt window
	A. Press Windows Key + X
	B. Click Run
	C. Type in cmd.exe and hit enter

#. ‘cd’ to the folder with get-pip.py in it

#. enter ::

    $ python get-pip.py

For more details, go `here <https://phoenixnap.com/kb/install-pip-windows>`_


Then use PIP to install Numpy ::

	$ pip3 install numpy

For more details, go `here <https://phoenixnap.com/kb/install-numpy>`_ 

``Resources:``

#. `Tutorial <https://numpy.org/learn/>`_
#. `Documentation <https://numpy.org/doc/stable/>`_


SciPy
'''''

SciPy is a Python library that adds functions for linear algebra, optimization, 
integration, and other scientific or engineering tasks.

``To Install:``

To use PIP to install SciPy::

	$ pip install scipy

For more details, go `here <https://www.scipy.org/install.html>`_

``Resources:``

#. `Getting Started <https://www.scipy.org/getting-started.html>`_
#. `Tutorial <https://docs.scipy.org/doc/scipy/reference/tutorial/index.html>`_
#. `Documentation <https://www.scipy.org/docs.html>`_


Cython
''''''

Cython is a compiler that helps in making C or C++ extensions for python.

``To Install:``

To download the newest release of Cython, visit https://cython.org/. 
Then unpack the file, enter the directory, and run::
	
	$ python setup.py install

To use PIP to install Cython::

	$ pip install Cython

For more details, go `here <https://cython.readthedocs.io/en/latest/src/quickstart/install.html>`_ .

For a list of Windows installers, go `here <https://www.lfd.uci.edu/~gohlke/pythonlibs/#cython>`_ .

``Resources:``

#. `Wiki <https://github.com/cython/cython/wiki>`_
#. `User's Guide <https://cython.readthedocs.io/en/latest/src/userguide/index.html>`_
#. `Cython <https://cython.org>`_


HDF5
''''

HDF5 (the Hierarchical Data Format version 5) is a format that supports large, 
complex data in a file directory like structure similar to how you might with your computer.

``To Install:``

Follow the instructions 
`here <https://support.hdfgroup.org/HDF5/faq/windows.html#:~:text=Download%20the%20HDF5%20pre%2Dcompiled,HDF5%20libraries%20and%20include%20files.>`_ 
to install HDF5.

To install from source code, 
follow the instructions `here <https://www.hdfgroup.org/downloads/hdf5/source-code/>`_ .

``Resources:``

#. `Examples <https://portal.hdfgroup.org/display/HDF5/HDF5+Examples>`_
#. `Learning HDF5 <https://portal.hdfgroup.org/display/HDF5/Learning+HDF5>`_
#. `Known Problems <https://portal.hdfgroup.org/display/support/HDF5%201.12.0#knownprob>`_


PyTables
''''''''

PyTables is a package for managing large hierarchical datasets.

``To Install:``

For a variety of installation instructions, 
follow the instructions `here <http://www.pytables.org/usersguide/installation.html>`_ .

``Resources:``

#. `FAQ <http://www.pytables.org/FAQ.html>`_
#. `Tutorial <http://www.pytables.org/usersguide/tutorials.html>`_
#. `Project Pointers <http://www.pytables.org/project_pointers.html>`_


BLAS
''''

BLAS (Basic Linear Algebra Subroutines) coordinates operations on vectors and .

``To Install:``

Follow the instruction methods `here <http://icl.cs.utk.edu/lapack-for-windows/lapack/>`_ 
to install BLAS.

Installation methods can be found `here <http://www.netlib.org/blas/#_software>`_ .

``Resources:``

#. `Documentation <http://www.netlib.org/blas/#_documentation>`_


LAPACK
''''''

LAPACK (Liner Algebra Package) is a software library for numerical liner algebra.

``To Install:``

Follow the instruction methods `here <https://icl.cs.utk.edu/lapack-for-windows/lapack/#libraries_intel>`_ 
to install LAPACK.

Installation methods can be found `here <http://www.netlib.org/blas/#_software>`_ .

``Resources:``

#. `FAQ <http://www.netlib.org/lapack/faq.html>`_
#. `User's Guide <http://www.netlib.org/lapack/lug/>`_


Numexpr
'''''''

Numexpr is a fast numerical evaluation tool for numpy, ensuring that 
expressions operating on arrays are faster and take up less memory.

``To Install:``

Based off of instructions found `here <https://numexpr.readthedocs.io/projects/NumExpr3/en/latest/user_guide.html>`_
, do the following:

#. Download `numexpr <https://github.com/pydata/numexpr>`_

#. Enter the directory where numexpr is located, and run ::

    $ python setup.py build
    $ python setup.py install

#. Then, enter a different directory and test numexpr ::

    $ python -c "import numexpr; numexpr.test()"

``Resources:``

#. `PyPi Project Homepage <https://pypi.org/project/numexpr/>`_
#. `Github Repository <https://github.com/pydata/numexpr>`_


--------
Website:
--------

Sphinx
''''''

A python based documentation generator that allows files to be written into HTML, LaTeX, 
ePub, Texinfo, pages, and plain text. Sphinx uses reStructuredText, which is a very 
straight-forward markup language.

``To Install:``

Either `install from source <https://www.sphinx-doc.org/en/master/usage/installation.html#installation-from-source>`_ 
, or `install with pip <https://www.sphinx-doc.org/en/master/usage/installation.html#install-pypi>`_ .

To install from the `github source <https://github.com/sphinx-doc/sphinx>`_ ::

#. Clone the Sphinx repository ::

    $ git clone https://github.com/sphinx-doc/sphinx

#. Enter the clone

#. Enter ::

    $ pip install .
    $ pip install git+https://github.com/sphinx-doc/sphinx

To install with pip:

#. Open Command Prompt and enter ::

    C:\> pip install -U sphinx

``Resources:``

#. `Sphinx <https://www.sphinx-doc.org/en/master/>`_
#. `Tutorial <http://matplotlib.sourceforge.net/sampledoc/>`_
#. `reStructuredText Cheat Sheet <https://docutils.sourceforge.io/docs/user/rst/cheatsheet.txt>`_


Sphinxcontrib-bibtex
''''''''''''''''''''

An extension allowing Sphinx to interact with BibTeX.

``To Install:``

#. Download the Sphinxcontrib-bibtex from `source <https://github.com/mcmtroffaes/sphinxcontrib-bibtex>`_ .

#. Un-zip the file (if necessary) and, in the command prompt, "cd" into the folder.

#. Enter ::

    python setup.py install


For more details about installing on windows with a setup.py file, go 
`here <https://docs.python.org/2/install/#standard-build-and-install>`_ .

``Resources:``

#. `Documentation <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/>`_ 
#. `Known Issues and Workarounds <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/usage.html#known-issues-and-workarounds>`_
#. `Example <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/quickstart.html#minimal-example>`_


PrettyTable
'''''''''''

PrettyTable is a python library that adds a lot of versatility to table creation.

``To Install:``

Visit the `Installation Wiki <https://code.google.com/archive/p/prettytable/wikis/Installation.wiki>`_ for instructions on how to 
install PrettyTable. 

``Resources:``

#. `Tutorial <https://code.google.com/archive/p/prettytable/wikis/Tutorial.wiki>`_


Cloud Sphinx
''''''''''''

Cloud is a Sphinx theme that PyNE uses to generate its 
HTML documentation (like this site).

``To Install:``

#. Download the `source code <https://foss.heptapod.net/doc-utils/cloud_sptheme/-/tree/branch/default>`_ .

#. Un-zip the file (if necessary) and, in the command prompt, "cd" into the folder.

#. Enter ::

    $ python setup.py install

For more details about installing on windows with a setup.py file, go 
`here <https://docs.python.org/2/install/#standard-build-and-install>`_ .

``Resources:``

#. `Documentation <https://cloud-sptheme.readthedocs.io/en/latest/>`_
#. `Source <https://foss.heptapod.net/doc-utils/cloud_sptheme>`_


Jupyter
'''''''

If you have downloaded Python through Anaconda Jupyter requirements should 
be satisfied, but it never hurts to make sure. 


``To Install:``

To download Jupyter Notebook, visit `Installing Jupyter-Notebook <https://jupyter.readthedocs.io/en/latest/install.html>`_ .

``Resources:``

#. `Additional Installation Information <https://jupyter.readthedocs.io/en/latest/install.html>`_
#. `Project Documentation <https://jupyter.readthedocs.io/en/latest/projects/doc-proj-categories.html#deployment>`_