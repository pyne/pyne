.. _macosdep:

===========================
MacOS Dependencies for PyNE
===========================

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

`Installing GCC <https://gcc.gnu.org/install/index.html>`_

If you have homebrew installed on your device, you can install GCC by entering ::

	$ brew install gcc

``Resources:``

#. `Documentation <https://gcc.gnu.org/onlinedocs/gfortran/#toc-Compiler-Characteristics-1>`_
#. `GCC FAQ <https://gcc.gnu.org/wiki/FAQ>`_
#. `GCC Glossary <https://gcc.gnu.org/wiki/GCC_glossary>`_


CMake
'''''

CMake is designed to build, test, and compile packages programs across platforms.

``To Install:``

To download CMake, you can visit `Downloads <https://cmake.org/download/>`_.

If you have homebrew installed on your device, you can install GCC by entering ::

	$ brew install gcc

``Resources:``

#. `Wiki <https://gitlab.kitware.com/cmake/community/-/wikis/home>`_
#. `Tutorial <https://cmake.org/cmake/help/latest/guide/tutorial/index.html>`_


Numpy
'''''

Numpy is a library in Python that adds support to matrices and arrays, 
as well as mathematical functions acting on those arrays.

``To Install:``

For the latest version using pip, in the command line enter ::

	$ pip install numpy

For the latest version using Conda, in the command line enter ::

	$ conda install numpy


``Resources:``

#. `Tutorial <https://numpy.org/learn/>`_
#. `Documentation <https://numpy.org/doc/stable/>`_


SciPy
'''''

SciPy is a Python library that adds functions for linear algebra, optimization, 
integration, and other scientific or engineering tasks.

``To Install:``

`Installation <https://www.scipy.org/install.html>`_

``Resources:``

#. `Getting Started <https://www.scipy.org/getting-started.html>`_
#. `Tutorial <https://docs.scipy.org/doc/scipy/reference/tutorial/index.html>`_
#. `Documentation <https://www.scipy.org/docs.html>`_


Cython
''''''

Cython is a compiler that helps in making C or C++ extensions for python.

``To Install:``

For the latest version, in the command line enter ::

	$ pip install cython

If you have homebrew installed on your device, you can install Cython by entering ::

	$ brew install cython

``Resources:``

#. `Wiki <https://github.com/cython/cython/wiki>`_
#. `User's Guide <https://cython.readthedocs.io/en/latest/src/userguide/index.html>`_
#. `Cython <https://cython.org>`_


HDF5
''''

HDF5 (the Hierarchical Data Format version 5) is a format that supports large, 
complex data in a file directory like structure similar to how you might with your computer.

``To Install:``

To install from source code, 
follow the instructions `here <https://www.hdfgroup.org/downloads/hdf5/source-code/>`_ .

For the latest version, in the command line enter ::

	$ brew install hdf5

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

For the latest version, in the command line enter ::

	$ brew install tables

``Resources:``

#. `FAQ <http://www.pytables.org/FAQ.html>`_
#. `Tutorial <http://www.pytables.org/usersguide/tutorials.html>`_
#. `Project Pointers <http://www.pytables.org/project_pointers.html>`_


LAPACK
''''''

LAPACK (Liner Algebra Package) is a software library for numerical liner algebra.

``To Install:``

Other installation methods can be found `here <http://www.netlib.org/lapack/#_software>`_

For the latest version, in the command line enter ::

	$ brew install lapack

``Resources:``

#. `FAQ <http://www.netlib.org/lapack/faq.html>`_
#. `User's Guide <http://www.netlib.org/lapack/lug/>`_


BLAS
''''

BLAS (Basic Linear Algebra Subroutines) coordinates operations on vectors and matrices.

``To Install:``

Other installation methods can be found `here <http://www.netlib.org/blas/#_software>`_

For the latest version, in the command line enter ::

	$ brew install openblas

``Resources:``

#. `Documentation <http://www.netlib.org/blas/#_documentation>`_


Numexpr
'''''''

Numexpr is a fast numerical evaluation tool for numpy, ensuring that 
expressions operating on arrays are faster and take up less memory.

``To Install:``

For the latest version, in the command line enter ::

	$ pip install numexpr

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

For the latest version, in the command line enter ::

	$ pip install sphinx

``Resources:``

#. `Sphinx <https://www.sphinx-doc.org/en/master/>`_
#. `Tutorial <http://matplotlib.sourceforge.net/sampledoc/>`_
#. `reStructuredText Cheat Sheet <https://docutils.sourceforge.io/docs/user/rst/cheatsheet.txt>`_


Sphinxcontrib-bibtex
''''''''''''''''''''

An extension allowing Sphinx to interact with BibTeX.

``To Install:``

For the latest version, in the command line enter ::

	$ pip install sphinxcontrib-bibtex

``Resources:``

#. `Documentation <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/>`_ 
#. `Known Issues and Workarounds <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/usage.html#known-issues-and-workarounds>`_
#. `Example <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/quickstart.html#minimal-example>`_


PrettyTable
'''''''''''

PrettyTable is a python library that adds a lot of versatility to table creation.

``To Install:``

For the latest version, in the command line enter ::

	$ pip install prettytable

``Resources:``

#. `Tutorial <https://code.google.com/archive/p/prettytable/wikis/Tutorial.wiki>`_


Cloud Sphinx
''''''''''''

Cloud is a Sphinx theme that PyNE uses to generate its 
HTML documentation (like this site).

``To Install:``

For the latest version, in the command line enter ::

	$ pip install cloud_sptheme

``Resources:``

#. `Documentation <https://cloud-sptheme.readthedocs.io/en/latest/>`_
#. `Source <https://foss.heptapod.net/doc-utils/cloud_sptheme>`_


Jupyter
'''''''

If you have downloaded Python through Anaconda Jupyter requirements should 
be satisfied, but it never hurts to make sure. 

You can check there version by entering ::

	$jupyter —-version

``To Install:``

If you are going to use Python 2:

For the latest version, in the command line enter ::

	$ pip install jupyter

If you are going to use Python 3:

For the latest version, in the command line enter ::

	$ pip3 install jupyter

``Resources:``

#. `Additional Installation Information <https://jupyter.readthedocs.io/en/latest/install.html>`_
#. `Project Documentation <https://jupyter.readthedocs.io/en/latest/projects/doc-proj-categories.html#deployment>`_