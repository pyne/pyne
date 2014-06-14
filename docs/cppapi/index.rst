.. _cppapi:

=================
C++ API
=================
The following files and libraries are part of the PyNE C++ interface:

.. toctree::
    :maxdepth: 1

    pyne
    nucname
    rxname
    data
    material
    enrichment
    extra_types
    h5wrap


-----------------
Using the C++ API
-----------------
While PyNE is a Python package, much of its core functionality exists purely in C/C++.
The python wrappers that expose these C-level utilities exist independently of the the 
C-code.  Therefore it is possible to use, compile, and link directly to PyNE core 
libraries from lower-level languages without having to touch Python at all.  
Additionally, this enables PyNE to be used from *within* other Python extension 
modules!

The API for PyNE functions on the C++ level is very similar to that which is exposed 
to Python.  The differences lie largely in that C functions are strongly typed, 
whereas PyNE's wrappers take heavy advantage of Python's duck typing.  For most use 
cases, the :ref:`usersguide`, the :ref:`pyapi`, and the header files should be 
sufficient to describe the C++ API.

PyNE creates a single shared object called ``libpyne.so`` (or equivalent) and is 
installed into Python's ``<prefix>/lib`` directory.  The 
headers for C/C++ are installed into the :file:`<prefix>/include/pyne` directory. 
Unfortunately, where :file:`<prefix>` is located changes based on the Python 
interpreter and arguments passed to ``setup.py`` at install time.
These often change from one system to another. However, these locations may be 
found from either the ``pyne`` or the ``pyne_config`` module.

.. code-block:: ipython 

    In [1]: from pyne import lib, includes

    In [2]: lib
    Out[2]: '/usr/lib64'

    In [3]: includes
    Out[3]: '/usr/include'

Therefore to link against pyne, add the following options to your linker:

.. code-block:: bash

    -lpyne -L$(python -c 'import pyne; print(pyne.lib)')

To include pyne from other C/C++ code, use the top-level pyne header:

.. code-block:: c++

    #include "pyne/pyne.h"

And be sure that the include directory is on your ``C_INLCUDE_PATH`` and 
``CPLUS_INLCUDE_PATH``.  You can do so from your compiler with:

.. code-block:: bash

    -I$(python -c 'import pyne; print(pyne.includes)')

.. _setup.py: https://github.com/pyne/pyne/blob/staging/setup.py


--------------------------------------------
Amalgamating PyNE into a Single Source File
--------------------------------------------
PyNE has a lot of great stuff in it! However, adding dependencies to C++ projects
can be annoying, frustrating, and error prone. It often seems easier to just rip 
out the functionality that you need and include it in your own project.  

*Good news!* PyNE offers a formal mechanism for combining some or all of its
C++ API into a single, redistributable source file and an accompanying header file.
This let's you use pyne in your projects without adding pyne as an external dependency.
This mechanism is known as *amalgamation*. 

In the top level pyne source code directory, there is an ``amalgamate.py`` script.
Simply run this script to combine all C++ source information into ``pyne.cpp`` and
``pyne.h`` files.  Run with no options to combine all commonly used C++ files.
Add options to modify the behavior.  Current options are:

.. code-block:: bash

    scopatz@ares ~/pyne $ ./amalgamate.py -h
    usage: amalgamate.py [-h] [-s SOURCE_PATH] [-i HEADER_PATH]
                         [-f FILES [FILES ...]]

    optional arguments:
      -h, --help            show this help message and exit
      -s SOURCE_PATH        Output *.cpp source path.
      -i HEADER_PATH        Output header path.
      -f FILES [FILES ...]  Files to amalgamate.

For example, to take only up through the rxname, amalgamate with:

.. code-block:: bash

    scopatz@ares ~/pyne $ ./amalgamate.py -s pyne.cc -i pyne.h -f license.txt \
    cpp/pyne.* cpp/extra_types.h cpp/h5wrap.h cpp/nucname.* cpp/rxname.*

`Cyclus <http://fuelcycle.org>`_ is an example of a project which uses an amalgamated
version of pyne.

