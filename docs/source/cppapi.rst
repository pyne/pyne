.. _cppapi:

=================
Using the C++ API
=================
While PyNE is a Python package, much of its core functionality exists purely in C/C++.
The python wrappers that expose these C-level utilities exist independently of the the 
C-code.  Therefore it is possible to use, compile, and link directly to PyNE core 
libraries from lower-level languages without having to touch Python at all.  Additionally, 
this enables PyNE to be used from *within* other Python extension modules!

The API for PyNE functions on the C++ level is very similar to that which is exposed to 
Python.  The differences lie largely in that C functions are strongly typed, whereas PyNE's
wrappers take heavy advantage of Python's duck typing.  For most use cases, the :ref:`usersguide`,
the :ref:`libref`, and the header files should be sufficient to describe the C++ API.

All PyNE shared objects are installed into the :file:`pyne/lib/` directory.  The headers
for both C/C++ and Cython are always installed into the :file:`pyne/includes/` directory.
Unfortunately, where :file:`pyne/` is located often changes from one system to another.
However, these locations may be found from the ``pyne_config`` module.

.. code-block:: ipython 

    In [1]: from pyne import pyne_config

    In [2]: pyne_config.lib
    Out[2]: '/usr/lib64/python2.7/site-packages/site-packages/pyne/lib'

    In [3]: pyne_config.includes
    Out[3]: '/usr/lib64/python2.7/site-packages/site-packages/pyne/includes'


The following table displays the C++ shared objects that are currently built (in Linux)
and which names to use when dynamically linking to them.  For an example of how to link 
please refer to PyNE's own `setup.py`_ file.  Additionally, feel free to contact the 
authors if you require additional assistance.

.. table:: Pure C++ Shared Object Libraries

    ======== =================== ===============
    Module   Filename            Link With
    ======== =================== ===============
    pyne     libpyne.so          'pyne'
    nucname  libpyne_nucname.so  'pyne_nucname'
    material libpyne_material.so 'pyne_material'
    ======== =================== ===============

.. _setup.py: https://github.com/pyne/pyne/blob/staging/setup.py
