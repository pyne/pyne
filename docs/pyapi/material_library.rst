.. _pyne_material_library:

==================================================
Material Libraries -- :mod:`pyne.material_library`
==================================================
This module contains the Material Library class, which is used to represent
collections of nuclearnmaterials throughout PyNE.

.. currentmodule:: pyne.material_library

All functionality may be found in the ``material_library`` package::

 from pyne import material_library

Material Libraries provide a way of assembling a collection of PyNE Material
objects, including the reading and writing of those collections in various
formats.  They behave in a dictionary-like way in which the key is a string
that is generally a semantically meaningful name, and the value is a Material
object.

The material_library class is presented below.  For more information please refer to :ref:`usersguide_material`.

*********************
MaterialLibrary Class
*********************
.. autoclass:: MaterialLibrary(lib=None, datapath="/materials", nucpath="/nucid")
   :members:
   :inherited-members:

************************************
MaterialLibrary Read/Write Functions
************************************

A key functionality of the `MaterialLibrary` class is the ability to read and
write material libraries in different formats.

.. autofunction:: from_hdf5

---------

.. autofunction:: write_hdf5

---------

.. autofunction:: from_json

---------

.. autofunction:: write_json


**************************************
MaterialLibrary Modification Functions
**************************************

Functions are provided to ensure robust behavior when working with changes to
`MaterialLibrary` collections.

.. autofunction:: add_material

---------

.. autofunction:: remove_material

---------

.. autofunction:: merge





