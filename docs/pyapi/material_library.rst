.. _pyne_material_library:

==================================================
Material Libraries -- :mod:`pyne.material_library`
==================================================
This module contains the Material Library class, which is used to represent
collections of nuclear materials throughout PyNE.

.. currentmodule:: pyne.material_library

All functionality may be found in the ``material_library`` package::

 from pyne import material_library

Material Libraries provide a way of assembling a collection of :ref:`pyne_material`
objects, including the reading and writing of those collections in various
formats.  They behave in a dictionary-like way in which the key is a string
that is generally a semantically meaningful name, and the value is a :ref:`pyne_material`
object. MaterialLibrary contains mainly a dictionary like object referencing a list of
:ref:`pyne_material` by the a string name key.

The MaterialLibrary class is presented below.  For more information please refer to :ref:`usersguide_material`.

HDF5 File Structure
-------------------

When using the `write_hdf5` method to write a material library in a group named
`my_mat_lib` (the user shall provide the name for this group in argument `datapath`), 
the default structure for the `HDF5` file is:
.. verbatim::
    /material_library/
    /----------------/my_mat_lib/
                     /----------/composition
                     /----------/nuclidelist
                     /----------/composition_metadata

Where, `/material_library` and `/material_library/my_mat_lib` are HDF5 groups. 

Previous HDF5 File Structure
-----------------------------

If the `datapath` or `/material_library` exist as a dataset in the file, 
then the old writing method will be used.

Old data structure looks like:
.. verbatim::
    /my_mat_lib/
    /----------/nucpath
    /my_mat_lib_metadata
    /nuclidelist

`my_mat_lib` (the `datapath` -- default `material`) is a HDF5 dataset containing the
array of material compositions, `nucpath` is a attribute containing the path to the
nuclide list (attached to the `datapath`).
`my_mat_lib_metadata` is a dataset containing an array of metadata of the materials.
`nuclidelist` is a dataset containing the list of nuclides composition the
materials.

`from_hdf5()` will detect the structure (old or new) of the file.

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



