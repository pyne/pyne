.. _pyne_material:

=================================
Materials -- :mod:`pyne.material`
=================================
This module contains the Material class, which is used to represent nuclear
materials throughout PyNE.

.. currentmodule:: pyne.material

All functionality may be found in the ``material`` package::

 from pyne import material

Materials are the primary container for mixtures of radionuclides. They map
nuclides to **mass fractions**, though they contain methods for converting
to/from atom fractions as well.  In many ways they take inspiration from numpy
arrays and python dictionaries.  Materials have two main attributes which define
them.

1. **comp**: a normalized composition mapping from nuclides (zzaaam-ints) to
   mass-fractions (floats).  
2. **mass**: the mass of the material.

The :ref:`pyne_material_library` class is available to efficiently manipulate
collections of materials.  The material class is presented below.  For more
information please refer to :ref:`usersguide_material`.

HDF5 File Structure
--------------------

When using the `write_hdf5` method to write a material in a group named
`my_mat` (the user shall provide the name for this group in argument `datapath`), 
the default structure for the `HDF5` file is:
.. verbatim::
    /material/
    /--------/my_mat/
             /------/composition
             /------/nuclidelist
             /------/composition_metadata

Where, `/material` and `/material/my_mat` are HDF5 groups, and 
`composition`, `nuclidelist` and `composition_metadata` are HDF5 datasets. 

The :ref:`pyne_material_library` class is available to efficiently manipulate
collections of materials.


Previous HDF5 File Structure
----------------------------

In files created with previous versions, it's possible that 
the `datapath` or `/material` already exist as a dataset in the file.
In such cases, the old writing method will be used.
Older data structure are deprecated but still available to written when providing a
 `nucpath` to the
`write_hdf5()` method, or when writing a material in a file with the old
data structure.
Old data structure looks like:
.. verbatim::
    /my_mat/
    /------/nucpath
    /my_mat_metadata
    /nuclidelist

`my_mat` (the `datapath` -- default `material`) is a HDF5 dataset containing the
material composition, `nucpath` is a attribute containing the path to the
nuclide list. The attribute is attached to the `datapath`.
`my_mat_metadata` is a dataset containing the metadata of the material.
`nuclidelist` is a dataset containing the list of nuclides composition the
material.

`from_hdf5()` will detect the structure (old or new) of the file (when using
`protocol1`).
   
The material class is presented below.  
========================================
*For more information please refer to :ref:`usersguide_material`.*
=======

**************
Material Class
**************
.. autoclass:: Material(comp, mass=-1.0, atoms_per_mol=-1.0, attrs=None)
   :members:
   :inherited-members:

*****************************
Material Generation Functions
*****************************
The following top-level module functions are used to generate materials from
various sources.

.. autofunction:: from_atom_frac

---------

.. autofunction:: from_hdf5

---------

.. autofunction:: from_text

****************
Material Library
****************
.. autoclass:: MaterialLibrary(lib=None, datapath="/materials", nucpath="/nucid")
   :members:
   :inherited-members:
