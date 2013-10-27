.. _pyne_material:

=================================
Materials -- :mod:`pyne.material`
=================================
This module contains the Material class, which is used to represent nuclear
materials throughout PyNE.

.. currentmodule:: pyne.material

All functionality may be found in the ``material`` package::

 from pyne import material

Materials are the primary container for radionuclides. They map nuclides to **mass weights**, 
though they contain methods for converting to/from atom fractions as well.
In many ways they take inspiration from numpy arrays and python dictionaries.  Materials
have two main attributes which define them.

1. **comp**: a normalized composition mapping from nuclides (zzaaam-ints) to mass-weights (floats).
2. **mass**: the mass of the material.

The material class is presented below.  For more information please refer to :ref:`usersguide_material`.

**************
Material Class
**************
.. autoclass:: Material(comp, mass=-1.0, atoms_per_mol=-1.0, attrs=None)
   :members:
   :inherited-members:

*****************************
Material Generation Functions
*****************************
The following top-level module functions are used to generate materials from various sources.

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

