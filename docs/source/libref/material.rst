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
.. autoclass:: Material(comp, mass=-1.0, name='', atoms_per_mol=-1.0)

    .. automethod:: molecular_weight(atoms_per_mol=-1.0)
    .. automethod:: norm_comp()
    .. automethod:: normalize()
    .. automethod:: mult_by_mass()
    .. automethod:: from_hdf5(filename, datapath, row=-1, protocol=1)
    .. automethod:: write_hdf5(filename, datapath="/material", nucpath="/nuc_zz", row=-0.0, chunksize=100)
    .. automethod:: from_text(filename)
    .. automethod:: write_text(filename)
    .. automethod:: to_atom_frac()
    .. automethod:: from_atom_frac(atom_fracs)
    .. automethod:: sub_mat(nuc_sequence, name="")
    .. automethod:: set_mat(nuc_sequence, value, name="")
    .. automethod:: del_mat(nuc_sequence, name="")
    .. automethod:: sub_range(lower=0, upper=10000000, name="")
    .. automethod:: set_range(lower=0, upper=10000000, value, name="")
    .. automethod:: del_range(lower=0, upper=10000000, name="")
    .. automethod:: sub_u(name="")
    .. automethod:: sub_pu(name="")
    .. automethod:: sub_fp(name="")
    .. automethod:: sub_lan(name="")
    .. automethod:: sub_act(name="")
    .. automethod:: sub_tru(name="")
    .. automethod:: sub_ma(name="")



*****************************
Material Generation Functions
*****************************
The following top-level module functions are used to generate materials from various sources.

.. autofunction:: from_atom_frac(atom_fracs, mass=-1.0, name='', atoms_per_mol=-1.0)

---------

.. autofunction:: from_hdf5(filename, datapath, row=-1, protocol=1)

---------

.. autofunction:: from_text(filename, mass=-1.0, name='', atoms_per_mol=-1.0)

