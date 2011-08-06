=================================
Materials -- :mod:`pyne.material`
=================================
This moudule contains the Material class, which is used to represent nuclear
materials throughout PyNE.

.. _material:

.. currentmodule:: pyne.material

All functionality may be found in the ``material`` package::

 from pyne.material import Material

Materials are the primary container for radionuclides. They map nuclides to mass weights.
In many ways they take inspiration from numpy arrays and python dictionaries.  Materials
have two main attributes which define them.

1. **comp**: a nomralized composition mapping from nuclides (zzaaam) to mass-weights (floats).
2. **mass**: the mass of the material.

By keeping the mass and the composition separate, operations that only affect one attribute
may be performed independent of the other.
    
**************
Material Class
**************
.. autoclass:: Material(comp, mass=-1, name='')

    .. automethod:: molecular_weight()
    .. automethod:: norm_comp()
    .. automethod:: normalize()
    .. automethod:: mult_by_mass()
    .. automethod:: load_from_hdf5(filename, groupname, row=-1)
    .. automethod:: load_from_text(filename)
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

