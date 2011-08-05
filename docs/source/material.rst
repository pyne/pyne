********************************************
Material Class -- :mod:`pyne.material`
********************************************
This moudule contains the Material class which is used to represent nuclear
materials throughout PyNE.

.. _material:

.. currentmodule:: pyne.material

All functionality may be found in the ``material`` package::

 from pyne.material import Material

.. autoclass:: Material(comp, [mass=-1, name=''])

    .. automethod:: atomic_weight()
    .. automethod:: norm_comp()
    .. automethod:: normalize()
    .. automethod:: mult_by_mass()
    .. automethod:: load_from_hdf5(filename, groupname, row=-1)
    .. automethod:: load_from_text(filename)
    .. automethod:: sub_mat(nuc_sequence, name="")
    .. automethod:: sub_u(name="")
    .. automethod:: sub_pu(name="")
    .. automethod:: sub_fp(name="")
    .. automethod:: sub_lan(name="")
    .. automethod:: sub_act(name="")
    .. automethod:: sub_tru(name="")
    .. automethod:: sub_ma(name="")
