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
may be performed independent of the other.  Additionally, most of the functionality is
implemented in a C++ class by the same name, so this interface is very fast and light-weight.
Materials may be initialized in a number of different ways.  For example, initializing from
dictionaries of compositions are shown below.

.. code-block:: ipython

    In [1]: from pyne.material import Material

    In [2]: leu = Material({'U238': 0.96, 'U235': 0.04}, 42, 'LEU')

    In [3]: leu
    Out[3]: pyne.material.Material({922350: 0.04, 922380: 0.96}, 42.0, 'LEU')    

    In [4]: nucvec = {10010:  1.0, 80160:  1.0, 691690: 1.0, 922350: 1.0,
       ...:           922380: 1.0, 942390: 1.0, 942410: 1.0, 952420: 1.0,
       ...:           962440: 1.0}

    In [5]: mat = Material(nucvec)

    In [6]: print mat
    Material: 
    mass = 9.0
    ----------
    H1     0.111111111111
    O16    0.111111111111
    TM169  0.111111111111
    U235   0.111111111111
    U238   0.111111111111
    PU239  0.111111111111
    PU241  0.111111111111
    AM242  0.111111111111
    CM244  0.111111111111

Materials may also be initilaized from plain text or HDF5 files (see :meth:`Material.load_from_text` and
:meth:`Material.load_from_hdf5`).  Once you have a Material instance, you can always obtain the unormalized
mass vector through :meth:`Material.mult_by_mass`.  Normalization routines to normalize the mass 
:meth:`Material.normalize` or the composition :meth:`Material.norm_comp` are also available.

.. code-block:: ipython

    In [7]: leu.mult_by_mass()
    Out[7]: {922350: 1.68, 922380: 40.32}

    In [8]: mat.normalize()

    In [9]: mat.mult_by_mass()
    Out[9]: {10010: 0.111111111111, 80160: 0.111111111111, 691690: 0.111111111111, 
       ...:  922350: 0.111111111111, 922380: 0.111111111111, 942390: 0.111111111111, 
       ...:  942410: 0.111111111111, 952420: 0.111111111111, 962440: 0.111111111111}

    In [10]: mat.mass
    Out[10]: 1.0

Furthermore, various arithemtic operations between Materials and numeric types are also defined.
Adding two Materials together will return a new Material whose values are the weighted union
of the two original. Multiplying a Material by 2, however, will simply double the mass.

.. code-block:: ipython

    In [11]: other_mat = mat * 2

    In [12]: other_mat
    Out[12]: pyne.material.Material({10010: 0.111111111111, 80160: 0.111111111111, 691690: 0.111111111111, 
       ...:                          922350: 0.111111111111, 922380: 0.111111111111, 942390: 0.111111111111, 
       ...:                          942410: 0.111111111111, 952420: 0.111111111111, 962440: 0.111111111111}, 
       ...:                          2.0, '')

    In [13]: other_mat.mass
    Out[13]: 2.0

    In [14]: weird_mat = leu + mat * 18

    In [15]: print weird_mat
    Material: 
    mass = 60.0
    ----------
    H1     0.0333333333333
    O16    0.0333333333333
    TM169  0.0333333333333
    U235   0.0613333333333
    U238   0.705333333333
    PU239  0.0333333333333
    PU241  0.0333333333333
    AM242  0.0333333333333
    CM244  0.0333333333333



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

