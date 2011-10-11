.. _usersguide_material:

=========
Materials
=========
This module contains the Material class, which is used to represent nuclear
materials throughout PyNE.

.. currentmodule:: pyne.material

All functionality may be found in the ``material`` package::

 from pyne.material import Material

Materials are the primary container for radionuclides. They map nuclides to **mass weights**,
though they contain methods for converting to/from atom fractions as well.
In many ways they take inspiration from numpy arrays and python dictionaries.  Materials
have two main attributes which define them.

1. **comp**: a normalized composition mapping from nuclides (zzaaam-ints) to mass-weights (floats).
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
    Out[3]: pyne.material.Material({922350: 0.04, 922380: 0.96}, 42.0, 'LEU', -1.0)

    In [4]: nucvec = {10010:  1.0, 80160:  1.0, 691690: 1.0, 922350: 1.0,
       ...:           922380: 1.0, 942390: 1.0, 942410: 1.0, 952420: 1.0,
       ...:           962440: 1.0}

    In [5]: mat = Material(nucvec)

    In [6]: print mat
    Material: 
    mass = 9.0
    atoms per molecule = -1.0
    -------------------------
    H1     0.111111111111
    O16    0.111111111111
    TM169  0.111111111111
    U235   0.111111111111
    U238   0.111111111111
    PU239  0.111111111111
    PU241  0.111111111111
    AM242  0.111111111111
    CM244  0.111111111111

Materials may also be initialized from plain text or HDF5 files (see :meth:`Material.from_text` and
:meth:`Material.from_hdf5`).  Once you have a Material instance, you can always obtain the unnormalized
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


Material Arithmetic
-------------------
Furthermore, various arithmetic operations between Materials and numeric types are also defined.
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
    atoms per molecule = -1.0
    -------------------------
    H1     0.0333333333333
    O16    0.0333333333333
    TM169  0.0333333333333
    U235   0.0613333333333
    U238   0.705333333333
    PU239  0.0333333333333
    PU241  0.0333333333333
    AM242  0.0333333333333
    CM244  0.0333333333333


Raw Attribute Access
--------------------
You may also change the attributes of a material directly without generating a new material instance.

.. code-block:: ipython

    In [16]: other_mat.mass = 10

    In [17]: other_mat.name = "Other Material"

    In [18]: other_mat.comp = {'H2': 3, 922350: 15.0}

    In [19]: print other_mat
    Material: Other Material
    mass = 10.0
    atoms per molecule = -1.0
    -------------------------
    H2     3.0
    U235   15.0

Of course when you do this you have to be careful because the composition and mass may now be out
of sync.  This may always be fixed with normalization.

.. code-block:: ipython

    In [20]: other_mat.norm_comp()

    In [21]: print other_mat
    Material: Other Material
    mass = 10.0
    atoms per molecule = -1.0
    -------------------------
    H2     0.166666666667
    U235   0.833333333333


Indexing & Slicing
------------------
Additionally (and very powerfully!), you may index into either the material or the composition 
to get, set, or remove sub-materials.  Generally speaking, the composition you may only index 
into by integer-key and only to retrieve the normalized value.  Indexing into the material allows the 
full range of operations and returns the unnormalized mass weight.  Moreover, indexing into
the material may be performed with integer-keys, string-keys, slices, or sequences of nuclides.

.. code-block:: ipython

    In [22]: leu.comp[922350]
    Out[22]: 0.04

    In [23]: leu['U235']
    Out[23]: 1.68

    In [24]: weird_mat['U':'Am']
    Out[24]: pyne.material.Material({922350: 0.0736, 922380: 0.8464, 942390: 0.04, 942410: 0.04}, 50.0, '', -1.0)

    In [25]: other_mat[:920000] = 42.0

    In [26]: print other_mat
    Material: Other Material
    mass = 50.3333333333
    atoms per molecule = -1.0
    -------------------------
    H2     0.834437086093
    U235   0.165562913907

    In [27]: del mat[962440, 'TM169', 'Zr90', 80160]

    In [28]: mat[:]
    Out[28]: pyne.material.Material({10010: 0.166666666667, 922350: 0.166666666667, 922380: 0.166666666667, 
       ...:                          942390: 0.166666666667, 942410: 0.166666666667, 952420: 0.166666666667}, 
       ...:                          0.666666666667, '', -1.0)

Other methods also exist for obtaining commonly used sub-materials, such as gathering the Uranium or 
Plutonium vector.  


Molecular Weights & Atom Fractions
----------------------------------
You may also calculate the molecular weight of a material via the :meth:`Material.molecular_weight` method.
This uses the :func:`pyne.data.nuc_weight` function to look up the atomic weight values of
the constituent nuclides.

.. code-block:: ipython

    In [29]: leu.molecular_weight()
    Out[29]: 237.9290388038301

Note that by default, materials are assumed to have one atom per molecule.  This is a poor
assumption for more complex materials.  For example, take water.  Without specifying the 
number of atoms per molecule, the molecular weight calculation will be off by a factor of 3.
This can be remedied by passing the correct number to the method.  If there is no other valid
number of molecules stored on the material, this will set the appropriate attribute on the 
class.

.. code-block:: ipython

    In [30]: h2o = Material({10010: 0.11191487328808077, 80160: 0.8880851267119192})

    In [31]: h2o.molecular_weight()
    Out[31]: 6.003521561343334

    In [32]: h2o.molecular_weight(3.0)
    Out[32]: 18.01056468403

    In [33]: h2o.atoms_per_mol
    Out[33]: 3.0

It is often also useful to be able to convert the current mass-weighted material to 
an atom fraction mapping.  This can be easily done via the :meth:`Material.to_atom_frac`
method.  Continuing with the water example, if the number of atoms per molecule is 
properly set then the atom fraction return is normalized to this amount.  Alternatively, 
if the atoms per molecule are set to its default state on the class, then a truly 
fractional number of atoms is returned.

.. code-block:: ipython

    In [34]: h2o.to_atom_frac()
    Out[34]: {10010: 2.0, 80160: 1.0}

    In [35]: h2o.atoms_per_mol = -1.0

    In [36]: h2o.to_atom_frac()
    Out[36]: {10010: 0.666666666667, 80160: 0.333333333333}

Additionally, you may wish to convert the an existing set of atom fractions to a 
new material stream.  This can be done with the :meth:`Material.from_atom_frac` method, 
which will clear out the current contents of the material's composition and replace
it with the mass-weighted values.  Note that 
when you initialize a material from atom fractions, the sum of all of the atom fractions
will be stored as the atoms per molecule on this class.  Additionally, if a mass is not 
already set on the material, the molecular weight will be used.

.. code-block:: ipython

    In [37]: h2o_atoms = {10010: 2.0, 'O16': 1.0}

    In [38]: h2o = Material()

    In [39]: h2o.from_atom_frac(h2o_atoms)

    In [40]: h2o.comp
    Out[40]: {10010: 0.111914873288, 80160: 0.888085126712}

    In [41]: h2o.atoms_per_mol
    Out[41]: 3.0

    In [42]: h2o.mass
    Out[42]: 18.01056468403

    In [43]: h2o.molecular_weight()
    Out[43]: 18.01056468403


Moreover, other materials may also be used to specify a new material from atom fractions.
This is a typical case for reactors where the fuel vector is convolved inside of another 
chemical form.  Below is an example of obtaining the Uranium-Oxide material from Oxygen
and low-enriched uranium.

.. code-block:: ipython

    In [44]: uox = Material()

    In [45]: uox.name = "UOX"

    In [46]: uox.from_atom_frac({leu: 1.0, 'O16': 2.0})

    In [47]: print uox
    Material: UOX
    mass = 269.918868043
    atoms per molecule = 3.0
    ------------------------
    O16    0.118516461895
    U235   0.0352593415242
    U238   0.846224196581

.. note:: Materials may be used as keys in a dictionary because they are hashable.


Further information on the Material class may be seen in the library reference :ref:`pyne_material`.

