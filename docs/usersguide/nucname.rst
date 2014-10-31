.. currentmodule:: pyne.nucname

.. _usersguide_nucname:

==========================
Nuclide Naming Conventions
==========================
One of the most frustrating aspects of nuclear data software is the large number
of different ways that people choose to nuclide names.  Functionally, there are 
three pieces of information that *should* be included in a radionuclide's name

1. **Z Number**: The number of protons.
2. **A Number**: The number of nucleons (neutrons + protons).
3. **Metastable Level**: The metastable state of the nucleus as defined by ENSDF.

Some common naming conventions exist.  The following are the ones currently 
supported by PyNE.  Functions to convert between forms may be seen in :ref:`name_cast`.

.. include:: ../nucnameforms.rst

If there are more conventions that you would like to see supported, please contact the :ref:`dev_team`.

--------------
Canonical Form
--------------
The ``zzzaaammmm`` integer form of nuclide names is the fundamental form of nuclide naming because
it accurately captures all of the needed information in the smallest amount of space.  Given that the 
Z-number may be up to three digits, A-numbers are always three digits, and the excitation level is
one digit, all possible nuclides are represented on the range ``0 <= zzzaaammmm < 10000000000``.  This 
falls well within 32 bit integers (but sadly outside of the smaller 16 bit ints).

On the other hand, ``name`` string representations may be anywhere from two characters (16 bits)
up to six characters (48 bits).  So in general, ``zzzaaammmm`` is smaller by 50%.  Other forms do 
not necessarily contain all of the required information (``MCNP``) or require additional storage 
space (``Serpent``).  It may seem pedantic to quibble over the number of bits per nuclide name, 
but these identifiers are used everywhere throughout nuclear code, so it behooves us to be as small
and fast as possible.

The other distinct advantage that integer forms have is that you can natively perform arithmetic
on them.  For example::

    # Am-242m
    nuc = 942420001

    # Z-number
    zz = nuc/10000

    # A-number
    aaa = (nuc/10)%1000

    # Meta-stable state
    m = nuc%10

Code internal to PyNE use ``zzzaaammmm``, and except for human readability, you should too!  
Natural elements are specified in this form by having zero A-number and excitation flags
(``zzzaaammmm('U') == 920000``).

Well-Defined vs Ambiguous Situations
....................................
In situations where the input naming convention is well-defined, it is *highly*
recommended that you use the direct ``<form>_to_id()`` functions (e.g. 
``mcnp_to_id()``) to convert from a nuclide in the given form to the id form 
representation. When a high level of quality assurance is required, it is 
advisable to require an specific input format to leverage the exactness of the 
direct-to-id functions.

However, in situations where arbitrary nuclide naming conventions are allowed, 
you must use the ``id()`` function. An example of such a situation is when accepting 
human input. This function attempts to resolve the underlying nuclide identifier. 
For most nuclides and most normal spellings, this resolution is straightforward. 
However, some nulcides are ambiguous between the various supported naming conventions.
For more information please refer to the 
:ref:`nucname theory manual <theorymanual_nucname>`.

---------------
Examples of Use
---------------

.. code-block:: ipython

    In [1]: from pyne import nucname

    In [2]: nucname.zzaaam('U235')
    Out[2]: 922350

    In [3]: nucname.zzaaam(10010)
    Out[3]: 10010

    In [4]: nucname.zzllaaam(942390)
    Out[4]: '95-Am-242m'

    In [5]: nucname.name(10010)
    Out[5]: 'H1'

    In [6]: nucname.serpent('AM242M')
    Out[6]: 'Am-242m'

    In [7]: nucname.name_zz['SR']
    Out[7]: 38

    In [8]: nucname.zz_name[57]
    Out[8]: 'LA'

    In [9]: nucname.alara('FE56')
    Out[9]: 'fe:56'



Further information on the naming module may be seen in the library reference :ref:`pyne_nucname`.
