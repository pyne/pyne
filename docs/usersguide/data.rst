.. _usersguide_data:

============
Atomic & Isotopic Data
============

.. currentmodule:: pyne.data

.. automodule:: pyne.data

For a complete specification for the classes in the ``cccc`` module, please
refer to the Library Reference entry for :ref:`pyne_cccc`.

***************************
Example Use of Data Class
***************************

The isotopic & elemental abundance data is usually loaded at the time the atomic_mass function is launched. This data now
exists in the C++ implementaion, behinds the scenes the nuc_data.h5 file is loaded and the appropriate data elements populated.


.. code-block:: ipython

   In [1]: from pyne.data import atomic_mass

   In [2]: print atomic_mass('2H')
   2.01410177812

If one chooses to, the nuclear data path can be set to some nonesense value and the data will be initialised from memory
instead.

.. code-block:: ipython

   In [1]: from pyne.pyne_config import pyne_conf
   
   In [2]: from pyne.data import atomic_mass
   
   In [3]: pyne_conf.NUC_DATA_PATH = b'some silly path that doesnt exist'

   In [4]: print atomic_mass('2H')
   2.01410177812

This is not particularly useful for Python users of PyNE, however the pure C++ users can now use the atomic and isotopic data from the amalgamated source and not need to carry the nuc_data.h5 with you.
