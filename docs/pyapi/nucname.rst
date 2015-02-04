.. _pyne_nucname:

============================================
Nuclide Naming Module -- :mod:`pyne.nucname`
============================================
This package is used to convert between various nuclide naming schemes.  
Currently the following naming conventions are supported:

.. _name_forms:

.. include:: ../nucnameforms.rst

.. currentmodule:: pyne.nucname

All functionality may be found in the ``nucname`` package::

 from pyne import nucname

This contains several zzaaam, zzzaaa, zzllaaam, name, MCNP, Groundstate and Serpent 
converter function as well as other helpful module attributes.

.. _name_cast:

-----------------------------------
Naming Convention Casting Functions
-----------------------------------

.. autofunction:: id(nuc)

-----

.. autofunction:: name(nuc)

-----

.. autofunction:: zzaaam(nuc)

-----

.. autofunction:: zzzaaa(nuc)

-----

.. autofunction:: zzllaaam(nuc)

-----

.. autofunction:: mcnp(nuc)

-----

.. autofunction:: serpent(nuc)

-----

.. autofunction:: nist(nuc)

-----

.. autofunction:: cinder(nuc)

-----

.. autofunction:: alara(nuc)

-----

.. autofunction:: groundstate(nuc)


-----------------------------------
Id Conversion Functions
-----------------------------------

.. autofunction:: zzaaam_to_id(nuc)

-----

.. autofunction:: zzzaaa_to_id(nuc)

---

.. autofunction:: zzllaaam_to_id(nuc)

---

.. autofunction:: mcnp_to_id(nuc)

-----

.. autofunction:: serpent_to_id(nuc)

-----

.. autofunction:: nist_to_id(nuc)

-----

.. autofunction:: cinder_to_id(nuc)

-----

.. autofunction:: alara_to_id(nuc)


-----------------------------------
Number Functions
-----------------------------------

.. autofunction:: znum(nuc)

-----

.. autofunction:: anum(nuc)

-----

.. autofunction:: snum(nuc)

-----------------------
Conversion Dictionaries
-----------------------

.. attribute:: name_zz

   Dictionary that is used to convert an elemental symbol (str) to its charge Z-number (int).
   For example::

      nucname.name_zz["HE"] = 2
      nucname.name_zz["U"]  = 92


.. attribute:: zz_name

   Dictionary that is used to convert a charge Z-number (int) to its elemental symbol (str).
   For example::

      nucname.name_zz[1]  = "H"
      nucname.name_zz[94] = "PU"



---------------------
Element Groups (name)
---------------------
Element groups for the Lanthanides, Actinides, Transuranics, Minor Actinides,
and Fission Products.

.. attribute:: nucname.LAN

.. attribute:: nucname.ACT

.. attribute:: nucname.TRU

.. attribute:: nucname.MA

.. attribute:: nucname.FP

The groups are defined as follows::

   nucname.LAN = set(['CE', 'DY', 'ER', 'EU', 'GD', 'HO', 'LA', 'LU', 'ND', 'PM', 'PR', 'SM', 'TB', 'TM', 'YB'])

   nucname.ACT = set(['AC', 'AM', 'BK', 'CF', 'CM', 'ES', 'FM', 'LR', 'MD', 'NO', 'NP', 'PA', 'PU', 'TH', 'U'])

   nucname.TRU = set(['AM', 'BH', 'BK', 'CF', 'CM', 'CN', 'DB', 'DS', 'ES', 'FL', 'FM', 'HS', 'LR', 'LV', 'MD', 
                      'MT', 'NO', 'NP', 'PU', 'RF', 'RG', 'SG'])

   nucname.MA  = set(['AM', 'BK', 'CF', 'CM', 'ES', 'FM', 'LR', 'MD', 'NO', 'NP'])

   nucname.FP  = set(['AG', 'AL', 'AR', 'AS', 'AT', 'AU', 'B',  'BA', 'BE', 'BI', 'BR', 'C',  'CA', 'CD', 'CE',
                            'CL', 'CO', 'CR', 'CS', 'CU', 'DY', 'ER', 'EU', 'F',  'FE', 'FR', 'GA', 'GD', 'GE',
                            'H',  'HE', 'HF', 'HG', 'HO', 'I',  'IN', 'IR', 'K',  'KR', 'LA', 'LI', 'LU', 'MG',
                            'MN', 'MO', 'N',  'NA', 'NB', 'ND', 'NE', 'NI', 'O',  'OS', 'P',  'PB', 'PD', 'PM',
                            'PO', 'PR', 'PT', 'RA', 'RB', 'RE', 'RH', 'RN', 'RU', 'S',  'SB', 'SC', 'SE', 'SI',
                            'SM', 'SN', 'SR', 'TA', 'TB', 'TC', 'TE', 'TI', 'TL', 'TM', 'V',  'W',  'XE',  'Y',
                            'YB', 'ZN', 'ZR'])


-------------------
Element Groups (zz)
-------------------
Element groups for the Lanthanides, Actinides, Transuranics, Minor Actinides,
and Fission Products.

.. attribute:: nucname.lan

.. attribute:: nucname.act

.. attribute:: nucname.tru

.. attribute:: nucname.ma

.. attribute:: nucname.fp

The groups are defined as follows::

   nucname.lan = set([57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71])

   nucname.act = set([89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103])

   nucname.tru = set([93,  94,  95,  96,  97,  98, 99, 100, 101, 102, 103, 104, 105,
                      106, 107, 108, 109, 110, 111, 112, 114, 116])

   nucname.ma  = set([93, 95, 96, 97, 98, 99, 100, 101, 102, 103])

   nucname.fp  = set([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
                     29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,
                     54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78,
                     79, 80, 81, 82, 83, 84, 85, 86, 87, 88])

