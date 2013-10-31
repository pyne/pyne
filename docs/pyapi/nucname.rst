.. _pyne_nucname:

============================================
Nuclide Naming Module -- :mod:`pyne.nucname`
============================================
This package is used to convert between various nuclide naming schemes.  
Currently the following naming conventions are supported:

.. _name_forms:

 #. **zzaaam**: This type places the charge of the nucleus out front, then has three
    digits for the atomic mass number, and ends with a metastable flag (0 = ground,
    1 = first excited state, 2 = second excited state, etc).  Uranium-235 here would
    be expressed as '922350'.
 #. **name**: This is the more common, human readable notation.  The chemical symbol
    (one or two characters long) is first, followed by the atomic weight.  Lastly if
    the nuclide is metastable, the letter *M* is concatenated to the end.  For
    example, 'H-1' and 'Am242M' are both valid.  Note that nucname will always
    return name form with the dash removed and all letters uppercase.
 #. **MCNP**: The MCNP format for entering nuclides is unfortunately
    non-standard.  In most ways it is similar to zzaaam form, except that it
    lacks the metastable flag.  For information on how metastable isotopes are
    named, please consult the MCNPX documentation for more information.
 #. **Serpent**: The serpent naming convention is similar to name form.
    However, only the first letter in the chemical symbol is uppercase, the
    dash is always present, and the the meta-stable flag is lowercase.  For
    instance, 'Am-242m' is the valid serpent notation for this nuclide.
 #. **NIST**: The NIST naming convention is also similar to the Serpent form.
    However, this convention contains no metastable information.  Moreover, the
    A-number comes before the element symbol.  For example, '242Am' is the
    valid NIST notation.
 #. **CINDER**: The CINDER format is similar to zzaaam form except that the
    placement of the Z- and A-numbers are swapped. Therefore, this format is
    effectively aaazzzm.  For example, '2420951' is the valid cinder notation
    for 'AM242M'.
 #. **ALARA**: In ALARA format, elements are denoted by the lower case atomic symbol. Isotopes are
    specified by appending a semicolon and A-number. For example, "fe" and "fe:56" represent
    elemental iron and iron-56 respectively. No metastable flag exists.

.. currentmodule:: pyne.nucname

All functionality may be found in the ``nucname`` package::

 from pyne import nucname

This contains several zzaaam, name, MCNP and Serpent converter function as
well as other helpful module attributes.

.. _name_cast:

-----------------------------------
Naming Convention Casting Functions
-----------------------------------

.. autofunction:: zzaaam(nuc)

-----

.. autofunction:: name(nuc)

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



-----------------------
Other nucname functions
-----------------------

.. autofunction:: current_form(nuc)

