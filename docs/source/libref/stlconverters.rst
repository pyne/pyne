.. _pyne_stlconverters:

============================================================
C++ Standard Library Converters -- :mod:`pyne.stlconverters`
============================================================
This module contains wrapper classes for commonly used constructs 
in the C++ standard library.  Becasue Cython does not yet do templating,
these classes must be declared and defined for every type.

.. currentmodule:: pyne.stlconverters

All functionality may be found in the ``stlconverters`` module::

 from pyne import stlconverters

This module is largely used by PyNE under the convers, in Cython and elsewhere.
However, these classes are of more general interest so feel free to use them in
your own code as well.

-----------
Set Proxies
-----------

.. autoclass:: SetInt

.. autoclass:: SetStr


-----------
Map Proxies
-----------

.. autoclass:: MapStrInt

.. autoclass:: MapIntStr

.. autoclass:: MapIntDouble

.. autoclass:: MapIntComplex
