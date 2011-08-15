.. _pyne_cccc:

================================
CCCC Formats -- :mod:`pyne.cccc`
================================

.. currentmodule:: pyne.cccc

The ``cccc`` module contains a number of classes for reading various cross
section, flux, geometry, and data files with specifications given by the
Committee for Computer Code Coordination. The following types of files can be
read using classes from this module: ISOTXS, DLAYXS, BRKOXS, RTFLUX, ATFLUX,
RZFLUX, MATXS, and SPECTR.

All functionality may be found in the ``cccc`` package::

    from pyne import cccc

The information below represents the complete specification of the classes in
the cccc module. For examples of usage, please refer to the User's Guide entry
for :ref:`usersguide_cccc`.

************
Isotxs Class
************

.. autoclass:: Isotxs(filename)
   :members:

************
Dlayxs Class
************

.. autoclass:: Dlayxs(filename)
   :members:

************
Brkoxs Class
************

.. autoclass:: Brkoxs(filename)

************
Rtflux Class
************

.. autoclass:: Rtflux(filename)

************
Atflux Class
************

.. autoclass:: Atflux(filename)

************
Rzflux Class
************

.. autoclass:: Rzflux(filename)

************
Matxs Class
************

.. autoclass:: Matxs(filename)
