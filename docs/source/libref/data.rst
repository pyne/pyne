.. _pyne_data:

======================================
Basic Nuclear Data -- :mod:`pyne.data`
======================================
This module provides a top-level interface for a variety of basic nuclear data needs.
This aims to provide quick access to very high fidelity data.  Usually values
are taken from the :file:`nuc_data.h5` library.

.. currentmodule:: pyne.data

All functionality may be found in the ``data`` module::

 from pyne import nucname

-------------
Atomic weight
-------------

.. autofunction:: nuc_weight(nuc)

.. data:: nuc_weight_map

    A mapping from zzaaam-nuclides to their mass weights.
    This is used by :func:`nuc_weight` under the hood.

