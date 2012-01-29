.. _pyne_data:

======================================
Basic Nuclear Data -- :mod:`pyne.data`
======================================
This module provides a top-level interface for a variety of basic nuclear data needs.
This aims to provide quick access to very high fidelity data.  Usually values
are taken from the :file:`nuc_data.h5` library.

.. currentmodule:: pyne.data

All functionality may be found in the ``data`` module::

 from pyne import data

-------------
Atomic Weight
-------------

.. autofunction:: nuc_weight(nuc)

.. data:: nuc_weight_map

    A mapping from zzaaam-nuclides to their mass weights.
    This is used by :func:`nuc_weight` under the hood.


----------
Decay Data
----------

.. autofunction:: half_life(nuc)

.. data:: half_life_map

    A mapping from zzaaam-nuclides to their half lives.
    This is used by :func:`half_life` under the hood.


.. autofunction:: decay_const(nuc)

.. data:: decay_const_map

    A mapping from zzaaam-nuclides to their decay constants.
    This is used by :func:`decay_const` under the hood.


--------------------------
Neutron Scattering Lengths
--------------------------

.. autofunction:: b(nuc)

.. data:: b_map

    A mapping from zzaaam-nuclides to their bound neuton scattering lengths.
    This is used by :func:`b` under the hood.

-----

.. autofunction:: b_coherent(nuc)

.. data:: b_coherent_map

    A mapping from zzaaam-nuclides to their bound coherent neuton scattering lengths.
    This is used by :func:`b_coherent` under the hood.

-----

.. autofunction:: b_incoherent(nuc)

.. data:: b_incoherent_map

    A mapping from zzaaam-nuclides to their bound incoherent neuton scattering lengths.
    This is used by :func:`b_coherent` under the hood.

