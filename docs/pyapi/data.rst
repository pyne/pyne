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

-----------
Atomic Mass
-----------

.. autofunction:: atomic_mass(nuc)

.. data:: atomic_mass_map

    A mapping from zzaaam-nuclides to their mass weights.
    This is used by :func:`atomic_mass` under the hood.

------------------------
Natural Abundance Ratios
------------------------

.. autofunction:: natural_abund(nuc)

.. data:: natural_abund_map

    A mapping from zzaaam-nuclides to their mass weights.
    This is used by :func:`natural_abund` under the hood.

-----------
Atomic Data
-----------

.. autofunction:: calculate_xray_data(nuc, k_conv, l_conv) 

----------
Decay Data
----------

.. autofunction:: half_life(nuc)

    Half life reported for a level in ENSDF level data.

.. autofunction:: decay_const(nuc)

.. autofunction:: branch_ratio(from_nuc, to_nuc, use_metastable=True)

    Branch ratio reported for a reaction in ENSDF level data. This is not
    necessarily equivalent to the value reported in a decay dataset by
    decay_branch_ratio but it should be.

.. autofunction:: state_energy(nuc, use_metastable=True)

.. autofunction:: decay_children(nuc, use_metastable=True)

.. autofunction:: id_from_level(nuc, level, special="")

.. autofunction:: metastable_id(nuc, level=1)

.. autofunction:: decay_half_life(from_nuc, to_nuc)

.. autofunction:: decay_half_life_byparent(parent)

.. autofunction:: decay_branch_ratio(from_nuc, to_nuc)

.. autofunction:: decay_branch_ratio_byparent(parent)

.. autofunction:: decay_photon_branch_ratio(from_nuc, to_nuc)

.. autofunction:: decay_photon_branch_ratio_byparent(parent)

.. autofunction:: decay_beta_branch_ratio(from_nuc, to_nuc)

.. autofunction:: decay_beta_branch_ratio_byparent(parent)

.. autofunction:: gamma_energy(parent)

.. autofunction:: gamma_photon_intensity(parent)

    This intensity is just the relative intensity for different photons in the
    same dataset. It needs to be multiplied by the output of
    :func:`decay_photon_branch_ratio` to give the percentage of parent decays.

.. autofunction:: gamma_conversion_intensity(parent)

.. autofunction:: gamma_total_intensity(parent)

.. autofunction:: gamma_from_to_byparent(parent)

.. autofunction:: gamma_from_to_byen(en, enerror=None)

.. autofunction:: gamma_parent(en, enerror=None)

.. autofunction:: gamma_xrays(parent)

.. autofunction:: alpha_energy(parent)

.. autofunction:: alpha_intensity(parent)

.. autofunction:: alpha_parent(en, enerror=None)

.. autofunction:: alpha_child_byen(en, enerror=None)

.. autofunction:: alpha_child_byparent(parent)

.. autofunction:: beta_endpoint_energy(parent)

.. autofunction:: beta_average_energy(parent)

.. autofunction:: beta_intensity(parent)

.. autofunction:: beta_parent(en, enerror=None)

.. autofunction:: beta_child_byen(en, enerror=None)

.. autofunction:: beta_child_byparent(parent)

.. autofunction:: ecbp_endpoint_energy(parent)

.. autofunction:: ecbp_average_energy(parent)

.. autofunction:: ec_intensity(parent)

.. autofunction:: beta_plus_intensity(parent)

.. autofunction:: ecbp_parent(en, enerror=None)

.. autofunction:: ecbp_child_byen(en, enerror=None)

.. autofunction:: ecbp_child_byparent(parent)

.. autofunction:: ecbp_xrays(parent)

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

