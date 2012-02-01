.. _pyne_xs_channels:

====================================================
Nuclear Reaction Channels -- :mod:`pyne.xs.channels`
====================================================
This module provides an easy interface to very quickly grab multi-group cross sections
from the cross section cache and collapse them to the appropriate group structure.
The more finely resolved group constants

.. currentmodule:: pyne.xs.channels

All functionality may be found in the ``channels`` module::

    from pyne.xs import channels

The following terminology applies for this module:

* **G**: The number of low-resolution energy groups, indexed by **g**.
* **N**: The number of high-resolution energy groups, indexed by **n**.
* **E_g**: The low-resolution group structure [MeV].  Must be monotonic.
* **E_n**: The high-resolution group structure [MeV].  Must be monotonic 
  in the same direction as **E_g**.
* **phi**: Neutron flux [n/cm^2/s]. 
* **E**: The incident neutron energy [MeV].
* **E'**: The exiting neutron energy [MeV].
* **theta**: The scattering angle [radians].
* **M_A**: The atomic mass [amu] of the target material.
* **T**: The temperature [kelvin] of the target material.


********
Channels
********
.. autofunction:: sigma_t(nuc, T=300.0, E_g=None, E_n=None, phi_n=None)

-----

.. autofunction:: sigma_f(nuc, E_g=None, E_n=None, phi_n=None)

-----

.. autofunction:: chi(nuc, E_g=None, E_n=None, phi_n=None, eres=101)

-----

.. autofunction:: sigma_s(nuc, T, E_g=None, E_n=None, phi_n=None)

-----

.. autofunction:: sigma_s_gh(nuc, T, E_g=None, E_n=None, phi_n=None)

-----

.. autofunction:: sigma_a(nuc, T, E_g=None, E_n=None, phi_n=None)

-----

.. autofunction:: sigma_a_reaction(nuc, rx, E_g=None, E_n=None, phi_n=None)

-----

.. autofunction:: metastable_ratio(nuc, rx, E_g=None, E_n=None, phi_n=None)
