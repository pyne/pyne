.. _pyne_xs_models:

=============================================
Cross Section Models -- :mod:`pyne.xs.models`
=============================================
This module provides functions to compute cross sections, and related quantities, from 
fundamental physical models.  From an empirical perspective, physical models are not as 
valuable as experimental data.  However, these functions exist to be called upon when
experimental data is not available.

.. currentmodule:: pyne.xs.models

All functionality may be found in the ``models`` module::

    from pyne.xs import models

Many of these functions may be called with either scalar or vector arguments.

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

Models API
-----------

.. automodule:: pyne.xs.models
    :members:
