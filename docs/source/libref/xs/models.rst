.. _pyne_xs_models:

=============================================
Cross Section Models -- :mod:`pyne.xs.models`
=============================================
This module provides functions to compute cross sections, and related quantities, from 
fundamental physical models.  From an empirical perspective, physical models are not a 
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


******************
Physical Constants
******************
.. data:: k

    Boltzmann's constant in [MeV/K].

.. data:: m_n

    Neutron mass in [amu].

*********************
Partial Energy Matrix
*********************
When collapsing a flux spectrum or multigroup cross-section from one more finely 
resolved group structure (**E_n**) to another coarser structure (**E_g**) it is 
useful to define a partial energy operator.  This describes the linear contribution
of each nth group to each gth group.  Calling this matrix **P_E**, it thereby solves
the equation:

.. math::

    \phi_g = P_E \cdot \phi_n

Note that for this transform to be valid both group structures must be monotonic in the 
same direction.  All elements of **P_E** live on the range [0,1].

.. autofunction:: partial_energy_matrix(E_g, E_n)

---------

.. autofunction:: partial_energy_matrix_mono(E_g, E_n, slope=-1)

************************
Group Collapse Functions
************************
.. autofunction:: phi_g(E_g, E_n, phi_n)

---------

.. autofunction:: group_collapse(sigma_n, phi_n, phi_g=None, partial_energies=None, E_g=None, E_n=None)


***************
Physical models
***************
.. autofunction:: chi(E)

---------

.. autofunction:: alpha(E_prime, E, theta, M_A=1.0, T=300.0)

---------

.. autofunction:: beta(E_prime, E, T=300.0)

---------

.. autofunction:: alpha_at_theta_0(E_prime, E, M_A=1.0, T=300.0)

---------

.. autofunction:: alpha_at_theta_pi(E_prime, E, M_A=1.0, T=300.0)

---------

.. autofunction:: one_over_gamma_squared(E)

---------

.. autofunction:: E_prime_min(E, M_A=1.0)

---------

.. autofunction:: sigma_s_const(b)

---------

.. autofunction:: sigma_s(E, b=1.0, M_A=1.0, T=300.0)
