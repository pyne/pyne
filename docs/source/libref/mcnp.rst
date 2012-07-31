.. _pyne_mcnp:

====================================================
MCNP Input and Output Interfaces -- :mod:`pyne.mcnp`
====================================================

.. currentmodule:: pyne.mcnp

.. automodule:: pyne.mcnp
  
There is a class for a variety of types of files that MCNP produces.  The
functionality of the module can be obtained by importing as such::

    from pyne import mcnp

The userguide page :ref:`usersguide_mcnp` includes an example of the usage of
the `mcnp.Inp` class.

**************
Inp Class
**************

The :ref:`pyne_mcnpcard` module is used heavily by `mcnp.Inp`, but the user
should not interact with the module on their own. `mcnpcard` contains a class
for each card in MCNP.

.. autoclass:: Inp
   :members:
   :inherited-members:

