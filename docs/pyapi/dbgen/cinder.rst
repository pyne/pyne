.. _pyne_dbgen_cinder:

=======================================
CINDER Data -- :mod:`pyne.dbgen.cinder`
=======================================
This module locates, parses, and adds CINDER cross section and fission product yield data to ``nuc_data.h5``.
Note that this module requires that the ``cinder.dat`` file exist within the ``DATAPATH`` directory.  
This often requires that MCNPX is installed.

.. currentmodule:: pyne.dbgen.cinder

All functionality may be found in the ``cinder`` module::

    from pyne.dbgen import cinder

CINDER Data API
===============

.. automodule:: pyne.dbgen.cinder
    :members:
