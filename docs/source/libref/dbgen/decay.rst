.. _pyne_dbgen_decay:

=================================================
Radioactive Decay Data -- :mod:`pyne.dbgen.decay`
=================================================
This module adds radioactive decay data via ENSDF to ``nuc_data.h5``.

.. currentmodule:: pyne.dbgen.decay

All functionality may be found in the ``decay`` module::

    from pyne.dbgen import decay

.. autofunction:: grab_ensdf_decay(build_dir="")

------------

.. autofunction:: parse_decay(build_dir="")

------------

.. autofunction:: make_atomic_decay_table(nuc_data, build_dir="")

------------

.. autofunction:: make_decay(nuc_data, build_dir)

