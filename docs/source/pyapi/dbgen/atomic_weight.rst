.. _pyne_dbgen_atomic_weight:

=================================================
Atomic Weights -- :mod:`pyne.dbgen.atomic_weight`
=================================================
This module adds high-fidelity atomic weights and natural abundance data to ``nuc_data.h5``.

.. currentmodule:: pyne.dbgen.atomic_weight

All functionality may be found in the ``atomic_weight`` module::

    from pyne.dbgen import atomic_weight

.. autofunction:: grab_kaeri_atomic_abund(build_dir="")

------------

.. autofunction:: parse_atomic_abund(build_dir="")

------------


.. autofunction:: make_atomic_weight_table(nuc_data, build_dir="")

------------

.. autofunction:: make_atomic_weight(nuc_data, build_dir)
