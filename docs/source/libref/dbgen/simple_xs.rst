.. _pyne_dbgen_simple_xs:

========================================================
Simple Cross Sections -- :mod:`pyne.dbgen.simple_xs`
========================================================
This module adds simple cross section data from KAERI to ``nuc_data.h5``.

.. currentmodule:: pyne.dbgen.simple_xs

All functionality may be found in the ``simple_xs`` module::

    from pyne.dbgen import simple_xs

.. autofunction:: grab_kaeri_simple_xs(build_dir="")

------------

.. autofunction:: get_xs_from_file(filename, eng, chan)

------------

.. autofunction:: parse_simple_xs(build_dir="")

------------

.. autofunction:: make_simple_xs_tables(nuc_data, build_dir="")

------------

.. autofunction:: make_simple_xs(nuc_data, build_dir)

