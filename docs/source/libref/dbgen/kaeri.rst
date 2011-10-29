.. _pyne_dbgen_kaeri:

========================================
KAERI Helpers -- :mod:`pyne.dbgen.kaeri`
========================================
This module provides tools for scraping nuclear data off of the KAERI website (http://atom.kaeri.re.kr).
These functions are used by other parts of dbgen.  Please use with respect.

.. currentmodule:: pyne.dbgen.kaeri

All functionality may be found in the ``kaeri`` module::

    from pyne.dbgen import kaeri

.. autofunction:: grab_kaeri_nuclide(nuc, build_dir="", n=None)

------------

.. autofunction:: parse_for_natural_isotopes(htmlfile)

------------

.. autofunction:: parse_for_all_isotopes(htmlfile)
