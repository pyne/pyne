.. _pyne_serpent:

**************************************
Serpent Support -- :mod:`pyne.serpent`
**************************************
Serpent_ is a continuous energy Monte Carlo reactor physics code.  Pyne
contains support for reading in serpents three types of output files, which are
all in Matlab's ``*.m`` format.  These files are all read in as Python
dictionaries and optionally written out to a corresponding ``*.py`` file.

.. currentmodule:: pyne.serpent

All functionality may be found in the ``serpent`` package::

 from pyne import serpent

-------
Parsers
-------

.. autofunction:: parse_res(resfile, write_py=False)

-----------

.. autofunction:: parse_dep(depfile, write_py=False, make_mats=True)

-----------

.. autofunction:: parse_det(detfile, write_py=False)

.. _Serpent: http://montecarlo.vtt.fi/
