.. _pyne_dbgen_scattering_lengths:

==================================================================
Neutron Scattering Lengths -- :mod:`pyne.dbgen.scattering_lengths`
==================================================================
This module provides a way to grab and store raw data for neutron scattering 
lengths.  This data comes from Neutron News, Vol. 3, No. 3, 1992, pp. 29-37 via 
a NIST webpage (http://www.ncnr.nist.gov/resources/n-lengths/list.html).  Please
contact Alan Munter, <alan.munter@nist.gov> for more information.

.. currentmodule:: pyne.dbgen.scattering_lengths

All functionality may be found in the ``scattering_lengths`` module::

    from pyne.dbgen import scattering_lengths

.. autofunction:: grab_scattering_lengths(build_dir="", file_out='scattering_lengths.html')

------------

.. autofunction:: nist_num(nist_data)

------------

.. autofunction:: parse_scattering_lengths(build_dir)

------------

.. autofunction:: make_scattering_lengths_table(nuc_data, build_dir="")

------------

.. autofunction:: make_scattering_lengths(nuc_data, build_dir)
