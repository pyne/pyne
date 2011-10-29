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

.. autofunction:: grab_cinder_dat(build_dir="")

------------

.. autofunction:: get_group_sizes(raw_data)

------------

.. autofunction:: make_mg_group_structure(nuc_data, build_dir="")

------------

.. autofunction:: make_mg_absorption(nuc_data, build_dir="")

------------

.. autofunction:: make_mg_fission(nuc_data, build_dir="")

------------

.. autofunction:: make_mg_gamma_decay(nuc_data, build_dir="")

------------

.. autofunction:: get_fp_sizes(raw_data)

------------

.. autofunction:: parse_neutron_fp_info(raw_data)

------------

.. autofunction:: make_neutron_fp_info(nuc_data, build_dir="")

------------

.. autofunction:: grab_photon_fp_info(raw_data)

------------

.. autofunction:: make_photon_fp_info(nuc_data, build_dir="")

------------

.. autofunction:: make_neutron_fp_yields(nuc_data, build_dir="")

------------

.. autofunction:: make_photon_fp_yields(nuc_data, build_dir)

------------

.. autofunction:: make_cinder(nuc_data, build_dir)
