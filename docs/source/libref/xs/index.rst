.. _pyne_xs:

=========================================
Cross Section Interface -- :mod:`pyne.xs`
=========================================
Pyne provides a top-level interface for computing (and caching) multigroup neutron cross sections.
These cross sections will be computed from a variety of available data sources (stored in nuc_data.h5).
In order of preference: 

#. CINDER 63-group cross sections,
#. A two-point fast/thermal interpolation (using 'simple_xs' data from KAERI),
#. or physical models.

Top-level functionality may be be found in this package's API module::

    from pyne.xs.api import *

In the future, this package should support generating multigroup cross sections from user-specified 
pointwise data sources (such as ENDF or ACE files).

**Cross Section Modules**

.. toctree::
    :maxdepth: 1

    models
    cache
    channels
