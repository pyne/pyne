.. _pyne_xs_cache:

===========================================
Cross Section Cache -- :mod:`pyne.xs.cache`
===========================================
This module contains the ``XSCache`` class, a singleton instance of this class ``xs_cache``, 
and this classes associated helper functions.  Largely, these helper functions are responsible
for pulling the appropriate data out of ``nuc_data.h5``.  

.. currentmodule:: pyne.xs.cache

All functionality may be found in the ``cache`` module::

 from pyne.xs import cache

The reason for using a cache for multi-group cross sections is to prevent excessive re-computation
and to be able to automatically invalidate all data when then group structure changes.

The following terminology applies for this module:

* **G**: The number of user specified energy groups, indexed by **g**.
* **N**: The number of high-resolution energy groups, indexed by **n**.
* **E_g**: The user specified group structure.  This is often at a 
  lower resolution than what is present in the data library.  Must be 
  monotonic.
* **E_n**: The high-resolution group structure.  This is often comes from 
  the CINDER data within ``nuc_data.h5``.  Must be monotonic in the same
  direction as **E_g**.  (Note: the default cinder group structure is 
  monotonically decreasing.)

The cache will allow you to circumvent reading data from cinder if you supply it 
with your own **E_n** and your own high-resolution cross sections.

**************
XSCache Class
**************
.. autoclass:: XSCache()

.. data:: xs_cache 

    A singleton instance of XSCache.  Rather than creating your own instance,
    you should just import and use this one::

        # API
        from pyne.xs.api import xs_cache

        # Directly
        from pyne.xs.cache import xs_cache


**********************
Cache Helper Functions
**********************
The following functions help the cache extract data from the ``nuc_data.h5`` file.

.. autofunction:: get_E_n()

---------

.. autofunction:: get_sigma_f_n(nuc)

---------

.. autofunction:: get_sigma_a_n(nuc)

---------

.. autofunction:: get_sigma_a_reaction_n(nuc, rx)

.. data:: ABSORPTION_RX

    Acceptable absorption reaction strings::

        set(["", 'np *', 'a  *', 'h  *', '2p *', '3n *', 'd  *', 'np/d',
             'na', '*', 'nd', 'g  *', '3n', 'np', 'nt', 't', 'nt *',
             'x  *', '4n *', 'na *', 'nd *', 't  *', 'a', 'c', '2p', 'd',
             'g', 'h', 'n', '4n', 'p', 'n  *', '2a', '2n *', 'x', '2n',
             'nh *', 'p  *'])

.. data:: ABSORPTION_RX_MAP 

    Aliases for absorption reactions::

        {'neutron': 'n',
         'gamma': 'g',
         'alpha': 'a',
         'proton': 'p',
         'trit': 't',
         'triton': 't',
         'deut': 'd',
         'deuteron': 'd',
        }

