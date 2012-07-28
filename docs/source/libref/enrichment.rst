.. _pyne_enrichment:

============================================
Enrichment -- :mod:`pyne.enrichment`
============================================
.. automodule:: pyne.enrichment


******************
Cascades
******************
.. autoclass:: Cascade
    :members: alpha, Mstar, j, k, N, M, x_feed_j, x_prod_j, x_tail_j, mat_feed, 
              mat_prod, mat_tail, l_t_per_feed, swu_per_feed, swu_per_prod

----------

.. autofunction:: default_uranium_cascade

*******************
Enrichment Solvers
*******************
.. autofunction:: multicomponent

------------

.. autofunction:: ltot_per_feed

****************
Helper Functions
****************
.. autofunction:: prod_per_feed

------------

.. autofunction:: tail_per_feed

------------

.. autofunction:: tail_per_prod

------------

.. autofunction:: alphastar_i
