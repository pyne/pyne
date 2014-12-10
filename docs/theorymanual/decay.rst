.. _theorymanual_decay:

===============================
Decay 
===============================

:Author: Anthony Scopatz

.. currentmodule:: pyne.decay

The Bateman equations governing radioactive decay are an important subexpression 
of generalized transmutation equations. In many cases, it is desirable to compute
decay on its own, outside of the presence of an neutron or photon field.  In this
case radioactive decay is a function solely on intrinsic physical parameters, 
namely half-lives. This document recasts the Bateman equations into a form that 
is better suited for computation than the traditional expression.

****************************************
Canonical Bateman Equations for Decay
****************************************
The canonical expression of the Bateman equations for a decay chain 
proceeding from a nuclide :math:`A` to a nuclide :math:`Z` at time 
:math:`t` following a specific path is as follows [1]_:

.. math::

    N_C(t) = \frac{N_1(0)}{\lambda_C} \cdot \gamma \cdot \sum_{i=1}^C \lambda_i c_{i} e^{-\lambda_i t}

The symbols in this expression have the following meaning:

================= =========================================================
symbol            meaning
================= =========================================================
:math:`C`         length of the decay chain
:math:`i`         index for ith species, on range [1, C]
:math:`j`         index for jth species, on range [1, C]
:math:`t`         time [seconds]
:math:`N_i(t)`    number density of the ith species at time t
:math:`t_{1/2,i}` half-life of the ith species
:math:`\lambda_i` decay constant of ith species, :math:`ln(2)/t_{1/2,i}`
:math:`\gamma`    The total branch ratio for this chain
================= =========================================================

Additionally, :math:`c_{i}` is defined as:

.. math::

    c_i = \prod_{j=1,i\ne j}^C \frac{\lambda_j}{\lambda_j - \lambda_i}

Furthermore, the total chain branch ratio is defined as the product of the 
branch ratio between any two species [2]_:

.. math:: 

    \gamma = \prod_{i=i}^{C-1} \gamma_{i \to i+1}

Minor modifications are needed for terminal species: the first nuclide of a 
decay chain and the ending stable species. By setting :math:`C=1`, the Bateman
equations can be reduced to simply:

.. math:: 

    N_C(t) = N_1(0) e^{-\lambda_1 t}

For stable species, the appropriate equation is derived by taking the limit
of when the decay constant of the stable nuclide (:math:`\lambda_C`) goes to 
zero.  Also notice that every :math:`c_i` contains exactly one :math:`\lambda_C`
in the numerator which cancels with the :math:`\lambda_C` in the denominator 
in front of the summation:

.. math::

    \lim_{\lambda_C \to 0} N_C(t) = N_1(0)  \gamma \left[e^{-0t} + \sum_{i=1}^{C-1} \lambda_i \left(\frac{1}{0 - \lambda_i} \prod_{j=1,i\ne j}^{C-1} \frac{\lambda_j}{\lambda_j - \lambda_i} \right) e^{-\lambda_i t} \right]

    N_C(t) = N_1(0)  \gamma \left[1.0 - \sum_{i=1}^{C-1} \left(\prod_{j=1,i\ne j}^{C-1} \frac{\lambda_j}{\lambda_j - \lambda_i} \right) e^{-\lambda_i t} \right]


*********************************************
Binary Reformulation of Bateman Equations
*********************************************
There are two main strategies can be used to construct a version of these equations that 
is better suited to computation, if not clarity. 

First, lets aim for minimizing the number of 
operations that must be performed to achieve the same result. This can be done 
by grouping constants together and pre-calculating them. This saves the computer from 
having to perform the same operations at run time.  It is possible to express the 
Bateman equations as a simple sum of exponentials

.. math::

    N_C(t) = N_1(0) \sum_{i=1}^C k_{i} e^{-\lambda_i t}

where the coefficients :math:`k_i` are defined as:

.. math:: 

    k_i = \frac{\gamma}{\lambda_C} \lambda_i c_i

If :math:`k_i` are computed at run time then the this expression results in much more
computational effort that than the original Bateman equations since :math:`\gamma/\lambda_C` 
are brought into the summation. However, when :math:`k_i` are pre-caluclated, 
many floating point operations are saved by avoiding explicitly computing :math:`c_i`.

The second strategy is to note that computers are much better at dealing with powers of
2 then then any other base, even :math:`e`. Thus the ``exp2(x)`` function, or :math:`2^x`,
is faster than the natural exponential function ``exp(x)``, :math:`e^x`.  As proof of this
the following are some simple timing results:

.. code-block:: python 

    In [1]: import numpy as np

    In [2]: r = np.random.random(1000) / np.random.random(1000)

    In [3]: %timeit np.exp(r)
    10000 loops, best of 3: 26.6 µs per loop

    In [4]: %timeit np.exp2(r)
    10000 loops, best of 3: 20.1 µs per loop

This is a savings of about 25%.  Since the core of the Bateman equations are exponentials, 
it is worthwhile to squeeze this algorithm as much as possible.  Luckily, the decay constant
provides an intrinsic mechanism to convert to base-2:

.. math::

    N_C(t) = N_1(0) \sum_{i=1}^C k_{i} e^{-\lambda_i t}

    N_C(t) = N_1(0) \sum_{i=1}^C k_{i} e^{\frac{-\ln(2)\cdot t}{t_{1/2,i}}}

    N_C(t) = N_1(0) \sum_{i=1}^C k_{i} 2^{\frac{-t}{t_{1/2,i}}}

This expression can be further collapsed by defining :math:`a` to be the precomputed 
exponent values:

.. math:: 

    a_i = \frac{-1}{t_{1/2,i}}

Thus, the final form of the binary representation of the Bateman equations are
as follows:

**General Formulation:** 

.. math:: 

    N_C(t) = N_1(0) \sum_{i=1}^C k_{i} 2^{a_i t}

**First Nuclide in Chain:**

.. math:: 

    N_C(t) = N_1(0) \cdot 2^{a_1 t}

**Stable Nuclide:**

.. math:: 

    N_C(t) = N_1(0) \left[1.0 + \sum_{i=1}^{C-1} \lim_{\lambda_C\to 0}(k_{i}) \cdot 2^{a_i t} \right]


With completely precomputed :math:`k`, :math:`a`, and the ``exp2()`` function, this 
formulation minimizes the number of floating point operations while completely 
preserving physics. No assumptions were made aside from the Bateman equations 
themselves in this proof.

Note that it is not possible to reduce the number of operations further.  This 
is because  :math:`k` and :math:`a` cannot be combined without adding further 
operations.

***************************************
Implementation Specific Approximations
***************************************
The above formulation holds generally for any decay chain.  However, certain 
approximations are used in practice to reduce the number of chains and terms 
that are calculated.

1. Decay chains coming from spontaneous fission are not tallied as they 
   lead to an explosion of the total number of chains while contributing to 
   extraordinarily rare branches.
2. Decay alphas are not treated as He-4 production.
3. For chains longer than length 2, any 
   term whose half-life is less than :math:`10^{-8}` of the sum of all 
   half-lives in the chain is dropped. This filtering prevents excessive
   calculation from species which do not significantly contribute to 
   end atom fraction. The threshold :math:`10^{-8}` was chosen as 
   because it is a reasonable naive estimate of floating point error after 
   many operations. If the filtering causes there to be less than 
   two terms in the summation, then the filtering is turned off and all
   terms are computed.
4. To prevent other sources of floating point error, a nuclide is determined 
   to be stable when :math:`\lambda_i < 10^{-16}`, rather than when 
   :math:`\lambda_i = 0.0`.
5. If a chain has any ``NaN`` decay constants, the chain in rejected.
6. If a chain has any infinite :math:`k`, the chain in rejected.

In principle, each of these statements is reasonable. However, they
may preclude desired behavior by users. In such a situation, these 
assumptions should be revisited.

**********************
Additional Information
**********************
For further discussion, please see:

* `the mailing list post <https://groups.google.com/d/topic/pyne-dev/CXmRfBSThDE/discussion>`_, 
* `the pull request, #614 <https://github.com/pyne/pyne/pull/614>`_, and
* `the benchmark study <http://nbviewer.ipython.org/github/pyne/sandbox/blob/master/origen-cmp.ipynb>`_.

Note that the benchmark study shows quite high agreement between this method
and ORIGEN v2.2.

**********
References
**********

.. [1] Jerzy Cetnar, General solution of Bateman equations for nuclear transmutations, 
       Annals of Nuclear Energy, Volume 33, Issue 7, May 2006, Pages 640-645, 
       http://dx.doi.org/10.1016/j.anucene.2006.02.004.

.. [2] Logan J. Harr. Precise Calculation of Complex Radioactive Decay Chains. M.Sc thesis
       Air Force Institute of Technology. 2007. http://www.dtic.mil/dtic/tr/fulltext/u2/a469273.pdf
