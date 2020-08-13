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
zero. Also notice that every :math:`c_i` contains exactly one :math:`\lambda_C`
in the numerator which cancels with the :math:`\lambda_C` in the denominator
in front of the summation:

.. math::

    \lim_{\lambda_C \to 0} N_C(t) = N_1(0)  \gamma \left[e^{-0t} + \sum_{i=1}^{C-1} \lambda_i \left(\frac{1}{0 - \lambda_i} \prod_{j=1,i\ne j}^{C-1} \frac{\lambda_j}{\lambda_j - \lambda_i} \right) e^{-\lambda_i t} \right]

    N_C(t) = N_1(0)  \gamma \left[1.0 - \sum_{i=1}^{C-1} \left(\prod_{j=1,i\ne j}^{C-1} \frac{\lambda_j}{\lambda_j - \lambda_i} \right) e^{-\lambda_i t} \right]

Now, certain chains have intermeadiate nuclides that are *almost* stable. For example, decaying
from Es-254 to Po-210 goes through U-238, which is very close to stable relative to all of the
other nuclides in the chain. This can trigger floating point precision issues, where certain
terms will underflow or overflow or generate NaNs. Obviously this is a situation to be avoided,
if at all possible. To handle this sitiuation, let's call :math:`p` the index of the nuclide
that is almost stable. We can then note that the Bateman equations can be reduced by the
observation that :math:`\lambda_p \ll \lambda_{i\ne p}` after we separate out the p-term
from the summation:

.. math::

   \frac{N_C(t)}{N_1(0)} = \frac{\gamma}{\lambda_C}\sum_{i\ne p}^C \left[\lambda_i \frac{\lambda_p}{\lambda_p - \lambda_i}
                                                        \left(\prod_{j\ne i,p}^C \frac{\lambda_j}{\lambda_j - \lambda_i}\right)
                                                        e^{-\lambda_i t}\right]
                           + \frac{\gamma}{\lambda_C} \lambda_p \left(\prod_{j\ne p}^C \frac{\lambda_j}{\lambda_j - \lambda_p} \right) e^{-\lambda_p t}

   \frac{N_C(t)}{N_1(0)} = \frac{\gamma}{\lambda_C}\sum_{i\ne p}^C \left[\lambda_i \frac{\lambda_p}{\lambda_p - \lambda_i}
                                                        \left(\prod_{j\ne i,p}^C \frac{\lambda_j}{\lambda_j - \lambda_i}\right)
                                                        e^{-\lambda_i t}\right]
                           + \frac{\gamma}{\lambda_C} \lambda_p \left(\prod_{j\ne p}^C \frac{\lambda_j}{\lambda_j - \lambda_p}\right) e^{-\lambda_p t}

   \frac{N_C(t)}{N_1(0)} = \frac{\gamma}{\lambda_C}\sum_{i\ne p}^C \left[\lambda_i \frac{\lambda_p}{- \lambda_i}
                                                        \left(\prod_{j\ne i,p}^C \frac{\lambda_j}{\lambda_j - \lambda_i}\right)
                                                        e^{-\lambda_i t}\right]
                           + \frac{\gamma}{\lambda_C} \lambda_p \left(\prod_{j\ne p}^C \frac{\lambda_j}{\lambda_j}\right) e^{-\lambda_p t}

   \frac{N_C(t)}{N_1(0)} = \frac{-\gamma\lambda_p}{\lambda_C}\sum_{i\ne p}^C \left[
                                                        \left(\prod_{j\ne i,p}^C \frac{\lambda_j}{\lambda_j - \lambda_i}\right)
                                                        e^{-\lambda_i t}\right]
                           + \frac{\gamma\lambda_p}{\lambda_C} e^{-\lambda_p t}

The above expression for intermediate nuclides that are almost stable is valid when the last
nuclide in the chain is unstable. When the last nuclide is stable, both the pth
(almost stable nuclide) and the Cth (last and stable nuclide) must be removed can be split off from
the summation and handled separately. As previously, then take :math:`\lambda_C \to 0` and :math:`\lambda_p \ll \lambda_{i\ne p,C}`.

.. math::

   \frac{N_C(t)}{N_1(0)} = \frac{\gamma}{\lambda_C}\sum_{i\ne p}^{C-1} \left[\lambda_i \frac{\lambda_C}{\lambda_C - \lambda_i} \frac{\lambda_p}{\lambda_p - \lambda_i}
                                                        \left(\prod_{j\ne i,p}^{C-1} \frac{\lambda_j}{\lambda_j - \lambda_i}\right)
                                                        e^{-\lambda_i t}\right]
                           + \frac{\gamma}{\lambda_C} \lambda_p \frac{\lambda_C}{\lambda_C - \lambda_p} \left(\prod_{j\ne p}^{C-1} \frac{\lambda_j}{\lambda_j - \lambda_p} \right) e^{-\lambda_p t}
                           + \frac{\gamma}{\lambda_C} \lambda_C \frac{\lambda_p}{\lambda_p - \lambda_C} \left(\prod_{j\ne p}^{C-1} \frac{\lambda_j}{\lambda_j - \lambda_C} \right) e^{-\lambda_C t}

   \frac{N_C(t)}{N_1(0)} = \gamma\sum_{i\ne p}^{C-1} \left[\frac{\lambda_i \lambda_p}{(\lambda_C - \lambda_i)(\lambda_p - \lambda_i)}
                                                        \left(\prod_{j\ne i,p}^{C-1} \frac{\lambda_j}{\lambda_j - \lambda_i}\right)
                                                        e^{-\lambda_i t}\right]
                           + \frac{\gamma\lambda_p}{\lambda_C - \lambda_p} \left(\prod_{j\ne p}^{C-1} \frac{\lambda_j}{\lambda_j} \right) e^{-\lambda_p t}
                           + \frac{\gamma\lambda_p}{\lambda_p - \lambda_C} \left(\prod_{j\ne p}^{C-1} \frac{\lambda_j}{\lambda_j} \right) e^{-\lambda_C t}

   \frac{N_C(t)}{N_1(0)} = -\gamma\sum_{i\ne p}^{C-1} \left[\left(\prod_{j\ne i}^{C-1} \frac{\lambda_j}{\lambda_j - \lambda_i}\right) e^{-\lambda_i t}\right]
                           + \frac{\gamma\lambda_p}{\lambda_C - \lambda_p} e^{-\lambda_p t}
                           + \frac{\gamma\lambda_p}{\lambda_p - \lambda_C} e^{-\lambda_C t}

   \frac{N_C(t)}{N_1(0)} = -\gamma\sum_{i\ne p}^{C-1} \left[\left(\prod_{j\ne i}^{C-1} \frac{\lambda_j}{\lambda_j - \lambda_i}\right) e^{-\lambda_i t}\right]
                           + \frac{\gamma\lambda_p}{\lambda_C - \lambda_p} \left(e^{-\lambda_p t} - e^{-\lambda_C t}\right)

   \frac{N_C(t)}{N_1(0)} = -\gamma\sum_{i\ne p}^{C-1} \left[\left(\prod_{j\ne i}^{C-1} \frac{\lambda_j}{\lambda_j - \lambda_i}\right) e^{-\lambda_i t}\right]
                           -\gamma e^{-\lambda_p t} + \gamma


Lastly, we must handle the degenerate case where two nuclides in a chain  have the same exact half-lives.
This unfortunate situation arrises out of the fundemental nuclear data. Let's call these the pth and qth
species. To prevent underflow, overflow, and NaNs, we must separate these nuclides out of the summation
and then take the limit as :math:`\lambda_q \to \lambda_p`.

.. math::

   \frac{N_C(t)}{N_1(0)} = \frac{\gamma}{\lambda_C}\sum_{i\ne p,q}^{C} \left[\lambda_i \left(\prod_{j\ne i}^{C} \frac{\lambda_j}{\lambda_j - \lambda_i}\right) e^{-\lambda_i t}\right]
                           + \frac{\gamma}{\lambda_C} \lambda_p \frac{\lambda_q}{\lambda_q - \lambda_p} \left(\prod_{j\ne p,q}^{C} \frac{\lambda_j}{\lambda_j - \lambda_p} \right) e^{-\lambda_p t}
                           + \frac{\gamma}{\lambda_C} \lambda_q \frac{\lambda_p}{\lambda_p - \lambda_q} \left(\prod_{j\ne p,q}^{C} \frac{\lambda_j}{\lambda_j - \lambda_q} \right) e^{-\lambda_q t}

   \frac{N_C(t)}{N_1(0)} = \frac{\gamma}{\lambda_C}\sum_{i\ne p,q}^{C} \left[\lambda_i \left(\prod_{j\ne i}^{C} \frac{\lambda_j}{\lambda_j - \lambda_i}\right) e^{-\lambda_i t}\right]
                           + \frac{\gamma\lambda_p^2}{\lambda_C} \left(\prod_{j\ne p,q}^{C} \frac{\lambda_j}{\lambda_j - \lambda_p} \right)
                             \lim_{\lambda_q\to\lambda_p}\frac{e^{-\lambda_p t} - e^{-\lambda_q t}}{\lambda_q - \lambda_p}

   \frac{N_C(t)}{N_1(0)} = \frac{\gamma}{\lambda_C}\sum_{i\ne p,q}^{C} \left[\lambda_i \left(\prod_{j\ne i}^{C} \frac{\lambda_j}{\lambda_j - \lambda_i}\right) e^{-\lambda_i t}\right]
                           + \frac{\gamma\lambda_p^2}{\lambda_C} \left(\prod_{j\ne p,q}^{C} \frac{\lambda_j}{\lambda_j - \lambda_p} \right) t e^{-\lambda_p t}


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


**Single Nuclide in Chain:**

.. math::

    k_i = 1

**Last Nuclide Unstable:**

.. math::

    k_i = \frac{\gamma}{\lambda_C} \lambda_i \prod_{j\ne i}^C \frac{\lambda_j}{\lambda_j - \lambda_i}


**Last Nuclide Stable:**

.. math::

    k_{i\ne C} = -\gamma \prod_{j=1,i\ne j}^{C-1} \frac{\lambda_j}{\lambda_j - \lambda_i}

    k_C = \gamma


**Last Nuclide Unstable and pth Almost Stable:**

.. math::

    k_{i\ne p} = -\frac{\gamma\lambda_p}{\lambda_C} \prod_{j\ne i,p}^C \frac{\lambda_j}{\lambda_j - \lambda_i}

    k_p = \frac{\gamma\lambda_p}{\lambda_C}


**Last Nuclide Stable and pth Almost Stable:**

.. math::

    k_{i\ne p,C} = -\gamma \prod_{j\ne i}^{C-1} \frac{\lambda_j}{\lambda_j - \lambda_i}

    k_p = -\gamma

    k_C = \gamma


**Half-life Degeneracy Between pth and qth:**

.. math::

    k_i = \frac{\gamma}{\lambda_C} \lambda_i \prod_{j\ne i}^C \frac{\lambda_j}{\lambda_j - \lambda_i}

    k_p = \frac{\gamma\lambda_p^2}{\lambda_C} t \prod_{j\ne p,q}^C \frac{\lambda_j}{\lambda_j - \lambda_p}

    k_q = 0



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

where the :math:`k_i` are as listed above.  However, for practical purposes, it is better to
compute the :math:`k_i` from half-lives rather than decay constants.  This is because they
provide less floating point error, fewer oppurtunities to underflow or overflow to NaN or infinity,
and a better mechanism for detecting stability. Thus, alternatively, the :math:`k_i` are computed
as:

**Single Nuclide in Chain:**

.. math::

    k_i = 1

**Last Nuclide Unstable:**

.. math::

    k_i = \gamma t_{1/2,i}^{C-2} t_{1/2,C} \prod_{j\ne i}^{C} \frac{1}{t_{1/2,i} - t_{1/2,j}}


**Last Nuclide Stable:**

.. math::

    k_i = -\gamma t_{1/2,i}^{C-2} \prod_{j\ne i}^{C-1} \frac{1}{t_{1/2,i} - t_{1/2,j}}

    k_C = \gamma


**Last Nuclide Unstable and pth Almost Stable:**

.. math::

    k_{i\ne p} = -\frac{\gamma t_{1/2,C}}{t_{1/2,p}} t_{1/2,i}^{C-2} \prod_{j\ne i,p}^C \frac{1}{t_{1/2,i} - t_{1/2,j}}

    k_p = \frac{\gamma t_{1/2,C}}{t_{1/2,p}}


**Last Nuclide Stable and pth Almost Stable:**

.. math::

    k_{i\ne p,C} = -\gamma t_{1/2,i}^{C-2} \prod_{j\ne i}^{C-1} \frac{1}{t_{1/2,i} - t_{1/2,j}}

    k_p = -\gamma

    k_C = \gamma


**Half-life Degeneracy Between pth and qth:**

.. math::

    k_i = \gamma t_{1/2,i}^{C-2} t_{1/2,C} \prod_{j\ne i}^{C} \frac{1}{t_{1/2,i} - t_{1/2,j}}

    k_p = \gamma\ln(2) t_{1/2,p}^{C-4} t_{1/2,C}  t \prod_{j\ne p,q}^C \frac{1}{t_{1/2,p} - t_{1/2,j}}

    k_q = 0

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

1. Decay chains coming from spontaneous fission are only optionally tallied as they
   lead to an explosion of the total number of chains while contributing to
   extraordinarily rare branches.
2. Decay alphas are not treated as He-4 production.
3. The :math:`k_i` and :math:`a_i` are filtered to reject terms where
   :math:`|k_i| / \max(|k_i|) < 10^{-16}`.
   This filtering prevents excessive
   calculation from species which do not significantly contribute to
   end atom fraction. The threshold :math:`10^{-16}` was chosen as
   because it is a reasonable naive estimate of floating point error after
   many operations. Note that we may filter only on the :math:`k_i` because
   :math:`2^{a_i t} \le 1`.  That is, the exponentional component can only
   reduce the magnitude of a term, not increase it.

In principle, each of these statements is reasonable. However, they
may preclude desired behavior by users. In such a situation, these
assumptions should be revisited.

**********************
Additional Information
**********************
For further discussion, please see:

* `the mailing list post <https://groups.google.com/d/topic/pyne-dev/CXmRfBSThDE/discussion>`_,
* `the pull request, #614 <https://github.com/pyne/pyne/pull/614>`_, and
* `the benchmark study <https://nbviewer.jupyter.org/github/pyne/sandbox/blob/master/origen-cmp.ipynb>`_.

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
