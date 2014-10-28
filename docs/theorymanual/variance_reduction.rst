.. _variance_reduction:

=============================
Variance Reduction
=============================

.. currentmodule::variancereduction

*****************************
CADIS Method
*****************************

The Consistent Adjoint-Driven Importance Sampling (CADIS) method [1] is a Monte
Carlo variance reduction method that utilizes a deterministic estimate of the
adjoint flux (the *importance*) to generate a biased source and weight windows
that optimize a Monte Carlo simulation relative to a detector response
function. Once major feature of the scheme is "consistancy", that is, weight
windows are choosen such that particles are always born within them.

In the CADIS method the response is defined as:

.. math::

    R = int dP q(P) \Phi^+(P),

where :math:`q` is a probability distribution function describing the source
strength as a function of the space and energy phase space variable :math:`P`.
and :math:`\Phi^+(P)` the adjoint flux relateive to the detector response
function being optimized. The CADIS method defines the biased source distribution as:

.. math::

    \hat{q}(P) = \frac{q(P) \Phi^+(P)}{R}.

The corresponding weight window lower bounds are defined by:

.. math::
    
    w_l(P) = \frac{R}{Phi^+(P)\frac{\beta + 1}{2}},

where :math:`\beta` is the ratio of the weight window upper bound to the weight
window lower bound (5 in MCNP5).

...................
PyNE implimentation
...................

The PyNE implimentation of the CADIS method is a mesh-based implimentation and
is designed to be used in conjuction with the mesh-based source sampling
capabilities in the source_sampling model. This means that above method, which
is continuous in phase space must be adapted for discretization of space (mesh
volume elements) and energy bins.

Source density (:math:`q'`) is the cannonical way of representing a mesh-based
source within PyNE. This means that the first step of the CADIS method within
PyNE 

.. math::
    
    \hat{q}(P) = 

...........
Assumptions
...........

Any assumptions (explicit or implicit) about this method that would impact use, conclusions, validity, etc.

......................
Additional Information
......................

Details about impact of this theory on method use.

...................
Sample Calculations
...................



.. math::
    
    A\vec{x} = \vec{b}



**********
References
**********

[1] Haghighat, A. and Wagner, J. C., "Monte Carlo Variance Reduction with 
    Deterministic Importance Functions," Progress in Nuclear Energy, Vol. 42,
    No. 1, pp. 25-53, 2003.


All of the justification for us listening to you.







:math:`A\vec{x} = \vec{b}`.

To have blocks of mathematical content, do this

.. math::
    
    A\vec{x} = \vec{b}
