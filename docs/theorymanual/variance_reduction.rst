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

    R = \int dP \, q(P) \, \Psi^+(P),

where :math:`q` is a probability distribution function describing the source
strength as a function of the phase space variable :math:`P`.  and
:math:`\Phi^+(P)` the adjoint flux relateive to the detector response function
being optimized. The CADIS method defines the biased source distribution as:

.. math::

    \hat{q}(P) = \frac{q(P) \, \Psi^+(P)}{R}.

The corresponding weight window lower bounds are defined by:

.. math::
    
    ww(P) = \frac{R}{\Psi^+(P) \, \frac{\beta + 1}{2}},

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

Source density (:math:`q'`, units: :math:`time^{-1}length{-3}`) is the
cannonical way of representing a mesh-based source within PyNE.This means that
the first step of the CADIS method within PyNE is create a :math:`q` PDF. The
total source strength :math:`q_tot` is first found, integrating the source density
over space (:math:`i`) and energy (:math:`j`):

.. math::
    
    q_{tot} = \sum_{i \in I} \sum_{j \in J} V_i \, q'_{i, j},

The :math:`q` PDF can then be defined by:

.. math::

    q_{i,j } = \frac{V_i \, q'_{i, j}}{q_tot} for i \in I, j \in J

The response can then be calculated by integrating the product of :math:`q` and the adjoint flux over all phase space:

.. math::
    
    R = \sum_{i \in I} \sum_{j \in J} \Psi_{i, j}^{+} \, \frac{V_i \, q'_{i, j}}{q_{tot}}

The weight window lower bound is then:

.. math::

    ww_{i, j} = \frac{R}{\Psi_{i, j}^{+} \, \frac{\beta + 1}{2}}.

These values tagged to the weight window output mesh and can be printed out as
an MCNP5 WWINP file.  The biased source strength is:

.. math::

   \hat{q}_{i, j} = \frac{\Psi_{i, j}^{+} \, q'_{i, j} \, V_i}{R \, q_{tot}}

However, the biased source strength is not the quantity on interest, because
the source_sampling module is expecting biased source densities. The biased
source densities that are tagged to the output mesh are:

.. math::

   \hat{q}'_{i, j} = \frac{\Psi_{i, j}^{+} \, q'_{i, j}}{R \, q_{tot}}.

...........
Assumptions
...........

The source density mesh and adjoint flux mesh must have the spacial bounds.

...................
Sample Calculations
...................

In this section the expected results for the the test_variancereduction.py unit
test "test_cadis_single_e". Consider 

+--------------------+-----------+
|                   |  aasdfasdf|
+--------------------+-----------+
|                  |   a       |
+--------------------+-----------+





.. math::
    
    q_{tot} = \sum_{i \in I} \sum_{j \in J} V_i \, q'_{i, j}
            = 8*2.9 + 2*2.26 +

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
