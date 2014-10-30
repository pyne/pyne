.. _variance_reduction:

=============================
Variance Reduction
=============================

.. currentmodule::variancereduction

************
CADIS Method
************

The Consistent Adjoint-Driven Importance Sampling (CADIS) method [1] is a Monte
Carlo variance reduction method that utilizes a deterministic estimate of the
adjoint flux (the *importance*) to generate a biased source and weight windows
that optimize a Monte Carlo simulation relative to a detector response
function. One major feature of the scheme is "consistency", that is, weight
windows are chosen such that particles are always born within them.

In the CADIS method the response is defined as:

.. math::

    R = \int \, q(P) \, \Psi^+(P) \, dP,

where :math:`q` is the probability distribution function describing the source
strength as a function of the phase space variable :math:`P` (which may
represent any combination of space, energy, and direction). :math:`\Phi^+(P)`
is the adjoint flux relative to the detector response function being
optimized. The CADIS method defines the biased source distribution as:

.. math::

    \hat{q}(P) = \frac{q(P) \, \Psi^+(P)}{R}.

The corresponding weight window lower bounds are defined by:

.. math::
    
    ww(P) = \frac{R}{\Psi^+(P) \, \frac{\beta + 1}{2}},

where :math:`\beta` is the ratio of the weight window upper bound to the weight
window lower bound (default of 5 in MCNP5).

...................
PyNE implementation
...................

The PyNE implementation of the CADIS method is a mesh-based implementation and
is designed to be used in conjunction with the mesh-based source sampling
capabilities in the source_sampling module. This means that the above method, which
is continuous in phase space must be adapted for discretization of space (mesh
volume elements) and energy (in energy bins).

Source density (:math:`q'`, units: :math:`time^{-1}length^{-3}`) is the
canonical quantity for representing a mesh-based source within PyNE. This
means that the first step of the CADIS method within PyNE is to create a :math:`q`
PDF from a source density mesh. The total source strength :math:`q_tot` is
first found by integrating the source density over space (:math:`i`) and energy
(:math:`j`):

.. math::
    
    q_{tot} = \sum_{i \in I} \, \sum_{j \in J} V_i \, q'_{i, j},

The :math:`q` PDF can then be defined by:

.. math::

    q_{i,j } = \frac{V_i \, q'_{i, j}}{q_tot} \, for i \in I, j \in J

The response can then be calculated by integrating the product of :math:`q` and the adjoint flux over all phase space:

.. math::
    
    R = \sum_{i \in I} \, \sum_{j \in J} \Psi_{i, j}^{+} \, \frac{V_i \, q'_{i, j}}{q_{tot}}

The weight window lower bound is then:

.. math::

    ww_{i, j} = \frac{R}{\Psi_{i, j}^{+} \, \frac{\beta + 1}{2}}.

These values are tagged to the weight window output mesh and can be printed out as
an MCNP5 WWINP file. In the event that the adjoint flux is 0 for some
:math:`(i, j)`, the :math:`ww_{i, j}` value is replaced with 0. MCNP5 will not
play the weight window game when a particle enters a region of phase space
where the weight window lower bound is 0.

The biased source strength is:

.. math::

   \hat{q}_{i, j} = \frac{\Psi_{i, j}^{+} \, q'_{i, j} \, V_i}{R \, q_{tot}}

However, the biased source strength is not the quantity of interest, because
the source_sampling module is expecting biased source densities. The biased
source densities that are tagged to the output mesh are:

.. math::

   \hat{q}'_{i, j} = \frac{\Psi_{i, j}^{+} \, q'_{i, j}}{R \, q_{tot}}.

...........
Assumptions
...........

The source density mesh and adjoint flux mesh must have the spatial bounds.

...................
Sample Calculations
...................

In this section the expected results for the the test_variancereduction.py unit
test "test_cadis_multiple_e" are calculated. Consider a 2D mesh with the
following properties.

+---------------------------+---------------------------+
| :math:`q' = [2.6, 2.5]`   | :math:`q' = [2.9, 0]`     |
|                           |                           |
| :math:`\Phi = [1.3, 1.4]` | :math:`\Phi = [1.7, 1.9]` |
|                           |                           |
| :math:`V = 2`             | :math:`V = 2`             |
+---------------------------+---------------------------+
| :math:`q' = [2.9, 2.8]`   | :math:`q' = [2.4, 2.2]`   |
|                           |                           |
| :math:`\Phi = [1.1, 1.2]` | :math:`\Phi = [0, 1.6]`   |
|                           |                           |
| :math:`V = 8`             | :math:`V = 8`             |
+---------------------------+---------------------------+

Here, the vector quantities represent values at two energy groups. First
calculate :math:`q_{tot}` and :math:`R`:

.. math::
    
    q_{tot} & = \sum_{i \in I} \, \sum_{j \in J} V_i \, q'_{i, j} \\
            & = 8 \cdot 2.9 + 8 \cdot 2.8 + 2 \cdot 2.6 + 2 \cdot 2.5 \\
            &   + 8 \cdot 2.4 + 8 \cdot 2.2 + 2 \cdot 2.9 + 2 \cdot 0  \\
            & = 98.4

.. math::
    
    R & = \sum_{i \in I} \, \sum_{j \in J} \Psi_{i, j}^{+} \, \frac{V_i \, q'_{i, j}}{q_{tot}} \\
      & = \frac{1}{98.4} (
               1.1 \cdot 8 \cdot 2.9 + 1.2 \cdot 8 \cdot 2.8 
              + 1.3 \cdot 2 \cdot 2.6 + 1.4 \cdot 2 \cdot 2.5 \\
              & \qquad \quad + 0 \cdot 8 \cdot 2.4 + 1.6 \cdot 8 \cdot 2.2 
              + 1.7 \cdot 2 \cdot 2.9 + 1.9 \cdot 2 \cdot 0 ) \\
      & = 1.0587398374


The expected results are:

.. math::

   \hat{q}' &= \frac{\Psi_{i, j}^{+} \, q'_{i, j}}{R \, q_{tot}} \\
            &= \frac{1}{98.4 \cdot 1.0587398374}
             [1.1 \cdot 2.9, 1.2 \cdot 2.8, 1.3 \cdot 2.6, 1.4 \cdot 2.5, \\
            & \qquad \qquad \qquad \qquad \qquad 0 \cdot 2.4, 1.6 \cdot 2.2, 1.7 \cdot 2.9, 1.9 \cdot 0] \\
            &= [0.0306200806, 0.0322518718, 0.0324438472, 0.0335956998, \\
            & \qquad 0.0, 0.0337876752, 0.0473219428, 0.0]

.. math::

    ww &= \frac{R}{\Psi_{i, j}^{+} \, \frac{\beta + 1}{2}} \\
       &= [ \frac{1.0587398374}{1.1 \cdot{3}}, \frac{1.0587398374}{1.2 \cdot{3}},
                 \frac{1.0587398374}{1.3 \cdot{3}}, \frac{1.0587398374}{1.4 \cdot{3}}, \\
       & \qquad  \frac{1.0587398374}{0 \cdot{3}}, \frac{1.0587398374}{1.6 \cdot{3}},
                 \frac{1.0587398374}{1.7 \cdot{3}}, \frac{1.0587398374}{1.9 \cdot{3}} ] \\
       &= [0.3208302538, 0.2940943993, 0.2714717532, 0.2520809137, \\
       & \qquad  0.0, 0.2205707995, 0.2075960465, 0.1857438311]

Notice that the value in the :math:`ww` vector that is a division by 0 has been replaced with 0.

**********
References
**********

[1] Haghighat, A. and Wagner, J. C., "Monte Carlo Variance Reduction with 
    Deterministic Importance Functions," Progress in Nuclear Energy, Vol. 42,
    No. 1, pp. 25-53, 2003.

