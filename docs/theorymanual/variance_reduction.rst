.. _variance_reduction:

=============================
Variance Reduction
=============================

:Author: Elliott Biondo, Kalin Kiesling

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

************
MAGIC Method
************

The Method of Automatic Generation of Importances by Calculation (MAGIC) is 
a global variance reduction technique in which an initial particle distribution, 
in the form of fluxes, populations, or weights is obtained and then used to 
generate mesh-based weight windows or importances. This method recognizes 
the initial particle distribution will be poor in some highly attenuated regions
but upon iteration of the MAGIC method, the solution will improve. Below are the
steps for the MAGIC method. [2]

1. Run MCNP, in analogue mode, to set up a flux meshtally. Multigroup cross 
   section data and a high energy cut-off, corresponding to a mean-free path no 
   greater than the mesh voxel size, should be used.

2. Process the resulting meshtally data by normalizing the flux to have a value 
   of 0.5 in the source (or highest) region. Use the normalized flux to create
   a new weight window file to be used for MCNP.

3. Modify the original MCNP input to use the generated weight window file and 
   run again.

4. If results are are sufficient, no further iterations are necessary. Else, 
   repeat starting from step 2 until desired flux results are obtained.
   
5. If a high energy cut-off was used, reduce the cut-off energy and repeat 
   iterations until the particle distribution is acceptable. The final iteration
   should be performed with the appropriate energy cut-off and cross section data.

...................
PyNE implementation
...................

The implementation of MAGIC in PyNE uses a PyNE meshtally object, which is the 
result of a meshtal file processed by PyNE's mcnp.Meshtally. Using the results 
of the meshtal file and a specified tolerance for relative error :math:`t` and 
null value :math:`\phi_o`, the flux will be normalized for each energy bin and 
then be used to generate a wwinp file to be used in a subsequent iteration. The 
steps are as follows:

1. Read meshtally and determine maximum flux :math:`\phi_m^k` for each enery bin :math:`k`.

2. Normalize flux :math:`\phi_i^k` in every mesh voxel :math:`i` for each energy 
   bin :math:`k` according to :math:`\phi_m^k` to obtain a new :math:`\phi_i^{'k}`. 
   If the relative error :math:`e_i^k` 
   for voxel :math:`i` in energy bin :math:`k` is larger than the tolerance 
   value :math:`t`, then set flux :math:`\phi_i^{'k}` to the null value :math:`\phi_o` 
   instead.
 
 
 .. math::  
        
        \text{If } e_i^k < t \text{ then, } \phi_i^{'k} = \frac{\phi_i^{k}}{2 \, \phi_m^k}
        
        \text{If } e_i^k > t \text{ then, } \phi_i^{'k} = \phi_o
        
        
3. Use new flux values to create a weight window tag on the provide meshtally 
   and use PyNE's Wwinp class to create a weight window mesh.
   
...................
Sample Calculations
...................

In this section, the expected results of the test_variancereduction.py unit test
"test_magic_multi_bins" are shown. In this test, a 3D 2x2 mesh is given. Each
voxel contains flux data corresponding to 2 energy bins. The mesh is described by
the following flux and relative error data.

+--------------------------+----------------------+--------------------------+-----------------------+
| :math:`\phi_1^{1} = 1.2` | :math:`e_1^1 = 0.11` | :math:`\phi_1^{2} = 3.3` | :math:`e_1^2 = 0.013` |
+--------------------------+----------------------+--------------------------+-----------------------+
| :math:`\phi_2^{1} = 1.6` | :math:`e_2^1 = 0.14` | :math:`\phi_2^{2} = 1.7` | :math:`e_2^2 = 0.19`  |
+--------------------------+----------------------+--------------------------+-----------------------+
| :math:`\phi_3^{1} = 1.5` | :math:`e_3^1 = 0.02` | :math:`\phi_3^{2} = 1.4` | :math:`e_3^2 = 0.16`  |
+--------------------------+----------------------+--------------------------+-----------------------+
| :math:`\phi_4^{1} = 2.6` | :math:`e_4^1 = 0.04` | :math:`\phi_4^{2} = 1.0` | :math:`e_4^2 = 0.09`  |
+--------------------------+----------------------+--------------------------+-----------------------+

First, the maximum flux for each energy bin is found. In this case the maximum 
for energy bin :math:`k = 1` occurs in voxel 4 :math:`\phi_4^{1} = 2.6` and in 
voxel 1 :math:`\phi_1^{2} = 3.3` for energy bin :math:`k = 2`. In the first 
energy bin, the flux values are normalized by :math:`\phi_m^1 = 2.6` and in the 
second :math:`\phi_m^2 = 3.3`. If the error tolerance is set :math:`t = 0.15` and the
null value set to :math:`\phi_o = 0.001`, then voxels 
2 and 3 in the second energy bin have errors larger than the tolerance and are 
therefore set to the null value while everything else is normalized. The following
is the result.

+------------------------------------------------------------------------------------+------------------------------------------------------------------------------------+
| :math:`\phi_1^{'1} = \frac{\phi_1^1}{2 \, \phi_m^1} = \frac{1.2}{2*2.6} = 0.23077` | :math:`\phi_1^{'2} = \frac{\phi_1^2}{2 \, \phi_m^2} = \frac{3.3}{2*3.3} = 0.5`     |
+------------------------------------------------------------------------------------+------------------------------------------------------------------------------------+
| :math:`\phi_2^{'1} = \frac{\phi_2^1}{2 \, \phi_m^1} = \frac{1.6}{2*2.6} = 0.30769` | :math:`\phi_2^{'2} = \phi_o = 0.001`                                               |
+------------------------------------------------------------------------------------+------------------------------------------------------------------------------------+
| :math:`\phi_3^{'1} = \frac{\phi_3^1}{2 \, \phi_m^1} = \frac{1.5}{2*2.6} = 0.28846` | :math:`\phi_3^{'2} = \phi_o = 0.001`                                               |
+------------------------------------------------------------------------------------+------------------------------------------------------------------------------------+
| :math:`\phi_4^{'1} = \frac{\phi_4^1}{2 \, \phi_m^1} = \frac{2.6}{2*2.6} = 0.5`     | :math:`\phi_4^{'2} = \frac{\phi_4^2}{2 \, \phi_m^2} = \frac{1.0}{2*3.3} = 0.12122` |
+------------------------------------------------------------------------------------+------------------------------------------------------------------------------------+

The values :math:`\phi_i^{'k}` are then set as the new weight window values.


**********
References
**********

[1] Haghighat, A. and Wagner, J. C., "Monte Carlo Variance Reduction with 
    Deterministic Importance Functions," Progress in Nuclear Energy, Vol. 42,
    No. 1, pp. 25-53, 2003.
    
[2] Davis, A. and Turner, A., "Comparison of global variance reduction 
    techniques for Monte Carlo radiation transport simulations of ITER," Fusion 
    Engineering and Design, Vol. 86, Issues 9-11, pp. 2698-2700, 2011.

