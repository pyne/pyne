.. _theorymanual_template:

===============================
Spatial Solvers
===============================

.. currentmodule:: pyne.spatialsolver

This page presents a brief summary of the underlying methods used in the 
:py:mod:`spatialsolver` module. For full details, consult "Development of a 
Quantitative Decision Metric for Selecting the Most Suitable Discretization 
Method for S_N Transport Problems", a disseration by Sebastian Schunert (2013).
Chapter 2 has a detailed description of the methods implemented in PyNE. All 
code in this module is derived from Dr. Schunert's PhD work.

The spatial solver module contains a suite of methods that solve the S_N 
formulation of the Boltzmann transport equation. What differentiates the 
methods are their treatment of the spatial aspect of the solution process. 
"The linear Boltzmann transport equation describes the evolution of the flux of
neutral particles, i.e. neutrons or photons, in a host medium. It can be 
obtained from the general Boltzmann transport equation by neglecting particle-
particle interactions, the dependence of the material properties of the host 
medium on the particle flux, and assuming that no electric force field is 
present" (p. 5). The S_N equations define the transport equation along specific
segments of the angular variable, omega, such that omega(n)=(mu,eta,xi)^T, with
n = 1,..,N. 

All of the methods discussed here are discontinous finite element methods (DFEM),
which can also be called nodal methods. 
"Common to all FEM schemes is that the solution of the PDE is approximated by a linear
combination of functions belonging to some finite dimensional trial function space. The un-
knowns of the FEM computation are the coefficients of the linear combination of trial functions,
also referred to as expansion coefficients." (P. 14)
Nodal methods share these properties:
* All function spaces are defined local to a mesh cell.
* Coupling between cells occurs only through their faces.
* Coupling between cells is only imposed in an integral sense.
* Increasing the order of the methods is achieved by increasing the local order of expansion. 

These methods are implemented and described below:
* :ref:`hodd`
* :ref:`ahotn`
* :ref:`ahotn-ll-ln`
* :ref:`ahotn-nefd`
* :ref:`dgfem-lagrange`
* :ref:`dgfem-ld`
* :ref:`dgfem-ll`
* :ref:`sct-step`

.. _hodd:
*****************************
HODD
*****************************

Higher Order Diamond Difference (HODD)

.. _ahotn:
*****************************
General AHOTN
*****************************

The AHOTN methods are a type of Transverse Moment Type Methods (TMB).  TMB methods are a nodal method 
derived for an arbitrary expansion with order lambda.  They constitute the “per mesh-cell system of 
equations from the spatial Legendre moments of the transport Eq. 2.23 augmented by closure/auxiliary 
relations derived via the transverse moment procedure, followed by an approximate direction-by-direction
analytical solution of the resulting 1D transport equation.” (1, p.33).  

TMB methods are particularly good at resulting in accurate solutions on “course spatial meshes, thus 
potentially increasing solution efficiency” (1,33).  This is needed, because most other methods developed
at the time such as DD failed at solving course spatial meshes, resulting in either extremely inaccurate
or negative solutions.  

The AHOTN (Arbitrarily High Order Transport Method of the Nodal type) methods are unique from most other 
TMB methods because they are developed to have a very compact weighted diamond difference (WDD) representation 
of the per-cell set of equations.  A full derivation of the AHOTN solutions to the SN equations can be found 
from pages 34 through 40 of [1].  

There are three solvers that use the AHOTN method to solve the neutron transport equation packaged with pyne.  
They are the following: Linear-Nodal (LN), Linear-Linear (LL), and the NEFD algorithm.  They are discussed 
further below.

.. _ahotn-ll-ln:
*****************************
AHOTN-LL/-LN
*****************************

Both the AHOTN Linear Linear and Linear Nodal methods use linear TMB approximations to solve the SN equations.  
The linear nature of this approximation allows for “achieving high accuracy at a reasonably short execution time, 
thus resulting in improved computational efficiency” [1,p.37].  The LL and LN methods also both “utilize moments 
of the balance equations, Eq. 2.23, satisfying mx + my + mz <= 1 augmented by three WDD equations per dimensions.  
Thus, the full set of LN and LL equations comprises four balance relations and twelve WDD equations.  Within this 
subsection the WDD equations for the LL and LN method are derived in three-dimensional Cartesian geometry” [1,p.38].

“The difference between the LN and LL method is that the LL method retains the bilinear leakage component while 
the LN neglects it. From an algorithmic (i.e. solution of the local equations within the mesh sweep) point of 
view the LN provides the least coupling among the set of equations while the LL has stronger coupling in the WDD 
relations than both LN and AHOTN-1 and AHOTN-1 has stronger coupling than LL and LN in the nodal balance 
equations” [1,p.40].

.. _ahotn-nefd:
*****************************
AHOTN-NEFD
*****************************

?????? Method:

"As the AHOTN method can be conveniently cast into a WDD form with all the AHOTN
specifics lumped into the spatial weights, a standard WDD solver can be used to solve the
per-cell AHOTN system of equations. Typically, the WDD relations Eq. 2.48 are solved for
the outflow face moments and substituted into the nodal balance relations Eq. 2.23 which are
then solved for the (Λ + 1) 3 unknown nodal flux moments (NEFD algorithm),"(1, p.37)

.. _dgfem:
*****************************
DGFEM General
*****************************

The Discontinuous Galerkin Finite Element Method (DGFEM) solver uses 
identical test and trial function spaces that are typically substituted into the weak 
form of the transport equation and tested against all members of the test space 
to obtain a per-cell system of equations.
Two families of DGFEM function spaces are mostly used in discretizing the spatial variable in 
the SN approximation of the transport equation: (1) the complete family and (2) the Lagrange family." (1,p27) 

"Assume that we formulate our function spaces such that we solve for point values of the flux, i.e. we use 
Lagrange polynomials as basis functions. Then, in two-dimensional triangular geometry and three-dimensional 
tetrahedral geometry the complete basis would require one flux value per corner point. The Lagrange basis 
would introduce more degrees of freedom that are not associated with the flux values in the corner points. 
In two-dimensional and three-dimensional Cartesian geometry the Lagrange family would result in one flux 
value per corner point. The complete basis would result in less degrees of freedom. For Λ = 1 for example, 
the Lagrange function space seems for more natural for Cartesian meshes, while the complete family appears 
to be a more natural choice for triangles/tetrahedra." (1,p.28)

When comparing these three included DGFEM solvers with a fixed expansion order, the Lagrange family
is more accurate, while the complete family executes faster.

.. _dgfem-lagrange:
*****************************
DGFEM-Lagrange
*****************************

DGFEM-Lagrange is part of the Lagrange family of solvers of the DGFEM type.  All Lagrange solvers are implemented 
on the Lagrange function space, rather than the complete function space.  For an explanation of both spaces see 
4.4.1 and 4.4.2 from [1].  "The DGFEM method for discretizing the SN equations was first suggested by Reed and 
Hill [36] for two-dimensional triangular cells using a basis of Lagrange polynomials: Each Lagrange basis function 
is associated with a support point at which its value is unity while it assumes a zero value at all other support 
points. The unknowns in Reed’s methods are then the flux values at the support points and the method’s order is 
related to the number of support points within a single cell." [1,p.26]

.. _dgfem-ld:
*****************************
DGFEM-LD
*****************************

Linearly Discontinuous Method:  "The linear discontinuous DGFEM method (LD) is the special case of the complete 
DGFEM method of order Λ = 1. It is special in that the local matrix T is of size 4 × 4 and therefore its inverse 
can be precomputed thus saving execution time. Following [24] we decided to implement the LD method distinctly 
from the arbitrary order complete DGFEM kernel in order to create a highly optimized method." (1, p.92)

.. _dgfem-ll:
********************************
DGFEM-LL
********************************

Part of the complete family?  Or is this the complete dgfem solver?

.. _sct-step:
*****************************
Sct-Step
*****************************

One of the problems with most spatial solvers is the inconsistent and sharp boundary conditions.  A new method 
called SCT-STEP was developed to try to fix this, by using a “Step approximation in all cells that are intersected 
by lines and planes of non-smoothness.” [1,p.4] The SCT-STEP method is an implementation of Duo’s SCT Step 
algorithm in three dimensional cartesian geometry.

“It is important to point out that this extension is highly non-trivial because of the tremendous increase in 
complexity of the tracking and cell-splitting algorithms involved. The new method is labeled SCT-Step method." (1,p.4)


Background on Duo’s two dimensional SCT algorithm:
"In two-dimensional geometry Duo[17] suggested tracking of the singular characteristic line through the mesh 
and applying a sub-cell approach in intersected cells to keep segments in these cells isolated from each 
other. For the solution of the subcell equations, Duo used the Step Characteristic method applied to each 
of the segments separately. For further details of the Duo’s SCT algorithm, references [17] and [1] may be 
consulted. The results found in these two references were that the SCT algorithm (1) restored convergence 
in the infinity norm for C 0 type problems and (2) improved accuracy and observed rate of convergence for 
C 0 and C 1 test problems. Encouraged by the success of Duo’s SCT algorithm we decided to implement a similar 
algorithm for three-dimensional Cartesian geometry." (1, p.105)

.. _advantages:
*************************************
Advantages & Disadvantages Discussion
*************************************

"The fastest executing methods are as expected the zeroth order Diamond Difference and
the Linear Discontinuous method. These methods are followed by the Linear-Linear and the
Linear Nodal methods which are about five and 9 times slower, respectively. The five fastest
methods are either constant or linear approximation (with reduced number of cross moments)
and none of them need to call an external linear solver subroutine either because the linear
system is presolved or because no linear system has to be solved.
With increasing the expansion order, the grind time increases dramatically which is mainly
driven by the linear solve time t s which makes up the fastest growing part of the grind time: the
LU decomposition’s execution time scales cubically with the number of degrees of freedom of the
linear system of equations, i.e. Λ 6 . It is therefore not surprising that among the arbitrary order
methods the DGFEM with complete function space requires the least execution time, HODD is
only marginally cheaper than AHOTN, and DGFEM with Lagrange function space and Λ > 1
surprisingly takes the longest execution time. The reason why DGLA-Λ with Λ > 1 features
much longer execution times than AHOTN or HODD of the same order is the significantly more
expensive solution of the linear system of equations. We conjecture that the structure of the
DGLA matrices causes the Lapack routine dgesv to execute slower." (1, p.103)

.. _assumptions:
***********
Solver Assumptions
***********

1.  All systems are steady state (no time dependence)
2.  All mediums are non-multiplying (sigmaf = 0)
3.  Isotropic scattering present.

.. _refs:
**********
References
**********

1. SCHUNERT, SEBASTIAN. Development of a Quantitative Decision Metric for Selecting the
Most Suitable Discretization Method for S N Transport Problems. (Under the direction of
Yousry Y. Azmy.)
36.  Add 36 from thesis!
