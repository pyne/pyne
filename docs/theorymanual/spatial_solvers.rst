.. _theorymanual_template:

===============================
Ahot Spatial Solvers
===============================

.. currentmodule:: pyne.spatialsolver

This page presents a brief summary of the underlying methods used in the :py:mod:`spatialsolver` module. For full details, consult "Development of a Quantitative Decision Metric for Selecting the Most Suitable Discretization Method for S_N Transport Problems", a disseration by Sebastian Schunert (2013). Chapter 2 has a detailed description of the methods implemented in PyNE. All code in this module is derived from Dr. Schunert's PhD work.

*****************************
Fundamental solver theory
*****************************

The spatial solver module contains a suite of methods that solve the S_N formulation of the Boltzmann transport equation. What differentiates the methods are their treatment of the spatial aspect of the solution process. "The linear Boltzmann transport equation describes the evolution of the flux of neutral particles, i.e. neutrons or photons, in a host medium. It can be obtained from the general Boltzmann transport equation by neglecting particle-particle interactions, the dependence of the material properties of the host medium on the particle flux, and assuming that no electric force field is present" (p. 5). The S_N equations define the transport equation along specific segments of the angular variable, omega, such that omega(n)=(mu,eta,xi)^T, with n = 1,..,N.

%The following are spatial discretization methods used to solve the S_N equations of the Neutron Transport Equation (nte). There are three major assumptions each makes: (1) all systems are steady state (no time dependence), (2) in a non-multiplying medium (sigmaf = 0), and (3) featuring isotropic scattering. The S_N equations define the transport equation along specific segments of the angular variable, omega, such that omega(n)=(mu,eta,xi)^T, with n = 1,..,N. A nodal method is then used to come up with a specific solution.  The nodal method is what differentiates each of these codes, as outlined below.

*****************************
General AHOTN theory
*****************************

The arbitrarily high-order transport method of the nodal type is a Transverse Moment Type Method.  It is derived for arbitrary expansion orders with the order Λ

*****************************
Ahotn-LL/-LN theory
*****************************

The Linear Linear and Linear Nodal Methods use linear TMB approximations to solve the SN equations.

"TMB approximation is a good compromise
for achieving high accuracy at a reasonably short execution time, thus resulting in improved
computational efficiency." (1, p.38)

*****************************
Ahotn-NEFD theory
*****************************

?????? Method:

"As the AHOTN method can be conveniently cast into a WDD form with all the AHOTN
specifics lumped into the spatial weights, a standard WDD solver can be used to solve the
per-cell AHOTN system of equations. Typically, the WDD relations Eq. 2.48 are solved for
the outflow face moments and substituted into the nodal balance relations Eq. 2.23 which are
then solved for the (Λ + 1) 3 unknown nodal flux moments (NEFD algorithm),"(1, p.37)

*****************************
DGFEM General Theory
*****************************

The DGFEM solvers use the Discontinuous Galerkin Finite Element Method (DGFEM) to solve the SN equations.  These use identical polynomial test and trial function spaces that are typically substituted into the weak form and tested against all members of the test space to obtain a per-cell system of equations." (1, p.25)  "In summary two families of DGFEM function spaces are mostly used in discretizing the spatial variable in the SN approximation of the transport equation: (1) the complete family and (2) the Lagrange family." (1,p27) 

"Assume that we formulate our function spaces such that we solve for point values of the flux, i.e. we use Lagrange polynomials as basis functions. Then, in two-dimensional triangular geometry and three-dimensional tetrahedral geometry the complete basis would require one flux value per corner point. The Lagrange basis would introduce more degrees of freedom that are not associated with the flux values in the corner points. In two-dimensional and three-dimensional Cartesian geometry the Lagrange family would result in one flux value per corner point. The complete basis would result in less degrees of freedom. For Λ = 1 for example, the Lagrange function space seems for more natural for Cartesian meshes, while the complete family appears to be a more natural choice for triangles/tetrahedra." (1,p.28)

When comparing these three included DGFEM solvers with a fixed expansion order, the Lagrange family
is more accurate, while the complete family excecutes faster.

*****************************
DGFEM-Lagrange theory
*****************************

DGFEM-Lagrange is part of the lagrange family of solvers of the DGFEM type.  All lagrange solvers are
implemented on the lagrange function space, rather than the complete function space.  For an explanation 
of both spaces see 4.4.1 and 4.4.2 from [1].  "The DGFEM method for discretizing the SN equations was first suggested by Reed and Hill [36] for two-dimensional triangular cells using a basis of Lagrange polynomials: Each Lagrange basis function is associated with a support point at which its value is unity while it assumes
a zero value at all other support points. The unknowns in Reed’s methods are then the flux
values at the support points and the method’s order is related to the number of support points
within a single cell."

*****************************
DGFEM-LD theory
*****************************

Linearly Discontinous Method:  "The linear discontinuous DGFEM method (LD) is the special case of the
complete DGFEM method of order Λ = 1. It is special in that the local matrix T is of size 4 × 4 and
therefore its inverse can be precomputed thus saving execution time. Following [24] we decided to implement
the LD method distinctly from the arbitrary order complete DGFEM kernel in order to create
a highly optimized method." (1, p.92)

********************************
DGFEM-LL theory (linear-linear?)
********************************

Part of the complete family?

*****************************
Sct-step theory
*****************************

The Sct-step solver is an implementation of Duo's Sct step algorithm in three dimensional cartesian geometry. [1, p.105]

SCT history and introduction to SCT-STEP:

"In addition, a novel method that explicitly tracks and eliminates lines and planes of non-
smoothness originating from “inconsistent” boundary conditions was developed and imple-
mented. This method uses the Step approximation in all cells that are intersected by lines
and planes of non-smoothness. It can be considered an extension of Duo’s Singular Character-
istic Tracking algorithm[17] to three spatial coordinates. However, it is important to point out
that this extension is highly non-trivial because of the tremendous increase in complexity of the
tracking and cell-splitting algorithms involved. The new method is labeled SCT-Step method." (1,p.4)

Maybe use this for the history and introduction instead:

"In two-dimensional geometry Duo[17] suggested tracking of the singular characteristic line
through the mesh and applying a sub-cell approach in intersected cells to keep segments in
these cells isolated from each other. For the solution of the subcell equations, Duo used the
Step Characteristic method applied to each of the segments separately. For further details of
the Duo’s SCT algorithm, references [17] and [1] may be consulted. The results found in these
two references were that the SCT algorithm (1) restored convergence in the infinity norm for
C 0 type problems and (2) improved accuracy and observed rate of convergence for C 0 and C 1
test problems. Encouraged by the success of Duo’s SCT algorithm we decided to implement a
similar algorithm for three-dimensional Cartesian geometry." (1, p.105)

SCT

*************************************
Advantages & Disadvantages discussion
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

***********
Solver Assumptions
***********

1.  All systems are steady state (no time dependence)
2.  All mediums are non-multiplying (sigmaf = 0)
3.  Isotropic scattering present.

**********************
Additional Information
**********************



**********
References
**********

1. SCHUNERT, SEBASTIAN. Development of a Quantitative Decision Metric for Selecting the
Most Suitable Discretization Method for S N Transport Problems. (Under the direction of
Yousry Y. Azmy.)
36.  Add 36 from thesis!
