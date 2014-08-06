.. _theorymanual_template:

===============================
Ahot Spatial Solvers
===============================

.. currentmodule:: pyne.spatialsolver

This page presents a detailed explanation of how :py:mod:`spatialsolver` works, including
underlying theory, any assocaited assumptions, and information about correct and
appropriate use of this physics. 


*****************************
Fundamental solver theory
*****************************

@TODO: add latex formatting for variables/equations, etc

"The linear Boltzmann transport equation describes the evolution of the flux of neutral particles,
i.e. neutrons or photons, in a host medium. It can be obtained from the general Boltzmann
transport equation by neglecting particle-particle interactions, the dependence of the material
properties of the host medium on the particle flux, and assuming that no electric force field is
present." [Thesis page 5 section 1.1].
The following are all spatial discretization methods used to solve the Sn equations of the Neutron Transport Equation (nte).  There are three major assumptions each makes: (1) all systems are steady state (no time dependence), (2) in a non-multiplying medium (sigmaf = 0), and (3) featuring isotropic scattering.  The Sn equations make the transport equation solvable by solving the transport equation along specific segments of omega, such that omega(n)=(mu,eta,xi)^T, with n = 1,..,N.  A nodal method is then used to come up with a specific solution.  The nodal method is what differentiates each of these codes, as outlined below.


*****************************
General AHOTN theory
*****************************

*****************************
Ahotn-LN theory 
*****************************

Linear Nodal Method:

Advantage: Very efficient when computing integral quantities [Sebastian thesis abstract]

**Merge with LL description?

*****************************
Ahotn-LL theory
*****************************

Linear Linear Method:

"The LN and LL methods both utilize moments of the balance equations, Eq. 2.23, satisfying
m x + m y + m z ≤ 1 augmented by three WDD equations per dimensions. Thus, the full set of
LN and LL equations comprises four balance relations and twelve WDD equations. Within this
subsection the WDD equations for the LL and LN method are derived in three-dimensional
Cartesian geometry. In particular, the WDD equations for the x-direction will be developed
for each of these two methods, but equivalent equations can be derived using the exact same
procedure for the y and z-direction.
Applying the operator Eq. 2.43 to the transport equation leads to an ordinary differential
h, i
equation for ψ m
<<INSERT THESIS LATEX HERE>
(2.51)
where χ km k ,m p (x) with k = y or k = z is given by:
<<INSERT THESIS LATEX HERE>
and +k = N, T and −k = S, B. In addition, the indices m k and m p are defined as follows:
if k = y then m k = m y and m p = m z , while for k = z we have m k = m z and m p = m y . For
making Eq. 2.51 amenable to a solution χ km k ,m p (x)’s dependence on x is approximated by:
<<INSERT THESIS LATEX HERE>
(2.53)where λ = 0 for the LN method and λ = 1 for the LL method." [thesis p.38-39]

*****************************
Ahotn-NEFD theory
*****************************

?????? Method:

"As the AHOTN method can be conveniently cast into a WDD form with all the AHOTN
specifics lumped into the spatial weights, a standard WDD solver can be used to solve the
per-cell AHOTN system of equations. Typically, the WDD relations Eq. 2.48 are solved for
the outflow face moments and substituted into the nodal balance relations Eq. 2.23 which are
then solved for the (Λ + 1) 3 unknown nodal flux moments (NEFD algorithm),"[Thesis p.37]

*****************************
DGFEM General Theory
*****************************
"The discontinuous Galerkin finite element method (DGFEM) uses identical polynomial test
and trial function spaces that are typically substituted into the weak form and tested against
all members of the test space to obtain a per-cell system of equations." [Thesis p.25]

There are three families of solvers under the general DGFEM label, linearly-discontinuous, lagrange and complete.  The Lagrange family is more accurate, while the complete family excecutes faster.

*****************************
DGFEM-LD theory
*****************************

Linearly Discontinous Method:

The linear discontinuous DGFEM method (LD) is the special case of the complete DGFEM
method of order Λ = 1. It is special in that the local matrix T is of size 4 × 4 and therefore its
inverse can be precomputed thus saving execution time. Following [24] we decided to implement
the LD method distinctly from the arbitrary order complete DGFEM kernel in order to create
a highly optimized method." [Thesis p.92]

"The first DGFEM methods which he refers to as linear discontinuous (LD) method uses the following approximation for the angular flux:
i,h
ψ m
p m (r) ,
ψ n i,h (r) =
(2.29)
m≤1
where m = m x + m y + m z ." [Thesis p.26]

Part of the complete family?

*****************************
DGFEM-LL theory (linear-linear?)
*****************************

Part of the complete family?

*****************************
DGFEM-Lagrange theory
*****************************

Part of the lagrange family?

*****************************
Sct-step theory
*****************************

"For the implementation of the SCT-Step method, the mesh sweep has to account for the
possibility of multiple distinct outflow segments on a single cell face as, for example, depicted in
Fig. 4.4. It is at the heart of the algorithm that the outflow averages are not mixed across the
boundaries imposed by the singular planes because mixing would defeat the initial purpose of
this algorithm: separating the solution slices illuminated by different boundary faces. Within
the described work, the SCT-Step algorithm was included into a standard sweep algorithm
ˆ n , for example the x-index runs fastest,
that sweeps the cells in a certain order dependent on Ω
followed by the y-index and finally the z-index runs slowest.
The subcell expressions, Eqs. 4.47, 4.48 and 4.49, are applied locally, i.e. whenever a
cell intersected by at least a singular planes is encountered, Eqs. 4.47 and 4.48 are used to
h, i
compute the segment’s outflow and cell averages. Then the segment volume-averaged flux ψ  ̄ n,m
is immediately collapsed into one cell-averaged scalar flux using Eq. 4.49.
This is in contrast to the treatment of the face-averaged fluxes, because they need to be
stored by illumination segment. Take for example a face that is intersected by the singular
characteristic as depicted on the left in Fig. 4.4: In this case, the cell downstream across this
109face from the one just solved, will be intersected by the singular characteristic, and the three
face-averaged segment fluxes will be needed as input for the three instances of Eq. 4.48 to be
solved. Similarly, for a face intersected by a single singular plane, two solution segments exist
and the two face-averaged segment fluxes will be required as input.
In the absence of reflective boundary conditions, at most n a = J · K + J + 1 angular face
information sets need to be stored (compared to at least I · J · K scalar fluxes), and therefore
the SCT-Step method does not increase the memory consumption significantly to store three
instead of one angular face flux per cell face. The memory consumption for solving the angular
face fluxes therefore increases to n a = 3n a .
The SCT-Step method does not superimpose an additional grid over the Cartesian mesh, it
merely splits a mesh cell appropriately into segments and collapses them immediately before the
cell’s solution is complete. It is important to contrast this to the idea of creating an unstructured
mesh specifically to isolate the illumination segments from each other e.g. as suggested in [11].
The disadvantage of this latter approach is that the mesh becomes dependent on the angular
ˆ n so that, for a practical algorithm, separate meshes need to be created for each
direction Ω
discrete ordinate and, in addition, restriction and prolongation operators need to be devised to
exchange information between these meshes.
The SCT-Step method is added to the selection of promising discretization methods because
it is expected to perform well in C 0 configurations where it is expected to restore cell-wise
convergence." [thesis p.109-110]


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
DGLA matrices causes the Lapack routine dgesv to execute slower." [thesis p.103]






Mathematical Details, some example math syntax:

You can use LaTeX format in rst. To do inline expressions (what $expression$ would
do in LaTeX) :math:`A\vec{x} = \vec{b}`.

To have blocks of mathematical content, do this

.. math::
    
    A\vec{x} = \vec{b}

Support is limited to a subset of LaTeX math by the conversion required for many output formats.

***********
Assumptions
***********

Any assumptions (explicit or implicit) about this method that would impact use, conclusions, validity, etc.

**********************
Additional Information
**********************

Details about impact of this theory on method use.

**********
References
**********

All of the justification for us listening to you.



