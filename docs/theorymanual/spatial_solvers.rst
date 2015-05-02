.. _theorymanual_spatialsolver:

===============================
Spatial Solvers
===============================

.. currentmodule:: pyne.spatialsolver

This page presents a brief summary of the implemented spatial discretization  methods used in the 
:py:mod:`spatialsolver` module. For full details, consult "Development of a 
Quantitative Decision Metric for Selecting the Most Suitable Discretization 
Method for SN Transport Problems", a dissertation by Sebastian Schunert [Schunert]_.
Chapter 2 contains a detailed description of the methods implemented in PyNE. All 
code in this module is derived from Dr. Schunert's PhD work.

The first order multigroup SN transport equations are a set of hyperbolic partial
differential equations that are derived from the neutron transport equation by discretization
of energy (multigroup) and direction (SN method). A single SN equations  approximates the 
neutron angular flux along a specific angular direction omega, such that omega(n)=(mu,
eta,xi)^T, with n = 1,..,N, for a specific energy group.
The spatial solver module contains a suite of discretization methods for the spatial dependence of the 
SN transport equations. 

All of the methods discussed here belong to a broader class of methods referred to as
nodal methods. Nodal methods are characterized by the following properties:

* All function spaces are defined local to a mesh cell.
* Coupling between cells occurs only through their faces.
* Coupling between cells is only imposed in an integral sense.
* Increasing the order of the methods is achieved by increasing the local order of expansion.
* As the test space typically contains a constant test function, nodal methods are conservative.

Nodal methods are closely related to discontinuous Finite Element Methods (DFEM) - a class of 
methods that was successfully deployed for the discretization of hyperbolic PDEs. 
Further, they are related to certain classes of Finite Volume Methods [Hesthaven]_. 

In contrast to Finite Difference Methods, the unknowns of DFEM methods are expansion coefficients 
of the flux shape within a mesh cell. Therefore, DFEM methods always solve for the approximation 
of the flux shape, not just for disconnected point values. In this regard, they are similar to 
continuous FEM methods (CFEMS). However, the flux values on the faces between two cells are not unique 
for DFEMs, while they are unique for CFEMs. The particular choice of functions for approximating the 
flux shape within a cell distinguishes the various methods from each other.  

These following methods are implemented and described below:

* :ref:`ahotn`
* :ref:`ahotn-nefd`
* :ref:`ahotn-ll-ln`
* :ref:`dgfem`
* :ref:`dgfem-lagrange`
* :ref:`dgfem-complete`
* :ref:`dgfem-ld`
* :ref:`sct-step`

.. _ahotn:

*****************************
General AHOTN
*****************************

The Arbitrarily High Order Transport Method of the Nodal type (AHOTN) is a class 
of methods that was developed based on physical intuition for providing accurate 
solutions to the SN transport equations on optically thick cells. As the AHOTN methods
are based on taking transverse moments of the SN equations, they are also referred to 
as transverse moment based (TMB) method.

TMB methods can derived for an arbitrary expansion order Λ denoting the expansion
order of the volumetric source term into Legendre polynomials. 
The set of equations for each angular direction and spatial mesh cell consists of
two types of equations:

* Volumetric moments of the SN equations with respect to Legendre polynomials (balance equations).
  These equations have volume as well as face unknowns with the total number being larger than the
  number of equations.
* Closure relations derived via the transverse moment formalism. The transverse moment formalism 
  is applied for the x, y and z direction separately and yields a ODE in the corresponding variable,
  respectively ([Schunert]_, p.33).  

TMB methods are particularly good at resulting in accurate solutions on coarse 
spatial meshes ([Schunert]_, p.33).  This is required, because many traditional methods,
such as diamond difference (DD) fail when applied to coarse spatial meshes, 
resulting in either extremely inaccurate or negative solutions.  

The AHOTN methods are unique compared to most other TMB methods because they are
developed to have a very compact weighted diamond difference (WDD) representation 
of the per-cell set of equations.  A full derivation of the AHOTN solutions to 
the SN equations can be found on pages 34-40 of [Schunert]_.  

There are three TMB solvers accessible in PyNE, discussed further below:

1. :ref:`NEFD <ahotn-nefd>`
2. :ref:`Linear-Nodal (LN) <ahotn-ll-ln>`
3. :ref:`Linear-Linear (LL) <ahotn-ll-ln>`

.. _ahotn-nefd:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AHOTN-NEFD
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The AHOTN method can be conveniently cast into a WDD form with all the AHOTN
specifics lumped into spatial weights - one for each angular direction, spatial dimension and spatial mesh cell. 
Thus, a standard weighted diamond difference solver (WWD) that is available in most first order SN transport codes can be used to solve the per-cell AHOTN system of equations. 
Typically, the WDD relations are solved for the outflow face moments and substituted into the nodal
balance relations, which are then solved for the (Λ + 1)^3 unknown nodal flux 
moments (this is called the NEFD algorithm) ([Schunert]_, p.37).

It is worth noting that in the limit of infinitely small cells, the AHOTN method
becomes identical to the higher order diamond difference (HODD) method 
([Schunert]_, p.37). 

The order of trial functions used in the AHOTN method is denoted by appending a
"-#". E.g., AHOTN-1 uses a first order expansion. 

.. _ahotn-ll-ln:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AHOTN-LL/-LN
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Both the AHOTN Linear Linear (LL) and Linear Nodal (LN) methods are linear TMB methods
that increase the level of approximations in order to  streamline the execution of the
methods. The linear nature of this approximation allows for 
achieving high accuracy in reasonably short execution time, thus resulting in 
improved computational efficiency" ([Schunert]_, p.37). 
The major difference to the AHOTN-1 method is that only those volume moments are retained
whose sum (over x,y, and z) is equal or less than one. Further, the leakage terms are appropriately
truncated to not introduce more unknowns. This reduces the number of equations to four balance relations and twelve WDD equations
in three-dimensional geometry ([Schunert]_, p.38).

The difference between the LN and LL method is that the LL method retains the 
bilinear leakage component while the LN neglects it. From an algorithmic (i.e.
solution of the local equations within the mesh sweep) point of view, the LN 
provides the least coupling among the set of equations. Conversely, the LL has
stronger coupling in the WDD relations than both LN and AHOTN-1, and AHOTN-1 
has stronger coupling than LL and LN in the nodal balance equations [[Schunert]_, p.40].

.. _dgfem:

*****************************
General DGFEM 
*****************************

The classical DFEM method is based on approximating the flux shape with polynomials within 
each cell. To derive a set of equations for the unknown expansion coefficients, the weak form 
of the SN transport equations is used. The weak form is an integral statement that is obtained
by multiplying the SN transport equations with a test function and integrating it over the extend
of a single mesh cell; finally using Green's theorem yields the weak form. 
The weak form is used by selecting identical polynomial test and trial functions sets and substituting
them into the weak form. As there are as many test functions as there are expansion coefficients, a closed
set of equations is obtained ([Schunert]_, p.25).

Two families of DGFEM function spaces are most commonly used in discretizing the
spatial variable in the SN approximation of the transport equation ([Schunert]_, p.27):

1. the *complete* family and 
2. the *Lagrange* family 

The *complete* expansion of order Λ uses all Legendre polynomials whose orders sum at most Λ (Note that it does not
matter which type of polynomials are used, instead of Legendre polynomials the set of monomials spanning the same space
could be used). In contrast the *Lagrange* uses all Legendre polynomials whose maximum moment in x, y, or z dimension
is less or equal than Λ. Thus, the *Lagrange* family of the same order comprises more members than the *complete* 
family.   

There are three DGFEM-type solvers accessible in PyNE, discussed further below:

1. :ref:`Lagrange <dgfem-lagrange>`
2. :ref:`Linear-Discontinuous (LD) <dgfem-complete>`
3. :ref:`Linear-Linear (LL) <dgfem-ld>`


When comparing the three included DGFEM solvers with a fixed expansion order, the *Lagrange* family
is more accurate, but the *complete* family executes more quickly.

.. _dgfem-lagrange:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
DGFEM-Lagrange
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The DGFEM-Lagrange family solver is implemented using Legendre polynomial as 
basis functions (remember that Lagrange in this context does not refer to using
Lagrange polynomials but rather to the members of the function space, or more accurately
their span). The DGFEM-Lagrange solver can handle arbitrary expansion orders. 
The number of unknowns per cell and angular direction is (Λ+1)^3.

.. _dgfem-complete:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
DGFEM-complete
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The same comments as for the DGFEM-Lagrange apply but the number of unknowns is
(Λ+3)*(Λ+2)*(Λ+1)/6. For Λ=1 a special version, the DGFEM-LD solver, was implemented. 

.. _dgfem-ld:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
DGFEM-LD
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The linear discontinuous DGFEM method is the special case of the complete DGFEM
method of order Λ = 1. It is special in that the local matrix implementing the 
equations is of size 4 × 4 and, therefore, its inverse can be precomputed --
saving execution time ([Schunert]_, p.92).


.. _sct-step:

*****************************
SCT-STEP
*****************************

One of the problems  most spatial solvers encounter when solving the SN equations is the
limited smoothness of the exact solution. Depending on the boundary conditions, realistic SN problems
either have discontinuous angular fluxes or discontinuous first derivatives. The lack of smoothness
stems from the non-smoothness of the boundary conditions at the corners of the domain (edges in three-dimensional geometry). 
This leads to a line (planes in 3D) of discontinuity within the domain which is referred to as Singular Characteristic (SC) 
(Singular Planes (SPs) in 3D).

In particular for problems featuring a discontinuous angular flux, standard spatial discretization methods
deliver inaccurate results that, depending on the chosen error norm, either do not converge to the exact solution (cell wise 
error norm) or converge at an extremely small convergence rate. 

A new method called Singular Characteristic Tracking Step solver (SCT-STEP) was developed by Schunert
and Azmy to try to rectify this, by using a "step approximation in all cells
that are intersected by lines and planes of non-smoothness" ([Schunert]_, p.4).
The SCT-STEP method is an extension of Duo’s SCT Step algorithm to three-
dimensional cartesian geometry.

The basic idea for the SCT-STEP is borrowed from Duo's method: In 
two-dimensional geometry Duo suggested tracking the singular characteristic line
through the mesh and applying a sub-cell approach in intersected cells to keep 
segments in these cells isolated from each other. For the solution of the subcell
equations, Duo used the Step Characteristic method applied to each of the 
segments separately ([Schunert]_, p.105).

The extension to three-dimensional geometry requires tracking the SC and SPs through the domain
and applying the step approximation in cells that are intersected either by the SC or a SP.

.. _advantages:

*************************************
Advantages & Disadvantages Discussion
*************************************

The fastest executing methods are, as expected, the zeroth order Diamond Difference and
the Linear Discontinuous method. These methods are followed by the Linear-Linear and the
Linear Nodal methods, which are about five and nine times slower, respectively. The five fastest
methods are either constant or linear approximation (with reduced number of cross moments),
and neither of these need to call an external linear solver subroutine, either because the linear
system is pre-solved or because no linear system needs to be solved.
With increasing the expansion order, the computation time increases dramatically. This is mainly
driven by the linear solve time, which makes up the fastest growing part of the computation time: the
LU decomposition’s execution time scales cubically with the number of degrees of freedom of the
linear system of equations, i.e. Λ^6 . It is therefore not surprising that among the arbitrary order
methods, the DGFEM with complete function space requires the least execution time, and HODD is
only marginally cheaper than AHOTN. DGFEM with Lagrange function space and Λ > 1, however,
surprisingly takes the longest execution time. The reason why DGLA-Λ with Λ > 1 features
much longer execution times than AHOTN or HODD of the same order is the significantly more
expensive solution of the linear system of equations. Schunert conjectures that the structure of the
DGLA matrices causes the Lapack routine dgesv to execute slower" ([Schunert]_, p.103).

The ultimate performance indicator is the computational efficiency which is the ability to 
obtain accurate results in as short a computational execution time as possible. While a detailed 
discussion of the various methods' computational efficiency is comprised in [Schunert]_, the following
summarizes these findings:

* If the problem features a discontinuous angular flux and the error is measured in the cell-wise infinity norm or 2-norm, the :ref:`sct-step`  method is the most 
  efficient method.
* If the error is measured in an integral error norm, i.e. computing region averaged fluxes or reaction rates, the :ref:`ahotn-ll-ln` are the most efficient methods.
* If the flux is continuous and the error is measured in a infinity norm or 2-norm, higher order methods perform better than lower order methods. For optically thick 
  problems with cell aspect ratios close to one, the :ref:`ahotn-nefd` method is most efficient. If more skewed aspect ratios are considered, the 
  :ref:`dgfem-complete` method is the most efficient. 

.. _assumptions:

*************************************
Solver Assumptions
*************************************

1.  All systems are steady state (no time dependence).
2.  All mediums are non-multiplying (sigmaf = 0). This is a limitation of the iteration structure around the 
    spatial solvers. The spatial solvers are not limited to non-multiplying media.
3.  Isotropic scattering present.

.. _refs:

*************************************
References
*************************************

.. [Schunert] SCHUNERT, SEBASTIAN. Development of a Quantitative Decision Metric for Selecting the Most Suitable Discretization Method for S N Transport Problems. (Under the direction of Yousry Y. Azmy.) http://repository.lib.ncsu.edu/ir/bitstream/1840.16/9048/1/etd.pdf
.. [Hesthaven] HESTHAVEN, J.S. and WARBURTON, T. Nodal Discontinuous Galerkin Methods. http://www.springer.com/us/book/9780387720654
