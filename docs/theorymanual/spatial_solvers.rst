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
The following are all spatial discretization methods used to solve the Sn equations of the Neutron Transport Equation (nte).  There are three major assumptions each makes: (1) all systems are steady state (no time dependence), (2) in a non-multiplying medium (sigmaf = 0), and (3) featuring isotropic scattering.  The Sn equations make the transport equation solvable by solving the transport equation along specific segments of omega, such that omega(n)=(mu,eta,xi)^T, with n = 1,..,N.  A nodal method is then used to come up with a specific solution.  The nodal method is what differentiates each of these codes, as outlined below.


*****************************
Ahotn-LN theory
*****************************

*****************************
Ahotn-LL theory
*****************************

*****************************
Ahotn-NEFD theory
*****************************

*****************************
DGFEM-LD theory
*****************************

*****************************
DGFEM-LL theory
*****************************

*****************************
DGFEM-Lagrange theory
*****************************

*****************************
Sct-step theory
*****************************

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



