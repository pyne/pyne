.. _pyne_spatialsolver:

==================================================
Spatialsolver Support -- :mod:`pyne.spatialsolver`
==================================================
Spatialsolver is a pyne module that contains seven neutron transport equation solvers.  Each
solver is its own unique nodal method.

The solvers included in this module are listed below.  The theory and methodology behind each
can be found in the pyne theory doccumentation.
 #. **AHOTN-LN**: Arbitrarily higher order transport method of the nodal type linear-nodal method
 #. **AHOTN-LL**:  Arbitrarily higher order transport method of the nodal type linear-linear method
 #. **AHOTN-NEFD**: Arbitrarily higher order transport method of the nodal type that makes use of the
    unknown nodal flux moments (NEFD algorithm).
 #. **DGFEM-LD**: The Discontinuous Galerkin Finite Element Method (DGFEM) with a linear discontinuous (LD)
    approximation for angular flux. (SEE PAGE 27 of thesis)
 #. **DGFEM-DENSE**: The Discontinuous Galerkin Finite Element Method (DGFEM) that use ??dense??? lagrange
    polynomials to "create a function space per dimension" [add citation thesis page 27].
Use Lagrange polynomials to create a function
space per dimension, i.e. in x, y and z direction.
 #. **DGFEM-LAGRANGE**:   The Discontinuous Galerkin Finite Element Method (DGFEM) that use lagrange
    polynomials to "create a function space per dimension" [add citation thesis page 27].
 #. **SCTSTEP**: SCT Step algorithm similar to Duo's SCT algorithm implemented in three dimensional Cartesian
    geometry.

As these are complicated solvers, they require a large amount of input data supplied by the user.  The
format we choose to take all this information in by is with a python dictionary.   Of the many key-pair values listed below, most are required, but some are optional.  The optional entries will be overriden by default values if not present/not specified. 

.. _dictionary_entries:

 #. **solver (xx)**: xxxxx

.. currentmodule:: pyne.spatialsolver

All functionality may be found in the ``spatialsolver`` package::

 from pyne import spatialsolver

Spatialsolver API
-----------

.. automodule:: pyne.spatialsolver
    :members:

.. _Spatialsolver: http://something.com/
