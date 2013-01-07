Direct Accelerated Geometry Monte Carlo (DAGMC) Toolkit
==========================================================

The Direct Acclerated Geometry Monte Carlo (DAGMC) Toolkit is an
interface (`DAGMC Source
<http://trac.mcs.anl.gov/projects/ITAPS/browser/MOAB/trunk/tools/dagmc>`_)
to the `MOAB mesh database
<http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB>`_ that provides the
methods necessary for ray tracing on a CAD-based geometric model.

This repository provides implementations of that interface for various
Monte Carlo radiation transport codes.  Depending on the distribution
limitations and software design of each code, there are many modes by
which we are delivering this capability.  There are 3 main
characteristics of a Monte Carlo physics package which define our
method of distribution.

Efforts are underway to make DAGMC available in the following physics
packages:
   * MCNP5: complete and in production use
   * Fluka: just beginning (12/2012)
   * Serpent: underway (12/2012)
   * OpenMC: planned for 2013
   * GEANT4: planned for 2013
   * MCNP6: planned for 2013

Geometry Interface
-------------------

In cases where the physics package has a cleanly defined geometry
interface (FLUKA), we are able to distribute a standalone collection of
methods that each user can compile and link with the physics package.

When the geometry interface is not cleanly defined, our modifications
include modifications to the original source code.  Therefore our
distribution mechanism depends the ability to integrate our
modifications into the main physics code development path.

Mainline Development Integration
----------------------------------

In cases where the authors of the physics package are willing to
integrate DAGMC as an option in their primary software distribution
(SHIFT), then the main distribution mechansim will be as part of that
software in its normal distribution channels.  This may either in be
in a shared software repository or in a regular release snapshot.

In cases where the primary authors prefer DAGMC to be distributed
separately as an enhancement for their software (MCNP5), the
distribution mechanism will be as a patch to their source code, since
we generally are not authorized to redistribute their code.

Building Documentation
-------------------------

The documentation for DAGMC has been configured such that the source is 
in the master branch, to be distributed with the source code.  The HTML 
version can be generated and pushed automatically to the gh-pages branch 
by running::

     make -f doc-src/Makefile gh-pages

in the top-level directory.
