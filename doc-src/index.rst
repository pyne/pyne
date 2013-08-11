.. DAGMC master file, created by
   sphinx-quickstart on Fri Aug 31 10:08:00 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

What is DAGMC?
======

The Direct Accelerated Geometry Monte Carlo [DAGMC] toolkit is a
component of the Mesh-Oriented datABase [MOAB_] that provides
fundamental functions for the ray-tracing and related geometry
operations of Monte Carlo radiation transport directly on complex 3-D
geometries such as those created in modern solid modeling software.

DAGMC is designed in a modular fashion with the expectation that it
can be integrated into a variety of Monte Carlo radiation tools. The
CNERG_ group at UW-Madison has been focusing its development on the
MCNP5_ software developed at `Los Alamos National Laboratory
<http://www.lanl.gov>`_ and distributed by the `Radiation Safety
Information Computing Center <http://rsicc.ornl.gov>`_.

We have prior experience integrating DAGMC with MCNPX, and ongoing
efforts to integrate DAGMC with other Monte Carlo physics packages
including: Fluka, Geant4, and Shift.

While we don't have a complete GUI, we currently rely on the Cubit_
software from Sandia.  It plays a role in our workflow that can
include importing CAD-files from other tools such as SolidWorks,
CATIA, etc.  A key technology for supporting different solid modeling
formats is CGM_.  In addition to defining the geometry, we rely on
Cubit for material assignment and can also support some other aspects
of input definition that are tied to the geometry (tallies and
variance reduction parameters, for example).  Some knowledge of MCNP
input files is necessary for other parameters such as material
definition, run control and source definition.

Both CGM_ and MOAB_ are developed by a team of collaborators at
Argonne, who are also working on improving some of the GUI tools
available for manipulating workflows like this.


.. toctree::
   :maxdepth: 1

   usersguide/index
   CNERG Support for DAGMC
   DagmcDevelopersGuide
   DagmcPublications
   Upcoming Features in DAGMC

.. _MOAB: http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB
.. _MCNP5: http://mcnp-green.lanl.gov/
.. _Cubit: http://cubit.sandia.gov
.. _CGM: http://trac.mcs.anl.gov/projects/ITAPS/wiki/CGM
.. _CNERG: http://cnerg.engr.wisc.edu

