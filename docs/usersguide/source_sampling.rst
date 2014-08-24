.. _usersguide_source_sampling:

==========================
Mesh-Based Source Sampling
==========================

.. currentmodule:: pyne.source_sampling

.. automodule:: pyne.source_sampling

The source sampling mode impliments mesh-based source sampling which can be
used as a component within Monte Carlo radiation transport codes. A MOAB mesh
tagged with energy-wise source densities (and optional biased source densities)
is supplied by the the user. The user can then supply six pseudorandom numbers
in order to generate a random sample of the particle birth parameters
(position, energy, statistic weight).  The source sampling module is written in
C++ and has fully-supported C++, Fortran, and Python interfaces, which
facilitates its use within physics codes. The source sampling module allows for
three sampling modes:

:analog:
  Particles are directly from a unmodified probability distribution function
  created from the source density mesh (i.e. positions/energies with high
  source density are sampled more often than those of low source density. 
:uniform:
  All mesh volume elements and energy bins are equiprobable and the statistical
  weight of particles is modified accordingly.
:user:
  In addition to source densities, the user supplies (on the same mesh) biased
  source densities. Particle birth parameters are then sampled on the basis of
  the biaed source densities, and the statisical weight of particles is
  modified accordingly

A complete description of the theory involved can be found in the 
source_sampling entry in the PyNE theory manual.


*************
C++ interface
*************

****************
Python interface
****************
The 

.. code-block:: python

.. code-block:: fortran

 program psuedo_mcnp
   implicit none
   double precision :: xxx, yyy, zzz, erg, wgt
   double precision, dimension(6) :: rands
   integer:: i, j, mode
   integer, parameter :: out_unit=2
   mode = 1
   call mcnp_sampling_setup(mode)
 
   open(unit=out_unit,file="samples.out", action="write", status="replace")
 
   do i=1,5000
     do j=1,6
       rands(j) = RAND()
     end do
       call mcnp_particle_birth(rands, xxx, yyy, zzz, erg, wgt)
       write(out_unit,*) xxx, yyy, zzz, erg, wgt
   end do
 end program psuedo_mcnp


For a complete specification for the classes in the ``cccc`` module, please
refer to the Library Reference entry for :ref:`pyne_source_sampling`.
