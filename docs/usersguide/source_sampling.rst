.. _usersguide_source_sampling:

==========================
Mesh-Based Source Sampling
==========================

.. currentmodule:: pyne.source_sampling

.. automodule:: pyne.source_sampling

The source sampling module implements mesh-based source sampling (both
Cartesian and tetrahedral) which can be used as a component within Monte Carlo
radiation transport codes. A MOAB mesh tagged with energy-wise source densities
(and optional biased source densities) is supplied by the user. The user
can then supply six pseudorandom numbers in order to generate a random sample
of the particle birth parameters (position, energy, statistical weight).  The
source sampling module is written in C++ and has fully-supported C++, Fortran,
and Python interfaces, which facilitates its use within physics codes. The
source sampling module allows for three sampling modes:

:analog:
  Particle birth parameters are sampled directly from a unmodified 
  probability distribution function created from the source density mesh 
  (i.e. positions/energies with high source density are sampled more often
  than those of low source density). 
:uniform:
  All mesh volume elements are sampled with equal probability. Energy bins are
  sampled in analog from the distribution within a given mesh volume element.
  Statistical weights of particles are modified accordingly.
:user:
  In addition to source densities, the user supplies (on the same mesh) biased
  source densities. Particle birth parameters are then sampled on the basis of
  the biased source densities, and the statistical weight of particles is
  modified accordingly. The biased source density tag has the same length as the
  source density tag. Alternatively, the tag may have a length of 1, in which
  case the bias is only applied spatially and energy groups are sampled in
  analog.

A complete description of the theory involved can be found in the 
source_sampling entry in the PyNE theory manual.


*************
C++ interface
*************

An object of the Sampler class is first instantiated using 1 of 2 constructors:

:analog and uniform constructor:
    Sampler(std::string filename,
            std::string src_tag_name,
            std::vector<double> e_bounds,
            bool uniform)
:user constructor:
    Sampler(std::string filename,
            std::string src_tag_name,
            std::vector<double> e_bounds,
            std::string bias_tag_name)

The "filename" is a MOAB mesh file (.h5m). The "src_tag_name" is the name
of the tag on the mesh that stores the source densities. The source density can
be specified for an arbitrary number of energy groups, stored as a MOAB vector
tag. The "e_bounds" parameter describes the upper and lower bounds of these 
energy groups. For example if the src_tag_name tag is a vector of length 175
(for 175 energy groups), e_bounds should be of length 176. The final parameter
differs for the two constructors. A bool can be supplied: false for analog
sampling or true for uniform sampling. Alternatively, a string can be supplied
to denote the name of a tag (on the same mesh) which supplies biased source
densities.

Once a Sampler object is created the Sampler.particle_birth() method is called.
This method takes a single argument: a vector of 6 pseudorandom number between
0 and 1. This method returns a vector of length 5 containing the sampled x position, y position,
z position, energy, and weight respectively.
eight
An example C++ program is supplied below. This program requires a mesh file
named "source.h5m" with a tag named "source_density" of length 1.

.. code-block:: cpp

  #include "stdlib.h"
  #include  <iostream>
  #include "pyne/source_sampling.h"
  
  int main(){
  
    std::string filename("source.h5m");
    std::string src_tag_name("source_density");
    std::vector<double> e_bounds;
    e_bounds.push_back(0); // 1 energy group, lower bound of 0 upper bound of 1
    e_bounds.push_back(1);
  
    pyne::Sampler s(filename, src_tag_name, e_bounds, true);
  
    std::vector<double> rands;
    int i;
    for(i=0; i<6; i++) rands.push_back((double)rand()/RAND_MAX);
  
    std::vector<double> samp = s.particle_birth(rands);
  
    std::cout<<"x: "<<samp[0]<<std::endl;
    std::cout<<"y: "<<samp[1]<<std::endl;
    std::cout<<"z: "<<samp[2]<<std::endl;
    std::cout<<"e: "<<samp[3]<<std::endl;
    std::cout<<"w: "<<samp[4]<<std::endl;
  
   return 0;
  } 

This program can be complied with:

.. code-block:: bash

  g++ test.cpp -o test -lMOAB -lpyne


****************
Python interface
****************

The Python interface mainly exists for the purpose of testing the Sampler class
with python.nose. It can be used in the same manner as the C++ class:

.. code-block:: python

 import numpy as np
 from random import uniform
 from pyne.source_sampling import Sampler
 
 s = Sampler("source.h5m", "source_density", np.array([0, 1]), True)
 samp = s.particle_birth([uniform(0, 1) for x in range(6)])
 
 print("x: {0}\ny: {1}\nz: {2}\ne: {3}\nw: {4}".format(
           samp[0], samp[1], samp[2], samp[3], samp[4]))


*****************
Fortran Interface
*****************

Because Fortran cannot store an instance of the Sampler class, to perform
source sampling from Fortran, a free-standing function "sampling_setup" is
called to create a global instance of the sampling class. This function takes
a single argument: an integer representing the problem mode (0: analog, 1:
uniform, 2: user). This function assumes the mesh file is "source.h5m" and
that the tag names are "source_density" and "biased_source_density". In 
addition, this function assumes that a file "e_bounds" is present which is
a plain text file containing the energy boundaries.

An example program using the Fortran interface is shown below:

.. code-block:: fortran

 program test
   implicit none
   double precision :: x, y, z, e, w
   double precision, dimension(6) :: rands
   integer:: i, j, mode

   mode = 1
   call sampling_setup(mode)
 
   do j=1,6
     rands(j) = RAND()
   end do

   call particle_birth(rands, x, y, z, e, w)
   print*, "x: ", x
   print*, "y: ", y
   print*, "x: ", z
   print*, "e: ", e
   print*, "w: ", w

 end program test

This program can be compiled like:

.. code-block:: bash

  gfortran test.F90 -lpyne -lstdc++ -o test

************************
Source Sampling in MCNP5
************************

Standard MCNP5 ships with an empty source subroutine "source.F90" which can
be completed by the user in order to implement any form of custom source
sampling. A source.F90 file has been written to allow for the use of PyNE
source sampling within MCNP5. This file is found in pyne/share/source.F90.
The simplest way to compile MCNP5 with the source subroutine is as follows:

  #. Obtain a copy of the MCNP5 source code.
  #. Navigate to the folder MCNP5/Source/src.
  #. Soft-link the following files into this folder:

     a. pyne/src/source_sampling.cpp
     b. pyne/src/source_sampling.h
     c. pyne/src/measure.cpp
     d. pyne/src/measure.h

  #. Remove the pre-existing empty source.F90 file.
  #. Soft-link pyne/src/source.F90.
  #. Open the file MCNP/Source/src/FILE.list.
  #. Edit line 78 to include the additional source files. It should look like "CXX_SRC := measure.cpp source_sampling.cpp".
  #. Compile MCNP5 using the standard build method.

Once MCNP5 is compiled, MCNP5 can be run normally. The file "source.h5m" must
be present in the working directory that MCNP5 is run from. This file should
contain source densities (on the "source_density" tag) and optionally biased
source densities (the "biased_source_density" tag). An "idum" card must be used
in the MCNP5 input file. This card should have two arguments. The first is the
sampling mode (0: analog, 1: uniform, 2: user). The second is the resample
limit for void rejection. For a given particle, if a source position is
selected in void (MCNP materal 0) the source position is resampled within the
selected mesh volume element until either a non-void position is found, or this
user-specified limit is researched.

For example, this "idum" card specifies uniform sampling with a resample limit
of 100:

.. code-block:: bash

  idum 1 100

