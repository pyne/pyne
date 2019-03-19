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
of the particle birth parameters (position, energy, statistical weight and cell
number). The source sampling module is written in C++ and has fully-supported
C++, Fortran, and Python interfaces, which facilitates its use within physics
codes. The source sampling module allows for six sampling modes:

:DEFAULT_ANALOG (mode 0):
  Particle birth parameters are sampled directly from a unmodified 
  probability distribution function created from the source density mesh voxels
  (i.e. positions/energies with high source density are sampled more often
  than those of low source density). 
:DEFAULT_UNIFORM (mode 1):
  All mesh volume elements are sampled with equal probability. Energy bins are
  sampled in analog from the distribution within a given mesh volume element.
  Statistical weights of particles are modified accordingly.
:DEFAULT_USER (mode 2):
  In addition to source densities, the user supplies (on the same mesh) biased
  source densities. Particle birth parameters are then sampled on the basis of
  the biased source densities, and the statistical weight of particles is
  modified accordingly. The biased source density tag has the same length as the
  source density tag. Alternatively, the tag may have a length of 1, in which
  case the bias is only applied spatially and energy groups are sampled in
  analog.
:SUBVOXEL_ANALOG (mode 3):
  Similar to DEFAULT_ANALOG, but the probability distribution function is
  created from the source density of mesh subvoxels.
:SUBVOXEL_UNIFORM (mode 4):
  Similar to DEFAULT_UNIFORM, but the probability distribution function is
  created from the source density of mesh subvoxels.
:SUBVOXEL_USER (mode 5):
  Similar to DEFAULT_USER, but the probability distribution function is
  created from the source density of mesh subvoxels. The tag may have a length
  of 1, number of photon energy groups, or maximum cells number in a voxel
  multiplied by the number of photon energy groups.

A complete description of the theory involved can be found in the 
source_sampling entry in the PyNE theory manual.


*************
C++ interface
*************

An object of the Sampler class is first instantiated using a constructor:

:constructor:
    Sampler(std::string filename,
            std::map<std::string, std::string> tag_names,
            std::vector<double> e_bounds,
            int mode)

The "filename" is a MOAB mesh file (.h5m). The "tag_names" is a map that stores
all the tag names need in the problem, such as the "src_tag_name" (required for
all modes), "bias_tag_name" (required in mode 2 and 5), "cell_number_tag_name"
(required for mode 3, 4, 5), "cell_fracs_tag_name" (required for mode 3, 4, 5).
The source density can be specified for an arbitrary number of energy groups,
stored as a MOAB vector tag. The "e_bounds" parameter describes the upper and
lower bounds of these energy groups. For example if the src_tag_name tag is a
vector of length 24 (for 24 energy groups), e_bounds should be of length 25.
The final parameter determins the source sampling mode.

Once a Sampler object is created the Sampler.particle_birth() method is called.
This method takes a single argument: a vector of 6 pseudorandom number between
0 and 1. This method returns a source particle containing the sampled x
position, y position, z position, energy, weight and cell number respectively.

An example C++ program is supplied below. This program requires a mesh file
named "source.h5m" with a tag named "source_density" of length 1.

.. code-block:: cpp

  #include "stdlib.h"
  #include  <iostream>
  #include "pyne/source_sampling.h"
  
  int main(){
  
    std::string filename("source.h5m");
    std::map<std::string, std::string> tag_names;
    tag_names.insert(std::pair<std::string, std::string>  ("src_tag_name",
    "source_density"));
    std::vector<double> e_bounds;
    e_bounds.push_back(0); // 1 energy group, lower bound of 0 upper bound of 1
    e_bounds.push_back(1);
  
    pyne::Sampler sampler(filename, tag_names, e_bounds, 0);
  
    std::vector<double> rands;
    int i;
    for(i=0; i<6; i++) rands.push_back((double)rand()/RAND_MAX);
  
    pyne::SourceParticle s = sampler.particle_birth(rands);
  
    std::cout<<"x: "<<s.get_x()<<std::endl;
    std::cout<<"y: "<<s.get_y()<<std::endl;
    std::cout<<"z: "<<s.get_z()<<std::endl;
    std::cout<<"e: "<<s.get_e()<<std::endl;
    std::cout<<"w: "<<s.get_w()<<std::endl;
    std::cout<<"c: "<<s.get_c()<<std::endl;
  
   return 0;
  } 

This program can be complied with:

.. code-block:: bash

  g++ test.cpp pyne/source_sampling.cpp pyne/measure.cpp -o test -lMOAB -lpyne


****************
Python interface
****************

The Python interface mainly exists for the purpose of testing the Sampler class
with python.nose. It can be used in the same manner as the C++ class:

.. code-block:: python

 import numpy as np
 from random import uniform
 from pyne.source_sampling import Sampler, SourceParticle
 
 tag_names = {"src_tag_name": "source_density"}
 sampler = Sampler("source.h5m", tag_names, np.array([0, 1]), 0)
 s = sampler.particle_birth([uniform(0, 1) for x in range(6)])
 
 print("x: {0}\ny: {1}\nz: {2}\ne: {3}\nw: {4}\nc: {5}".format(
           s.x, s.y, s.z, s.e, s.w, s.c))


*****************
Fortran Interface
*****************

Because Fortran cannot store an instance of the Sampler class, to perform
source sampling from Fortran, a free-standing function "sampling_setup" is
called to create a global instance of the sampling class. This function takes
a single argument: an integer representing the problem mode (0: DEFAULT_USER, 1:
DEFAULT_UNIFORM, 2: DEFAULT_USER, 3: SUBVOXEL_ANALOG, 4: SUBVOXEL_UNIFORM,
5: SUBVOXEL_USER). This function assumes the mesh file is "source.h5m" and
that the tag names are "source_density", "biased_source_density",
"cell_number_tag_name" and "cell_fracs_tag_name". In addition, this function
assumes that a file "e_bounds" is present which is a plain text file containing
the energy boundaries.

An example program using the Fortran interface is shown below:

.. code-block:: fortran

 program test
   implicit none
   double precision :: x, y, z, e, w
   double precision, dimension(6) :: rands
   integer:: i, j, mode, c

   mode = 1
   call sampling_setup(mode)
 
   do j=1,6
     rands(j) = RAND()
   end do

   call particle_birth(rands, x, y, z, e, w, c)
   print*, "x: ", x
   print*, "y: ", y
   print*, "x: ", z
   print*, "e: ", e
   print*, "w: ", w
   print*, "c: ", c

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
  #. Soft-link pyne/share/source.F90.
  #. Open the file MCNP/Source/src/FILE.list.
  #. Edit line 78 to include the additional source files. It should look like "CXX_SRC := measure.cpp source_sampling.cpp".
  #. Compile MCNP5 using the standard build method.

Once MCNP5 is compiled, MCNP5 can be run normally. The file "source.h5m" and
"e_bounds" must be present in the working directory that MCNP5 is run from.
The file "source.h5m" should contain source densities (on the "source_density"
tag) and optionally biased source densities (the "biased_source_density" tag).
The file "e_bounds" should contain the energy boundaries of the photon energy
groups used in the activation calculations. An "idum" card must be used
in the MCNP5 input file. This card should have three arguments. The first is the
sampling mode (0: DEFAULT_ANALOG, 1: DEFAULT_UNIFORM, 2: DEFAULT_USER,
3: SUBVOXEL_ANALOG, 4: SUBVOXEL_UNIFORM, 5: SUBVOXEL_USER). The second is the
resample limit for void and cell rejections. For a given particle, if a source
position is selected in void (MCNP material 0) or in a cell that disagrees with the
cell number, the source position is resampled within the selected mesh volume
element until either a correct position is found, or this user-specified limit
is researched. The third argument should specify the particle type: 1 for
neutrons, 2 for photons.

For example, this "idum" card specifies uniform sampling with a resample limit
of 100 with source particles specified as photons:

.. code-block:: bash

  idum 1 100 2

