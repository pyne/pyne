Monte Carlo Code-Specific Steps for FluDAG
+++++++++++++++++++++++++++++++++++++++++++++

There are three varieties of code-specific steps:

1. defining attributes of the geometry using groups in CUBIT
2. defining DAGMC runtime parameters using input file syntax
3. changes to the command-line


Geometry Metadata
''''''''''''''''''

In FluDAG, the geometry file can be used to define material 
assignments, and eventually we would like to add the capability to 
define boundary conditions and tally locations.
 
Assigning Materials & Densities
..................................

The generic workflow description includes details on
:ref:`grouping-basics` , but a specific naming convention is required
for FluDAG. To define materials, the FLUKA material name must be 
provided in the group name. The format for the group
name is as follows: :: M_[mat_name]

For example, suppose we wish to add volumes 1 through 5 to a group
that defineds the material to be iron, then the following command 
would be used.
::
    group "M_IRON" add volume 1 to 5
    
Compounds are also supported by FluDAG, for example, if we wish to have volume 6 
belong to a group whose material name is STAINLESS then we can can use 
::
    group "M_STAINLESS" add volume 6

Be aware that there are several predefined material names in Fluka, and they
are appropriately treated by FluDAG. 
    
*Note: All volumes must belong to a group, if they do not have any information
FluDAG will not assign material information.

The implicit complement is automatically assigned the value 1 + the id of the 
highest numberd volume.

Defining the Graveyard
..............................
* **Defining the "graveyard": vacuum boundaries**

A typical usage of Monte Carlo codes  include a volume that extends 
to infinity with an importance of 0 that bounds the active volumes of interest.
Since solid models cannot include any infinite volumes, it is
necessary to place a finite volume of importance 0 to define the
problem boundary. You will need to surround the entire geometry with a
shell of finite thickness, known as the "graveyard".  Any geometric
shape can be used for this; however, a cubic shell is usually preferred.  This
shell volume will represent the outside world and will be the volume
where all of the particles are terminated.

To create this volume create two volumes in CUBIT with the same shape,
same center, and one slightly larger than the other.  Subtract the
smaller from the larger.  The remaining volume is the graveyard.

Like the material definitions and boundary conditions discussed in the
previous section. The graveyard is defined by assigning it a specific
group name the following keyword:
::
    group "M_BLCKHOLE" add volume X
   
Consider a geometry with 99 volumes that all fit within a cube
centered at the origin with side-length 99 cm.  To create a graveyard
for this problem in CUBIT, you could issue the following commands:
::
    cubit_prompt> create brick x 100
    cubit_prompt> create brick x 105
    cubit_prompt> subtract vol 100 from vol 101
    cubit_prompt> group "M_BLCKHOLE" add vol 102


When FLuDAG is run the all particles that enter volumes in group "M_BLCKHOLE" 
will be killed, this is effectively the same as the concept of importance 
in MCNP.


Scoring Assignments
..................
We do not currently support scoring assignments through group names. The user must manually
add these to the Fluka input deck.

The proposed naming scheme would be the following, 
::
     group "[tally_type]_[particle_name]" add volume <list>
     
For example
::
     group "usrtrack_neutron" add volume 1 2 5 6
     group "usrbdx_proton" add volume 1 2 4 9


Preparing the FluDAG Input File
''''''''''''''''''''''''''''''''''''
The FluDAG (Fluka) input file will look almost identical to the originating
Fluka input file. The exception will be the removal of all data between
the cards GEOBEGIN and GEOEND, i.e. all native Fluka input data. The last entry 
on the line of GEOBEGIN should be FLUGG. 

For example the most simple valid Fluka geometry is as follows, 
::
     GEOBEGIN                                                              COMBNAME
         0    0          
     SPH S1         0.0 0.0 0.0 50.0
     CELL1        5 +S1
     CELL2        5 +S1
     GEOEND

To run this geometry with FluDAG, remove all data between GEOBEGIN and GEOEND, and 
switch the last entry to FLUGG, 
::
     GEOBEGIN                                                              FLUGG
     GEOEND


Running FluDAG
'''''''''''''''''''
Running FluDAG bears some similarity to running FLUGG, first create the CAD geometry of
the problem you wish to run. In order to produce the material assignment data from the 
CAD geometry we mus t


Running DAG-MCNP5 is identical to running the standard MCNP5, but a
few new keywords have been added to the command-line to specify the
necessary files.

:``gcad=<geom_file>``: (required) The ``geom_file`` is the geometry
                       file that contains your geometric model, either
                       in the ACIS (\*.sat) format or the MOAB (\*.h5m)
                       format.  If this entry is not present,
                       DAG-MCNP5 will assume that it is running in
                       standard MCNP5 mode.  This runtime parameter is
                       described in more detail above.

:``ftol=<faceting_tolerance>``: (optional) [default: 1e-3] This is a
                               real number that provides guidance to
                               the faceting engine regarding the
                               maximum distance between a facet and
                               the surface it is representing.  It is
                               only used when reading an ACIS (\*.sat)
                               ``geom_file``.  When reading a MOAB
                               (\*.h5m) file, the facets have already
                               been generated and this setting is
                               ignored.  This runtime parameter is
                               described in more detail above.

:``fcad=<facet_file>: (optional) The ``facet_file`` is written by
                           DAG-MCNP5 in the MOAB (\*.h5m) format.  When
                           an ACIS file is read by DAG-MCNP5, a number
                           of pre-processing and initialization steps
                           are necessary.  Since these can be time
                           consuming, the user has the option to
                           create a ``facet_file`` the first time that
                           they use a geometry and then use that
                           ``facet_file`` with the ``gcad`` keyword in
                           subsequent uses.  This runtime parameter is
                           described in more detail above.


:``lcad=<log_file>``: (optional) The ``log_file`` is a skeleton of an
                           MCNP file for the cells and surfaces in
                           your geometry.  This file is created by
                           DAG-MCNP5 to communicate the material
                           assignments, boundary conditions, and
                           tallies that you defined in your geometry.
                           If you give a name other than the default
                           (``lcad``) for this file on the command-line,
                           that file will be used instead of the one
                           generated automatically by DAG-MCNP5.  This
                           is useful to make small changes to your
                           material assignments and/or importances,
                           but **can not** be used to change the
                           geometry.  It is up to the user to ensure
                           that the ``log_file`` being used
                           corresponds to the geometry file in
                           question.  This runtime parameter is unique
                           to the DAG-MCNP5 implementation of DAGMC.
