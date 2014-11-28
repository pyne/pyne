Monte Carlo Code-Specific Steps for FluDAG
+++++++++++++++++++++++++++++++++++++++++++++

There are several varieties of code-specific steps:

1. defining attributes of the geometry using groups in CUBIT
2. producing material assignments in FLUKA input format from the h5m file, with the help of FluDAG
3. preparing the FLUKA input file for running with DAGMC
4. inserting the material assignments into the FLUKA input deck


Geometry Metadata
''''''''''''''''''

In FluDAG, the geometry file can be used to define material 
assignments, and eventually we would like to add the capability to 
define boundary conditions and tally locations.
 
Assigning Materials & Densities
..................................

The generic workflow description includes details on
:ref:`grouping-basics`, but a specific naming convention is required
for FluDAG. 
To define 
materials, the FLUKA material name must be 
provided in the group name. The format for the group
name is as follows:
::
    M_[material_name]

For example, suppose we wish to add volumes 1 through 5 to a group
that defines the material to be iron.  The following command 
would be used.
::
    group "M_IRON" add volume 1 to 5
    
This will produce in the input file,
::
    ASSIGNMA        IRON         1
    ASSIGNMA        IRON         2
    ASSIGNMA        IRON         3
    ASSIGNMA        IRON         4
    ASSIGNMA        IRON         5
    
Compounds are also supported by FluDAG, for example, if we wish to have volume 6 
belong to a group whose material name is STAINLESS then we can can use 
::
    group "M_STAINLESS" add volume 6

This will produce in the input file:
::
    MATERIAL                                        26                    STAINLES  
    *...+....1....+....2....+....3....+....4....+....5....+....6....+....7...
    ASSIGNMA    STAINLES         6

*Note: Material names longer than 8 characters are truncated to the first 8 
characters. 

Be aware that there are several predefined material names in FLUKA, and they
are appropriately treated by FluDAG. 
    
*Note: All volumes must belong to a group, if they do not have any information
FluDAG will not assign material information.

The implicit complement is automatically assigned the value 1 + the id of the 
highest numberd volume.

Defining the Graveyard
..............................

A typical usage of Monte Carlo codes include a volume that extends 
to infinity with an importance of 0 that bounds the active volumes of interest.
Since solid models cannot include any infinite volumes, it is
necessary to place a finite volume of importance 0 to define the
problem boundary.  You will need to surround the entire geometry with a
shell of finite thickness, known as the "graveyard".  Any geometric
shape can be used for this, however a cubic shell is usually preferred.  This
shell volume will represent the outside world and will be the volume
where all of the particles are terminated.

To create this volume create two volumes in CUBIT with the same shape,
same center, and one slightly larger than the other.  Subtract the
smaller from the larger.  The remaining volume is the graveyard.

As with the material definitions discussed in the previous section the 
graveyard is defined by assigning the volume a keyword
group name,
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


When FLuDAG is run all the particles that enter volumes in group "M_BLCKHOLE" 
will be killed.  This is effectively the same as the concept of importance 
in MCNP.


Scoring Assignments
..................
We do not currently support scoring assignments through group names. The user must manually
add these to the FLUKA input deck.

The proposed naming scheme would be the following, 
::
     group "[tally_type]_[particle_name]" add volume <list>
     
For example
::
     group "usrtrack_neutron" add volume 1 2 5 6
     group "usrbdx_proton" add volume 1 2 4 9


Preparing the FluDAG Input File
''''''''''''''''''''''''''''''''''''
The FluDAG (FLUKA) input file will look almost identical to the originating
Fluka input file. The exception will be the removal of all data between
the cards GEOBEGIN and GEOEND, i.e. all native Fluka input data. The last entry 
on the line of GEOBEGIN should be FLUGG. 

For example the most simple valid FLUKA geometry is as follows, 
::
     GEOBEGIN                                                              COMBNAME
         0    0          
     SPH S1         0.0 0.0 0.0 50.0
     CELL1        5 +S1
     CELL2        5 -S1
     GEOEND

To run this geometry with FluDAG, remove all data between GEOBEGIN and GEOEND, and 
switch the last entry to FLUGG, 
::
     GEOBEGIN                                                              FLUGG
     GEOEND


Running FluDAG
'''''''''''''''''''
Running FluDAG bears some similarity to running FLUGG: the first step is to create the CAD 
geometry of the problem you wish to run. In order to produce the material assignment 
data from the CAD geometry we must first facet the file:
::
     dagmc_preproc -f <facet_tol> <cad_file.sat> -o <name.h5m>
     
This will facet the geometry file to a tolerance of <facet_tol> and produce a faceted file
called <name.h5m>. From that facet file we can produce the material "snippet" file
::
     /path/to/fludag/executable/mainfludag <name.h5m>
     
Will load the named h5m file and produce the material assignments information. 
This information should then be pasted into the FLUKA input file and any adjustments
that need to be made should be made, for example adding the density of non standard 
materials, or adding your scoring information. **Please note that the user must always 
include the additional material and compound information themselves and take
responsibility to ensure that the FLUKA material index number does not overlap with one
produced by FluDAG.**

The FluDAG calculation is now ok to run, 
::
     $FLUPRO/flutil/rfluka -e <path/to/fludag/executable/mainfludag> \
     -d <path/to/h5m/file/name.h5m> \
     ++{standard fluka options}++ <fludag_input_file>

