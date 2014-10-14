Monte Carlo Code-Specific Steps for DAG-Tripoli4
+++++++++++++++++++++++++++++++

There are three varieties of code-specific steps:

1. defining attributes of the geometry using groups in CUBIT
2. defining DAGMC runtime parameters using input file syntax
3. changes to the command-line


Geometry Metadata
--------

The current version of DAG-Tripoli4 allows the definition of material
compositions and boundary conditions in the geometry.

Assigning Materials & Densities
'''''''''''''''''''''''''''''''''

The generic workflow description includes details on
:ref:`grouping-basics`, but a specific naming convention is required
for DAG-Tripoli4. To define materials, the Tripoli composition name
must be provided in the group name. The format for the group name is
as follows:
::
     comp_[compname]

``[compname]`` is replaced by the composition name that will be
specified in the Tripoli4 input file.  ``[density]`` is replaced by
either the atomic density or the mass density.  Like the MCNP cell
cards, positive values are atomic densities in [atoms/barn-cm] and
negative values are mass densities in [g/cc].

For example, suppose UO2 is composition ``Uoxide`` in the problem and
volumes 4 through 18 consist of this material.  To assign materials to
these volumes, the following command would be used:
::
     group "comp_Uoxide" add vol 4 to 18

*Note: If a volume is not assigned to a specific group, when run in
DAGMC it will be treated as a void; the material for that cell will be
zero. This can actually become a fairly useful debugging tool to
identify volumes that were not assigned to their appropriate group.*

If you would like to assign a material to the explicit complement, you
can use the same mechanism, but add the string ``_comp`` to the end of
the group name.  Since DAGMC only recognizes those groups that contain
an entity, it is necessary to add a volume to this group, recognizing
that this material will NOT be assigned to that volume.  (_It is often
convenient to assign the graveyard volume (see below) to the implicit
complement material group to minimize confusion._) For example, if you
would like the explicit complement to be modeled with the composition
named "air":
::
     create group "comp_air_comp"

Defining Boundary Conditions
'''''''''''''''''''''''''''''''

There are two general classes of boundary condition supported by
DAG-Tripoli4. a vacuum boundary and reflecting surfaces, and they are
implemented in different ways.

* **Defining the "graveyard": vacuum boundaries**

A vacuum boundary condition is typically defined in Tripoli4 by simply
having a surface with not defined volume on its other side.  Since
DAGMC's implicit complement is also defined this way, it is not
possible to use this convention in DAG-Tripoli4.  Instead, volumes are
created on the vacuum side of those same surfaces, but those volumes
are placed in a special group, known as the "graveyard" to change
their behavior.  Any geometric shape can be used for this; however, a
cubic shell is often preferred.  This shell volume will represent the
outside world and will be the cell where all of the particles are
terminated (thus it will have an importance of zero).

[need more on implicit complement in context of Tripoli]

To create this volume create two volumes in CUBIT with the same shape,
same center, and one slightly larger than the other.  Subtract the
smaller from the larger.  The remaining volume is the graveyard.

Like the material definitions and boundary conditions discussed in the
previous section. The graveyard is defined by assigning it a specific
group name, one of the following keywords:
::
    graveyard
    outside.world
    rest.of.world

Consider a geometry with 99 volumes that all fit within a cube
centered at the origin with side-length 99 cm.  To create a graveyard
for this problem in CUBIT, you could issue the following commands:
::
    cubit_prompt> create brick x 100
    cubit_prompt> create brick x 105
    cubit_prompt> subtract vol 100 from vol 101
    cubit_prompt> group "graveyard" add vol 102

When DAG-Tripoli4 is run, the surfaces of volume 102 (or any other
volumes included in the group) will be defined as having only one
volume - the one on the OTHER SIDE relative to voluem 102. (_Note:
this assumes that the two ``create brick`` commands generate volumes
numbered 100 and 101, respectively, and that the Boolean subtraction
results in a new volume number 102.

If you have boundary conditions (reflecting, white, or periodic) it is
not required that you surround them with the bounding volume, but is
not incorrect to do so.  Only areas where particles should escape need
to be enclosed.  However, it is often easiest to simply create a
single graveyard that covers all directions and volumes of the system.

* **Surface boundary conditions: reflection**

Surface boundary conditions are similarly enforced by specifying a
group name. This type of attribute (surface boundary condition) is
only required if reflective or white boundary conditions are used in
the problem.  If not, this section may be skipped.  *Note that
periodic boundary conditions are not yet supported.*

Specifying reflecting and white boundary conditions are fairly
straightforward.  The group names for reflecting and white are
respectively:
::
    spec.reflect
    white.reflect

Suppose surfaces 10 and 11 are reflecting boundary conditions.  To
specify these as reflecting surfaces, the following group would be
created:
::
    group "spec.reflect" add surf 10 11

DAGMC Runtime Parameters
'''''''''''''''''''''''''''

The DAGMC-Tripoli input file is formatted just like any other Tripoli
input file but using the ``DAGMC_GEOMETRY`` block to indicate the
geometry.  This block has the following parameters:

+---------------------------------------+----------------------------------+
|<geometry_filename>                    |  This must be the first parameter|
+---------------------------------------+----------------------------------+
|facet_tol <double faceting tolerance>  | (optional: default=0.001)        |
+---------------------------------------+----------------------------------+
|facet_file <string faceting filename>  | (optional)                       |
+---------------------------------------+----------------------------------+
|check_src_cell <"off"|"false"|"no">    | (optional: default=on)           |
+---------------------------------------+----------------------------------+
|usecad <"on"|"true"|"yes">             | (optional: default=off)          |
+---------------------------------------+----------------------------------+
|distlimit <"on"|"true"|"yes">          | (optional: default=off)          |
+---------------------------------------+----------------------------------+
|tolerance <double ray firing tolerance>| (optional: default=1e-8)         |
+---------------------------------------+----------------------------------+


These parameters are described in :ref:`the workflow description
<additional_parameters>`.  In addition to many other Tripoli input
blocks, it is important to define the material compositions that have
been assigned in the previous step.

Running DAGMC-Tripoli
'''''''''''''''''''''''

Running DAGMC-Tripoli is identical to running the standard Tripoli.

