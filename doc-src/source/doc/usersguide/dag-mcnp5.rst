Monte Carlo Code-Specific Steps for DAG-MCNP5
+++++++++++++++++++++++++++++++++++++++++++++

There are three varieties of code-specific steps:

1. defining attributes of the geometry using groups in CUBIT
2. defining DAGMC runtime parameters using input file syntax
3. changes to the command-line


Geometry Metadata
''''''''''''''''''

In DAG-MCNP5, the geometry file can be used to define material and
density assignments, boundary conditions, and tally locations.

Assigning Materials & Densities
..................................

The generic workflow description includes details on
:ref:`grouping-basics` , but a specific naming convention is required
for DAG-MCNP5. To define materials, both the MCNP material ID and
density must be provided in the group name. The format for the group
name is as follows: :: mat_[matid]_rho_[density]

``[matid]`` is replaced by the material ID that will be specified in
the MCNP input file.  ``[density]`` is replaced by either the atomic
density or the mass density.  Like the MCNP cell cards, positive
values are atomic densities in [atoms/barn-cm] and negative values are
mass densities in [g/cc].

For example, suppose UO2 is material 7 in the problem, has an atomic
density of 0.0223 and volumes 4 through 18 consist of this material.
To assign materials to these volumes, the following command would be
used:
::
     group "mat_7_rho_0.0223" add vol 4 to 18

*Note: If a volume is not assigned to a specific group, when run in
DAGMC it will be treated as a void; the material for that cell will
be zero. This can actually become a fairly useful debugging tool to
identify volumes that were not assigned to their appropriate group.*

If you would like to assign a material to the explicit complement, you
can use the same mechanism, but add the string ``_comp`` to the end of
the group name.  Since DAGMC only recognizes those groups that contain
an entity, it is necessary to add a volume to this group, recognizing
that this material will NOT be assigned to that volume.  (*It is often
convenient to assign the graveyard volume (see below) to the implicit
complement material group to minimize confusion.*) For example, if you
would like the explicit complement to be modeled as material 9 with
density 1 g/cc:
::
     create group "mat_9_rho_-1_comp"

Defining Boundary Conditions
..............................

There are two general classes of boundary condition supported by
DAG-MCNP5. a vacuum boundary and reflecting surfaces, and they are
implemented in different ways.

* **Defining the "graveyard": vacuum boundaries**

A typical usage of MCNP5 includes a volume that extends to infinity
with an importance of 0 that bounds the active volumes of interest.
Since solid models cannot include any infinite volumes, it is
necessary to place a finite volume of importance 0 to define the
problem boundary. You will need to surround the entire geometry with a
shell of finite thickness, known as the "graveyard".  Any geometric
shape can be used for this; however, a cubic shell is preferred.  This
shell volume will represent the outside world and will be the cell
where all of the particles are terminated (thus it will have an
importance of zero).

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


When DAG-MCNP5 is run, the importance of volume 102 (or any other
volumes included in the group) will be set to zero. (_Note: this
assumes that the two ``create brick`` commands generate volumes
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

Tally Assignments
..................

It is also possible, although not required, to specify tallies in the
geometry.  The general form for adding this meta-data is to create a
group of volumes or surfaces and encode the meta-data in the names of
those groups.

The user has the option of specifying tallies in the geometry
directly.  It is still possible to specify tallies in the MCNP input
file, however, the user has to make sure that the tally indices are
not duplicated lest a fatal error will occur.  Tallies are specified
as group names in the following format:
::
      tally_[CUBIT tally ID].[tally type keyword].[particles]

The ``[CUBIT tally ID]`` field is an integer from 0 to 99.  Different
tally types may have the same CUBIT ID and are still consistent.  The
tally number in MCNP is 10 times the CUBIT ID plus the tally type
index (e.g. 4 for cell flux tallies).

The ``[tally type keyword]`` is one of the following for each type of
tally:

+----------+------------------+
|Tally Type|tally type keyword|
+----------+------------------+
|f1        |surf.current      |
+----------+------------------+
|f2        |surf.flux         |
+----------+------------------+
|f4        |cell.flux         |
+----------+------------------+
|f6        |cell.heating      |
+----------+------------------+
|f7        |cell.fission      |
+----------+------------------+
|f8        |pulse.height      |
+----------+------------------+

Also \*tallies (the tally result times the incident particle energy)
are possible by placing an "e" before the tally type.  So to make a
\*f2 tally, the keyword would be ``esurf_flux``.  Pulse height (f8) tallies
have the option to include charge as well.  This is done by placing a
"q" before the keyword as in ``qpulse_height``.

The ``[particles]`` tag is a string stating which particles will be
tallied.  To tally both photons and neutrons, set the tag to "np".
The default is neutrons only.  Should this be tag be omitted, only
neutrons will be tallied.

Some CUBIT commands to do tallies:
::
    group "tally_0.surf.current" add surf 1 to 4
    group "tally_0.cell.flux.p" add vol 7
    group "tally_1.ecell.heating.np" add vol 2 6
    group "tally_6.cell.heating.n" add vol 2 6
    group "tally_7.cell.flux.p" add vol 1 to 3
    group "tally_12.pulse.height.p" add vol 10 to 14
    group "tally_14.qpulse.height.p" add vol 10 to 14

The above are equivalent to following MCNP definitions:
::
    f1:n 1 2 3 4 T
    f4:p 7 T
    *f16:n,p 2 6 T
    f66:n 2 6 T
    f74:p 1 2 3 T
    f128:p 10 11 12 13 14 T
    +f148:p 10 11 12 13 14 T

*(Note: the current convention is to always add a tally bin for the
total across all cells/volumes.)*

Preparing the DAG-MCNP5 Input File
''''''''''''''''''''''''''''''''''''

The DAG-MCNP5 input file contains only the data cards section of a
standard MCNP5 input file.  There are no cell or surface cards
included in the input file.

In addition to many other MCNP5 data cards, it is important to define
the materials that have been assigned in step 2.D.i.a above and any
tally modifiers, as desired, for the tallies defined in step 2.D.i.a
above.

A new data card has been added to DAG-MCNP5 to define parameters for
the DAGMC geometry capability.  These parameters are described in
:ref:`additional_parameters`.
::
    Form: dagmc  keyword1=value   keyword2=value
           check_src_cell: behavior of CEL variable in SDEF card
                           on  [default] standard interpretation for 
                                         CEL variable: source rejection
                           off           no cell rejection - assume that 
                                         sampled position is in cell CEL
        overlap_thickness: allows particle tracking through small overlaps
                           {real} [default=0.0]
                   usecad: toggle usage of solid model geometry
                           off [default] ray-tracing limited to facets
                           on            ray-tracing performed on solid model 
                                         geometry surfaces
                distlimit: toggle usage of flight distance sampled from 
                           physics to accelerate ray-tracing search
                           off [default] do not use physics flight distance
                           on            do use physics flight distance


Running DAG-MCNP5
'''''''''''''''''''

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
