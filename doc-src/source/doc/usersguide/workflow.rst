Using DAGMC - The DAGMC Workflow
-----------------------

The following steps are necessary to perform an analysis with most
DAGMC versions of Monte Carlo radiation transport codes. (see below
for more details on each step)

1. create and/or prepare a solid model using the CAD/solid-modeling tool of your choice
2. pre-processing of the solid model geometry in CUBIT
    * assign materials and densities
    * define boundary conditions
    * imprint & merge geometry
    * export that model in the ACIS format
3.  prepare native MC code input file
    * understanding additional ``dagmc`` parameters
4. Run MC code, possibly with additional command line parameters

Preparing Solid Models
++++++++++++++++++++++

In theory, solid models can be prepared in any modeling software
system (e.g. SolidWorks, Pro/E, Catia, etc).  What is most important
about the choice of solid modeling system is the ability to export to
a format that can be imported by CUBIT, in particular:

* ACIS (*.sat)
* STEP (*.stp, *.STEP, etc)

There are a number of concepts to consider while preparing a solid
model; however, the most important considerations are small gaps and
overlaps that might exist in the model. These gaps and overlaps can
lead to rapid failure when running a DAGMC-based analysis. The
following steps are provided to help make a more robust model *before*
running your DAGMC-based analysis.

Be aware: obtaining a robust model may be an iterative and time
consuming process. In some cases, the validity of the model will
require running a DAGMC-based analysis and assessing whether or not
the model yielded expected results or a small enough number of lost
particles. If the results did not meet expectations, changes to the
model may be in order.

Knowing the model
"""""""""""""""""

The first consideration to address is where the solid model originated
and for what purpose. In many instances, models constructed for
manufacturing purposes will have tolerances that are undesirable for
particle transport applications. For example, a gap might exist
between fuel pellets and the cladding wall for a PWR fuel rod. While
this is perfectly acceptable for an individual manufacturing the rod,
the gap could potentially cause present problems in a DAGMC-based
analysis, depending on how it is modeled.

Knowing who created the model and to what purpose provides a starting
point for preparing the model. If it was made with particle transport
in mind, then very little work may be needed; but as with the example
above, some models may require changes to accommodate the needs of a
DAGMC-based analysis.

Identifying weaknesses in the model
"""""""""""""""""""""""""""""""""""""

When assessing a model that is to be used for particle transport two
primary concerns must addressed. These concerns are:

    * Gaps
    * Overlaps

Gaps occur when the surfaces of two volumes/parts that should be in
contact are set apart from each instead of having coincident
surfaces. The size of the gap is generally unimportant, for most solid
modeling programs, a gap is a gap. The desired result is to have all
surfaces of volumes/parts to be coincident. If coincidence is not
achieved, particles may become lost when entering the region between
the surfaces.

Overlaps are found where two or more volumes/parts encroach upon the
same space. As with gaps, the magnitude of the overlapping volume is
usually unimportant.  When a particle enters a region of overlap, it
may not correctly determine which volume/part it resides in. If this
occurs, the particle may become lost.

Identifying gaps and overlaps may be difficult and time consuming;
however, some 3D modeling programs like SolidWorks have built in tools
to identify these occurrences. Rely on the modeling program to
identify these errors (the gaps and overlaps) and use the steps in the
next section to change, reduce and remove their effect on the model.

Modifying your model
"""""""""""""""""""""""

Once the gaps and overlaps in the model have been identified, the
three following methods may be used to change, reduce and remove their
effect on the model.

* Create "voided" geometries
* Modify volume/part dimensions
* Remove superfluous details

Each method is discussed in detail below:

As with the fuel rod example mentioned above, some geometries that are
'gaps' are also important. Instead of removing the gap entirely (by
changing the dimensions of the cladding or the fuel to force
coincidence), a new volume/part could be modeled that coincided with
the outer diameter of the fuel AND the inner diameter of the
cladding. Now a "voided" geometry occupies the previously unaccounted
for region. By specifying these "voided" geometries in a DAGMC-based
analysis, the physical importance of the region can be retained while
accomodating the requirement of having coincident surfaces.

Another method to resolve gaps and overlaps is to simply change the
dimensions of the volume/part (eg: making a dimension several cm
bigger or smaller to ensure coincidence surfaces). In many instances
this method could compromise the physics of the solution and is then
undesirable. However, in other instances, this solution is very
logical. One particularly significant example is if different volumes
were modeled with different unit systems. For example, one volume/part
might have been model in [in] while its neighbor was modeled in [cm];
while the surfaces may be nearly coincidence, rounding errors might
prevent coincidence from occurring. A simple change to one dimension
may hardly change the volume/part's characteristics yet result in
coincidence.

Finally, superfluous details may prevent a volume/part from coinciding
with its neighbors properly. A potential solution is to simply remove
the superfluous detail to simplfy the model and ensure the desired
surfaces are coincident. Some volumes/parts will inherently hurt the
model's effectiveness either due to its complex features or small
dimensions. A volume/part's effect on the model cannot truly be
assessed until a DAGMC-based analysis is run. This final method is
usually implemented in an attempt to reduce the number of lost particles
while maintaining the most important characteristics of the system.

*Note: Of all steps, the removal of superfluous details is the most
 subjective and heavily dependent on the model's intended
 application.*

Assessing your model
""""""""""""""""""""

Lost particles are undesirable; lost particles usually indicate
weaknesses and failures within the geometry. While the goal of the
DAGMC project is to guarantee that there will never be lost particles,
they can occur even on robust geometries.  It is up to the
user/analyst to determine what lost particle rate they consider
acceptable.  The UW-Madison group usually considers lost particle
rates that are less than 1/50,000 to be a threshold for most problems.
It is important to understand whether particles are being lost from an
important region of your phase space.

[Insert note on the implicit complement here]

Pre-processing Solid Models using CUBIT
+++++++++++++++++++++++++++++++++++++++++

*Note: For large models, the steps described below can be very tedious
and time consuming.  To accelerate the process, an automated approach
is available for models that have been properly prepared in the native
solid modeling software.  This AutomatedCubitConversion process is
described elsewhere, but reading the information below will provided
the knowledge-base needed to understand the automation process.*

This section focuses on steps that are independent of the MC code used
for analysis.  Additional steps for `DAG-MCNP5 <#S2Di>`_ and
`DAG-Tripoli4 <#S2Dii>`_ may be based on the instructions given here,
but are provided in separate parts of section 2.D below.

Importing the Solid Model
"""""""""""""""""""""""""""

The first step in CUBIT is to import the generated solid
model. Depending on the complexity of the model, this step can take
several seconds up to a half an hour. As an initial user, it is
recommend to start with simple models and geometries to obtain a
better understanding of CUBIT.

Imprint and Merge
"""""""""""""""""

For a DAGMC-based analysis to work properly, all of the surfaces must
be imprinted and merged.  Imprinting creates a common surface
interface between touching volumes.  Merging then takes the touching
surfaces and makes them into one surface.

To imprint, issue the following command:
::
     imprint body all

Should the imprint be successful, then the next step is to merge the
geometry. Sometimes it may be important to specify a merge tolerance.
To set the tolerance and merge, issue the following commands:
:: 
    merge tol 5e-7
    merge all

This process can be very time consuming. For large models of several
thousand volumes, the imprint and merge steps can take up to three
hours. However, small geometries (on the order of 100 volumes) the
process is rather quick.

.. _grouping-basics:

Grouping Volumes and Surfaces
"""""""""""""""""""""""""""""

A DAGMC-based analysis allows a number of attributes of the geometry
to be defined within the geometry file. These characteristics
generally relate to the physical behavior of the volume, for example
its material definition or boundary conditions.

Before the discussion of specific attributes, the practice of
"grouping" needs to be explained. A group is essentially a collection
of volumes or surfaces that share a common attribute; the practical
usage of "grouping" will be explained in the next section.

The general format for creating/adding volumes to a group is:
::
    group "group.name" add vol/surf ...

For example, to create a group called "moderator" containing volumes
5, 6, 7, and 10, the following command would be used:
::
    group "moderator" add vol 5 to 8 10

Another example, shows that groups don't have to be just contain
volumes, but can contain surfaces too. Below the group
"shield.boundary" is created with surfaces 16 and 37:
::
    group "shield.boundary" add surf 16 37

Due to the importance of using the ``group`` command reading the CUBIT
manual section on its full usage is highly recommended.

Finishing Up and Final Notes
""""""""""""""""""""""""""""

Before exporting, it is vital to set attributes on.  This saves the
absolute volume and surface IDs as well as any group specifications.
Failing to do this will result in fatal errors.  Make sure to type the
following:
::
     set attribute on

Finally, export the file as an ACIS file with a .sat extension.  If
you are using a version of CUBIT newer that v12.x, be sure to set the
ACIS geometry level back to version 19:
::
     set geom version 1900

For the remainder of this documentation, the geometry file will be
referred to as "geom.sat". Also, as noted before, the CUBIT conversion
process can be automated as described on the follow webpage:
AutomatedCubitConversion.

.. _additional_parameters:

Additional Input Parameters for DAGMC
++++++++++++++++++++++++++++++++++++++

DAGMC introduces a number of new input variables that control the
behavior of DAGMC when it interacts with your geometry.  This section
describes the conceptual behavior of those parameters and later
sections will describe their specific implementation and syntax within
each code.

Geometry File (basic)
"""""""""""""""""""""

* required
* Default: none

This file contains the geometry that has been created and
pre-processed using the workflow described above.  This file can be
either an ACIS geometry (usually with a ``.sat`` file extension) or a
MOAB facet file (usually with a ``.h5m`` file extension).

Faceting Tolerance (basic)
""""""""""""""""""""""""""

* optional
* Default: 0.001

One of the first steps performed by DAGMC is to generate a faceted
representation of your solid model.  To ensure a faithful
representation of the model, the facets are constrained such that all
points on each facet are within a tolerance of the nearest points on
the exact surface representation.  A smaller tolerance results in a
more faithful representation of the surface at the penalty of
producing more facets.  The user can control the faceting tolerance
using when they invoke their simulation, either on the command line or
in the input file, depending on the MC code being used for the
analysis.  This option only has an effect with the geometry file is a
solid model and not when it is a facet file.

Facet File (basic)
""""""""""""""""""

* optional
* Default: none

For some models, the initial processing can be time consuming.  When
reading a solid model geometry, this option causes a processed file to
be written that can be used on subsequent analyses.  This file will be
processed with the facet tolerance as defined above.  This facet
tolerance cannot be changed when the file is reused.

Overlap Thickness (advanced)
"""""""""""""""""""""""""""""

* optional
* Default: 0.0

Often CAD geometry has small overlaps between adjacent volumes due to
file translation or imprecise modeling. The particle tracking
algorithm can accommodate small overlaps if the maximum thickness of
overlaps is approximately known.

Source Cell Treatment (intermediate)
""""""""""""""""""""""""""""""""""""

* optional
* Default: on (same behavior as native code)

The implementation of this option is specific to the Monte Carlo code
being used.  Please refer to the documentation for your underlying
Monte Carlo code.

Use CAD geometry (advanced)
"""""""""""""""""""""""""""""

* optional
* Default: off

When this option is turned on, the ray-firing process finds the
intersection with the CAD-based solid model itself, and not just with
the faceted representation of that model.  The facet-based ray-firing
solution is used as an initial guess to solve for the point on the
actual CAD surface where the ray-surface intersection takes place.
This option is only available when the DAGMC toolkit has been linked
to the ACIS geometry libraries directly and not when it has been
linked via CUBIT.

Use Distance Limit (experimental)
"""""""""""""""""""""""""""""""""

* optional
* Default: off

This option allows a previously determined distance to the next
collision to be used to accelerate the search for ray-surface
intersections.  Any candidate ray-surface intersection that is at a
distance beyond the known distance to collision will be rejected,
including bounding box tests in the OBB tree search.

