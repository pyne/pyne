Introduction to FluDAG
######################

The Direct Accelerated Geometry for Monte Carlo code (DAGMC) is intended to work with any
Monte Carlo code as a faster form of geometry calls.
The Monte Carlo code “FLUKA” has a mechanism for making external geometry calls. This
mechanism has been implemented for calls into the geometry part of the Geant4 Monte Carlo code.
The FLUKA-Geant4 combination is known ad “FLUGG”. Here we follow much of the FLUGG
methodology in order to combine FLUKA with DAGMC into a combined system called “FluDAG”.
The following sections start by mirroring existing FLUGG documentation and then become specific to
FluDAG.

FLUKA Geometry Input and Initialization
***************************************

FLUKA Geometry
==============

  *  Detector elements are called regions, or “unions of zones”
  *  Zones result from the intersection of bodies that are elementary geometric shapes
  *  Each point can belong to one and only one region
  *  Regions can be replicated through rotations and/or translation to create a “lattice”. This facility
     is not used by FluDAG.
  *  Lattices are used only for scoring/bookkeeping

FLUKA Geometry Initialization
=============================

  *  Regions are assigned properties (medium, biasing, etc)
  *  Regions are assigned a number
  *  Any geometry package interfaced with FLUKA must provide a region numeration in order to be
     able to assign all the information necessary to the simulation at each region.
  *  To implement geometry initialization external to FLUKA the function
     JOMIWR(ngeom, lingeo, lougeo, irtrue) must be implemented.

     ◦ For FLUGG, this function is defined in the Wrapper/src directory file WrapInit.cc as
       jomiwr(G4int& nge, const G4int& lin, const G4int& lou, G4int& flukaReg)

     ◦ This function, with the identical signature, must be defined for FluDAG as well.

     ◦ A table describing the arguments and their use for FLUKA, FLUGG, and FluDAG is shown below.


Table 1: JOMIWR Arguments 
-------------------------

+------------------+-------------+-----------------------------+------------------------+
| FLUKA doc name   | FLUGG name  | FLUKA doc                   | FluDAG Use             |
+==================+=============+=============================+========================+           
| ngeo             | nge         | last memory location        |      none              |
|                  |             | used in blank common        |                        | 
+------------------+-------------+-----------------------------+------------------------+
| lingeo           | lin         | geometry input unit         |      none              |
+------------------+-------------+-----------------------------+------------------------+
| lougeo           | lou         | region output unit          |      none              |
+------------------+-------------+-----------------------------+------------------------+
| irtrue           | flukaReg    | number of regions for       | assigned number of     |
|                  |             | this problem                | volumes, from          | 
|                  |             |                             | DAG->num_entities(3)   |
+------------------+-------------+-----------------------------+------------------------+

Source: FLUKA documentation “geomcr.dvi” of May 17, 1999;
FLUGG wrapper file WrapInit.cc; and FluDAG code implementation

FLUKA Geometry Calls
********************
"The second task of geometry is to allow navigation inside the detector model, i.e. the geometrical
propagation of particle tracks in the detector, finding their intersection with the volume boundaries,
etc.”

    *  All variables are in double precision

    *  FLUKA physics always assumes that a step crossing a physical boundary (a boundary between
       different regions) is shortened by the geometry routines at the boundary distance.

    *  During charged particle tracking FLUKA can, and usually does, ask the geometry for steps
       longer than the one actually ‘wished’ (requested?).
       → the geometry should not make the assumption that the next step will start at the end point
         of the previous one, or at least should not rely on this assumption

1  Tracking call to G1WR
=========================
*G1wr calculates the distance traveled in the present zone/region and the number of the next zone/region
to be entered by the particle.*

    1. This routine is called from the FLUKA tracking routine GEOFAR (or GEOMTR for navigation
        in the magnetic field).
    2. Error Conditions handled by the FLUKA code

       *   Error Code IRPRIM = -33: a real error occurred. The starting point lies inside the region
           whose number is designated by IR but the geometry is unable to find which region,
           IRPRIM, will be entered after exiting IR.

       *   Error Code IRPRIM = -3: the starting point does not lie inside IR. This situation can occur
           very frequently during charged particle tracking, e.g.

           ◦ when a particle reaching a boundary is reflected back or

           ◦ by the multiple scattering angle applied after the step, or

           ◦ by the initial deflection of the new step in the new region [1]_.
           The geometry is supposed to be able to detect and flag these situations, particularly when
           occurring ON a boundary, which is the typical case.
           The physics at the interfaces is kept as meaningful as possible, taking into account

           ◦ initial or final step deflections

           ◦ materials of the old and new region

           ◦ grazing angle, etc.

           *Interface effects dominate electron transport and are a fundamental ingredient of the
           response of detectors made by alternating layers.*

    3. Until the subsequent steps are shorter than DSNEAR, FLUKA does not call the geometry at all,
       keeping all the information computed from the previous step.
    4. Charged particles moving in a magnetic field do not follow straight trajectories, so their curved
       path is split in to chords. When the path intersects volume boundaries the arc segment is chosen
       in such a way to minimize the distance between the chord and the real path.
    5. G1WR does boundary rounding, which depends on the distance between the particle and the
       boundary. A particle is taken to be at the right of th eboundary if its direction vector is heading
       right when within the boundary rounding, regardless if the particle is really ‘right’ or not, and
       likewise for the left side.

.. [1] FLUKA has a complex multiple scattering algorithm with two deflections, one at the beginning and one 
       at the end of a step.

1.1 DAG Functions used to implement g1wr()
------------------------------------------

ray_fire
	vol     	= volume to fire the ray at
 	point   	= point from which to start the ray
	dir		= unit direction vector of the ray
	next_surf      	= output:  next surface intersected by the ray; 0 if none found
	next_surf_dist 	= output: distane to the next surface.  0 => undefined, don't use

next_vol 
        Get the volume on the other side of a surface
	surface 	= next_surf from ray_fire
	old_volume 	= a volume on one side of the surface
	new_volume	= output:  the volume on the other side of the surface

Table 2:  G1WR arguments:  FLUKA documentation, FLUGG WrapG1.cc, FluDAG
-----------------------------------------------------------------------

+-------+------------+----------------+--------------------------+---------------------+
|       | FLUKA name | FLUGG name     |       FLUKA doc          |     FluDAG Use      |
+=======+============+================+==========================+=====================+
|       |     X      | G4double& pSx  | Cartesian coordinates    | Passed to ray_fire  |
|       +------------+----------------+                          |                     |
| Input |     Y      | G4double& pSy  |                          |                     |
|       +------------+----------------+                          |                     |
|       |     Z      | G4double& pSz  |                          |                     |
|       +------------+----------------+--------------------------+---------------------+
|       |    WB(3)   | G4double* pV   | Direction cosines vector | Passed to ray_fire  |
|       +------------+----------------+--------------------------+---------------------+
|       |     IR     | G4int& oldReg  | Region number            | Use in              |
|       |            |                |                          | DAG->entity_by_index|
|       |            |                |                          | to get vol entity   |
|       +------------+----------------+--------------------------+---------------------+
|       |  IRLTTC    | const          | Lattice number           | Ignore              |
|       |            | G4int& oldLttc |                          |                     |
|       +------------+----------------+--------------------------+---------------------+
|       |   DIST0    | G4double&      | Step length suggested by | Updated             |
|       |            | propStep       | the physics              |                     |
|       +------------+----------------+--------------------------+---------------------+
|       |    NASC    | G4int& nascFlag| Number of the surface hit| Not used            |
|       |            |                | by the particle at the   |                     |
|       |            |                | preceding step (-1 if the|                     |
|       |            |                | particle changed its     |                     |
|       |            |                | direction or if it is a  |                     |
|       |            |                | new particle)            |                     |
+-------+------------+----------------+--------------------------+---------------------+
|       |     S      | G4double&      | Step length approved     | = next_surf_dist    |
|       |            | retStep        |                          | from ray_fire()     |
|       +------------+----------------+--------------------------+---------------------+
| Output|  IRPRIM    | G4int& newReg  | Region number of the     | Assigned index of   |
|       |            |                | particle after S         | newvol gotten from  |
|       |            |                | - or - error code        | next_vol()          |
|       +------------+----------------+--------------------------+---------------------+
|       |   NASC     | G4int&         | Number of the surface hit| Not used            |
|       |            | nascFlag       | by the particle (1 for   |                     |
|       |            |                | normal tracking, 0       |                     |
|       |            |                | otherwise)               |                     |
|       +------------+----------------+--------------------------+---------------------+
|       |  DSNEAR    | G4double& saf  | Param with value <= to   | Not used            |
|       |            |                | the minimum of the       |                     |
|       |            |                | distances between the    |                     |
|       |            |                | particle and each bdry   |                     |
|       |            |                | (0 if complex).          |                     |
|       |            |                | s -= s*3.0e-9;           |                     |
|       |            |                | saf = s/cm.              |                     |
|       +------------+----------------+--------------------------+---------------------+
|       | IRLTNW     | G4int& newLttc | New lattice number       | Not used            |
+-------+------------+----------------+--------------------------+---------------------+
|Optiona| IRLTGG     | G4int& LttcFlag| Number of crossed        | Not used            |
|l if   |            |                | lattices                 |                     |
|IRLTNW +------------+----------------+--------------------------+---------------------+
|=IRLTTC| SLTCHN(0:20| G4double* sLt  | Array containing the step| Not used            |
|       | 00)        |                | in each crossed lattice  |                     |
|       +------------+----------------+--------------------------+---------------------+
|       | JRLTGG(0:20| G4int* jrLt    | Id of crossed lattices   | Not used            |
|       | 00)        |                |                          |                     |
+-------+------------+----------------+--------------------------+---------------------+

2  LK calls
===========
*Look for region number*

2.1 LKWR
--------

Table 3:  LKWR arguments:  FLUKA documentation, FLUGG WrapLookZ.cc, FluDAG
---------------------------------------------------------------------------
fluka_funcs:lkwr() is called by fluka_funcs:g1wr() to determine the particle volume

+-------+------------+----------------+--------------------------+---------------------+
|       | FLUKA name | FLUGG name     |       FLUKA doc          |     FluDAG Use      |
+=======+============+================+==========================+=====================+
|       |     X      | G4double& pSx  | Cartesian coordinates    | Test if this point  |
|       +------------+----------------+ of the particle          | is in a given       |
| Input |     Y      | G4double& pSy  |                          | volume              |
|       +------------+----------------+                          |                     |
|       |     Z      | G4double& pSz  |                          |                     |
|       +------------+----------------+--------------------------+---------------------+
|       |    WB(3)   | G4double* pV   | Direction cosines vector | Not used            |
|       +------------+----------------+--------------------------+---------------------+
|       |     IR     | G4int& oldReg  | Current region number    | Assigned to newReg  |
|       |            |                | of the particle (can be  | if the particle is  |
|       |            |                | dummy)                   | on a boundary       |
|       +------------+----------------+--------------------------+---------------------+
|       |  IRLTTC    | const          | Current lattice number of| Ignore              |
|       |            | G4int& oldLttc | the particle-can be dummy|                     |
+-------+------------+----------------+--------------------------+---------------------+
| Output|  NMEDG     | G4int& nextReg | Region corresponding to  | Assigned the index  |
|       |            |                | XYZ                      | of the volume where |
|       |            |                |                          | the point is found  |
|       +------------+----------------+--------------------------+---------------------+
|       |  IRPRIM    | G4int&         | Error code: -3 for       | Set in fluka_funcs: |
|       |            | flagErr        | problems, 0 otherwise    | lkwr()              |
|       +------------+----------------+--------------------------+---------------------+
|       | IRLTNW     | G4int& newLttc | New lattice number       | Not used            |
+-------+------------+----------------+--------------------------+---------------------+

DAG Functions used to implement lkwr()
______________________________________

DAG->num_entities(3)
      Get the number of entities in the geometry with dimension = 3, i.e. the number of volumes

DAG->entity_by_index(3,i)
      volume = ith volume entity

DAG->point_in_volume(volume, xyz, is_inside)
      volume = volume to test
      xyz = double vector representing location to test for volume containment
      is_inside = output: 0 if xyz is outside volume; 1 if inside; -1 if on the boundary

Psuedo-code for FluDAG implementation of lkwr()
_______________________________________________

lkwr(xyz, oldRegion, Region, flagErr)
{
    int num_vols = DAG->num_entities(3)
    for (i=1 to i=num_vols)
    {
         volume = DAG->entity_by_index(3,i)
         DAG->point_in_volume(volume, xyz, is_inside)
         if (is_inside == -1)
         {
               Region = oldRegion;
               flaggErr = oldRegion;
               return;
         }
          if (is_inside == 1)
          {
                Region = i;
                flagErr = i // flagErr must be set like this for unknown reasons
                return;
          }
      }
 }

2.2 LKFXWR
----------
This tracking routine is similar in purpose to lkwr() and identical in signature. It is intended to fix
“particular conditions” and is not used by FluDAG.

2.3 LKMGWR
----------
This routine is called from FLUKA’s GEOMAG and localizes a particle for magnetic field tracking. As
with LKFXWR its signature is identical to that of lkwr().

2.4 LKDBWR
----------
This routine is for geometry debugging and is not implemented in FLUGG.

3 NRML 
======
*Compute the normal to a boundary call*

Compute normal unit-vector in global coordinates.
Fluka requires a normal vector exiting from the final position (of the particle) volume, that is: entering
in the volume of the initial position. Geant4 always computes the normal vector exiting from the
volume.

In GetLocalExitNormal() call the volume is the pre-step volume (so the G4 normal vector sign is
opposite of the fluka-required normal). If IsLocalExitNormalValid=false, the normal value is
computed from init-step volume (in this case the sign must change), or from end-step volume (the sign
is the same). The normal vector is always computed on the boundary between volumes, in global
coordinates (to take rotation of parameterised volumes in the hierarchy into consideration).

So: nrmlwr returns an inwards pointing unit normal of the shape for the surface closest to the point
returned by the navigator (last step end-point).


Table 4:  nrml()  arguments:  FLUKA documentation, FLUGG WrapLookZ.cc, FluDAG
------------------------------------------------------------------------------

+-------+------------+----------------+--------------------------+---------------------+
|       | FLUKA name | FLUGG name     |       FLUKA doc          |     FluDAG Use      |
+=======+============+================+==========================+=====================+
|       |     X      | G4double& pSx  | Cartesian coordinates    | Particle location   |
|       +------------+----------------+ of the particle          |                     | 
| Input |     Y      | G4double& pSy  |                          |                     |
|       +------------+----------------+                          |                     |
|       |     Z      | G4double& pSz  |                          |                     |
|       +------------+----------------+--------------------------+---------------------+
|       |    TX      | G4double& pVx  | Direction cosines vector | Not used            |
|       +------------+----------------+                          |                     |
|       |    TY      | G4double& pVy  |                          |                     |
|       +------------+----------------+                          |                     |
|       |    TZ      | G4double& pVz  |                          |                     |
+-------+------------+----------------+--------------------------+---------------------+
| Global|            |                |                          |  next_surf          |
+-------+------------+----------------+--------------------------+---------------------+
| Output|  UN        | G4idouble*     | Normal vector            | Returned            |
+-------+------------+----------------+--------------------------+---------------------+
