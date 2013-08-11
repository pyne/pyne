Understanding Ray History State in DAG-MCNP5
============================================

Documenting Calls to the track() Method in MCNP5
---------------------------------------------------

The method used in MCNP5 to find the distance to the next boundary is
track().  DAGMC completely circumvents any call to track() and uses
its own version.

As a first step to understanding the generic set of use cases for the
ray-to-boundary method, and in particular, the different states of the
physics code during these uses cases, this document will itemize these
cases for MCNP5.

1. (hstory.F90): non-forced collision neutral particle transport
2. (hstory.F90): forced collision neutral particle transport
3. (tally.F90): tracking to tally segmenting surface (not relevant for DAGMC)
4. (transm.F90): optical thickness of material for point detector/dxtran
5. (electr.F90): electron substep transport  (also in the updated electron_history.F90)

Mapping the possible Ray History States
-----------------------------------------

For each new particle history, the ray history is reset to be empty.

The only mechanism by which a ray history can grow is during ray firing:
MBErrorCode result = DAG->ray_fire(vol, point, dir, next_surf, next_surf_dist, &history, 
                                   (use_dist_limit ? dist_limit : 0 ) )

Theere are 2 circumstances that reset a ray history:
* particle history termination
* a change in direction other than reflection
    * note that during reflection, the previous direction is changed
      to make it look like a boundary crossing/streaming event
      (dagmc_surf_reflection)

During reflection, we reset the history to the last interaction,
erasing the tail except for the last facet.

It is also possible to rollback the history removing the most recent interaction.

Other History State Variables
-------------------------------

* visited_surface (bool)
    * turned on when a surface interaction is handled either in surface crossing or reflection
    * turned off immediately following the next ray fire
* last_uvw (unit direction vector)
    * set during reflection to spoof a surface crossing
    * set following every ray fire to record direction of last ray_fire
    * note: not reset for each history because is not used until after
      the first ray_fire in each history
* last_nps (int)
    * 

