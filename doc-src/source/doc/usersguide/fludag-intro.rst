___________________
FluDAG Intro
___________________
FluDAG is the combination of DAGMC with FLUKA, a Monte Carlo code created in Italy and used at CERN and elsewhere.
FLUKA provides a way of initiating calls to user-written functions in two ways

a.  through key words and records in its input format 
b.  by passing an integer flag to the main FLUKA routine.

Causing FLUKA to use DAGMC geometry is a matter of method b.  When the FLUKA main routine is called 
with an integer flag, a preset  collection of routines defined internally to FLUKA are not called.  
Instead their counterparts, as defined in external libraries that can be created and linked
in by the user, are called.


FluDAG User Features
====================
+-----------------+---------+--------------+------------+----------+
|                 | Fluka   |   Planned    |  Progress  |  Done    |
+-----------------+---------+--------------+------------+----------+
| Particles       |   X     |              |     X      |          |
+-----------------+---------+--------------+------------+----------+
|   neutrons      |   X     |      X       |            |          |
+-----------------+---------+--------------+------------+----------+
|   photons       |   X     |      X       |            |          |
+-----------------+---------+--------------+------------+----------+
|   charged       |   X     |              |     X      |          |
+-----------------+---------+--------------+------------+----------+
|   neutrino      |   X     |              |            |          |
+-----------------+---------+--------------+------------+----------+
| Trackng         |   X     |              |     X      |          |
+-----------------+---------+--------------+------------+----------+
|   calc normal   |   X     |      X       |            |          |
+-----------------+---------+--------------+------------+----------+
|   raytrace      |   X     |      X       |            |          |
+-----------------+---------+--------------+------------+----------+
|   history       |   X     |      X       |            |          |
+-----------------+---------+--------------+------------+----------+
| Input           |   X     |              |     X      |          |
+-----------------+---------+--------------+------------+----------+
|    .sat file    |         |              |     X      |          |
+-----------------+---------+--------------+------------+----------+
|   mat assign    |   X     |              |     X      |          |
+-----------------+---------+--------------+------------+----------+
|   low mat       |   X     |      X       |            |          |
+-----------------+---------+--------------+------------+----------+
