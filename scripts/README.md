Contributed Scripts for Workflow Improvements
================================================

This directory contains a number of scripts that have been contributed
to help facilitate the workflow with DAGMC.

Automated Imprint/Merge & Graveyard
-------------------------------------

The script `finish_dagmc_geom.bash` also relies on the Cubit journal
file `finish_dagmc_geom.jou`, to perform the following steps
automatically on an ACIS file.  For a given geometry file, say
`geom.sat`, this do the following steps:

* imprint & merge all volumes/surfaces
* export an ACIS file without graveyard for visualization: `geom_ng.sat`
* add a graveyard volume with an inner surface that is a brick 10%
  larger than the geometry bounding box and an outer surface 1% larger
  than the inner surface
* assign that volume to the `graveyard` group
* export an ACIS file with graveyard for transport: `geom_g.sat`
* convert the non-graveyard geometry to STL: `geom.stl`
* convert the graveyard geometry to H5M using default settings: `geom.h5m`



