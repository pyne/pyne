#! /bin/bash

EXTRA_ARGS=$@

rm -rf CMakeCache.txt

cmake \
-D DAGMC_MESHTALLY_SOURCE=$HOME/DAGMC/MCNP5/dagmc/ \
-D MOAB_HOME=$HOME/dagmc_bld/MOAB \
$EXTRA_ARGS \
..
