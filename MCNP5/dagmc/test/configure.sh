#! /bin/bash

EXTRA_ARGS=$@

rm -rf CMakeCache.txt

cmake \
-D DAGMC_MESHTALLY_SOURCE=$HOME/research/alt_tallies/DAGMC/MCNP5/dagmc/ \
-D MOAB_HOME=$HOME/opt/moab-4.6.0 \
$EXTRA_ARGS \
..
