#! /bin/bash

EXTRA_ARGS=$@

rm -rf CMakeCache.txt

cmake \
-D MOAB_HOME=$HOME/dagmc_bld/MOAB \
-D DAGMC_FLUDAG_SOURCE=$HOME/DAGMC/FluDAG/src \
$EXTRA_ARGS \
../src
