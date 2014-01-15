#! /bin/bash

EXTRA_ARGS=$@

rm -rf CMakeCache.txt

cmake \
-D DAGMC_FLUDAG_SOURCE=$HOME/DAGMC/FluDAG/src/ \
-D MOAB_HOME=$HOME/dagmc_bld/MOAB   \
-D GTEST_HOME=$HOME/DAGMC/gtest \
$EXTRA_ARGS \
..
