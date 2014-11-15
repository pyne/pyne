#! /bin/bash

EXTRA_ARGS=$@

rm -rf CMakeCache.txt

cmake \
-D DAGMC_TALLY_SOURCE=$HOME/DAGMC/MCNP5/dagmc/ \
-D MOAB_HOME=$HOME/dagmc_bld/MOAB \
-D GTEST_HOME=$HOME/DAGMC/gtest \
$EXTRA_ARGS \
..
