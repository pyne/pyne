#! /bin/bash

EXTRA_ARGS=$@

rm -rf CMakeCache.txt

cmake \
-D MOAB_HOME=$HOME/dagmc_bld/MOAB \
$EXTRA_ARGS \
../src
