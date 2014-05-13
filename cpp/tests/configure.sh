#! /bin/bash

EXTRA_ARGS=$@

rm -rf CMakeCache.txt

cmake \
-D HDF5_HOME=$HOME/dagmc_bld/HDF5 \
-D PYNE_HOME=$HOME/.local/lib/python2.7/site-packages/pyne \
$EXTRA_ARGS \
..
