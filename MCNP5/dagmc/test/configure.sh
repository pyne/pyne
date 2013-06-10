#! /bin/bash

EXTRA_ARGS=$@

rm -rf CMakeCache.txt

cmake \
-D DAGMC_MESHTALLY_SOURCE=$HOME/DAGMC/MCNP5/dagmc/test/ \
$EXTRA_ARGS \
..
