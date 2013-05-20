#! /bin/bash

EXTRA_ARGS=$@

rm -rf CMakeCache.txt

cmake \
-D DAGMC_MESHTALLY_SOURCE=$HOME/research/alt_tallies/DAGMC/MCNP5/dagmc \
$EXTRA_ARGS \
..
