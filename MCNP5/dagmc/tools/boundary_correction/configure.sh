#! /bin/bash

EXTRA_ARGS=$@

rm -rf CMakeCache.txt

cmake \
-D MOAB_DIR=$HOME/opt/moab/ \
$EXTRA_ARGS \
..
