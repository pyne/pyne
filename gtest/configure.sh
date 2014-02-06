#! /bin/bash

EXTRA_ARGS=$@

rm -rf CMakeCache.txt

cmake \
$EXTRA_ARGS \
../gtest-1.7.0
