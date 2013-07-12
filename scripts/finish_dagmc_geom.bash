#!/bin/bash
#
# This script will accomplish the following tasks for an ACIS geometry "geom.sat"
#   * imprint and merge the geometry
#   * export a copy of the geometry with no graveyard as "geom_ng.sat"
#   * add a graveyard that is 10% larger than the bounding box of the geometry
#   * export the graveyard geometry as "geom_g.sat"
#   * convert the no graveyard geometry to an STL file for visualization
#   * convert the graveyard geometry to an H5M file for transport

fname_input=$1
bname=`basename $1 .sat`
fname_output=${bname}_g.sat
fname_no_graveyard=${bname}_ng.sat
fname_stl=${bname}.stl
fname_h5m=${bname}.h5m

cubit -nographics -batch -noecho -nojournal -information off finish_dagmc_geom.jou fname_input=\"$fname_input\" \
         fname_output=\"$fname_output\" fname_no_graveyard=\"$fname_no_graveyard\"

./mbconvert $fname_no_graveyard $fname_stl
./dagmc_preproc $fname_output -o $fname_h5m