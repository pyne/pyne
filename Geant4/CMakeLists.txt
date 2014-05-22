######################################################################
#
# CMakeLists.txt for dagmc Geant4
#
# $Id: CMakeLists.txt 6346 2014-01-27 03:33:21Z cary $
#
######################################################################

# Add in directories in layering order, lowest to highest
set(SCI_INDENT "${SCI_INDENT}  ")
message(STATUS "${SCI_INDENT}Adding dagsolid subdir.")
add_subdirectory(dagsolid)

