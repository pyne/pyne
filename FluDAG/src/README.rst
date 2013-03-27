CMake Notes
===========
CMakeLists.txt is set up to create targets named "readvol" and "mainfludag.  
You will need a local.cmake file in the same directory as the src directory.
local.cmake should define the following three variables:

# Locally defined vars
set(HOME /filespace/example/...)
set(FLUPRO ${HOME}/path/to/FLUKA)
set(MOAB_HOME ${HOME}/path/to/MOAB)

CMakeLists.txt has an INCLUDE directive that loads local.cmake.

To use cmake, go to an empty directory outside of src (e.g. the bin that is 
parallel to src) and type, in the case of bin,
>cmake ../src

This will create a makefile that will allow you to type
>make readvol

-or-

>make mainfludag

When complete, there will be a readvol or mainfludag target in the bin directory, 
along with many files created by the cmake system.

Note 1:  The cmake system creates Makefiles.  This is why the old make system
'Makefile' has been renamed to 'makefile'.

Note 2:  cmake doesn't care where you are when you run it, you just need to refer 
to the directory that contains the top level CMakeLists.txt.

Make Notes 
==========
NOTE:  MAKE IS NOT CURRENTLY IN USE
-------------------------------------

Makefile should be modified for local conditions.
There are three ways to do this:
1.  Export definitions in your startup script, i.e.
    via bash:
    export MOAB=${HOME}/dagmc_bld/MOAB
    export FLUDAG=${HOME}/DAGMC/FluDAG

2.  Modify makefile.local, which is include'd by
    Makefile, to define local MOAB and FLUDAG
    locations.  This will overwrite the variables
    in the bash script.

3.  Modify the Makefile directory directly by defining 
    local MOAB and FLUDAG locations.   

Workflow in src to make and run mainfludag: 
> make clean 
> make obj
> make
If it doesn't exist
> mkdir myrun
> cd myrun
> cp ../input/* .
> cp ../mainfludag .
> source runme
This creates a fluka_xxx subdirectory that will be deleted if there are no errors.

