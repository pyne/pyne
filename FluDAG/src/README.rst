CMake Notes
===========
CMakeLists.txt is set up to create targets named "readvol" and "mainfludag.  
Three localization variables are needed.

* HOME      - the path to your home directory
* FLUPRO    - the path to the FLUKA directory
* MOAB_HOME - the path to the MOAB directory

Any or all of these may be defined in the shell initialization script, e.g.
.bashrc, or they can be defined with the call to cmake, as in

>cmake -DMOAB_HOME=$HOME/path/to/MOAB ../src

If _at least one_ of the three variables is not defined, the file 
local.cmake is INCLUDE'd.  

Sample local.cmake
------------------
# Locally defined vars
set(HOME /path/to/home_directory)
set(FLUPRO ${HOME}/path/to/FLUKA)
set(MOAB_HOME ${HOME}/path/to/MOAB)



To use cmake, go to an empty directory outside of src 
> cd ..
> mkdir bin
> cd bin
> cmake ../src

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

