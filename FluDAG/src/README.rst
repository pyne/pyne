CMake Notes
===========
CMakeLists.txt is set up to create the target named "mainfludag.  
Three localization variables are needed.

* HOME      - the path to your home directory
* FLUPRO    - the path to the FLUKA directory
* MOAB_HOME - the path to the MOAB directory

To compile fludag, 
- create a sister directory to the current src directory and move to it.
> mkdir ../bld
> cd ../bld
> cp ../configure.sh .

- modify the local copy of configure.sh for your environment and execute it.
  This command sets some local variables and runs cmake on the src directory.
> ./configure.sh

- Create the fludag executable, "mainfludag":
> make

-or-

> make mainfludag

When complete, there will be a mainfludag target in the bld directory, 
along with possible other targets and many files created by the cmake system.

Note:  cmake doesn't care where you are when you run it, you just need to refer 
to the directory that contains the top level CMakeLists.txt.

Workflow in src to run mainfludag: 
> mkdir myrun
> cd myrun
> cp ../input/* .
> cp ../mainfludag .
> source runme
This creates a fluka_xxx subdirectory that will be deleted if there are no errors.

ctest Notes
===========
A google test framework has been added in ./test (DAGMC/FluDAG/src/test)
See the README.rst in ./test.
