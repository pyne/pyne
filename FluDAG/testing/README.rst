FluDAG Testing
--------------
--------------

Testing for the FluDAG implementation consists of two major components:

1. Interface testing provides test cases for individual geometry
   queries, as they pass through the interface defined for FLUKA in
   the FluGG interface.

2. Integrated testing provides test cases for full physics analysis
   using the FluDAG interface.

To create and run the interface tests you need the google test library,
do the stepse listed below.

a.  Set up the localization variables, see Localization below

b.  In the current (testing) directory make a build subdirectory and move to it:
testing> mkdir build
testing> cd build

c.  In the build directory call cmake, referencing the iface directory:
build> cmake ../iface

d.  Several files and subdirectories will be created in the build directory.
Build the tests by typing
build> make issuetest

e.  Run the test by typing
build> ./issuetest

All tests should pass.

Google Test Library
____________________
If you or the system you are working on already has Google Test installed and libgtest.a built, 
you can skip this step.  The CMakeLists.txt includes a directive that will find it.

Otherwise, you can build libgtest from the this repository.

Build libgtest
------------------
1. In FluDag/testing/gtest create a lib directory and move to it 
gtest> mkdir lib
gtest> cd lib

2.  Create the Makefile with the gtest distribution's CMakeLists.txt (currently version 1.6.0)
lib> cmake ../gtest-1.6.0

3.  Build the library
lib> make gtest

When done you will have FluDAG/testing/gtest/lib/libgtest.a.

Localization
_______________
Regarding the cmake file that controls the build, as with the source build
three localization variables are required.

*  If $FLUPRO is defined as an environment variables, it will be read and included 
in the CMakeLists.txt file.  When installing FLUKA, it is typical for a $FLUPRO 
environment variable to be defined.

*  $MOAB_HOME is unlikely to be defined as an environment variable, and CMakeLists.txt
expects the user to either enter a -DMOAB_HOME="..." definition, or define it in 
local.cmake, which is INCLUDE'd by CMakeLists.txt.

*  $TOPDIR likewise should be defined on the command line when calling cmake,
or its definition should be placed in local.cmake

To summarize, three variables are locally defined for cmake.
 - $FLUPRO    - environment variable, -Ddefine, or local.cmake
 - $MOAB_HOME - -Ddefine or local.cmake
 - $TOPDIR    - -Ddefine or local.cmake
