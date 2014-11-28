FluDAG Testing
--------------
--------------

Testing for the FluDAG implementation consists of two major components:

1. Interface testing provides test cases for individual geometry
   queries, as they pass through the interface defined for FLUKA in
   the FluGG interface.

2. Integrated testing provides test cases for full physics analysis
   using the FluDAG interface.

List of Tests
_______________
issuetest
nrml

To create and run the interface tests you need the google test library installed
in your system.  The cmake configuration will find it.

To create the tests follow the instructions below.

a.  In the current (testing) directory make a build subdirectory and move to it:

testing> mkdir build
testing> cd build

Note: the name of the build subdirectory is arbitrary.

b.  In the build directory call cmake, referencing the iface directory:

build> cmake ../iface DTOPDIR=directory/containing/DAGMC/FluDAG -DMOAB_HOME=path/to/MOAB 

Note 1: if you don't have $FLUPRO as an environment variable (normally it is 
        defined as part of the FLUKA install) you must add -DFLUPRO=path/to/fluka-libs
Note 2: The variables defined with -D switches can also be defined in a file
        in the iface directory named 'local.cmake'

d.  Several files and subdirectories will be created in the build directory.
You can then build the tests by typing the following command where 'testname' is one of the 
names in the List of Tests above

build> make 'testname'

e.  Run the test by typing
build> ./'testname'

All tests should pass.

Google Test Library
____________________
If you or the system you are working on already has Google Test installed and libgtest.a built, 
you can skip this step.  The CMakeLists.txt includes a directive that will find it.

Otherwise, you can build libgtest from the this repository.

Build libgtest if necessary
------------------
Note:  this step will become unnecessary when the gtest lib is automatically built by the cmake facility.

1. In FluDag/testing/gtest create a lib directory and move to it 
gtest> mkdir lib
gtest> cd lib

2.  Create the Makefile with the gtest distribution's CMakeLists.txt (currently version 1.6.0)
lib> cmake ../gtest-1.6.0

3.  Build the library
lib> make gtest

When done you will have FluDAG/testing/gtest/lib/libgtest.a.

Notes on Localization 
_____________________
Regarding the cmake file that controls the build three localization variables are required.

*  If $FLUPRO is defined as an environment variables, it will be read and included 
in the CMakeLists.txt file.  When installing FLUKA, it is typical for a $FLUPRO 
environment variable to be defined.

*  $MOAB_HOME is unlikely to be defined as an environment variable, and CMakeLists.txt
expects the user to either enter a -DMOAB_HOME="..." definition, or define it in 
local.cmake, which is INCLUDE'd by CMakeLists.txt.

*  $TOPDIR is the directory where the DAGMC repository is installed.  It likewise can be 
defined on the command line when calling cmake, or its definition should be placed in local.cmake.

To summarize, three variables are locally defined for cmake.
 - $FLUPRO    - environment variable, -Ddefine, or local.cmake
 - $MOAB_HOME - -Ddefine or local.cmake
 - $TOPDIR    - -Ddefine or local.cmake

Debugging
_________

>cmake path/to/CMakeLists.txt -DCMAKE_BUILD_TYPE=Debug
