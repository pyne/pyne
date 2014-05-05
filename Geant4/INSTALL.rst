DagSolid Install instructions
===============================

To install DagSolid in the base DagSolid directory, i.e. that containing src, tests etc. Create a new directory, "bld" for example and move to it

    mkdir bld
    cd bld

Now issue a CMake command to build the files required to make DagSolid

    cmake ../. -DMOAB_DIR=<path to moab> -DGEANT_DIR=<path to geant> -DDAGSOLID_DIR=<path to DagSolid>
    
Now build and install
   
    make
    make install

This last step install the DagSolid Library (libdagsolid.so) into the lib directory. If you would like to run test, cd into tests and then run

    cd test
    ./dagsolid_unit_tests

You should also add the DagSolid /lib directory to your LD_LIBRARY_PATH
