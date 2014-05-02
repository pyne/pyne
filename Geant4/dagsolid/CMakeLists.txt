PROJECT (DagSolid)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.4)

# disable threading
add_definitions(-DGTEST_HAS_PTHREAD=0)

# Run cmake on the CMakeLists.txt in these subdirectories
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(test)

