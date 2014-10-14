PROJECT (FluDAG)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.4)

add_definitions(-DGTEST_HAS_PTHREAD=0)

# Run cmake on the CMakeLists.txt in these subdirectories
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(tests)

# ToDo:  Maybe Add a POST_BUILD custom command that moves the executables
