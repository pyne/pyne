# CMAKE build script for KDE boundary correction tool

PROJECT(KDEBoundaryCorrection)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Add MOAB_DIR to list of include directories
INCLUDE_DIRECTORIES("${MOAB_DIR}/include")

# Define the MOAB libraries that are needed for this project
SET(MOAB_LIBRARIES
    ${MOAB_DIR}/lib/libMOAB.so 
    ${MOAB_DIR}/lib/libdagmc.so
)

# Compile boundary correction script executable
ADD_EXECUTABLE(boundary.exe ../boundary.cpp)
TARGET_LINK_LIBRARIES(boundary.exe ${MOAB_LIBRARIES})
