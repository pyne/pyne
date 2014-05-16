PROJECT (PyneDag)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

IF (NOT DEFINED HOME) 
	set (HOME "$ENV{HOME}")
ENDIF (NOT DEFINED HOME) 

SET(hdf5_include ${HDF5_HOME}/include)
SET(pyne_include ${PYNE_HOME}/include)

SET(hdf5_lib_dir ${HDF5_HOME}/lib)
SET(pyne_lib_dir ${PYNE_HOME}/lib)
SET(libs
	${hdf5_lib_dir}/libhdf5.so
	${pyne_lib_dir}/libpyne.so
	${pyne_lib_dir}/libpyne_material.so
	${pyne_lib_dir}/libpyne_nucname.so
	${pyne_lib_dir}/libpyne_data.so
)


# Where to look for includes
INCLUDE_DIRECTORIES(${hdf5_include} ${pyne_include})

ADD_EXECUTABLE(test_material test_createMaterial.cpp)
TARGET_LINK_LIBRARIES(test_material ${libs})

# Create the executable that can call flukam_
ADD_EXECUTABLE(pyne_mat pyne_mat.cpp)
TARGET_LINK_LIBRARIES(pyne_mat ${libs})
