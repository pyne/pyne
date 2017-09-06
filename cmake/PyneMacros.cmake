# set platform preprocessor macro
macro(pyne_set_platform)
  set(PYNE_PLATFORM "__${CMAKE_SYSTEM_NAME}__")
  if(APPLE)
    set(PYNE_PLATFORM "__APPLE__")
  elseif(WIN32)
    if(MSVC)
        set(PYNE_PLATFORM "__WIN_MSVC__")
    elseif(CMAKE_COMPILER_IS_GNUC OR CMAKE_COMPILER_IS_GNUCXX)
        set(PYNE_PLATFORM "__WIN_GNUC__")
    endif(MSVC)
  else(APPLE)
    set(PYNE_PLATFORM "__LINUX__")
  endif(APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D${PYNE_PLATFORM}")
  message("-- Pyne platform defined as: ${PYNE_PLATFORM}")
endmacro()

# Fortran settings
# FFLAGS depend on the compiler
macro(pyne_setup_fortran)
  # languages
  enable_language(Fortran)

  # Augment the Fortran implicit link libraries
  message(STATUS "CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES: "
          ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})
  if (APPLE)
    message(STATUS "CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES Before Fix: "
            ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
    set(LIBGCC_S)
    # The previous method found the gcc_s library by version,
    # find_library(LIBGCC_S_PATH gcc_s.${gcc_s_ver}
    # but this always found the library regardless of version, and then
    # the name as set was different from what was found.  This new way
    # ensures that the name is consistent with what was found.
    set(LIBGCC_S_PATH)
    foreach (gcc_s_ver 3 2 1)
      find_library(LIBGCC_S_PATH gcc_s.${gcc_s_ver}
        PATHS ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}
        NODEFAULTPATH
      )
      if (LIBGCC_S_PATH)
        break()
      endif()
    endforeach()
    message(STATUS "LIBGCC_S_PATH = ${LIBGCC_S_PATH}.")
    if(LIBGCC_S_PATH)
      get_filename_component(LIBGCC_S ${LIBGCC_S_PATH} NAME)
      # Pull off leading lib and trailing .dylib.
      string(REGEX REPLACE "^lib" "" LIBGCC_S ${LIBGCC_S})
      string(REGEX REPLACE "\\.dylib$" "" LIBGCC_S ${LIBGCC_S})
    endif()
    message(STATUS "LIBGCC_S = ${LIBGCC_S}.")
    if(LIBGCC_S)
      set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES
          ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} ${LIBGCC_S})
    endif()
  endif()
  message(STATUS "CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES: "
          ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})

  get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

  if(Fortran_COMPILER_NAME MATCHES "gfortran.*")
    # gfortran
    set(CMAKE_Fortran_FLAGS_RELEASE
        "-funroll-all-loops -c -fpic -fdefault-real-8 -fdefault-double-8")
    set(CMAKE_Fortran_FLAGS_DEBUG
        "-c -fpic -fdefault-real-8 -fdefault-double-8")
  elseif(Fortran_COMPILER_NAME MATCHES "ifort.*")
    # ifort (untested)
    set(CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O2 -r8")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g -r8")
  elseif (Fortran_COMPILER_NAME MATCHES "g77")
    # g77
    set(CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O2 -m32")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
  else(Fortran_COMPILER_NAME MATCHES "gfortran.*")
    message("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message("Fortran compiler: " ${Fortran_COMPILER_NAME})
    message ("No optimized Fortran compiler flags are known, we just try -fpic...")
    set(CMAKE_Fortran_FLAGS_RELEASE "-fpic")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-fpic")
  endif(Fortran_COMPILER_NAME MATCHES "gfortran.*")
endmacro()


#  add lib to pyne list
macro(add_lib_to_pyne _name _source)
  # add the library
  add_library(${_name} ${_source})
  # add it to the list of pyne libraries
  set(PYNE_LIBRARIES ${PYNE_LIBRARIES} ${_name})
endmacro()


# Print pyne logo
macro(pyne_print_logo)
  set(cat_prog cat)
  if(WIN32)
    set(cat_prog type)
  endif(WIN32)
  execute_process(COMMAND ${cat_prog}
                  ${PROJECT_SOURCE_DIR}/cmake/logo.txt
                  OUTPUT_VARIABLE variable)
  message("${variable}")
endmacro()


# determine if spatial solver module should be built
macro(pyne_set_build_spatial_solver)
  IF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    IF(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.6" AND
       NOT APPLE AND NOT "$ENV{CONDA_BUILD}")
      SET(BUILD_SPATIAL_SOLVER true)
    ELSE()
      SET(BUILD_SPATIAL_SOLVER false)
    ENDIF()
  ENDIF()
endmacro()


# set build type
macro(pyne_set_build_type)
  # Default to release build type
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "The build type" FORCE)
  endif(NOT CMAKE_BUILD_TYPE)

  # quiets fortify_source warnings when not compiling with optimizations
  # in linux distros where compilers were compiled with fortify_source enabled by
  # default (e.g. Arch linux).
  STRING(TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE)
  IF(NOT ${BUILD_TYPE} STREQUAL "release")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0")
  ENDIF()
  MESSAGE("-- Build type: ${CMAKE_BUILD_TYPE}")
endmacro()


# Setup the RPATH correctly
macro(pyne_configure_rpath)
  # use, i.e. don't skip the full RPATH for the build tree
  SET(CMAKE_SKIP_BUILD_RPATH FALSE)

  # when building, don't use the install RPATH already
  # (but later on when installing)
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # the RPATH to be used when installing, but only if it's not a system directory
  LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
       "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
  IF("${isSystemDir}" STREQUAL "-1")
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    GET_FILENAME_COMPONENT(cxxCompilerRoot ${CMAKE_CXX_COMPILER} DIRECTORY)
    GET_FILENAME_COMPONENT(cxxCompilerRoot ${cxxCompilerRoot} DIRECTORY)
    IF(NOT "${CMAKE_INSTALL_RPATH}" STREQUAL "${cxxCompilerRoot}")
      SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}:${cxxCompilerRoot}/lib")
    ENDIF (NOT "${CMAKE_INSTALL_RPATH}" STREQUAL "${cxxCompilerRoot}")
  ENDIF("${isSystemDir}" STREQUAL "-1")
  MESSAGE("-- CMAKE_INSTALL_RPATH: ${CMAKE_INSTALL_RPATH}")
endmacro()

