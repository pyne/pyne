INCLUDE(DownloadAndExtract)

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

# C++ settings
macro(pyne_setup_cxx)
  set(CMAKE_CXX_STANDARD 11)
endmacro()

macro(pyne_set_asm_platform)
  # first set OS
  if (WIN32)
    set(_plat "win")
  elseif(APPLE)
    set(_plat "apple")
  else()
    set(_plat "linux")
  endif()
  # next set compiler
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(_plat "${_plat}-gnu")
  elseif (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR
    ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang"))
    set(_plat "${_plat}-clang")
  else()
    set(_plat "${_plat}-NOTFOUND")
  endif()
  set(PYNE_ASM_PLATFORM "${_plat}")
endmacro()

# Fortran settings
# FFLAGS depend on the compiler
macro(pyne_setup_fortran)
    enable_language(Fortran)
    if(NOT CMAKE_Fortran_COMPILER)
        message(FATAL_ERROR "No Fortran compiler found")
    endif()
    message(STATUS "Using Fortran compiler: ${CMAKE_Fortran_COMPILER}")
endmacro(pyne_setup_fortran)


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
  SET(BUILD_SPATIAL_SOLVER false)
  IF ( ENABLE_SPATIAL_SOLVERS )
    MESSAGE("-- Checking whether to build spatial solvers")
    MESSAGE("-- -- Checking CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
    IF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      MESSAGE("-- -- -- Checking CMAKE_CXX_COMPILER_VERSION: ${CMAKE_CXX_COMPILER_VERSION}")
      MESSAGE("-- -- -- Checking if APPLE: ${APPLE}")
      IF(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.6" AND
        NOT APPLE )
        SET(BUILD_SPATIAL_SOLVER true)
      ELSE()
        SET(BUILD_SPATIAL_SOLVER false)
      ENDIF()
    ENDIF()
  ENDIF( ENABLE_SPATIAL_SOLVERS)
  MESSAGE("-- Build spatial solvers: ${BUILD_SPATIAL_SOLVER}")
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

macro(pyne_set_fast_compile)
  if(NOT DEFINED PYNE_FAST_COMPILE)
    set(PYNE_FAST_COMPILE TRUE)
  endif()
  message(STATUS "PyNE Fast Compile: ${PYNE_FAST_COMPILE}")
endmacro()

macro(pyne_download_platform)
if(PYNE_FAST_COMPILE)
  # Download bateman solver from PyNE data
  download_platform("http://raw.githubusercontent.com/pyne/data/master" "decay"
                      ".cpp" ".s")

  # Download CRAM solver from PyNE data
  download_platform("http://raw.githubusercontent.com/pyne/data/master" "cram"
                         ".c" ".s")
endif()
endmacro()

# fast compile with assembly, if available.
macro(fast_compile _srcname _gnuflags _clangflags _otherflags)
  get_filename_component(_base "${_srcname}" NAME_WE)  # get the base name, without the extension
  # get the assembly file name
  if (PYNE_ASM_PLATFORM)
    set(_asmname "${_base}-${PYNE_ASM_PLATFORM}.s")
  else()
    set(_asmname "${_base}-NOTFOUND")
  endif()

  # pick the filename to compile, either source or assembly
  if(NOT PYNE_FAST_COMPILE)
    message(STATUS "Not fast compiling ${_srcname} since PyNE fast compile is disabled.")
    set(_filename "${_srcname}")
  elseif(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${_asmname}")
    message(STATUS "Not fast compiling ${_srcname} since ${CMAKE_CURRENT_SOURCE_DIR}/${_asmname} does not exist.")
    set(_filename "${_srcname}")
  else()
    message(STATUS "Compiling ${_srcname} fast from assembly ${_asmname}")
    set(_filename "${_asmname}")
  endif()
    set(PYNE_SRCS "${_filename}" "${PYNE_SRCS}")

  # set some compile flags for the selected file
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set_source_files_properties("${_filename}" PROPERTIES COMPILE_FLAGS "${_clangflags}")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set_source_files_properties("${_filename}" PROPERTIES COMPILE_FLAGS "${_gnuflags}")
  else()
    set_source_files_properties("${_filename}" PROPERTIES COMPILE_FLAGS "${_otherflags}")
  endif()
endmacro()
