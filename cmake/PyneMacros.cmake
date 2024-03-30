include(DownloadAndExtract)

# Configure scikit-build
macro(pyne_configure_skbuild)
  if(SKBUILD)
    # Scikit-build installs files to ${SKBUILD_PLATLIB_DIR}/install-directory
    # So, set bin directory to root environment (install prefix/bin)
    set(CMAKE_INSTALL_BINDIR ${SKBUILD_SCRIPTS_DIR})
  endif()
endmacro()


# Configure RPATH
macro(pyne_configure_rpath)
  if(APPLE)
    set(RPATH "@loader_path")
  elseif(UNIX)
    set(RPATH "$ORIGIN")
  else()
    set(RPATH OFF) # Windows
  endif()
endmacro()

# Set platform preprocessor macro
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
  message(STATUS "PyNE platform is defined as: ${PYNE_PLATFORM}")
endmacro()

# Set the platform for the assembler
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
  message(STATUS "PyNE assembler platform is defined as: ${PYNE_ASM_PLATFORM}")
  
  # Enable assembly
  enable_language(ASM)
endmacro()

# Fortran settings
# FFLAGS depend on the compiler
macro(pyne_setup_fortran)

  # Enable Fortran
  enable_language(Fortran)

  # Augment the Fortran implicit link libraries
  message(STATUS "CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES: "
          "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}")
  if (APPLE)
    message(STATUS "CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES Before Fix: "
            "${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")
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
        ${DEPS_LIB_HINTS}
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
          "${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")

  get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

  if(Fortran_COMPILER_NAME MATCHES "gfortran.*")
    # gfortran
    set(CMAKE_Fortran_FLAGS_RELEASE
        "-funroll-all-loops -fpic -fdefault-real-8 -fdefault-double-8")
    set(CMAKE_Fortran_FLAGS_DEBUG
        "-fpic -fdefault-real-8 -fdefault-double-8")
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

  # add -fallow-argument-mismatch to fix build with gfortran 10+
  # https://github.com/pyne/pyne/issues/1416
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
endmacro()

# Print pyne logo
macro(print_pyne_logo)
  set(cat_prog cat)
  if(WIN32)
    set(cat_prog type)
  endif(WIN32)
  execute_process(COMMAND ${cat_prog}
                  ${PROJECT_SOURCE_DIR}/cmake/logo.txt
                  OUTPUT_VARIABLE PYNE_LOGO)
  message("${PYNE_LOGO}")
endmacro()


# Set if spatial solver module should be built
macro(pyne_set_build_spatial_solver)
  set(BUILD_SPATIAL_SOLVER OFF)
  if ( ENABLE_SPATIAL_SOLVERS )
    message(STATUS "Checking whether to build spatial solvers")
    message(STATUS "Checking CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      message(STATUS "Checking CMAKE_CXX_COMPILER_VERSION: ${CMAKE_CXX_COMPILER_VERSION}")
      message(STATUS "Checking if APPLE: ${APPLE}")
      if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.6" AND NOT APPLE)
        message(STATUS "Building spatial solvers")
        set(BUILD_SPATIAL_SOLVER ON)
      else()
        message(WARNING "Not building spatial solvers.")
        set(BUILD_SPATIAL_SOLVER OFF)
      endif()
    endif()
  endif()
  message(STATUS "Build spatial solvers: ${BUILD_SPATIAL_SOLVER}")
endmacro()


# Set build type
macro(pyne_set_build_type)
  # Default to release build type
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "The build type" FORCE)
  endif(NOT CMAKE_BUILD_TYPE)

  # quiets fortify_source warnings when not compiling with optimizations
  # in linux distros where compilers were compiled with fortify_source enabled by
  # default (e.g. Arch linux).
  string(TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE)
  if(NOT ${BUILD_TYPE} STREQUAL "release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0")
  endif()
  message( STATUS "Build type: ${CMAKE_BUILD_TYPE}")
endmacro()

# Download source files
macro(pyne_download_files)
  # Download bateman solver from PyNE data
  download_src("http://raw.githubusercontent.com/pyne/data/master" "decay" ".cpp")
  # Download CRAM solver from PyNE data
  download_src("http://raw.githubusercontent.com/pyne/data/master" "cram" ".c")

  if(PYNE_FAST_COMPILE)
    if(NOT WIN32)
      # Download bateman solver from PyNE data
      download_platform_specific("http://raw.githubusercontent.com/pyne/data/master" "decay" ".s")
      # Download CRAM solver from PyNE data
      download_platform_specific("http://raw.githubusercontent.com/pyne/data/master" "cram" ".s")
    else()
      message(WARNING "Not downloading compiled files on Windows.")
      set(PYNE_FAST_COMPILE OFF)
    endif()
  endif()
endmacro()

# fast compile with assembly, if available.
macro(fast_compile _srcname _gnuflags _clangflags _otherflags)
  get_filename_component(_base "${_srcname}" NAME_WE)  # get the base name, without the extension
  # get the assembly file name
  if(PYNE_ASM_PLATFORM)
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
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set_source_files_properties("${_filename}" PROPERTIES COMPILE_FLAGS "${_clangflags}")
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set_source_files_properties("${_filename}" PROPERTIES COMPILE_FLAGS "${_gnuflags}")
  else()
    set_source_files_properties("${_filename}" PROPERTIES COMPILE_FLAGS "${_otherflags}")
  endif()
endmacro()
