# set platform preprocessor macro
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

# Fortran settings
# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -c -fpic -fdefault-real-8 -fdefault-double-8")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -fdefault-real-8 -fdefault-double-8")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3 -r8")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g -r8")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

macro( add_lib_to_pyne _name _source )
  # add the library
  add_library(${_name}             ${_source})
  # add it to the list of pyne libraries
  set(PYNE_LIBRARIES ${PYNE_LIBRARIES} ${_name})
endmacro()

macro( print_logo )
  set(cat_prog cat)
  if(WIN32)
    set(cat_prog type)
  endif(WIN32)
  execute_process(COMMAND ${cat_prog} ${PROJECT_SOURCE_DIR}/cmake/logo.txt OUTPUT_VARIABLE variable)
  message("${variable}")
endmacro()
