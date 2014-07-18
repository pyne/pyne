#
# This macro builds the f2py module
#   - target_name
#   - 
#
macro (BuildF2pyModule target_name modulename module_pyf_name filelist1)

 set(filelist ${filelist1}  ${ARGN})
 set(filename temp_script.py)
 # Copy all the files
 EXECUTE_PROCESS(COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/${module_pyf_name} ${CMAKE_CURRENT_BINARY_DIR} )
 FOREACH( f ${filelist})
  EXECUTE_PROCESS(COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/${f} ${CMAKE_CURRENT_BINARY_DIR} )
 ENDFOREACH(f)
 # write the script that will build the f2py extension
 SET(filename ${CMAKE_CURRENT_BINARY_DIR}/${filename} )
 FILE(WRITE ${filename} "import sys\n")
 FILE(APPEND ${filename} "from numpy.f2py import main\n")
 FILE(APPEND ${filename} "sys.argv = [''] +'-c --f77exec=${CMAKE_Fortran_COMPILER} --f90exec=${CMAKE_Fortran_COMPILER} -m ${modulename} ${modulename}.pyf ${filelist} -llapack'.split()\n")
 FILE(APPEND ${filename} "main()\n")

 # We had the normal target of the module
 add_custom_target(${target_name} ALL DEPENDS ${modulename}.so) 

 # TODO : to be corrected with the filelist is more than one file.
 # ... and a special target to build vertex.so, that depends on the sources files
# add_custom_command(OUTPUT  ${modulename}.so 
#  COMMAND echo See `pwd`/f2pyBuild.log for logs
#  COMMAND ${TRIQS_PYTHON_INTERPRETER} temp_script.py > f2pyBuild.log 2>&1
#  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${filelist} ${CMAKE_CURRENT_SOURCE_DIR}/${module_pyf_name} 
# Need a DEPENDS that hooks up spatialsolver.py with this?
  )

 #FILE(RELATIVE_PATH rel ${CMAKE_SOURCE_DIR}/Modules ${CMAKE_CURRENT_SOURCE_DIR}/)
 #install (FILES ${CMAKE_CURRENT_BINARY_DIR}/${modulename}.so
 # DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/${ExecutableName}/triqs/${rel}/..)

endmacro (BuildF2pyModule)
