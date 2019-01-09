# Macros for downloading a tar file and extracting it if a sample file
# doesn't already exist.

macro(download_and_extract _url _checkfile)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${_checkfile}")
    # don't need to do anything
    message(STATUS "${_checkfile} exists, no need to download or extract!")
  else()
    get_filename_component(_base "${_url}" NAME)
    # download the file if we need to
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${_base}")
      message(STATUS "Downloading ${_url} -> ${_base}")
      file(DOWNLOAD "${_url}" "${CMAKE_CURRENT_SOURCE_DIR}/${_base}"
           SHOW_PROGRESS STATUS _rtn)
      list(GET _rtn 0 _rtncode)
      if(NOT 0 EQUAL _rtncode)
        message(FATAL_ERROR ${_rtn})
      endif()
    endif()
    # extract the file
    message(STATUS "Extracting ${_base}")
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xvf "${_base}"
                    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
  endif()
endmacro()


macro(download_platform _base_url _base_name _ext_src _ext_plat)
  # first download the src file
  set(_url "${_base_url}/${_base_name}.tar.gz")
  set(_checkfile "${_base_name}${_ext_src}")
  download_and_extract("${_url}" "${_checkfile}")
  # now download the platform specific file
  set(_url "${_base_url}/${_base_name}-${PYNE_ASM_PLATFORM}.tar.gz")
  set(_checkfile "${_base_name}-${PYNE_ASM_PLATFORM}${_ext_plat}")
  download_and_extract("${_url}" "${_checkfile}")
endmacro()

