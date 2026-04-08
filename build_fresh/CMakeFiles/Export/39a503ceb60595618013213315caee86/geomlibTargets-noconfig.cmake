#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "geomlib::geomlib" for configuration ""
set_property(TARGET geomlib::geomlib APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(geomlib::geomlib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libgeomlib.a"
  )

list(APPEND _cmake_import_check_targets geomlib::geomlib )
list(APPEND _cmake_import_check_files_for_geomlib::geomlib "${_IMPORT_PREFIX}/lib/libgeomlib.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
