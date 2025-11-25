#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "geomlib::geomlib" for configuration "Release"
set_property(TARGET geomlib::geomlib APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(geomlib::geomlib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libgeomlib.a"
  )

list(APPEND _cmake_import_check_targets geomlib::geomlib )
list(APPEND _cmake_import_check_files_for_geomlib::geomlib "${_IMPORT_PREFIX}/lib/libgeomlib.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
