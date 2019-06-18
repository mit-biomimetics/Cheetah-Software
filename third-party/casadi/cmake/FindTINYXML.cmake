# Get package info using pkg-config
find_package(PkgConfig)
pkg_search_module(TINYXML tinyxml)

include(canonicalize_paths)
canonicalize_paths(TINYXML_LIBRARY_DIRS)

if(NOT TINYXML_INCLUDE_DIR)
find_path(TINYXML_INCLUDE_DIR
    tinyxml.h
#  HINTS $ENV{SUNDIALS}/include
  )
endif()

# Set standard flags
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TINYXML DEFAULT_MSG TINYXML_LIBRARIES TINYXML_INCLUDE_DIR)
