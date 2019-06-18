# Get package info using pkg-config
find_package(PkgConfig)
pkg_search_module(BONMIN bonmin)

include(canonicalize_paths)
canonicalize_paths(BONMIN_LIBRARY_DIRS)

# add osx frameworks to BONMIN_LIBRARIES
if(BONMIN_LIBRARIES)
  if(APPLE)
    # turn "-framework;foo;-framework;bar;other" into "-framework foo;-framework bar;other"
    string(REPLACE "-framework;" "-framework " BONMIN_LDFLAGS_OTHER "${BONMIN_LDFLAGS_OTHER}")
    # take "-framework foo;-framework bar;other" and add only frameworks to BONMIN_LIBRARIES
    foreach(arg ${BONMIN_LDFLAGS_OTHER})
      if("${arg}" MATCHES "-framework .+")
        set(BONMIN_LIBRARIES "${BONMIN_LIBRARIES};${arg}")
      endif("${arg}" MATCHES "-framework .+")
    endforeach(arg ${BONMIN_LDFLAGS_OTHER})
  endif(APPLE)
endif(BONMIN_LIBRARIES)

# Set standard flags
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BONMIN DEFAULT_MSG BONMIN_LIBRARIES BONMIN_INCLUDE_DIRS)
