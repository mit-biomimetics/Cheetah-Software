# Get package info using pkg-config
find_package(PkgConfig)
pkg_search_module(CLP clp)

include(canonicalize_paths)
canonicalize_paths(CLP_LIBRARY_DIRS)

# add osx frameworks to CLP_LIBRARIES
if(CLP_LIBRARIES)
  if(APPLE)
    # turn "-framework;foo;-framework;bar;other" into "-framework foo;-framework bar;other"
    string(REPLACE "-framework;" "-framework " CLP_LDFLAGS_OTHER "${CLP_LDFLAGS_OTHER}")
    # take "-framework foo;-framework bar;other" and add only frameworks to CLP_LIBRARIES
    foreach(arg ${CLP_LDFLAGS_OTHER})
      if("${arg}" MATCHES "-framework .+")
        set(CLP_LIBRARIES "${CLP_LIBRARIES};${arg}")
      endif("${arg}" MATCHES "-framework .+")
    endforeach(arg ${CLP_LDFLAGS_OTHER})
  endif(APPLE)
endif(CLP_LIBRARIES)

# Set standard flags
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CLP DEFAULT_MSG CLP_LIBRARIES CLP_INCLUDE_DIRS)
