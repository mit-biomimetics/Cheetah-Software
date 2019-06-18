# Get package info using pkg-config
find_package(PkgConfig)
pkg_search_module(IPOPT ipopt)

include(canonicalize_paths)
canonicalize_paths(IPOPT_LIBRARY_DIRS)

# add osx frameworks to IPOPT_LIBRARIES
if(IPOPT_LIBRARIES)
  if(APPLE)
    # turn "-framework;foo;-framework;bar;other" into "-framework foo;-framework bar;other"
    string(REPLACE "-framework;" "-framework " IPOPT_LDFLAGS_OTHER "${IPOPT_LDFLAGS_OTHER}")
    # take "-framework foo;-framework bar;other" and add only frameworks to IPOPT_LIBRARIES
    foreach(arg ${IPOPT_LDFLAGS_OTHER})
      if("${arg}" MATCHES "-framework .+")
        set(IPOPT_LIBRARIES "${IPOPT_LIBRARIES};${arg}")
      endif("${arg}" MATCHES "-framework .+")
    endforeach(arg ${IPOPT_LDFLAGS_OTHER})
  endif(APPLE)
endif(IPOPT_LIBRARIES)

# Callback support
if(IPOPT_INCLUDEDIR)
  if(EXISTS ${IPOPT_INCLUDEDIR}/IpIpoptData.hpp AND EXISTS ${IPOPT_INCLUDEDIR}/IpOrigIpoptNLP.hpp AND EXISTS ${IPOPT_INCLUDEDIR}/IpTNLPAdapter.hpp  AND EXISTS ${IPOPT_INCLUDEDIR}/IpDenseVector.hpp AND EXISTS ${IPOPT_INCLUDEDIR}/IpExpansionMatrix.hpp)
    set(WITH_IPOPT_CALLBACK TRUE)
  else()
   list(APPEND IPOPT_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/external_packages/ipopt/3.10")
   set(WITH_IPOPT_CALLBACK TRUE)
   # message(STATUS "Detected an IPOPT configuration without development headers. Build will proceed, but without callback functionality. To enable it, see https://github.com/casadi/casadi/wiki/enableIpoptCallback")
  endif()
endif(IPOPT_INCLUDEDIR)

# Set standard flags
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(IPOPT DEFAULT_MSG IPOPT_LIBRARIES IPOPT_INCLUDE_DIRS)
