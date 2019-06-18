message(STATUS "Looking for ecos")

find_path(ECOS_CORE_INCLUDE_DIR
    ecos.h
  HINTS $ENV{ECOS}/include /usr/local/include/ecos ~/local/include/ecos
)

find_path(ECOS_AMD_INCLUDE_DIR
    amd.h
  HINTS $ENV{ECOS}/external/amd/include /usr/local/include/amd ~local/include/amd /usr/local/include/ecos ~local/include/ecos
)

find_path(ECOS_LDL_INCLUDE_DIR
    ldl.h
  HINTS $ENV{ECOS}/external/ldl/include /usr/local/include/ldl ~local/include/ldl /usr/local/include/ecos ~local/include/ecos /usr/local/include/ecos/ldl ~local/include/ecos/ldl 
)

find_path(ECOS_SSC_INCLUDE_DIR
    SuiteSparse_config.h
  HINTS $ENV{ECOS}/external/SuiteSparse_config /usr/local/include/SuiteSparse_config ~local/include/SuiteSparse_config /usr/local/include/ecos ~local/include/ecos /usr/local/include/ecos/SuiteSparse_config ~local/include/ecos/SuiteSparse_config 
)

if(ECOS_CORE_INCLUDE_DIR AND ECOS_AMD_INCLUDE_DIR AND ECOS_LDL_INCLUDE_DIR AND ECOS_SSC_INCLUDE_DIR)
  set(ECOS_INCLUDE_DIR "${ECOS_CORE_INCLUDE_DIR};${ECOS_AMD_INCLUDE_DIR};${ECOS_LDL_INCLUDE_DIR};${ECOS_SSC_INCLUDE_DIR}")
endif()

if(ECOS_INCLUDE_DIR)
  message(STATUS "Found ecos include directory: ${ECOS_INCLUDE_DIR}")
else()
  message(STATUS "Could not find ecos include dir")
endif()

find_library(ECOS_LIBRARY
  NAMES ecos
  PATHS ~/local/lib /usr/local/lib
  $ENV{ECOS} $ENV{ECOS}/lib
)

if(ECOS_LIBRARY)
  set(FOUND_ECOS_LIBS TRUE)
  set(ECOS_LIBRARIES ${ECOS_LIBRARY})
endif()

if(ECOS_LIBRARIES)
  set(ECOS_LIBRARIES ${ECOS_LIBRARIES})
  message(STATUS "Found ecos libraries ${ECOS_LIBRARIES}")
  set(ECOS_FOUND_LIBS TRUE)
else()
  set(ECOS_FOUND_LIBS FALSE)
  message(STATUS "Could not find ecos libraries ${ECOS_LIBRARIES}")
endif()

if(ECOS_INCLUDE_DIR AND ECOS_FOUND_LIBS)
  set(ECOS_FOUND TRUE)
else()
  set(ECOS_FOUND FALSE)
  message(STATUS "ECOS: Cound not find ecos. Try setting ECOS env var.")
endif()

# Standard ecos build depends on rt for profiling
if(WIN32)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
else()
  find_library(RT rt)
  set(ECOS_LIBRARIES ${ECOS_LIBRARIES} "${RT}")
endif()
