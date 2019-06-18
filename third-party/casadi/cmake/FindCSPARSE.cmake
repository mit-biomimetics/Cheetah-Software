message(STATUS "Looking for csparse")

find_path(CSPARSE_INCLUDE_DIR
    cs.h
  HINTS $ENV{CSPARSE}/include /usr/include/suitesparse
  )

if(CSPARSE_INCLUDE_DIR)
  message(STATUS "Found csparse include directory: ${CSPARSE_INCLUDE_DIR}")
else()
  message(STATUS "Could not fund csparse include dir")
endif()

# libraries
set(CSPARSE_LIBS_LIST
  csparse
  cxsparse
  )

set(CSPARSE_LIBRARIES)
foreach(LIB ${CSPARSE_LIBS_LIST})
  find_library(CSPARSE_LIB_${LIB}
    NAMES ${LIB}
    HINTS $ENV{CSPARSE}/lib)
  if(CSPARSE_LIB_${LIB})
#    message(STATUS "Found ${LIB}: ${CSPARSE_LIB_${LIB}}")
    set(CSPARSE_LIBRARIES ${CSPARSE_LIBRARIES} ${CSPARSE_LIB_${LIB}})
#  else()
#    message(STATUS "Could not find lib${LIB}")
  endif()
endforeach()

if(CSPARSE_LIBRARIES)
  set(CSPARSE_LIBRARIES ${CSPARSE_LIBRARIES})
  message(STATUS "Found csparse libraries ${CSPARSE_LIBRARIES}")
  set(CSPARSE_FOUND_LIBS TRUE)
else()
  message(STATUS "Could not find csparse libraries")
endif()

if(CSPARSE_INCLUDE_DIR AND CSPARSE_FOUND_LIBS)
  set(CSPARSE_FOUND TRUE)
else()
  message(STATUS "CSPARSE: Cound not find csparse. Try setting CSPARSE env var.")
endif()
