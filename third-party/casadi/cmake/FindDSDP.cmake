message(STATUS "Looking for dsdp")

find_path(DSDP_INCLUDE_DIR
    dsdp5.h
  HINTS $ENV{DSDP}/include /usr/include/dsdp
  )

if(DSDP_INCLUDE_DIR)
  message(STATUS "Found dsdp include directory: ${DSDP_INCLUDE_DIR}")
else()
  message(STATUS "Could not find dsdp include dir")
endif()

# libraries
set(DSDP_LIBS_LIST
  dsdp
  )

set(DSDP_LIBRARIES)
foreach(LIB ${DSDP_LIBS_LIST})
  find_library(DSDP_LIB_${LIB}
    NAMES ${LIB}
    HINTS $ENV{DSDP}/lib)
  if(DSDP_LIB_${LIB})
#    message(STATUS "Found ${LIB}: ${DSDP_LIB_${LIB}}")
    set(DSDP_LIBRARIES ${DSDP_LIBRARIES} ${DSDP_LIB_${LIB}})
#  else()
#    message(STATUS "Could not find lib${LIB}")
  endif()
endforeach()

if(DSDP_LIBRARIES)
  set(DSDP_LIBRARIES ${DSDP_LIBRARIES})
  message(STATUS "Found dsdp libraries ${DSDP_LIBRARIES}")
  set(DSDP_FOUND_LIBS TRUE)
else()
  set(DSDP_FOUND_LIBS FALSE)
  message(STATUS "Could not find dsdp libraries ${DSDP_LIBRARIES}")
endif()

if(DSDP_INCLUDE_DIR AND DSDP_FOUND_LIBS)
  set(DSDP_FOUND TRUE)
else()
  set(DSDP_FOUND FALSE)
  message(STATUS "DSDP: Cound not find dsdp. Try setting DSDP env var.")
endif()
