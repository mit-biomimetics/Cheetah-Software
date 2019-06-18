message(STATUS "Looking for blasfeo")

find_path(BLASFEO_INCLUDE_DIR
    blasfeo_common.h
  HINTS $ENV{BLASFEO}/include  /opt/blasfeo/include
  )

if(BLASFEO_INCLUDE_DIR)
  message(STATUS "Found blasfeo include directory: ${BLASFEO_INCLUDE_DIR}")
else()
  message(STATUS "Could not find blasfeo include dir")
endif()

find_library(BLASFEO_LIB
    NAMES blasfeo
    HINTS $ENV{BLASFEO}/lib /opt/blasfeo/lib
)

if(BLASFEO_LIB)
  set(BLASFEO_LIBRARIES ${BLASFEO_LIB})
  message(STATUS "Found blasfeo libraries ${BLASFEO_LIBRARIES}")
  set(BLASFEO_FOUND_LIBS TRUE)
else()
  message(STATUS "Could not find blasfeo libraries")
endif()

if(BLASFEO_INCLUDE_DIR AND BLASFEO_FOUND_LIBS)
  set(BLASFEO_FOUND TRUE)
else()
  message(STATUS "BLASFEO: Cound not find blasfeo. Try setting BLASFEO env var.")
endif()
