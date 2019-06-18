# TRY TO FIND THE INCLUDE DIRECTORY
find_path(SQIC_INCLUDE_DIR
  sqic.mod
  PATHS $ENV{SQIC_INCLUDE_DIR})

if(NOT SQIC_INCLUDE_DIR)
  message(STATUS "SQIC: Cound not find include dir. Try setting SQIC_INCLUDE_DIR env var.")
endif()

find_library(SQIC_LIB_SNBLAS
  snblas
  PATHS $ENV{SQIC_LIBRARY_DIR})

if(NOT SQIC_LIB_SNBLAS)
  message(STATUS "SQIC: Cound not find library snblas. Try setting SQIC_LIBRARY_DIR env var.")
endif()

find_library(SQIC_LIB_SQIC
  sqic
  PATHS $ENV{SQIC_LIBRARY_DIR})

if(NOT SQIC_LIB_SQIC)
  message(STATUS "SQIC: Cound not find library sqic. Try setting SQIC_LIBRARY_DIR env var.")
endif()

if(SQIC_INCLUDE_DIR AND SQIC_LIB_SNBLAS AND SQIC_LIB_SQIC)
  set(SQIC_LIBRARIES ${SQIC_LIB_SNBLAS} ${SQIC_LIB_SQIC})
  set(SQIC_FOUND TRUE)
endif()
