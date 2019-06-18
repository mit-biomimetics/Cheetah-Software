# TRY TO FIND THE WSMP LIBRARY

find_library(WSMP_LIB
    names wsmp64 wsmp
    HINTS $ENV{WSMP})

if(WSMP_LIB)
  message(STATUS "Found WSMP: ${WSMP_LIB}")
  set(WSMP_LIBRARIES ${WSMP_LIB})
  set(WSMP_FOUND TRUE)
else()
  set(WSMP_FOUND FALSE)
  message(STATUS "Could not find WSMP; looking in environmental variable WSMP ($ENV{WSMP})")
endif()
