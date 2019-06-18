# TRY TO FIND THE MA27 LIBRARY

find_library(MA27_LIB
    NAMES MA27
    HINTS /usr/local/lib/)

if(MA27_LIB)
  message(STATUS "Found MA27: ${MA27_LIB}")
  set(MA27_LIBRARIES ${MA27_LIB})
  set(MA27_FOUND TRUE)
else(MA27_LIB)
  set(MA27_FOUND FALSE)
  message(STATUS "Could not find MA27")
endif(MA27_LIB)
