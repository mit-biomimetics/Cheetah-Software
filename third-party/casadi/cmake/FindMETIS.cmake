# TRY TO FIND THE METIS LIBRARY

find_library(METIS_LIB
    names coinmetis metis
    HINTS $ENV{METIS})

if(METIS_LIB)
  message(STATUS "Found METIS: ${METIS_LIB}")
  set(METIS_LIBRARIES ${METIS_LIB})
  set(METIS_FOUND TRUE)
else(METIS_LIB)
  set(METIS_FOUND FALSE)
  message(STATUS "Could not find METIS; looking in environmental variable METIS ($ENV{METIS})")
endif(METIS_LIB)
