# TRY TO FIND THE MA57 LIBRARY

find_library(MA57_LIB
    NAMES MA57
    HINTS /usr/local/lib/)

if(MA57_LIB)
  message(STATUS "Found MA57: ${MA57_LIB}")
  set(MA57_LIBRARIES ${MA57_LIB})
  set(MA57_FOUND TRUE)
else(MA57_LIB)
  set(MA57_FOUND FALSE)
  message(STATUS "Could not find MA57")
endif(MA57_LIB)
