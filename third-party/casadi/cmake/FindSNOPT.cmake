# Library
find_library(SNOPT_LIBRARY
snopt7
HINTS $ENV{SNOPT} $ENV{SNOPT}/lib
)

if(SNOPT_LIBRARY)
  set(SNOPT_LIBRARIES ${SNOPT_LIBRARY})
  message(STATUS "Found SNOPT library: ${SNOPT_LIBRARY}")
else()
  message(STATUS "Could not find SNOPT libs")
endif()

if(SNOPT_LIBRARIES)
  set(SNOPT_FOUND TRUE)
endif()
