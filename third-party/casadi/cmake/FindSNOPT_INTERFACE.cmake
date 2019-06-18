# Library
find_library(SNOPT_INTERFACE_LIBRARY
snopt7_c
HINTS $ENV{SNOPT} $ENV{SNOPT}/lib $ENV{SNOPT_INTERFACE} $ENV{SNOPT_INTERFACE}/lib
)

if(SNOPT_INTERFACE_LIBRARY)
  set(SNOPT_INTERFACE_LIBRARIES ${SNOPT_INTERFACE_LIBRARY})
  message(STATUS "Found SNOPT interface library: ${SNOPT_INTERFACE_LIBRARY}")
else()
  message(STATUS "Could not find SNOPT interface libs")
endif()

if(SNOPT_INTERFACE_LIBRARIES)
  set(SNOPT_INTERFACE_FOUND TRUE)
endif()

find_path(SNOPT_INTERFACE_INCLUDE_DIR
snopt_cwrap.h
HINTS $ENV{SNOPT}/include $ENV{SNOPT_INTERFACE}/include
)

if(SNOPT_INTERFACE_INCLUDE_DIR)
   set(SNOPT_INTERFACE_FOUND_INCLUDE TRUE)
   message(STATUS "Found SNOPT interface include dir: ${SNOPT_INTERFACE_INCLUDE_DIR}")
else()
   message(STATUS "Could not find SNOPT interface include dir")
endif()
