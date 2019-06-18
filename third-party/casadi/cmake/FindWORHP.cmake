# TRY TO FIND THE INCLUDE DIRECTORY
find_file(WORHP_LICENSE_KEY
worhp.lic
HINTS $ENV{WORHP} $ENV{WORHP_LICENSE_FILE}
)

if(WORHP_LICENSE_KEY)
  message(STATUS "Worhp license detected: ${WORHP_LICENSE_KEY}")
else()
  message(STATUS "Could not find worhp license key, looking in $ENV{WORHP}")
endif()

find_path(WORHP_INCLUDE_DIR
worhp.h
HINTS $ENV{WORHP}/include/worhp/ $ENV{WORHP})

if(WORHP_INCLUDE_DIR)
  message(STATUS "Worhp include dir detected: ${WORHP_INCLUDE_DIR}")
else()
  message(STATUS "Could not find worhp include dir, looking in $ENV{WORHP}/include/worhp/")
endif()

find_library(WORHP_LIBRARY
  worhp
  HINTS $ENV{WORHP}/lib/ $ENV{WORHP})


if(WORHP_LIBRARY)
  message(STATUS "Worhp library detected: ${WORHP_LIBRARY}")
else()
  message(STATUS "Could not find worhp library, looking in $ENV{WORHP}/lib/")
endif()

set(WORHP_LIBRARIES ${WORHP_LIBRARY})


if(WORHP_LIBRARY AND WORHP_INCLUDE_DIR)
  message(STATUS "Worhp interface ready with these libraries: ${WORHP_LIBRARIES} ")
  set(WORHP_FOUND TRUE)
else()
  message(STATUS "Will not compile worhp interface")
endif()
