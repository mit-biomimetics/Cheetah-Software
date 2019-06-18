message(STATUS "Looking for mosek")

find_path(MOSEK_INCLUDE_DIR
    mosek.h
  HINTS ~/mosek/7/tools/platform/linux64x86/h $ENV{MOSEK}/h $ENV{MOSEK}/7/tools/platform/linux64x86/h
)

if(MOSEK_INCLUDE_DIR)
  message(STATUS "Found mosek include directory: ${MOSEK_INCLUDE_DIR}")
else()
  message(STATUS "Could not find mosek include dir")
endif()

find_library(MOSEK_LIBRARY
  NAMES mosek64
  PATHS ~/mosek/7/tools/platform/linux64x86/bin
  $ENV{MOSEK}/bin $ENV{MOSEK}/7/tools/platform/linux64x86/bin ~/local/lib /usr/local/lib
)

if(MOSEK_LIBRARY)
  set(FOUND_MOSEK_LIBS TRUE)
  set(MOSEK_LIBRARIES ${MOSEK_LIBRARY})
endif()

if(MOSEK_LIBRARIES)
  set(MOSEK_LIBRARIES ${MOSEK_LIBRARIES})
  message(STATUS "Found mosek libraries ${MOSEK_LIBRARIES}")
  set(MOSEK_FOUND_LIBS TRUE)
else()
  set(MOSEK_FOUND_LIBS FALSE)
  message(STATUS "Could not find mosek libraries ${MOSEK_LIBRARIES}")
endif()

if(MOSEK_INCLUDE_DIR AND MOSEK_FOUND_LIBS)
  set(MOSEK_FOUND TRUE)
else()
  set(MOSEK_FOUND FALSE)
  message(STATUS "MOSEK: Cound not find mosek. Try setting MOSEK env var.")
endif()


