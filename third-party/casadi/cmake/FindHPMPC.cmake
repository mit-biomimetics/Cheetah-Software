message(STATUS "Looking for hpmpc")

find_path(HPMPC_INCLUDE_DIR
    mpc_solvers.h
  HINTS $ENV{HPMPC}/include /opt/hpmpc/include
  )

if(HPMPC_INCLUDE_DIR)
  #set(HPMPC_INCLUDE_DIR ${HPMPC_INCLUDE_DIR} ${B_INCLUDE_DIR})
  message(STATUS "Found hpmpc include directory: ${HPMPC_INCLUDE_DIR}")
else()
  message(STATUS "Could not find hpmpc include dir")
endif()

find_library(HPMPC_LIB
    NAMES hpmpc
    HINTS $ENV{CSPARSE}/lib /opt/hpmpc/lib
)

if(HPMPC_LIB)
  set(HPMPC_LIBRARIES ${HPMPC_LIB} ${BLASFEO_LIBRARIES})
  message(STATUS "Found hpmpc libraries ${HPMPC_LIBRARIES}")
  set(HPMPC_FOUND_LIBS TRUE)
else()
  message(STATUS "Could not find hpmpc libraries")
endif()

if(HPMPC_INCLUDE_DIR AND HPMPC_FOUND_LIBS)
  set(HPMPC_FOUND TRUE)
else()
  message(STATUS "HPMPC: Cound not find hpmpc. Try setting HPMPC env var.")
endif()


