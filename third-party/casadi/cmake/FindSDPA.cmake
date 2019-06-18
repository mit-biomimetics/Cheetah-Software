find_path(SDPA_INCLUDE_DIR
sdpa_call.h
~/local/sdpa.7.3.1/include/
~/local/include/
/usr/local/include/
/usr/local/sdpa.7.3.1/include/
)

find_path(MUMPS_INCLUDE_DIR
dmumps_root.h
${SDPA_INCLUDE_DIR}/../mumps/build/include/
~/local/include/
/usr/local/include/)

if(SDPA_INCLUDE_DIR AND MUMPS_INCLUDE_DIR)
  set(FOUND_SDPA_INCLUDES TRUE)
  set(SDPA_INCLUDE_DIRS ${SDPA_INCLUDE_DIR} ${MUMPS_INCLUDE_DIR})
  message(STATUS "Found SDPA include dirs: ${SDPA_INCLUDE_DIRS}")
  message(STATUS "Found MUMPS include dirs: ${MUMPS_INCLUDE_DIR}")
else()
  message(STATUS "Could not find SDPA include dirs")
endif()

find_library(SDPA_LIBRARY
  NAMES sdpa
  PATHS ${SDPA_LIB_DIR} ${SDPA_INCLUDE_DIR}/../lib ~/local/lib /usr/local/lib)

find_library(MUMPS_LIBRARY
  NAMES dmumps
  PATHS ${MUMPS_INCLUDE_DIR}/../lib ~/local/lib /usr/local/lib)

find_library(MUMPS_COMMON_LIBRARY
  NAMES mumps_common
  PATHS ${MUMPS_INCLUDE_DIR}/../lib ~/local/lib /usr/local/lib)

find_library(MUMPS_PORD_LIBRARY
  NAMES pord
  PATHS ${MUMPS_INCLUDE_DIR}/../lib ~/local/lib /usr/local/lib)

find_library(MUMPS_MPISEQ_LIBRARY
  NAMES mpiseq
  PATHS ${MUMPS_INCLUDE_DIR}/../libseq ~/local/lib /usr/local/lib)

if(SDPA_LIBRARY AND MUMPS_LIBRARY AND MUMPS_COMMON_LIBRARY AND MUMPS_PORD_LIBRARY AND MUMPS_MPISEQ_LIBRARY)
  set(FOUND_SDPA_LIBS TRUE)
  set(SDPA_LIBRARIES ${SDPA_LIBRARY} ${MUMPS_LIBRARY} ${MUMPS_COMMON_LIBRARY} ${MUMPS_PORD_LIBRARY} ${MUMPS_MPISEQ_LIBRARY} lapack pthread)
  message(STATUS "Found SDPA libs: ${SDPA_LIBRARIES}")
else()
  message(STATUS "Could not find SDPA libs")
endif()

if(FOUND_SDPA_INCLUDES AND FOUND_SDPA_LIBS)
  set(FOUND_SDPA TRUE)
endif()
