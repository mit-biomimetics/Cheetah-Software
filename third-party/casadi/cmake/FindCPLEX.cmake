# From martinsch/pgmlink and modified
# This module finds cplex.
#
# User can give CPLEX enviromental variable as a hint
#
# It sets the following variables:
#  CPLEX_FOUND              - Set to false, or undefined, if cplex isn't found.
#  CPLEX_INCLUDE_DIRS       - include directory
#  CPLEX_LIBRARIES          - library files
# if(WIN32)
#   execute_process(COMMAND cmd /C set CPLEX_STUDIO_DIR OUTPUT_VARIABLE CPLEX_STUDIO_DIR_VAR ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
  
#   if(NOT CPLEX_STUDIO_DIR_VAR)
#     MESSAGE(FATAL_ERROR "Unable to find CPLEX: environment variable CPLEX_STUDIO_DIR<VERSION> not set.")
#   endif()
  
#   STRING(REGEX REPLACE "^CPLEX_STUDIO_DIR" "" CPLEX_STUDIO_DIR_VAR ${CPLEX_STUDIO_DIR_VAR})
#   STRING(REGEX MATCH "^[0-9]+" CPLEX_WIN_VERSION ${CPLEX_STUDIO_DIR_VAR})
#   STRING(REGEX REPLACE "^[0-9]+=" "" CPLEX_STUDIO_DIR_VAR ${CPLEX_STUDIO_DIR_VAR})
#   file(TO_CMAKE_PATH "${CPLEX_STUDIO_DIR_VAR}" CPLEX_ROOT_DIR_GUESS) 
  
#   set(CPLEX_WIN_VERSION ${CPLEX_WIN_VERSION} CACHE STRING "CPLEX version to be used.")
#   set(CPLEX_ROOT_DIR "${CPLEX_ROOT_DIR_GUESS}" CACHE PATH "CPLEX root directory.")
  
#   MESSAGE(STATUS "Found CLPEX version ${CPLEX_WIN_VERSION} at '${CPLEX_ROOT_DIR}'")
  
#   STRING(REGEX REPLACE "/VC/bin/.*" "" VISUAL_STUDIO_PATH ${CMAKE_C_COMPILER})
#   STRING(REGEX MATCH "Studio [0-9]+" CPLEX_WIN_VS_VERSION ${VISUAL_STUDIO_PATH})
#   STRING(REGEX REPLACE "Studio " "" CPLEX_WIN_VS_VERSION ${CPLEX_WIN_VS_VERSION})
  
#   if(${CPLEX_WIN_VS_VERSION} STREQUAL "9")
#     set(CPLEX_WIN_VS_VERSION 2008)
#   elseif(${CPLEX_WIN_VS_VERSION} STREQUAL "10")
#     set(CPLEX_WIN_VS_VERSION 2010)
#   elseif(${CPLEX_WIN_VS_VERSION} STREQUAL "11")
#     set(CPLEX_WIN_VS_VERSION 2012)
#   else()
#     MESSAGE(FATAL_ERROR "CPLEX: unknown Visual Studio version at '${VISUAL_STUDIO_PATH}'.")
#   endif()
  
#   set(CPLEX_WIN_VS_VERSION ${CPLEX_WIN_VS_VERSION} CACHE STRING "Visual Studio Version")
  
#   if("${CMAKE_C_COMPILER}" MATCHES "amd64")
#     set(CPLEX_WIN_BITNESS x64)
#   else()
#     set(CPLEX_WIN_BITNESS x86)
#   endif()
  
#   set(CPLEX_WIN_BITNESS ${CPLEX_WIN_BITNESS} CACHE STRING "On Windows: x86 or x64 (32bit resp. 64bit)")

#   MESSAGE(STATUS "CPLEX: using Visual Studio ${CPLEX_WIN_VS_VERSION} ${CPLEX_WIN_BITNESS} at '${VISUAL_STUDIO_PATH}'")

#   if(NOT CPLEX_WIN_LINKAGE)
#     set(CPLEX_WIN_LINKAGE mda CACHE STRING "CPLEX linkage variant on Windows. One of these: mda (dll, release), mdd (dll, debug), mta (static, release), mtd (static, debug)")
#   endif(NOT CPLEX_WIN_LINKAGE)

#   # now, generate platform string
#   set(CPLEX_WIN_PLATFORM "${CPLEX_WIN_BITNESS}_windows_vs${CPLEX_WIN_VS_VERSION}/stat_${CPLEX_WIN_LINKAGE}")

# else()

  if( "$ENV{CPLEX}" STREQUAL "" )
    set(CPLEX_ROOT_DIR "" CACHE PATH "CPLEX root directory.")
  else()
    set(CPLEX_ROOT_DIR $ENV{CPLEX})
  endif()
  set(CPLEX_WIN_PLATFORM "")  

FIND_PATH(CPLEX_INCLUDE_DIR
  ilcplex/cplex.h
  HINTS ${CPLEX_ROOT_DIR}/cplex/include ${CPLEX_ROOT_DIR}/include
  PATHS ENV C_INCLUDE_PATH
        ENV C_PLUS_INCLUDE_PATH
        ENV INCLUDE_PATH
  )


if(WITH_CPLEX_SHARED)
  FIND_LIBRARY(CPLEX_LIBRARY
  NAMES cplex${CPLEX_VERSION} cplex
  HINTS ${CPLEX_ROOT_DIR}/cplex/bin/${CPLEX_WIN_PLATFORM} #windows
        ${CPLEX_ROOT_DIR}/cplex/bin/x86-64_debian4.0_4.1 #unix
        ${CPLEX_ROOT_DIR}/cplex/bin/x86-64_sles10_4.1 #unix 
        ${CPLEX_ROOT_DIR}/cplex/bin/x86-64_linux #unix 
        ${CPLEX_ROOT_DIR}/cplex/bin/x86-64_osx #osx 
        ${CPLEX_ROOT_DIR}/cplex/bin/x86-64_darwin #osx 
        ${CPLEX_ROOT_DIR}/cplex/bin/x64_win64 #windows
        ${CPLEX_ROOT_DIR}/bin/${CPLEX_WIN_PLATFORM} #windows
        ${CPLEX_ROOT_DIR}/bin/x86-64_debian4.0_4.1 #unix
        ${CPLEX_ROOT_DIR}/bin/x86-64_sles10_4.1 #unix 
        ${CPLEX_ROOT_DIR}/bin/x86-64_linux #unix 
        ${CPLEX_ROOT_DIR}/bin/x86-64_osx #osx 
        ${CPLEX_ROOT_DIR}/bin/x86-64_darwin #osx 
        ${CPLEX_ROOT_DIR}/bin/x64_win64 #windows
  PATHS ENV LIBRARY_PATH #unix
        ENV LD_LIBRARY_PATH #unix
  )
  message(STATUS "CPLEX Library: ${CPLEX_LIBRARY}")
else()

  FIND_LIBRARY(CPLEX_LIBRARY
    NAMES cplex${CPLEX_WIN_VERSION} cplex
    HINTS ${CPLEX_ROOT_DIR}/cplex/lib/${CPLEX_WIN_PLATFORM} #windows
          ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_debian4.0_4.1/static_pic #unix
          ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_sles10_4.1/static_pic #unix 
          ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_linux/static_pic #unix 
          ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_osx/static_pic #osx 
          ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_darwin/static_pic #osx 
    PATHS ENV LIBRARY_PATH #unix
          ENV LD_LIBRARY_PATH #unix
    )
  message(STATUS "CPLEX Library: ${CPLEX_LIBRARY}")
endif()

INCLUDE(FindPackageHandleStandardArgs)

if(WITH_CPLEX_SHARED)

  FIND_PACKAGE_HANDLE_STANDARD_ARGS(CPLEX DEFAULT_MSG 
   CPLEX_LIBRARY CPLEX_INCLUDE_DIR)

  IF(CPLEX_FOUND)
    SET(CPLEX_INCLUDE_DIRS ${CPLEX_INCLUDE_DIR})
    SET(CPLEX_LIBRARIES ${CPLEX_LIBRARY} )
    IF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
      SET(CPLEX_LIBRARIES "${CPLEX_LIBRARIES};m;pthread")
    ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  ENDIF(CPLEX_FOUND)

  MARK_AS_ADVANCED(CPLEX_LIBRARY CPLEX_INCLUDE_DIR)

else()
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(CPLEX DEFAULT_MSG 
   CPLEX_LIBRARY CPLEX_INCLUDE_DIR)

  IF(CPLEX_FOUND)
    SET(CPLEX_INCLUDE_DIRS ${CPLEX_INCLUDE_DIR})
    SET(CPLEX_LIBRARIES ${CPLEX_LIBRARY} )
    IF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
      SET(CPLEX_LIBRARIES "${CPLEX_LIBRARIES};m;pthread")
    ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  ENDIF(CPLEX_FOUND)

  MARK_AS_ADVANCED(CPLEX_LIBRARY CPLEX_INCLUDE_DIR)

endif()

