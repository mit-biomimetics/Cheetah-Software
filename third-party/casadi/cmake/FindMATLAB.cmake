# - this module looks for Matlab
# Defines:
#  MATLAB_INCLUDE_DIR: include path for mex.h, engine.h
#  MATLAB_LIBRARIES:   required libraries: libmex, etc
#  MATLAB_MEX_LIBRARY: path to libmex.lib
#  MATLAB_MX_LIBRARY:  path to libmx.lib
#  MATLAB_MAT_LIBRARY:  path to libmat.lib
#  MATLAB_ENG_LIBRARY: path to libeng.lib
#  MATLAB_ROOT: path to Matlab's root directory
#  MATLAB_MEX_EXE: path to mex program
#  MATLAB_MEX_EXT: mex extension
#  MATLAB_CFLAGS, MATLAB_CLINKER_FLAGS, MATLAB_CLIBS (and same for CXX and Fortran)
#
# Note: You cannot use MATLAB_CLIBS etc to "target_link_libraries" as that gets confused 
# by the flag that specifies where the matlab libraries are (at least on Windows for Visual Studio).
# This variable also contains system libraries etc so it's probably not a good idea to use it in 
# your CMake file. You should probably use
#
#  include_directories(${MATLAB_INCLUDE_DIR})
#  add_definitions(${MATLAB_CXXFLAGS})
#  target_link_libraries(yourmexfile ${MATLAB_LIBRARIES} )
#
# Reason for MATLAB_CXXFLAGS: on linux, mex files need to be compiled with -fPIC, but that means all 
# linked libraries need to be compiled with -fPIC as well.

# This is a derivative work of file FindMatlab.cmake released with
# CMake v2.8, because the original seems to be a bit outdated
#
# This file uses ideas from Gerardus (Author: Ramon Casero <rcasero@...>, Tom Doel)
# and a post by Kent Williams
# http://www.cmake.org/pipermail/cmake/2013-December/056593.html
#
# Heavily modified by Kris Thielemans for WIN32 and recent versions of matlab

#=============================================================================
# Copyright 2005-2009 Kitware, Inc.
# Copyright 2014 University College London

# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

# If MATLAB_ROOT was defined in the environment, use it.
if (NOT MATLAB_ROOT AND NOT $ENV{MATLAB_ROOT} STREQUAL "")
  set(MATLAB_ROOT $ENV{MATLAB_ROOT} CACHE PATH "set this if CMake does not find it automatically")
  file(TO_CMAKE_PATH "${MATLAB_ROOT}" MATLAB_ROOT)
endif()


if (NOT MATLAB_ROOT)

  # get path to the Matlab executable
  find_program(MATLAB_EXE_PATH matlab
        PATHS /usr/local/bin)
  if (MATLAB_EXE_PATH)
    # found it, so let's find its absolute path and go one level up to find MATLAB_ROOT
    message(STATUS "MATLAB executable found: ${MATLAB_EXE_PATH}")
    # remove symbolic links
    get_filename_component(MATLAB_EXE_PATH "${MATLAB_EXE_PATH}" REALPATH )
    # find directory of executable
    get_filename_component(my_MATLAB_ROOT "${MATLAB_EXE_PATH}" PATH )
    # find root dir
    get_filename_component(my_MATLAB_ROOT "${my_MATLAB_ROOT}" PATH )
    # store it in the cache
    set(MATLAB_ROOT "${my_MATLAB_ROOT}" CACHE PATH "Location of MATLAB files")
  else()
    # First set MATLAB_ROOT to an empty string but as a cached variable
    # this avoids CMake creating a local variable with the same name
    set(MATLAB_ROOT "" CACHE PATH "Location of MATLAB files")

    # now go and look in more places
    if(WIN32)
      # Search for a version of Matlab available, starting from the most modern one to older versions
      foreach(MATVER "8.5" "8.4" "8.3" "8.2" "8.1" "8.0" "7.16" "7.15" "7.14" "7.13" "7.12" "7.11" "7.10" "7.9" "7.8" "7.7" "7.6" "7.5" "7.4")
        if((NOT DEFINED MATLAB_ROOT)
            OR ("${MATLAB_ROOT}" STREQUAL "")
            OR ("${MATLAB_ROOT}" STREQUAL "/registry"))
          get_filename_component(MATLAB_ROOT "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\${MATVER};MATLABROOT]" ABSOLUTE)
        endif()
      endforeach()
      if("${MATLAB_ROOT}" STREQUAL "/registry")
        set(MATLAB_ROOT "")
      endif()

    elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")   # Check if this is a Mac
      # we look in the applications folder
      # Search for a version of Matlab available, starting from the most modern one to older versions
      foreach(MATVER "R2015b" "R2015a" "R2014b" "R2014a" "R2013b" "R2013a" "R2012b" "R2012a" "R2011b" "R2011a" "R2010b" "R2010a" "R2009b" "R2009a" "R2008b")
        if((NOT DEFINED MATLAB_ROOT) OR ("${MATLAB_ROOT}" STREQUAL ""))
          if(EXISTS /Applications/MATLAB_${MATVER}.app)
            set(MATLAB_ROOT /Applications/MATLAB_${MATVER}.app)
          endif()
        endif()
      endforeach()
    else()
      # could go and look in some other default places for other OSses
    endif()
  endif()
endif()

if ("${MATLAB_ROOT}" STREQUAL "")
  message(STATUS "MATLAB not found. Set MATLAB_ROOT")
  # TODO should really skip rest of configuration as it will all fail anyway.
else()
  message(STATUS "MATLAB_ROOT set to ${MATLAB_ROOT}")
endif()

# Find out where MATLAB libraries are
IF (NOT MATLAB_LIBRARIES)
  if (WIN32)
    # Directory name depending on whether the Windows architecture is 32
    # bit or 64 bit
    if(CMAKE_SIZEOF_VOID_P MATCHES "4")
      set(WINDIR "win32")
    elseif(CMAKE_SIZEOF_VOID_P MATCHES "8")
      set(WINDIR "win64")
    else()
      message(FATAL_ERROR
        "CMAKE_SIZEOF_VOID_P (${CMAKE_SIZEOF_VOID_P}) doesn't indicate a valid platform")
    endif()

    # Folder where the MEX libraries are, depending on the Windows compiler
    if(${CMAKE_GENERATOR} MATCHES "Visual Studio 6")
      set(MATLAB_LIBRARIES_DIR "${MATLAB_ROOT}/extern/lib/${WINDIR}/microsoft/msvc60")
    elseif(${CMAKE_GENERATOR} MATCHES "Visual Studio 7")
      # Assume people are generally using Visual Studio 7.1,
      # if using 7.0 need to link to: ../extern/lib/${WINDIR}/microsoft/msvc70
      set(MATLAB_LIBRARIES_DIR "${MATLAB_ROOT}/extern/lib/${WINDIR}/microsoft/msvc71")
      # set(MATLAB_LIBRARIES_DIR "${MATLAB_ROOT}/extern/lib/${WINDIR}/microsoft/msvc70")
    elseif(${CMAKE_GENERATOR} MATCHES "Borland")
      # Assume people are generally using Borland 5.4,
      # if using 7.0 need to link to: ../extern/lib/${WINDIR}/microsoft/msvc70
      set(MATLAB_LIBRARIES_DIR "${MATLAB_ROOT}/extern/lib/${WINDIR}/microsoft/bcc54")
      # set(MATLAB_LIBRARIES_DIR "${MATLAB_ROOT}/extern/lib/${WINDIR}/microsoft/bcc50")
      # set(MATLAB_LIBRARIES_DIR "${MATLAB_ROOT}/extern/lib/${WINDIR}/microsoft/bcc51")
    elseif(${CMAKE_GENERATOR} MATCHES "Visual Studio*")
      # If the compiler is Visual Studio, but not any of the specific
      # versions above, we try our luck with the microsoft directory
      set(MATLAB_LIBRARIES_DIR "${MATLAB_ROOT}/extern/lib/${WINDIR}/microsoft/")
    endif()

  else(WIN32)

    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      set(LIBRARY_EXTENSION .dylib)
    else()
      set(LIBRARY_EXTENSION .so)
    endif()
    
    # There seems to be no way to let cmake do a search in subdirectories, so use "find"
    execute_process(
      COMMAND find "${MATLAB_ROOT}/bin" -name libmex${LIBRARY_EXTENSION}
      COMMAND xargs echo -n
      OUTPUT_VARIABLE mex_lib
      )
    if (mex_lib)
      # find directory of executable
      get_filename_component(MATLAB_LIBRARIES_DIR ${mex_lib} PATH )
    endif()
  endif(WIN32)

  if (NOT MATLAB_LIBRARIES_DIR)
    set(MATLAB_LIBRARIES_DIR ${MATLAB_ROOT})
  endif()

  message(STATUS "Will look for MATLAB libraries in ${MATLAB_LIBRARIES_DIR}")

endif() # MATLAB_LIBRARIES

if (WIN32)
  set(BATEXT .bat)
  set(LIBPREFIX lib)
endif()


# This is common to all platforms:

# Get paths to the Matlab MEX libraries
find_library(MATLAB_MEX_LIBRARY
    ${LIBPREFIX}mex
    HINTS ${MATLAB_LIBRARIES_DIR}
    )
find_library(MATLAB_MX_LIBRARY
    ${LIBPREFIX}mx
    HINTS ${MATLAB_LIBRARIES_DIR}
    )
find_library(MATLAB_MAT_LIBRARY
    ${LIBPREFIX}mat
    HINTS ${MATLAB_LIBRARIES_DIR}
    )
find_library(MATLAB_ENG_LIBRARY
    ${LIBPREFIX}eng
    HINTS ${MATLAB_LIBRARIES_DIR}
    )
find_library(MATLAB_UT_LIBRARY
    ${LIBPREFIX}ut
    HINTS ${MATLAB_LIBRARIES_DIR}
    )

# Get path to the include directory
find_path(MATLAB_INCLUDE_DIR
    "mex.h"
    HINTS  "${MATLAB_ROOT}/extern/include"  "${MATLAB_ROOT}/include"
    )

find_program( MATLAB_MEX_PATH mex${BATEXT}
             HINTS ${MATLAB_ROOT}/bin
             DOC "The mex program path"
            )

find_program( MATLAB_MEXEXT_PATH mexext${BATEXT}
             HINTS ${MATLAB_ROOT}/bin
             DOC "The mexext program path"
            )

# find mex extension
execute_process(
    COMMAND ${MATLAB_MEXEXT_PATH}
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE MATLAB_MEX_EXT
    )

set(MATLAB_LIBRARIES
  ${MATLAB_MEX_LIBRARY} ${MATLAB_MX_LIBRARY} ${MATLAB_ENG_LIBRARY} ${MATLAB_UT_LIBRARY}
  CACHE PATH "Libraries to link mex files"
)

######################### find flags using "mex -v"
# This is based on Kent Williams's post at
# http://www.cmake.org/pipermail/cmake/2013-December/056593.html

# mex -v outputs all the settings used for building MEX files, so 
# we can use it to grab the important variables needed

# This sets MATLAB_CFLAGS, MATLAB_CLINKER_FLAGS, MATLAB_CLIBS (and same for CXX and Fortran)

macro(MATLAB_GETFLAGS FILENAME)
 execute_process(COMMAND ${MATLAB_MEX_PATH} -v -n ${FILENAME}
  OUTPUT_VARIABLE mexOut
  ERROR_VARIABLE mexErr)

# parse mex output line by line by turning file into CMake list of lines
string(REGEX REPLACE "\r?\n" ";" _mexOut "${mexOut}")
foreach(line ${_mexOut})  
  if("${line}" MATCHES "[\t ]+DEFINES *:") # on Linux
    string(REGEX REPLACE "[\t ]+DEFINES *: *" "" mexDefines "${line}")
  elseif("${line}" MATCHES "[\t ]+COMPDEFINES *:") # on Windows
    string(REGEX REPLACE "[\t ]+COMPDEFINES *: *" "" mexDefines "${line}")
  elseif("${line}" MATCHES "[\t ]+LINKFLAGS *:")
    string(REGEX REPLACE "[\t ]+LINKFLAGS *: *" "" mexLdFlags "${line}")
    # get rid of /implib statement used on some older windows-matlab versions (refers to temp file)
    string(REGEX REPLACE "/implib:\".*\"" "" mexLdFlags "${mexLdFlags}")
  elseif("${line}" MATCHES "[\t ]+LINKLIBS *:")
    string(REGEX REPLACE "[\t ]+LINKLIBS *: *" "" mexLdLibs "${line}")
  elseif("${line}" MATCHES "[\t ]+LINKEXPORT *:")
    string(REGEX REPLACE "[\t ]+LINKEXPORT *: *" "" mexLdExport "${line}")
  elseif("${line}" MATCHES "[\t ]+CXXFLAGS *[:=]")
    string(REGEX REPLACE "[\t ]+CXXFLAGS *[:=] *" "" mexCxxFlags "${line}")
    #message(STATUS "mexcxx ${mexCxxFlags}")
  elseif("${line}" MATCHES "[\t ]+CFLAGS *[:=]")
    string(REGEX REPLACE "[\t ]+CFLAGS *[:=] *" "" mexCFlags "${line}")
  elseif("${line}" MATCHES "[\t ]+FFLAGS *[:=]")
    string(REGEX REPLACE "[\t ]+FFLAGS *[:=] *" "" mexFFlags "${line}")
  # pre-2014a flags
  elseif("${line}" MATCHES "[\t ]+CXXLIBS *[:=]")
    string(REGEX REPLACE "[\t ]+CXXLIBS *[:=] *" "" mexCxxLibs "${line}")
  elseif("${line}" MATCHES "[\t ]+CLIBS *[:=]")
    string(REGEX REPLACE "[\t ]+CLIBS *[:=] *" "" mexCLibs "${line}")
  elseif("${line}" MATCHES "[\t ]+FLIBS *[:=]")
    string(REGEX REPLACE "[\t ]+FLIBS *[:=] *" "" mexFLibs "${line}")
  elseif("${line}" MATCHES "[\t ]+LDFLAGS *[:=]")
    string(REGEX REPLACE "[\t ]+LDFLAGS *[:=] *" "" mexLdFlags "${line}")
  elseif("${line}" MATCHES "[\t ]+LDCXXFLAGS *[:=]")
    string(REGEX REPLACE "[\t ]+LDCXXFLAGS *[:=] *" "" mexLdCxxFlags "${line}")
  elseif("${line}" MATCHES "[\t ]+LDCFLAGS *[:=]")
    string(REGEX REPLACE "[\t ]+LDCFLAGS *[:=] *" "" mexLdCFlags "${line}")
  elseif("${line}" MATCHES "[\t ]+LDFLAGS *[:=]")
    string(REGEX REPLACE "[\t ]+LDFLAGS *[:=] *" "" mexLdFFlags "${line}")
  endif()
endforeach()
endmacro()

# Since 2014a or so, mex can only be used for one type of source and no longer
# reports all flags. We therefore need to run mex with different file types.
#### C
MATLAB_GETFLAGS(${PROJECT_SOURCE_DIR}/src/cmake/FindMATLAB_mextest.c)
set(MATLAB_CFLAGS "${mexDefines} ${mexCFlags}" CACHE STRING "Flags to compile C MATLAB Mex files (or libraries that link with them)")
set(MATLAB_CLINKER_FLAGS "${mexLdFlags} ${mexLdCFlags} ${mexLdExport}" CACHE STRING "Flags to link MATLAB C Mex files")
set(MATLAB_CLIBS "${mexLdLibs} ${mexCLibs}" CACHE STRING "Flags with libraries to link MATLAB C Mex files")

#### C++
MATLAB_GETFLAGS(${PROJECT_SOURCE_DIR}/src/cmake/FindMATLAB_mextest.cxx)
set(MATLAB_CXXFLAGS "${mexDefines} ${mexCxxFlags}" CACHE STRING "Flags to compile C++ MATLAB Mex files (or libraries that link with them)")
set(MATLAB_CXXLINKER_FLAGS "${mexLdFlags} ${mexLdCxxFlags} ${mexLdExport}" CACHE STRING "Flags to link MATLAB C++ Mex files")
set(MATLAB_CXXLIBS "${mexLdLibs} ${mexCxxLibs}" CACHE STRING "Flags with libraries to link MATLAB C++ Mex files")

#### Fortran
MATLAB_GETFLAGS(${PROJECT_SOURCE_DIR}/src/cmake/FindMATLAB_mextest.f)
set(MATLAB_FFLAGS "${mexDefines} ${mexFFlags}" CACHE STRING "Flags to compile Fortran MATLAB Mex files (or libraries that link with them)")
set(MATLAB_FLINKER_FLAGS "${mexLdFlags} ${mexLdFFlags} ${mexLdExport}" CACHE STRING "Flags to link MATLAB Fortran Mex files")
set(MATLAB_FLIBS "${mexLdLibs} ${mexFLibs}" CACHE STRING "Flags with libraries to link MATLAB Fortran Mex files")

# handle the QUIETLY and REQUIRED arguments and set MATLAB_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MATLAB "MATLAB not found. If you do have it, set MATLAB_ROOT and reconfigure" 
  MATLAB_ROOT MATLAB_INCLUDE_DIR  MATLAB_LIBRARIES
  MATLAB_MEX_PATH
  MATLAB_MEXEXT_PATH
  MATLAB_MEX_EXT
)

mark_as_advanced(
  MATLAB_EXE_PATH
  MATLAB_LIBRARIES
  MATLAB_MEX_LIBRARY
  MATLAB_MX_LIBRARY
  MATLAB_ENG_LIBRARY
  MATLAB_INCLUDE_DIR
  MATLAB_FOUND
  MATLAB_MAT_LIBRARY
  MATLAB_MEX_PATH
  MATLAB_MEXEXT_PATH
  MATLAB_MEX_EXT
  MATLAB_CFLAGS
  MATLAB_CXXFLAGS
  MATLAB_FFLAGS
  MATLAB_CLINKER_FLAGS
  MATLAB_CXXLINKER_FLAGS
  MATLAB_FLINKER_FLAGS
  MATLAB_CLIBS
  MATLAB_CXXLIBS
  MATLAB_FLIBS
)
