##
##     This file is part of qpOASES.
##
##     qpOASES -- An Implementation of the Online Active Set Strategy.
##     Copyright (C) 2007-2017 by Hans Joachim Ferreau, Andreas Potschka,
##     Christian Kirches et al. All rights reserved.
##
##     qpOASES is free software; you can redistribute it and/or
##     modify it under the terms of the GNU Lesser General Public
##     License as published by the Free Software Foundation; either
##     version 2.1 of the License, or (at your option) any later version.
##
##     qpOASES is distributed in the hope that it will be useful,
##     but WITHOUT ANY WARRANTY; without even the implied warranty of
##     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##     See the GNU Lesser General Public License for more details.
##
##     You should have received a copy of the GNU Lesser General Public
##     License along with qpOASES; if not, write to the Free Software
##     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
##



##
##     Filename:  CMakeLists.txt
##     Author:    Hans Joachim Ferreau (thanks to Milan Vukov)
##     Version:   3.2
##     Date:      2007-2017
##

cmake_minimum_required(VERSION 2.6)


PROJECT(qpOASES CXX)
SET(PACKAGE_NAME "qpOASES")
SET(PACKAGE_VERSION "3.2.0")
SET(PACKAGE_SO_VERSION "3.2")
SET(PACKAGE_DESCRIPTION "An implementation of the online active set strategy")
SET(PACKAGE_AUTHOR "Hans Joachim Ferreau, Andreas Potschka, Christian Kirches et al.")
SET(PACKAGE_MAINTAINER "Hans Joachim Ferreau, Andreas Potschka, Christian Kirches et al.")
SET(PACKAGE_URL "https://projects.coin-or.org/qpOASES")

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/libs)

#IF( NOT CMAKE_VERBOSE_MAKEFILE )
    #SET( CMAKE_VERBOSE_MAKEFILE OFF )
#ENDIF( NOT CMAKE_VERBOSE_MAKEFILE )

############################################################
#################### compiler flags ########################
############################################################

# warnings are disabled! (qpOASES has many warnings)
set(CMAKE_CXX_FLAGS "-O3 -no-pie -ggdb -w")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__NO_COPYRIGHT__ -Wno-error")
IF ( UNIX )
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wfloat-equal -Wshadow -DLINUX")
    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_DEBUG} -O3 -finline-functions")
ELSEIF( WINDOWS )
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -nologo -EHsc -DWIN32")
ENDIF()

SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D__DEBUG__")

############################################################
######################## rpath #############################
############################################################
# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/casadi")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
   set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/casadi")
endif("${isSystemDir}" STREQUAL "-1")



############################################################
#################### build and install #####################
############################################################
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)

# compile qpOASES libraries
FILE(GLOB SRC src/*.cpp)

# library
ADD_LIBRARY(qpOASES SHARED ${SRC})
INSTALL(TARGETS qpOASES
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION lib
)
SET_TARGET_PROPERTIES(
    qpOASES
    PROPERTIES
    SOVERSION ${PACKAGE_SO_VERSION}
    )

# headers
INSTALL(FILES include/qpOASES.hpp
  DESTINATION include)
INSTALL(DIRECTORY include/qpOASES
  DESTINATION include
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN "*.ipp"
  PATTERN ".svn" EXCLUDE)


##
##   end of file
##
