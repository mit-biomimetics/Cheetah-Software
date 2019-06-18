##
##	This file is part of qpOASES.
##
##	qpOASES -- An Implementation of the Online Active Set Strategy.
##	Copyright (C) 2007-2017 by Hans Joachim Ferreau, Andreas Potschka,
##	Christian Kirches et al. All rights reserved.
##
##	qpOASES is free software; you can redistribute it and/or
##	modify it under the terms of the GNU Lesser General Public
##	License as published by the Free Software Foundation; either
##	version 2.1 of the License, or (at your option) any later version.
##
##	qpOASES is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##	See the GNU Lesser General Public License for more details.
##
##	You should have received a copy of the GNU Lesser General Public
##	License along with qpOASES; if not, write to the Free Software
##	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
##



##
##	Filename:  make_osx.mk
##	Author:    Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
##	Version:   3.2
##	Date:      2007-2017
##

################################################################################
# user configuration

# include directories, relative
IDIR =   ${TOP}/include
SRCDIR = ${TOP}/src
BINDIR = ${TOP}/bin

# MacOSX SDK
SYSROOT = /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk
SDK = -isysroot ${SYSROOT} -stdlib=libc++

# Matlab include directory (ADAPT TO YOUR LOCAL SETTINGS!)
MATLAB_IDIR   = /Applications/MATLAB_R2013a.app/extern/include/
MATLAB_LIBDIR =

# system or replacement BLAS/LAPACK
REPLACE_LINALG = 1

ifeq ($(REPLACE_LINALG), 1)
	LIB_BLAS =   ${SRCDIR}/BLASReplacement.o
	LIB_LAPACK = ${SRCDIR}/LAPACKReplacement.o
	LA_DEPENDS = ${LIB_LAPACK} ${LIB_BLAS}
else
	LIB_BLAS =   -framework Accelerate
	LIB_LAPACK =
	LA_DEPENDS =
endif

# choice of sparse solver: NONE, MA27, or MA57
# If choice is not 'NONE', BLAS and LAPACK replacements must not be used
USE_SOLVER = NONE

ifeq ($(USE_SOLVER), MA57)
	LIB_SOLVER = /usr/local/lib/libhsl_ma57.a /usr/local/lib/libfakemetis.a
	DEF_SOLVER = SOLVER_MA57
else ifeq ($(USE_SOLVER), MA27)
	LIB_SOLVER = /usr/local/lib/libhsl_ma27.a
	DEF_SOLVER = SOLVER_MA27
else
	LIB_SOLVER =
	DEF_SOLVER = SOLVER_NONE
endif

################################################################################
# do not touch this

CPP = clang++
CC  = clang
AR  = ar
RM  = rm
F77 = gfortran
ECHO = echo
CD = cd

# file extensions
OBJEXT = o
LIBEXT = a
DLLEXT = dylib
EXE =
MEXOCTEXT = mex
DEF_TARGET = -o $@
SHARED = -dynamiclib ${SDK} -lgcc_s.10.5 -ldylib1.o

# 32 or 64 depending on target platform
BITS = $(shell getconf LONG_BIT)

# decide on MEX interface extension
ifeq ($(BITS), 32)
	MEXEXT = mexglx
else
	MEXEXT = mexa64
endif

CPPFLAGS = ${SDK} -Wall -pedantic -Wshadow -Wfloat-equal -Wconversion -Wsign-conversion -O3 -fPIC -DLINUX
#          -g -D__DEBUG__ -D__NO_COPYRIGHT__ -D__SUPPRESSANYOUTPUT__ -D__USE_SINGLE_PRECISION__

# libraries to link against when building qpOASES .so files
LINK_LIBRARIES = ${LIB_LAPACK} ${LIB_BLAS} -lm ${LIB_SOLVER}
LINK_LIBRARIES_WRAPPER = -lm ${LIB_SOLVER}

# how to link against the qpOASES shared library
QPOASES_LINK = -L${BINDIR}  -lqpOASES -L${SYSROOT}/usr/lib/System -lgcc_s.10.5 -lcrt1.o
QPOASES_LINK_WRAPPER = -L${BINDIR} -lqpOASES_wrapper

# link dependencies when creating executables
LINK_DEPENDS = ${LA_DEPENDS} ${BINDIR}/libqpOASES.${LIBEXT} ${BINDIR}/libqpOASES.${DLLEXT}
LINK_DEPENDS_WRAPPER = ${BINDIR}/libqpOASES_wrapper.${LIBEXT} ${BINDIR}/libqpOASES_wrapper.${DLLEXT}


##
##	end of file
##
