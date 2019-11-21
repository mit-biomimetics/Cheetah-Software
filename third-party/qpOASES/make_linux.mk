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
##	Filename:  make_linux.mk
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

# Matlab include directory (ADAPT TO YOUR LOCAL SETTINGS!)
#MATLAB_IDIR   = ${HOME}/Programs/matlab/extern/include/
MATLAB_IDIR   = /usr/local/matlab/extern/include/
MATLAB_LIBDIR = /usr/local/matlab/bin/glnxa64/

# system or replacement BLAS/LAPACK
REPLACE_LINALG = 1

ifeq ($(REPLACE_LINALG), 1)
	LIB_BLAS =   ${SRCDIR}/BLASReplacement.o
	LIB_LAPACK = ${SRCDIR}/LAPACKReplacement.o
else
	LIB_BLAS =   /usr/lib/libblas.so.3gf
	LIB_LAPACK = /usr/lib/liblapack.so.3gf
#	LIB_BLAS = ${MATLAB_LIBDIR}/libmwblas.so
#	LIB_LAPACK = ${MATLAB_LIBDIR}/libmwlapack.so
endif

# choice of sparse solver: NONE, MA27, or MA57
# If choice is not 'NONE', BLAS and LAPACK replacements must not be used
USE_SOLVER = NONE
#USE_SOLVER = MA57

ifeq ($(USE_SOLVER), MA57)
	LIB_SOLVER = ${MATLAB_LIBDIR}/libmwma57.so
	DEF_SOLVER = SOLVER_MA57
	LINKHSL = -Wl,-rpath=${MATLAB_LIBDIR}
else ifeq ($(USE_SOLVER), MA27)
	LIB_SOLVER = /usr/local/lib/libhsl_ma27.a
	DEF_SOLVER = SOLVER_MA27
	LINKHSL =
else
	LIB_SOLVER =
	DEF_SOLVER = SOLVER_NONE
	LINKHSL =
endif

################################################################################
# do not touch this

CPP = g++
CC  = gcc
AR  = ar
RM  = rm
F77 = gfortran
ECHO = echo
CD = cd
CP = cp

# file extensions
OBJEXT = o
LIBEXT = a
DLLEXT = so
EXE =
MEXOCTEXT = mex
DEF_TARGET = -o $@
SHARED = -shared

# 32 or 64 depending on target platform
BITS = $(shell getconf LONG_BIT)

# decide on MEX interface extension
ifeq ($(BITS), 32)
	MEXEXT = mexglx
else
	MEXEXT = mexa64
endif



CPPFLAGS = -Wall -pedantic -Wshadow -Wfloat-equal -O3 -Wconversion -Wsign-conversion -fPIC -DLINUX -D__USE_LONG_INTEGERS__ -D__USE_LONG_FINTS__ -D${DEF_SOLVER} -D__NO_COPYRIGHT__
#          -g -D__DEBUG__ -D__NO_COPYRIGHT__ -D__SUPPRESSANYOUTPUT__ -D__USE_SINGLE_PRECISION__

# libraries to link against when building qpOASES .so files
LINK_LIBRARIES = ${LIB_LAPACK} ${LIB_BLAS} -lm ${LIB_SOLVER}
LINK_LIBRARIES_WRAPPER = -lm ${LIB_SOLVER} -lstdc++

# how to link against the qpOASES shared library
QPOASES_LINK = -L${BINDIR} -Wl,-rpath=${BINDIR} ${LINKHSL} -lqpOASES
QPOASES_LINK_WRAPPER = -L${BINDIR} -Wl,-rpath=${BINDIR} ${LINKHSL} -lqpOASES_wrapper

# link dependencies when creating executables
LINK_DEPENDS = ${LIB_LAPACK} ${LIB_BLAS} ${BINDIR}/libqpOASES.${LIBEXT} ${BINDIR}/libqpOASES.${DLLEXT}
LINK_DEPENDS_WRAPPER = ${BINDIR}/libqpOASES_wrapper.${LIBEXT} ${BINDIR}/libqpOASES_wrapper.${DLLEXT}


##
##	end of file
##
