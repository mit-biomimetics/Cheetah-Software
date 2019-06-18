##
##	This file is part of qpOASES.
##
##	qpOASES -- An Implementation of the Online Active Set Strategy.
##	Copyright (C) 2007-2015 by Hans Joachim Ferreau, Andreas Potschka,
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
##	Filename:  src/Makefile
##	Author:    Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
##	Version:   3.2
##	Date:      2007-2015
##

include ../make.mk

##
##	flags
##

IFLAGS      =  -I. \
               -I${IDIR}

QPOASES_OBJECTS = \
	SQProblem.${OBJEXT} \
	QProblem.${OBJEXT} \
	QProblemB.${OBJEXT} \
	SQProblemSchur.${OBJEXT} \
	Bounds.${OBJEXT} \
	Constraints.${OBJEXT} \
	SubjectTo.${OBJEXT} \
	Indexlist.${OBJEXT} \
	Flipper.${OBJEXT} \
	Utils.${OBJEXT} \
	Options.${OBJEXT} \
	Matrices.${OBJEXT} \
	MessageHandling.${OBJEXT} \
	SparseSolver.${OBJEXT}


QPOASES_EXTRAS_OBJECTS = \
	SolutionAnalysis.${OBJEXT} \
	OQPinterface.${OBJEXT}

QPOASES_DEPENDS = \
	${IDIR}/qpOASES.hpp \
	${IDIR}/qpOASES/SQProblem.hpp \
	${IDIR}/qpOASES/QProblem.hpp \
	${IDIR}/qpOASES/Flipper.hpp \
	${IDIR}/qpOASES/QProblemB.hpp \
	${IDIR}/qpOASES/Bounds.hpp \
	${IDIR}/qpOASES/Constraints.hpp \
	${IDIR}/qpOASES/SubjectTo.hpp \
	${IDIR}/qpOASES/Indexlist.hpp \
	${IDIR}/qpOASES/Utils.hpp \
	${IDIR}/qpOASES/Constants.hpp \
	${IDIR}/qpOASES/Types.hpp \
	${IDIR}/qpOASES/Options.hpp \
	${IDIR}/qpOASES/Matrices.hpp \
	${IDIR}/qpOASES/MessageHandling.hpp \
	${IDIR}/qpOASES/UnitTesting.hpp


##
##	targets
##

all: ${LINK_DEPENDS}


${BINDIR}/libqpOASES.${LIBEXT}: ${QPOASES_OBJECTS} ${QPOASES_EXTRAS_OBJECTS}
	@${ECHO} "Creating static lib" $@
	@${AR} r $@ $^

${BINDIR}/libqpOASES.${DLLEXT}: ${QPOASES_OBJECTS} ${QPOASES_EXTRAS_OBJECTS}
	@${ECHO} "Creating shared lib" $@
	@${CPP} ${DEF_TARGET} ${SHARED} $^ ${LINK_LIBRARIES}

clean:
	@${ECHO} "Cleaning up (src)"
	@${RM} -f *.${OBJEXT} *.${LIBEXT} *.${DLLEXT}

clobber: clean


%.${OBJEXT}: %.cpp ${QPOASES_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} $<


##
##	end of file
##
