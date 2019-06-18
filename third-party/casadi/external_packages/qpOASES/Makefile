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
##	Filename:  Makefile
##	Author:    Hans Joachim Ferreau
##	Version:   3.2
##	Date:      2007-2015
##

include make.mk

##
##	targets
##


all: src examples
#src_aw testing

src:
	@cd $@; ${MAKE} -s 

#src_aw:
#	@cd $@; ${MAKE} -s 

examples: src
	@cd $@; ${MAKE} -s

doc:
	@cd $@; ${MAKE} -s 

testing: src
	@cd testing/cpp; ${MAKE} -s

test: testing
	@cd testing/cpp; ${MAKE} -s runTests

debugging:
	@cd $@; ${MAKE} -s 

clean:
	@cd src               && ${MAKE} -s clean
	@cd examples          && ${MAKE} -s clean
	@cd bin               && ${RM} -f *.* *{EXE}
	@cd testing/cpp       && ${MAKE} -s clean

#	&& cd src_aw            && ${MAKE} -s clean && cd .. \
#	&& cd debugging         && ${MAKE} -s clean && cd .. \


clobber: clean

scilab:
	@echo Compiling Scilab interface...
	@cd ./interfaces/scilab/; ${MAKE} -s

python: all
	cd ./interfaces/python/ && python setup.py build_ext --inplace

pythoninstall: all
	cd ./interfaces/python/ && python setup.py install

c_wrapper:
	@echo Compiling C interface...
	@cd ./interfaces/c/; ${MAKE} -s

.PHONY : all src examples doc testing debugging clean clobber scilab python phythoninstall c_wrapper


##
##   end of file
##
