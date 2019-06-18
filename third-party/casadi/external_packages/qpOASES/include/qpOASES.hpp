/*
 *  This file is part of qpOASES.
 *
 *  qpOASES -- An Implementation of the Online Active Set Strategy.
 *  Copyright (C) 2007-2015 by Hans Joachim Ferreau, Andreas Potschka,
 *  Christian Kirches et al. All rights reserved.
 *
 *  qpOASES is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  qpOASES is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with qpOASES; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *  \file include/qpOASES.hpp
 *  \author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *  \version 3.2
 *  \date 2007-2015
 */


#if defined(__SINGLE_OBJECT__) || defined(__C_WRAPPER__)

#include <MessageHandling.cpp>
#include <Utils.cpp>
#include <Indexlist.cpp>
#include <SubjectTo.cpp>
#include <Bounds.cpp>
#include <Constraints.cpp>
#include <BLASReplacement.cpp>
#include <LAPACKReplacement.cpp>
#include <Matrices.cpp>
#include <Options.cpp>
#include <QProblemB.cpp>
#include <Flipper.cpp>
#include <QProblem.cpp>
#include <SQProblem.cpp>
#include <SparseSolver.cpp>
#include <SQProblemSchur.cpp>

#ifndef __C_WRAPPER__
#include <OQPinterface.cpp>
#include <SolutionAnalysis.cpp>
#endif

#else /* default compilation mode */

#include <qpOASES/QProblemB.hpp>
#include <qpOASES/QProblem.hpp>
#include <qpOASES/SQProblem.hpp>
#include <qpOASES/SQProblemSchur.hpp>
#include <qpOASES/extras/OQPinterface.hpp>
#include <qpOASES/extras/SolutionAnalysis.hpp>

#endif
