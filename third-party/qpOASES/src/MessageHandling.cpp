/*
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2017 by Hans Joachim Ferreau, Andreas Potschka,
 *	Christian Kirches et al. All rights reserved.
 *
 *	qpOASES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpOASES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpOASES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file src/MessageHandling.cpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.2
 *	\date 2007-2017
 *
 *	Implementation of the MessageHandling class including global return values.
 *
 */


#include <stdio.h>

#ifdef __MATLAB__
  #include "mex.h"
#endif

#include <qpOASES/MessageHandling.hpp>
#include <qpOASES/Utils.hpp>


BEGIN_NAMESPACE_QPOASES

/** Default file to display messages. */
FILE* stdFile = stdout;



#ifndef __SUPPRESSANYOUTPUT__

/** Defines pairs of global return values and messages. */
MessageHandling::ReturnValueList returnValueList[] =
{
/* miscellaneous */
{ SUCCESSFUL_RETURN, "Successful return", VS_VISIBLE },
{ RET_DIV_BY_ZERO, "Division by zero", VS_VISIBLE },
{ RET_INDEX_OUT_OF_BOUNDS, "Index out of bounds", VS_VISIBLE },
{ RET_INVALID_ARGUMENTS, "At least one of the arguments is invalid", VS_VISIBLE },
{ RET_ERROR_UNDEFINED, "Error number undefined", VS_VISIBLE },
{ RET_WARNING_UNDEFINED, "Warning number undefined", VS_VISIBLE },
{ RET_INFO_UNDEFINED, "Info number undefined", VS_VISIBLE },
{ RET_EWI_UNDEFINED, "Error/warning/info number undefined", VS_VISIBLE },
{ RET_AVAILABLE_WITH_LINUX_ONLY, "This function is available under Linux only", VS_HIDDEN },
{ RET_UNKNOWN_BUG, "The error occurred is not yet known", VS_VISIBLE },
{ RET_PRINTLEVEL_CHANGED, "Print level changed", VS_VISIBLE },
{ RET_NOT_YET_IMPLEMENTED, "Requested function is not yet implemented.", VS_VISIBLE },
/* Indexlist */
{ RET_INDEXLIST_MUST_BE_REORDERD, "Index list has to be reordered", VS_VISIBLE },
{ RET_INDEXLIST_EXCEEDS_MAX_LENGTH, "Index list exceeds its maximal physical length", VS_VISIBLE },
{ RET_INDEXLIST_CORRUPTED, "Index list corrupted", VS_VISIBLE },
{ RET_INDEXLIST_OUTOFBOUNDS, "Physical index is out of bounds", VS_VISIBLE },
{ RET_INDEXLIST_ADD_FAILED, "Adding indices from another index set failed", VS_VISIBLE },
{ RET_INDEXLIST_INTERSECT_FAILED, "Intersection with another index set failed", VS_VISIBLE },
/* SubjectTo / Bounds / Constraints */
{ RET_INDEX_ALREADY_OF_DESIRED_STATUS, "Index is already of desired status", VS_VISIBLE },
{ RET_ADDINDEX_FAILED, "Adding index to index set failed", VS_VISIBLE },
{ RET_REMOVEINDEX_FAILED, "Removing index from index set failed", VS_VISIBLE },
{ RET_SWAPINDEX_FAILED, "Cannot swap between different indexsets", VS_VISIBLE },
{ RET_NOTHING_TO_DO, "Nothing to do", VS_VISIBLE },
{ RET_SETUP_BOUND_FAILED, "Setting up bound index failed", VS_VISIBLE },
{ RET_SETUP_CONSTRAINT_FAILED, "Setting up constraint index failed", VS_VISIBLE },
{ RET_MOVING_BOUND_FAILED, "Moving bound between index sets failed", VS_VISIBLE },
{ RET_MOVING_CONSTRAINT_FAILED, "Moving constraint between index sets failed", VS_VISIBLE },
{ RET_SHIFTING_FAILED, "Shifting of bounds/constraints failed", VS_VISIBLE },
{ RET_ROTATING_FAILED, "Rotating of bounds/constraints failed", VS_VISIBLE },
/* QProblem */
{ RET_QPOBJECT_NOT_SETUP, "The QP object has not been setup correctly, use another constructor", VS_VISIBLE },
{ RET_QP_ALREADY_INITIALISED, "QProblem has already been initialised", VS_VISIBLE },
{ RET_NO_INIT_WITH_STANDARD_SOLVER, "Initialisation via extern QP solver is not yet implemented", VS_VISIBLE },
{ RET_RESET_FAILED, "Reset failed", VS_VISIBLE },
{ RET_INIT_FAILED, "Initialisation failed", VS_VISIBLE },
{ RET_INIT_FAILED_TQ, "Initialisation failed due to TQ factorisation", VS_VISIBLE },
{ RET_INIT_FAILED_CHOLESKY, "Initialisation failed due to Cholesky decomposition", VS_VISIBLE },
{ RET_INIT_FAILED_HOTSTART, "Initialisation failed! QP could not be solved!", VS_VISIBLE },
{ RET_INIT_FAILED_INFEASIBILITY, "Initial QP could not be solved due to infeasibility!", VS_VISIBLE },
{ RET_INIT_FAILED_UNBOUNDEDNESS, "Initial QP could not be solved due to unboundedness!", VS_VISIBLE },
{ RET_INIT_FAILED_REGULARISATION, "Initialisation failed as Hessian matrix could not be regularised", VS_VISIBLE },
{ RET_INIT_SUCCESSFUL, "Initialisation done", VS_VISIBLE },
{ RET_OBTAINING_WORKINGSET_FAILED, "Failed to obtain working set for auxiliary QP", VS_VISIBLE },
{ RET_SETUP_WORKINGSET_FAILED, "Failed to setup working set for auxiliary QP", VS_VISIBLE },
{ RET_SETUP_AUXILIARYQP_FAILED, "Failed to setup auxiliary QP for initialised homotopy", VS_VISIBLE },
{ RET_NO_CHOLESKY_WITH_INITIAL_GUESS, "Externally computed Cholesky factor cannot be combined with an initial guess", VS_VISIBLE },
{ RET_NO_EXTERN_SOLVER, "No extern QP solver available", VS_VISIBLE },
{ RET_QP_UNBOUNDED, "QP is unbounded", VS_VISIBLE },
{ RET_QP_INFEASIBLE, "QP is infeasible", VS_VISIBLE },
{ RET_QP_NOT_SOLVED, "Problems occurred while solving QP with standard solver", VS_VISIBLE },
{ RET_QP_SOLVED, "QP successfully solved", VS_VISIBLE },
{ RET_UNABLE_TO_SOLVE_QP, "Problems occurred while solving QP", VS_VISIBLE },
{ RET_INITIALISATION_STARTED, "Starting problem initialisation...", VS_VISIBLE },
{ RET_HOTSTART_FAILED, "Unable to perform homotopy due to internal error", VS_VISIBLE },
{ RET_HOTSTART_FAILED_TO_INIT, "Unable to initialise problem", VS_VISIBLE },
{ RET_HOTSTART_FAILED_AS_QP_NOT_INITIALISED, "Unable to perform homotopy as previous QP is not solved", VS_VISIBLE },
{ RET_ITERATION_STARTED, "Iteration", VS_VISIBLE },
{ RET_SHIFT_DETERMINATION_FAILED, "Determination of shift of the QP data failed", VS_VISIBLE },
{ RET_STEPDIRECTION_DETERMINATION_FAILED, "Determination of step direction failed", VS_VISIBLE },
{ RET_STEPLENGTH_DETERMINATION_FAILED, "Determination of step direction failed", VS_VISIBLE },
{ RET_OPTIMAL_SOLUTION_FOUND, "Optimal solution of neighbouring QP found", VS_VISIBLE },
{ RET_HOMOTOPY_STEP_FAILED, "Unable to perform homotopy step", VS_VISIBLE },
{ RET_HOTSTART_STOPPED_INFEASIBILITY, "Premature homotopy termination because QP is infeasible", VS_VISIBLE },
{ RET_HOTSTART_STOPPED_UNBOUNDEDNESS, "Premature homotopy termination because QP is unbounded", VS_VISIBLE },
{ RET_WORKINGSET_UPDATE_FAILED, "Unable to update working sets according to initial guesses", VS_VISIBLE },
{ RET_MAX_NWSR_REACHED, "Maximum number of working set recalculations performed", VS_VISIBLE },
{ RET_CONSTRAINTS_NOT_SPECIFIED, "Problem does comprise constraints! You have to specify new constraints' bounds", VS_VISIBLE },
{ RET_INVALID_FACTORISATION_FLAG, "Invalid factorisation flag", VS_VISIBLE },
{ RET_UNABLE_TO_SAVE_QPDATA, "Unable to save QP data", VS_VISIBLE },
{ RET_STEPDIRECTION_FAILED_TQ, "Abnormal termination due to TQ factorisation", VS_VISIBLE },
{ RET_STEPDIRECTION_FAILED_CHOLESKY, "Abnormal termination due to Cholesky factorisation", VS_VISIBLE },
{ RET_CYCLING_DETECTED, "Cycling detected", VS_VISIBLE },
{ RET_CYCLING_NOT_RESOLVED, "Cycling cannot be resolved, QP is probably infeasible", VS_VISIBLE },
{ RET_CYCLING_RESOLVED, "Cycling probably resolved", VS_VISIBLE },
{ RET_STEPSIZE, "", VS_VISIBLE },
{ RET_STEPSIZE_NONPOSITIVE, "", VS_VISIBLE },
{ RET_SETUPSUBJECTTOTYPE_FAILED, "Setup of SubjectToTypes failed", VS_VISIBLE },
{ RET_ADDCONSTRAINT_FAILED, "Addition of constraint to working set failed", VS_VISIBLE },
{ RET_ADDCONSTRAINT_FAILED_INFEASIBILITY, "Addition of constraint to working set failed", VS_VISIBLE },
{ RET_ADDBOUND_FAILED, "Addition of bound to working set failed", VS_VISIBLE },
{ RET_ADDBOUND_FAILED_INFEASIBILITY, "Addition of bound to working set failed", VS_VISIBLE },
{ RET_REMOVECONSTRAINT_FAILED, "Removal of constraint from working set failed", VS_VISIBLE },
{ RET_REMOVEBOUND_FAILED, "Removal of bound from working set failed", VS_VISIBLE },
{ RET_REMOVE_FROM_ACTIVESET, "Removing from active set:", VS_VISIBLE },
{ RET_ADD_TO_ACTIVESET, "Adding to active set:", VS_VISIBLE },
{ RET_REMOVE_FROM_ACTIVESET_FAILED, "Removing from active set failed", VS_VISIBLE },
{ RET_ADD_TO_ACTIVESET_FAILED, "Adding to active set failed", VS_VISIBLE },
{ RET_CONSTRAINT_ALREADY_ACTIVE, "Constraint is already active", VS_VISIBLE },
{ RET_ALL_CONSTRAINTS_ACTIVE, "All constraints are active, no further constraint can be added", VS_VISIBLE },
{ RET_LINEARLY_DEPENDENT, "New bound/constraint is linearly dependent", VS_VISIBLE },
{ RET_LINEARLY_INDEPENDENT, "New bound/constraint is linearly independent", VS_VISIBLE },
{ RET_LI_RESOLVED, "Linear independence of active constraint matrix successfully resolved", VS_VISIBLE },
{ RET_ENSURELI_FAILED, "Failed to ensure linear independence of active constraint matrix", VS_VISIBLE },
{ RET_ENSURELI_FAILED_TQ, "Abnormal termination due to TQ factorisation", VS_VISIBLE },
{ RET_ENSURELI_FAILED_NOINDEX, "QP is infeasible", VS_VISIBLE },
{ RET_ENSURELI_FAILED_CYCLING, "QP is infeasible", VS_VISIBLE },
{ RET_BOUND_ALREADY_ACTIVE, "Bound is already active", VS_VISIBLE },
{ RET_ALL_BOUNDS_ACTIVE, "All bounds are active, no further bound can be added", VS_VISIBLE },
{ RET_CONSTRAINT_NOT_ACTIVE, "Constraint is not active", VS_VISIBLE },
{ RET_BOUND_NOT_ACTIVE, "Bound is not active", VS_VISIBLE },
{ RET_HESSIAN_NOT_SPD, "Projected Hessian matrix not positive definite", VS_VISIBLE },
{ RET_HESSIAN_INDEFINITE, "Hessian matrix is indefinite", VS_VISIBLE },
{ RET_MATRIX_SHIFT_FAILED, "Unable to update matrices or to transform vectors", VS_VISIBLE },
{ RET_MATRIX_FACTORISATION_FAILED, "Unable to calculate new matrix factorisations", VS_VISIBLE },
{ RET_PRINT_ITERATION_FAILED, "Unable to print information on current iteration", VS_VISIBLE },
{ RET_NO_GLOBAL_MESSAGE_OUTPUTFILE, "No global message output file initialised", VS_VISIBLE },
{ RET_DISABLECONSTRAINTS_FAILED, "Unable to disable constraints", VS_VISIBLE },
{ RET_ENABLECONSTRAINTS_FAILED, "Unable to enable constraints", VS_VISIBLE },
{ RET_ALREADY_ENABLED, "Bound or constraint is already enabled", VS_VISIBLE },
{ RET_ALREADY_DISABLED, "Bound or constraint is already disabled", VS_VISIBLE },
{ RET_NO_HESSIAN_SPECIFIED, "No Hessian matrix has been specified", VS_VISIBLE },
{ RET_USING_REGULARISATION, "Using regularisation as Hessian matrix is not positive definite", VS_VISIBLE },
{ RET_EPS_MUST_BE_POSITVE, "Eps for regularisation must be sufficiently positive", VS_VISIBLE },
{ RET_REGSTEPS_MUST_BE_POSITVE, "Maximum number of regularisation steps must be non-negative", VS_VISIBLE },
{ RET_HESSIAN_ALREADY_REGULARISED, "Hessian has been already regularised", VS_VISIBLE },
{ RET_CANNOT_REGULARISE_IDENTITY, "Identity Hessian matrix cannot be regularised", VS_VISIBLE },
{ RET_CANNOT_REGULARISE_SPARSE, "Sparse matrix cannot be regularised as diagonal entry is missing", VS_VISIBLE },
{ RET_NO_REGSTEP_NWSR, "No additional regularisation step could be performed due to limits", VS_VISIBLE },
{ RET_FEWER_REGSTEPS_NWSR, "Fewer additional regularisation steps have been performed due to limits", VS_VISIBLE },
{ RET_CHOLESKY_OF_ZERO_HESSIAN, "Cholesky decomposition of (unregularised) zero Hessian matrix", VS_VISIBLE },
{ RET_ZERO_HESSIAN_ASSUMED, "Zero Hessian matrix assumed as null pointer passed without specifying hessianType", VS_VISIBLE },
{ RET_CONSTRAINTS_ARE_NOT_SCALED, "(should not be thrown, no longer in use)", VS_VISIBLE },
{ RET_INITIAL_BOUNDS_STATUS_NYI, "(should not be thrown, no longer in use)", VS_VISIBLE },
{ RET_ERROR_IN_CONSTRAINTPRODUCT, "Error in user-defined constraint product function", VS_VISIBLE },
{ RET_FIX_BOUNDS_FOR_LP, "All initial bounds must be fixed when solving an (unregularised) LP", VS_VISIBLE },
{ RET_USE_REGULARISATION_FOR_LP, "Set options.enableRegularisation=BT_TRUE for solving LPs", VS_VISIBLE },
/* SQProblem */
{ RET_UPDATEMATRICES_FAILED, "Unable to update QP matrices", VS_VISIBLE },
{ RET_UPDATEMATRICES_FAILED_AS_QP_NOT_SOLVED, "Unable to update matrices as previous QP is not solved", VS_VISIBLE },
/* Utils */
{ RET_UNABLE_TO_OPEN_FILE, "Unable to open file", VS_VISIBLE },
{ RET_UNABLE_TO_WRITE_FILE, "Unable to write into file", VS_VISIBLE },
{ RET_UNABLE_TO_READ_FILE, "Unable to read from file", VS_VISIBLE },
{ RET_FILEDATA_INCONSISTENT, "File contains inconsistent data", VS_VISIBLE },
/* Options */
{ RET_OPTIONS_ADJUSTED,	"Options needed to be adjusted for consistency reasons", VS_VISIBLE },
/* SolutionAnalysis */
{ RET_UNABLE_TO_ANALYSE_QPROBLEM, "Unable to analyse (S)QProblem(B) object", VS_VISIBLE },
/* Benchmark */
{ RET_NWSR_SET_TO_ONE, "Maximum number of working set changes was set to 1", VS_VISIBLE },
{ RET_UNABLE_TO_READ_BENCHMARK, "Unable to read benchmark data", VS_VISIBLE },
{ RET_BENCHMARK_ABORTED, "Benchmark aborted", VS_VISIBLE },
{ RET_INITIAL_QP_SOLVED, "Initial QP solved", VS_VISIBLE },
{ RET_QP_SOLUTION_STARTED, "Solving QP no.", VS_VISIBLE },
{ RET_BENCHMARK_SUCCESSFUL, "Benchmark terminated successfully", VS_VISIBLE },
/* Sparse matrices */
{ RET_NO_DIAGONAL_AVAILABLE, "Sparse matrix does not have entries on full diagonal", VS_VISIBLE },
{ RET_DIAGONAL_NOT_INITIALISED, "Diagonal data of sparse matrix has not been initialised", VS_VISIBLE },
/* Dropping of infeasible constraints */
{ RET_ENSURELI_DROPPED, "Linear independence resolved by dropping blocking constraint", VS_VISIBLE },
/* Schur complement computations */
{ RET_KKT_MATRIX_SINGULAR, "KKT matrix is singular", VS_VISIBLE },
{ RET_QR_FACTORISATION_FAILED, "QR factorization of Schur complement failed", VS_VISIBLE },
{ RET_INERTIA_CORRECTION_FAILED, "Inertia correction of KKT matrix failed", VS_VISIBLE },
{ RET_NO_SPARSE_SOLVER, "No Sparse Solver installed", VS_VISIBLE },
/* Simple exitflags */
{ RET_SIMPLE_STATUS_P1, "QP problem could not be solved within given number of iterations", VS_VISIBLE },
{ RET_SIMPLE_STATUS_P0, "QP problem solved", VS_VISIBLE },
{ RET_SIMPLE_STATUS_M1, "QP problem could not be solved due to an internal error", VS_VISIBLE },
{ RET_SIMPLE_STATUS_M2, "QP problem is infeasible (and thus could not be solved)", VS_VISIBLE },
{ RET_SIMPLE_STATUS_M3, "QP problem is unbounded (and thus could not be solved)", VS_VISIBLE },
/* IMPORTANT: Terminal list element! */
{ TERMINAL_LIST_ELEMENT, "", VS_HIDDEN }
};

#else /* __SUPPRESSANYOUTPUT__ */

MessageHandling::ReturnValueList returnValueList[1]; /* Do not use messages for embedded platforms! */

#endif /* __SUPPRESSANYOUTPUT__ */



/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/

/*
 *	M e s s a g e H a n d l i n g
 */
MessageHandling::MessageHandling( )
{
	errorVisibility   = VS_VISIBLE;
	warningVisibility = VS_VISIBLE;
	infoVisibility    = VS_VISIBLE;

	outputFile = stdFile;
	errorCount = 0;
}

/*
 *	M e s s a g e H a n d l i n g
 */
MessageHandling::MessageHandling( FILE* _outputFile )
{
	errorVisibility   = VS_VISIBLE;
	warningVisibility = VS_HIDDEN;
	infoVisibility    = VS_HIDDEN;

	outputFile = _outputFile;
	errorCount = 0;
}

/*
 *	M e s s a g e H a n d l i n g
 */
MessageHandling::MessageHandling(	VisibilityStatus _errorVisibility,
									VisibilityStatus _warningVisibility,
		 							VisibilityStatus _infoVisibility
									)
{
	errorVisibility   = _errorVisibility;
	warningVisibility = _warningVisibility;
	infoVisibility    = _infoVisibility;

	outputFile = stdFile;
	errorCount = 0;
}

/*
 *	M e s s a g e H a n d l i n g
 */
MessageHandling::MessageHandling( 	FILE* _outputFile,
									VisibilityStatus _errorVisibility,
									VisibilityStatus _warningVisibility,
		 							VisibilityStatus _infoVisibility
									)
{
	errorVisibility   = _errorVisibility;
	warningVisibility = _warningVisibility;
	infoVisibility    = _infoVisibility;

	outputFile = _outputFile;
	errorCount = 0;
}



/*
 *	M e s s a g e H a n d l i n g
 */
MessageHandling::MessageHandling( const MessageHandling& rhs )
{
	errorVisibility   = rhs.errorVisibility;
	warningVisibility = rhs.warningVisibility;
	infoVisibility    = rhs.infoVisibility;

	outputFile = rhs.outputFile;
	errorCount = rhs.errorCount;
}


/*
 *	~ M e s s a g e H a n d l i n g
 */
MessageHandling::~MessageHandling( )
{
	#ifndef __SUPPRESSANYOUTPUT__
	if ( ( outputFile != 0 ) && ( outputFile != stdout ) && ( outputFile != stderr ) )
		fclose( outputFile );
 	#endif /* __SUPPRESSANYOUTPUT__ */
}


/*
 *	o p e r a t o r =
 */
MessageHandling& MessageHandling::operator=( const MessageHandling& rhs )
{
	if ( this != &rhs )
	{
		errorVisibility   = rhs.errorVisibility;
		warningVisibility = rhs.warningVisibility;
		infoVisibility    = rhs.infoVisibility;

		outputFile = rhs.outputFile;
		errorCount = rhs.errorCount;
	}

	return *this;
}


/*
 *	t h r o w E r r o r
 */
returnValue MessageHandling::throwError(	returnValue Enumber,
											const char* additionaltext,
											const char* functionname,
											const char* filename,
											const unsigned long linenumber,
											VisibilityStatus localVisibilityStatus
											)
{
	/* consistency check */
	if ( Enumber <= SUCCESSFUL_RETURN )
		return throwError( RET_ERROR_UNDEFINED,0,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );

	/* Call to common throwMessage function if error shall be displayed. */
	if ( errorVisibility == VS_VISIBLE )
		return throwMessage( Enumber,additionaltext,functionname,filename,linenumber,localVisibilityStatus,"ERROR" );
	else
		return Enumber;
}


/*
 *	t h r o w W a r n i n g
 */
returnValue MessageHandling::throwWarning(	returnValue Wnumber,
											const char* additionaltext,
											const char* functionname,
											const char* filename,
											const unsigned long linenumber,
											VisibilityStatus localVisibilityStatus
											)
{
	/* consistency check */
  	if ( Wnumber <= SUCCESSFUL_RETURN )
		return throwError( RET_WARNING_UNDEFINED,0,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );

	/* Call to common throwMessage function if warning shall be displayed. */
	if ( warningVisibility == VS_VISIBLE )
		return throwMessage( Wnumber,additionaltext,functionname,filename,linenumber,localVisibilityStatus,"WARNING" );
  	else
  		return Wnumber;
}


/*
 *	t h r o w I n f o
 */
returnValue MessageHandling::throwInfo(	returnValue Inumber,
										const char* additionaltext,
										const char* functionname,
										const char* filename,
										const unsigned long linenumber,
										VisibilityStatus localVisibilityStatus
										)
{
	/* consistency check */
	if ( Inumber < SUCCESSFUL_RETURN )
		return throwError( RET_INFO_UNDEFINED,0,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );

	/* Call to common throwMessage function if info shall be displayed. */
	if ( infoVisibility == VS_VISIBLE )
		return throwMessage( Inumber,additionaltext,functionname,filename,linenumber,localVisibilityStatus,"INFO" );
	else
		return Inumber;
}


/*
 *	r e s e t
 */
returnValue MessageHandling::reset( )
{
	setErrorVisibilityStatus(   VS_VISIBLE );
	setWarningVisibilityStatus( VS_HIDDEN );
	setInfoVisibilityStatus(    VS_HIDDEN );

	setOutputFile( stdFile );
	setErrorCount( 0 );

	return SUCCESSFUL_RETURN;
}


/*
 *	l i s t A l l M e s s a g e s
 */
returnValue MessageHandling::listAllMessages( )
{
	#ifndef __SUPPRESSANYOUTPUT__
	int_t keypos = 0;
	char myPrintfString[MAX_STRING_LENGTH];

	/* Run through whole returnValueList and print each item. */
	while ( returnValueList[keypos].key != TERMINAL_LIST_ELEMENT )
	{
		snprintf( myPrintfString,MAX_STRING_LENGTH," %d - %s \n",(int)keypos,returnValueList[keypos].data );
		myPrintf( myPrintfString );

		++keypos;
	}
	#endif /* __SUPPRESSANYOUTPUT__ */

	return SUCCESSFUL_RETURN;
}



/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/


/*
 *	t h r o w M e s s a g e
 */
returnValue MessageHandling::throwMessage(	returnValue RETnumber,
											const char* additionaltext,
											const char* functionname,
											const char* filename,
											const unsigned long linenumber,
											VisibilityStatus localVisibilityStatus,
											const char* RETstring
											)
{
	#ifndef __SUPPRESSANYOUTPUT__

	int_t keypos = 0;
	char myPrintfString[MAX_STRING_LENGTH];

	/* 1) Determine number of whitespace for output. */
	char whitespaces[MAX_STRING_LENGTH];
	int_t numberOfWhitespaces = (errorCount-1)*2;

	if ( numberOfWhitespaces < 0 )
		numberOfWhitespaces = 0;

	if ( numberOfWhitespaces > 40 )
		numberOfWhitespaces = 40;

	if ( numberOfWhitespaces >= (int_t)MAX_STRING_LENGTH )
		numberOfWhitespaces = (int_t)MAX_STRING_LENGTH-1;

	memset( whitespaces, ' ', (size_t) numberOfWhitespaces );
	whitespaces[numberOfWhitespaces] = '\0';

	/* 2) Find error/warning/info in list. */
	while ( returnValueList[keypos].key != TERMINAL_LIST_ELEMENT )
	{
		if ( returnValueList[keypos].key == RETnumber )
			break;
		else
			++keypos;
	}

	if ( returnValueList[keypos].key == TERMINAL_LIST_ELEMENT )
	{
		throwError( RET_EWI_UNDEFINED,0,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
		return RETnumber;
	}

	/* 3) Print error/warning/info. */
	if ( ( returnValueList[keypos].globalVisibilityStatus == VS_VISIBLE ) && ( localVisibilityStatus == VS_VISIBLE ) )
	{
		if ( errorCount < 0 )
		{
			myPrintf( "\n" );
			errorCount = 0;
		}

		if ( errorCount > 0 )
		{
			snprintf( myPrintfString,MAX_STRING_LENGTH,"%s->", whitespaces );
			myPrintf( myPrintfString );
		}

		if ( additionaltext == 0 )
		{
			#ifdef __DEBUG__
			snprintf(	myPrintfString,MAX_STRING_LENGTH,"%s (%s, %s:%d): \t%s\n",
						RETstring,functionname,filename,(int_t)linenumber,returnValueList[keypos].data
						);
			#else
			snprintf(	myPrintfString,MAX_STRING_LENGTH,"%s:  %s\n",
						RETstring,returnValueList[keypos].data
						);
			#endif
			myPrintf( myPrintfString );
		}
		else
		{
			#ifdef __DEBUG__
			snprintf(	myPrintfString,MAX_STRING_LENGTH,"%s (%s, %s:%d): \t%s %s\n",
						RETstring,functionname,filename,(int_t)linenumber,returnValueList[keypos].data,additionaltext
						);
			#else
			snprintf(	myPrintfString,MAX_STRING_LENGTH,"%s:  %s %s\n",
						RETstring,returnValueList[keypos].data,additionaltext
						);
			#endif
			myPrintf( myPrintfString );
		}

		/* take care of proper indention for subsequent error messages */
		if ( RETstring[0] == 'E' )
		{
			++errorCount;
		}
		else
		{
			if ( errorCount > 0 )
				myPrintf( "\n" );
			errorCount = 0;
		}
	}

	#endif /* __SUPPRESSANYOUTPUT__ */

	return RETnumber;
}

/****************************************************************************/
/* S T A T I C  P U B L I C */
/****************************************************************************/

const char* MessageHandling::getErrorCodeMessage(	const returnValue _returnValue
													)
{
	#ifndef __SUPPRESSANYOUTPUT__

	int_t keypos = 0;
	
	/* 2) Find error/warning/info in list. */
	while ( returnValueList[keypos].key != TERMINAL_LIST_ELEMENT )
	{
		if ( returnValueList[keypos].key == _returnValue )
			break;
		else
			++keypos;
	}

	if ( returnValueList[keypos].key == TERMINAL_LIST_ELEMENT )
	{
		return "Unknown error code";
	}

	return (returnValueList[keypos].data != 0) ? returnValueList[keypos].data : "No message for this error code";
    
	#else /* __SUPPRESSANYOUTPUT__ */

	return "No message for this error code";

	#endif /* __SUPPRESSANYOUTPUT__ */
}


/*****************************************************************************
 *  G L O B A L  M E S S A G E  H A N D L E R                                *
 *****************************************************************************/


/** Global message handler for all qpOASES modules.*/
#if defined(__DSPACE__) || defined(__XPCTARGET__)
static MessageHandling globalMessageHandler( stdFile,VS_VISIBLE,VS_HIDDEN,VS_HIDDEN );
#endif


/*
 *	g e t G l o b a l M e s s a g e H a n d l e r
 */
MessageHandling* getGlobalMessageHandler( )
{
	#ifndef __DSPACE__
    #ifndef __XPCTARGET__
	static MessageHandling globalMessageHandler( stdFile,VS_VISIBLE,VS_VISIBLE,VS_VISIBLE );
	#endif /* __DSPACE__ */
    #endif /* __XPCTARGET__ */

	return &globalMessageHandler;
}


END_NAMESPACE_QPOASES


/*
 *	end of file
 */
