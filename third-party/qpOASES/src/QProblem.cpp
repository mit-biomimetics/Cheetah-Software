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
 *	\file src/QProblem.cpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.2
 *	\date 2007-2017
 *
 *	Implementation of the QProblem class which is able to use the newly
 *	developed online active set strategy for parametric quadratic programming.
 */


#include <qpOASES/QProblem.hpp>
#include <qpOASES/LapackBlasReplacement.hpp>


BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/


/*
 *	Q P r o b l e m
 */
QProblem::QProblem( ) : QProblemB( )
{
	freeConstraintMatrix = BT_FALSE;
	A = 0;

	lbA = 0;
	ubA = 0;

	sizeT = 0;
	T = 0;
	Q = 0;

	Ax = 0;
	Ax_l = 0;
	Ax_u = 0;

	constraintProduct = 0;

	tempA = 0;
	ZFR_delta_xFRz = 0;
	delta_xFRy = 0;
	delta_xFRz = 0;
	tempB = 0;
	delta_yAC_TMP = 0;
	tempC = 0;
}


/*
 *	Q P r o b l e m
 */
QProblem::QProblem( int_t _nV, int_t _nC, HessianType _hessianType, BooleanType allocDenseMats ) 
	: QProblemB( _nV,_hessianType,allocDenseMats )
{
	int_t i;

	/* consistency checks */
	if ( _nV <= 0 )
	{
		_nV = 1;
		THROWERROR( RET_INVALID_ARGUMENTS );
	}

	if ( _nC < 0 )
	{
		_nC = 0;
		THROWERROR( RET_INVALID_ARGUMENTS );
	}

	if ( _nC > 0 )
	{
		freeConstraintMatrix = BT_FALSE;
		A = 0;

		lbA = new real_t[_nC];
		for( i=0; i<_nC; ++i ) lbA[i] = 0.0;

		ubA = new real_t[_nC];
		for( i=0; i<_nC; ++i ) ubA[i] = 0.0;
	}
	else
	{
		/* prevent segmentation faults in case nC == 0
		 * (avoiding checks for A!=0 around all calls to A->... */
		freeConstraintMatrix = BT_TRUE;
		A = new DenseMatrix( );

		lbA = 0;
		ubA = 0;
	}

	constraints.init( _nC );

	delete[] y; /* y of no constraints version too short! */
	y = new real_t[_nV+_nC];
	for( i=0; i<_nV+_nC; ++i ) y[i] = 0.0;

	if (allocDenseMats == BT_TRUE)
	{
		sizeT = getMin( _nV,_nC );
		T = new real_t[sizeT*sizeT];
		Q = new real_t[_nV*_nV];
	}
	else
	{
		sizeT = 0;
		T = 0;
		Q = 0;
	}

	if ( _nC > 0 )
	{
		Ax = new real_t[_nC];
		Ax_l = new real_t[_nC];
		Ax_u = new real_t[_nC];
	}
	else
	{
		Ax = 0;
		Ax_l = 0;
		Ax_u = 0;
	}

	constraintProduct = 0;

	tempA = new real_t[_nV];			/* nFR */
	ZFR_delta_xFRz = new real_t[_nV];	/* nFR */
	delta_xFRz = new real_t[_nV];		/* nZ */

	if ( _nC > 0 )
	{
		tempB = new real_t[_nC];			/* nAC */
		delta_xFRy = new real_t[_nC];		/* nAC */
		delta_yAC_TMP = new real_t[_nC];	/* nAC */
		tempC = new real_t[_nC];			/* nAC */
	}
	else
	{
		tempB = 0;
		delta_xFRy = 0;
		delta_yAC_TMP = 0;
		tempC = 0;
	}

	flipper.init( (uint_t)_nV,(uint_t)_nC );
}


/*
 *	Q P r o b l e m
 */
QProblem::QProblem( const QProblem& rhs ) : QProblemB( rhs )
{
	freeConstraintMatrix = BT_FALSE;
	A = 0;

	copy( rhs );
}


/*
 *	~ Q P r o b l e m
 */
QProblem::~QProblem( )
{
	clear( );
}


/*
 *	o p e r a t o r =
 */
QProblem& QProblem::operator=( const QProblem& rhs )
{
	if ( this != &rhs )
	{
		clear( );
		QProblemB::operator=( rhs );
		copy( rhs );
	}

	return *this;
}


/*
 *	r e s e t
 */
returnValue QProblem::reset( )
{
	int_t i;
	int_t nV = getNV( );
	int_t nC = getNC( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );


	/* 1) Reset bounds, Cholesky decomposition and status flags. */
	if ( QProblemB::reset( ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_RESET_FAILED );

	/* 2) Reset constraints. */
	constraints.init( nC );

	/* 3) Reset TQ factorisation. */
	if ( T!=0 )
		for( i=0; i<sizeT*sizeT; ++i )
			T[i] = 0.0;

	if ( Q!=0 )
		for( i=0; i<nV*nV; ++i )
			Q[i] = 0.0;

	/* 4) Reset constraint product pointer. */
	constraintProduct = 0;

	/* 5) Reset flipper object */
	flipper.init( (uint_t)nV,(uint_t)nC );

	return SUCCESSFUL_RETURN;
}


/*
 *	i n i t
 */
returnValue QProblem::init(	SymmetricMatrix *_H, const real_t* const _g, Matrix *_A,
							const real_t* const _lb, const real_t* const _ub,
							const real_t* const _lbA, const real_t* const _ubA,
							int_t& nWSR, real_t* const cputime,
							const real_t* const xOpt, const real_t* const yOpt,
							const Bounds* const guessedBounds, const Constraints* const guessedConstraints,
							const real_t* const _R
							)
{
	int_t i;
	int_t nV = getNV( );
	int_t nC = getNC( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Consistency checks. */
	if ( isInitialised( ) == BT_TRUE )
	{
		THROWWARNING( RET_QP_ALREADY_INITIALISED );
		reset( );
	}

	if ( guessedBounds != 0 )
	{
		for( i=0; i<nV; ++i )
		{
			if ( guessedBounds->getStatus( i ) == ST_UNDEFINED )
				return THROWERROR( RET_INVALID_ARGUMENTS );
		}
	}

	if ( guessedConstraints != 0 )
	{
		for( i=0; i<nC; ++i )
			if ( guessedConstraints->getStatus( i ) == ST_UNDEFINED )
				return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	/* exclude these possibilities in order to avoid inconsistencies */
	if ( ( xOpt == 0 ) && ( yOpt != 0 ) && ( ( guessedBounds != 0 ) || ( guessedConstraints != 0 ) ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	if ( ( _R != 0 ) && ( ( xOpt != 0 ) || ( yOpt != 0 ) || ( guessedBounds != 0 ) || ( guessedConstraints != 0 ) ) )
		return THROWERROR( RET_NO_CHOLESKY_WITH_INITIAL_GUESS );

	/* 2) Setup QP data. */
	if ( setupQPdata( _H,_g,_A,_lb,_ub,_lbA,_ubA ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 3) Call to main initialisation routine. */
	return solveInitialQP( xOpt,yOpt,guessedBounds,guessedConstraints,_R, nWSR,cputime );
}


/*
 *	i n i t
 */
returnValue QProblem::init(	const real_t* const _H, const real_t* const _g, const real_t* const _A,
							const real_t* const _lb, const real_t* const _ub,
							const real_t* const _lbA, const real_t* const _ubA,
							int_t& nWSR, real_t* const cputime,
							const real_t* const xOpt, const real_t* const yOpt,
							const Bounds* const guessedBounds, const Constraints* const guessedConstraints,
							const real_t* const _R
							)
{
	int_t i;
	int_t nV = getNV( );
	int_t nC = getNC( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Consistency checks. */
	if ( isInitialised( ) == BT_TRUE )
	{
		THROWWARNING( RET_QP_ALREADY_INITIALISED );
		reset( );
	}

	if ( guessedBounds != 0 )
	{
		for( i=0; i<nV; ++i )
		{
			if ( guessedBounds->getStatus( i ) == ST_UNDEFINED )
				return THROWERROR( RET_INVALID_ARGUMENTS );
		}
	}

	if ( guessedConstraints != 0 )
	{
		for( i=0; i<nC; ++i )
			if ( guessedConstraints->getStatus( i ) == ST_UNDEFINED )
				return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	/* exclude these possibilities in order to avoid inconsistencies */
	if ( ( xOpt == 0 ) && ( yOpt != 0 ) && ( ( guessedBounds != 0 ) || ( guessedConstraints != 0 ) ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	if ( ( _R != 0 ) && ( ( xOpt != 0 ) || ( yOpt != 0 ) || ( guessedBounds != 0 ) || ( guessedConstraints != 0 ) ) )
		return THROWERROR( RET_NO_CHOLESKY_WITH_INITIAL_GUESS );

	/* 2) Setup QP data. */
	if ( setupQPdata( _H,_g,_A,_lb,_ub,_lbA,_ubA ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 3) Call to main initialisation routine. */
	return solveInitialQP( xOpt,yOpt,guessedBounds,guessedConstraints,_R, nWSR,cputime );
}


/*
 *	i n i t
 */
returnValue QProblem::init(	const char* const H_file, const char* const g_file, const char* const A_file,
							const char* const lb_file, const char* const ub_file,
							const char* const lbA_file, const char* const ubA_file,
							int_t& nWSR, real_t* const cputime,
							const real_t* const xOpt, const real_t* const yOpt,
							const Bounds* const guessedBounds, const Constraints* const guessedConstraints,
							const char* const R_file
							)
{
	int_t i;
	int_t nV = getNV( );
	int_t nC = getNC( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Consistency checks. */
	if ( isInitialised( ) == BT_TRUE )
	{
		THROWWARNING( RET_QP_ALREADY_INITIALISED );
		reset( );
	}

	if ( guessedBounds != 0 )
	{
		for( i=0; i<nV; ++i )
		{
			if ( guessedBounds->getStatus( i ) == ST_UNDEFINED )
				return THROWERROR( RET_INVALID_ARGUMENTS );
		}
	}

	if ( guessedConstraints != 0 )
	{
		for( i=0; i<nC; ++i )
			if ( guessedConstraints->getStatus( i ) == ST_UNDEFINED )
				return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	/* exclude these possibilities in order to avoid inconsistencies */
	if ( ( xOpt == 0 ) && ( yOpt != 0 ) && ( ( guessedBounds != 0 ) || ( guessedConstraints != 0 ) ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	if ( ( R_file != 0 ) && ( ( xOpt != 0 ) || ( yOpt != 0 ) || ( guessedBounds != 0 ) || ( guessedConstraints != 0 ) ) )
		return THROWERROR( RET_NO_CHOLESKY_WITH_INITIAL_GUESS );

	/* 2) Setup QP data from files. */
	if ( setupQPdataFromFile( H_file,g_file,A_file,lb_file,ub_file,lbA_file,ubA_file ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_UNABLE_TO_READ_FILE );

	if ( R_file == 0 )
	{
		/* 3) Call to main initialisation routine. */
		return solveInitialQP( xOpt,yOpt,guessedBounds,guessedConstraints,0, nWSR,cputime );
	}
	else
	{
		/* Also read Cholesky factor from file and store it directly into R [thus... */
		returnValue returnvalue = readFromFile( R, nV,nV, R_file );
		if ( returnvalue != SUCCESSFUL_RETURN )
			return THROWWARNING( returnvalue );

		/* 3) Call to main initialisation routine. ...passing R here!] */
		return solveInitialQP( xOpt,yOpt,guessedBounds,guessedConstraints,R, nWSR,cputime );
	}
}



/*
 *	h o t s t a r t
 */
returnValue QProblem::hotstart(	const real_t* const g_new,
								const real_t* const lb_new, const real_t* const ub_new,
								const real_t* const lbA_new, const real_t* const ubA_new,
								int_t& nWSR, real_t* const cputime,
								const Bounds* const guessedBounds, const Constraints* const guessedConstraints
								)
{
	int_t i, nActiveFar;
	int_t nV = getNV ();
	int_t nC = getNC ();
	real_t starttime = 0.0;
	real_t auxTime = 0.0;

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );


	/* Possibly update working sets according to guesses for working sets of bounds and constraints. */
	if ( ( guessedBounds != 0 ) || ( guessedConstraints != 0 ) )
	{
		if ( cputime != 0 )
			starttime = getCPUtime( );

		const Bounds*      actualGuessedBounds      = ( guessedBounds != 0 )      ? guessedBounds      : &bounds;
		const Constraints* actualGuessedConstraints = ( guessedConstraints != 0 ) ? guessedConstraints : &constraints;

		if ( setupAuxiliaryQP( actualGuessedBounds,actualGuessedConstraints ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		status = QPS_AUXILIARYQPSOLVED;

		/* Allow only remaining CPU time for usual hotstart. */
		if ( cputime != 0 )
		{
			auxTime = getCPUtime( ) - starttime;
			*cputime -= auxTime;
		}
	}

	returnValue returnvalue = SUCCESSFUL_RETURN;

	/* Simple check for consistency of bounds and constraints. */
	if ( areBoundsConsistent(lb_new, ub_new, lbA_new, ubA_new) != SUCCESSFUL_RETURN )
		return setInfeasibilityFlag(returnvalue,BT_TRUE);

	++count;

	int_t nWSR_max = nWSR;
	int_t nWSR_performed = 0;

	real_t cputime_remaining = INFTY, *pcputime_rem;
	real_t cputime_needed = 0.0;

	real_t farbound = options.initialFarBounds;

	/* writeQpDataIntoMatFile( "qpData.mat" ); */
	/* writeQpWorkspaceIntoMatFile( "qpWorkspace.mat" ); */

	if ( haveCholesky == BT_FALSE )
	{
		returnvalue = setupInitialCholesky( );
		if (returnvalue != SUCCESSFUL_RETURN)
			return THROWERROR(returnvalue);
	}

	BooleanType isFirstCall = BT_TRUE;

	if ( options.enableFarBounds == BT_FALSE )
	{
		/* Automatically call standard solveQP if regularisation is not active. */
		returnvalue = solveRegularisedQP(	g_new,lb_new,ub_new,lbA_new,ubA_new,
											nWSR,cputime,0,
											isFirstCall
											);
	}
	else
	{
		real_t *ub_new_far = new real_t[nV];
		real_t *lb_new_far = new real_t[nV];
		real_t *ubA_new_far = new real_t[nC];
		real_t *lbA_new_far = new real_t[nC];

		/* possibly extend initial far bounds to largest bound/constraint data */
		if (ub_new)
			for (i = 0; i < nV; i++)
				if ((ub_new[i] < INFTY) && (ub_new[i] > farbound)) farbound = ub_new[i];
		if (lb_new)
			for (i = 0; i < nV; i++)
				if ((lb_new[i] > -INFTY) && (lb_new[i] < -farbound)) farbound = -lb_new[i];
		if (ubA_new)
			for (i = 0; i < nC; i++)
				if ((ubA_new[i] < INFTY) && (ubA_new[i] > farbound)) farbound = ubA_new[i];
		if (lbA_new)
			for (i = 0; i < nC; i++)
				if ((lbA_new[i] > -INFTY) && (lbA_new[i] < -farbound)) farbound = -lbA_new[i];

		updateFarBounds(	farbound,nV+nC,
							lb_new,lb_new_far, ub_new,ub_new_far,
							lbA_new,lbA_new_far, ubA_new,ubA_new_far
							);

		for ( ;; )
		{
			nWSR = nWSR_max;
			if ( cputime != 0 )
			{
				cputime_remaining = *cputime - cputime_needed;
				pcputime_rem = &cputime_remaining;
			}
			else
				pcputime_rem = 0;

			/* Automatically call standard solveQP if regularisation is not active. */
			returnvalue = solveRegularisedQP(	g_new,lb_new_far,ub_new_far,lbA_new_far,ubA_new_far,
												nWSR,pcputime_rem,nWSR_performed,
												isFirstCall
												);

			nWSR_performed  = nWSR;
			cputime_needed += cputime_remaining;
			isFirstCall     = BT_FALSE;

			/* Check for active far-bounds and move them away */
			nActiveFar = 0;
			farbound *= options.growFarBounds;

			if ( infeasible == BT_TRUE )
			{
				if ( farbound >= INFTY )
				{
					returnvalue = RET_HOTSTART_STOPPED_INFEASIBILITY;
					break; // goto farewell;
				}

				updateFarBounds(	farbound,nV+nC,
									lb_new,lb_new_far, ub_new,ub_new_far,
									lbA_new,lbA_new_far, ubA_new,ubA_new_far
									);
			}
			else if ( status == QPS_SOLVED )
			{
				real_t tol = farbound/options.growFarBounds * options.boundTolerance;

				for ( i=0; i<nV; ++i )
				{
					if ( ( ( lb_new == 0 ) || ( lb_new_far[i] > lb_new[i] ) ) && ( getAbs ( lb_new_far[i] - x[i] ) < tol ) )
						++nActiveFar;
					if ( ( ( ub_new == 0 ) || ( ub_new_far[i] < ub_new[i] ) ) && ( getAbs ( ub_new_far[i] - x[i] ) < tol ) )
						++nActiveFar;
				}
				for ( i=0; i<nC; ++i )
				{
					if ( ( ( lbA_new == 0 ) || ( lbA_new_far[i] > lbA_new[i] ) ) && ( getAbs ( lbA_new_far[i] - Ax[i] ) < tol ) )
						++nActiveFar;
					if ( ( ( ubA_new == 0 ) || ( ubA_new_far[i] < ubA_new[i] ) ) && ( getAbs ( ubA_new_far[i] - Ax[i] ) < tol ) )
						++nActiveFar;
				}

				if ( nActiveFar == 0 )
					break;

				status = QPS_HOMOTOPYQPSOLVED;

				if ( farbound >= INFTY )
				{
					unbounded = BT_TRUE;
					returnvalue = RET_HOTSTART_STOPPED_UNBOUNDEDNESS;
					goto farewell;
				}

				updateFarBounds(	farbound,nV+nC,
									lb_new,lb_new_far, ub_new,ub_new_far,
									lbA_new,lbA_new_far, ubA_new,ubA_new_far
									);
			}
			else
			{
				/* some other error when solving QP */
				break;
			}

			/* advance ramp offset to avoid Ramping cycles */
			rampOffset++;
		}

		farewell:
			/* add time to setup auxiliary QP */
			if ( cputime != 0 )
				*cputime = cputime_needed + auxTime;
			delete[] lbA_new_far; delete[] ubA_new_far;
			delete[] lb_new_far; delete[] ub_new_far;
	}

	return ( returnvalue != SUCCESSFUL_RETURN ) ? THROWERROR( returnvalue ) : returnvalue;
}


/*
 *	h o t s t a r t
 */
returnValue QProblem::hotstart(	const char* const g_file,
								const char* const lb_file, const char* const ub_file,
								const char* const lbA_file, const char* const ubA_file,
								int_t& nWSR, real_t* const cputime,
								const Bounds* const guessedBounds, const Constraints* const guessedConstraints
								)
{
	int_t nV = getNV( );
	int_t nC = getNC( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* consistency check */
	if ( g_file == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );


	/* 1) Allocate memory (if bounds exist). */
	real_t* g_new   = new real_t[nV];
	real_t* lb_new  = ( lb_file != 0 )  ? new real_t[nV] : 0;
	real_t* ub_new  = ( ub_file != 0 )  ? new real_t[nV] : 0;
	real_t* lbA_new = ( lbA_file != 0 ) ? new real_t[nC] : 0;
	real_t* ubA_new = ( ubA_file != 0 ) ? new real_t[nC] : 0;


	/* 2) Load new QP vectors from file. */
	returnValue returnvalue;
	returnvalue = loadQPvectorsFromFile(	g_file,lb_file,ub_file,lbA_file,ubA_file,
											g_new,lb_new,ub_new,lbA_new,ubA_new
											);
	if ( returnvalue != SUCCESSFUL_RETURN )
	{
		if ( ubA_file != 0 )
			delete[] ubA_new;
		if ( lbA_file != 0 )
			delete[] lbA_new;
		if ( ub_file != 0 )
			delete[] ub_new;
		if ( lb_file != 0 )
			delete[] lb_new;
		delete[] g_new;

		return THROWERROR( RET_UNABLE_TO_READ_FILE );
	}


	/* 3) Actually perform hotstart. */
	returnvalue = hotstart(	g_new,lb_new,ub_new,lbA_new,ubA_new,
							nWSR,cputime,
							guessedBounds,guessedConstraints
							);


	/* 4) Free memory. */
	if ( ubA_file != 0 )
		delete[] ubA_new;
	if ( lbA_file != 0 )
		delete[] lbA_new;
	if ( ub_file != 0 )
		delete[] ub_new;
	if ( lb_file != 0 )
		delete[] lb_new;
	delete[] g_new;

	return returnvalue;
}


/*
 * s o l v e C u r r e n t E Q P
 */
returnValue QProblem::solveCurrentEQP(	const int_t n_rhs,
										const real_t* g_in,
										const real_t* lb_in,
										const real_t* ub_in,
										const real_t* lbA_in,
										const real_t* ubA_in,
										real_t* x_out,
										real_t* y_out
										)
{
	if ( ( x_out == 0 ) || ( y_out == 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	returnValue returnvalue = SUCCESSFUL_RETURN;
	int_t ii, jj;
	int_t nV  = getNV( );
	int_t nC  = getNC( );
	int_t nFR = getNFR( );
	int_t nFX = getNFX( );
	int_t nAC = getNAC( );

	real_t *delta_xFX = new real_t[nFX];
	real_t *delta_xFR = new real_t[nFR];
	real_t *delta_yAC = new real_t[nAC];
	real_t *delta_yFX = new real_t[nFX];

	/* 1) Determine index arrays. */
	int_t* FR_idx;
	int_t* FX_idx;
	int_t* AC_idx;

	bounds.getFree( )->getNumberArray( &FR_idx );
	bounds.getFixed( )->getNumberArray( &FX_idx );
	constraints.getActive( )->getNumberArray( &AC_idx );

	for ( ii = 0 ; ii < (nV+nC)*n_rhs; ++ii )
		y_out[ii] = 0.0;

	for ( ii = 0 ; ii < n_rhs; ++ii )
	{
		returnvalue = determineStepDirection(
			g_in, lbA_in, ubA_in, lb_in, ub_in, BT_FALSE, BT_FALSE,
			delta_xFX, delta_xFR, delta_yAC, delta_yFX );

		for ( jj = 0; jj < nFX; ++jj )
			x_out[FX_idx[jj]] = delta_xFX[jj];
		for ( jj = 0; jj < nFR; ++jj )
			x_out[FR_idx[jj]] = delta_xFR[jj];
		for ( jj = 0; jj < nFX; ++jj )
			y_out[FX_idx[jj]] = delta_yFX[jj];
		for ( jj = 0; jj < nAC; ++jj )
			y_out[nV+AC_idx[jj]] = delta_yAC[jj];

		g_in += nV;
		lb_in += nV;
		ub_in += nV;
		lbA_in += nC;
		ubA_in += nC;
		x_out += nV;
		y_out += nV+nC;
	}


	delete[] delta_yFX;
	delete[] delta_yAC;
	delete[] delta_xFR;
	delete[] delta_xFX;

	return returnvalue;
}



/*
 *	g e t W o r k i n g S e t
 */
returnValue QProblem::getWorkingSet( real_t* workingSet )
{
	int_t nV = this->getNV();

	if ( workingSet == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* At which limit are the bounds active? */
	getWorkingSetBounds( workingSet );

	/* At which limit are the contraints active? */
	getWorkingSetConstraints( &(workingSet[nV]) );

	return SUCCESSFUL_RETURN;
}


/*
 *	g e t W o r k i n g S e t B o u n d s
 */
returnValue QProblem::getWorkingSetBounds( real_t* workingSetB )
{
	return QProblemB::getWorkingSetBounds( workingSetB );
}


/*
 *	g e t W o r k i n g S e t C o n s t r a i n t s
 */
returnValue QProblem::getWorkingSetConstraints( real_t* workingSetC )
{
	int_t i;
	int_t nC = this->getNC();

	if ( workingSetC == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	for ( i=0; i<nC; ++i )
	{
		switch ( constraints.getStatus(i) )
		{
			case ST_LOWER: workingSetC[i] = -1.0; break;
			case ST_UPPER: workingSetC[i] = +1.0; break;
			default:       workingSetC[i] =  0.0; break;
		}
	}

	return SUCCESSFUL_RETURN;
}



/*
 *	g e t N Z
 */
int_t QProblem::getNZ( ) const
{
	/* nZ = nFR - nAC */
	return getNFR( ) - getNAC( );
}


/*
 *	g e t D u a l S o l u t i o n
 */
returnValue QProblem::getDualSolution( real_t* const yOpt ) const
{
	int_t i;

	for( i=0; i<getNV( )+getNC( ); ++i )
		yOpt[i] = y[i];

	/* return optimal dual solution vector
	 * only if current QP has been solved */
	if ( ( getStatus( ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( getStatus( ) == QPS_HOMOTOPYQPSOLVED )  ||
		 ( getStatus( ) == QPS_SOLVED ) )
	{
		return SUCCESSFUL_RETURN;
	}
	else
	{
		return RET_QP_NOT_SOLVED;
	}
}



/*
 *	s e t C o n s t r a i n t P r o d u c t
 */
returnValue QProblem::setConstraintProduct( ConstraintProduct* const _constraintProduct )
{
	constraintProduct = _constraintProduct;

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t P r o p e r t i e s
 */
returnValue QProblem::printProperties( )
{
	#ifndef __SUPPRESSANYOUTPUT__

	/* Do not print properties if print level is set to none! */
	if ( options.printLevel == PL_NONE )
		return SUCCESSFUL_RETURN;

	char myPrintfString[MAX_STRING_LENGTH];

	myPrintf( "\n#################   qpOASES  --  QP PROPERTIES   #################\n" );
	myPrintf( "\n" );

	/* 1) Variables properties. */
	snprintf( myPrintfString,MAX_STRING_LENGTH,  "Number of Variables: %4.1d\n",(int)getNV( ) );
	myPrintf( myPrintfString );

	if ( bounds.hasNoLower( ) == BT_TRUE )
			myPrintf( "Variables are not bounded from below.\n" );
		else
			myPrintf( "Variables are bounded from below.\n" );

	if ( bounds.hasNoUpper( ) == BT_TRUE )
			myPrintf( "Variables are not bounded from above.\n" );
		else
			myPrintf( "Variables are bounded from above.\n" );

	myPrintf( "\n" );


	/* 2) Constraints properties. */
	snprintf( myPrintfString,MAX_STRING_LENGTH,  "Total number of Constraints:      %4.1d\n",(int)getNC( ) );
	myPrintf( myPrintfString );

	snprintf( myPrintfString,MAX_STRING_LENGTH,  "Number of Equality Constraints:   %4.1d\n",(int)getNEC( ) );
	myPrintf( myPrintfString );

	snprintf( myPrintfString,MAX_STRING_LENGTH,  "Number of Inequality Constraints: %4.1d\n",(int)(getNC( )-getNEC( )) );
	myPrintf( myPrintfString );

	if ( getNC( ) > 0 )
	{
		if ( constraints.hasNoLower( ) == BT_TRUE )
				myPrintf( "Constraints are not bounded from below.\n" );
			else
				myPrintf( "Constraints are bounded from below.\n" );

		if ( constraints.hasNoUpper( ) == BT_TRUE )
				myPrintf( "Constraints are not bounded from above.\n" );
			else
				myPrintf( "Constraints are bounded from above.\n" );
	}

	myPrintf( "\n" );


	/* 3) Further properties. */
	switch ( hessianType )
	{
		case HST_ZERO:
			myPrintf( "Hessian is zero matrix (i.e. actually an LP is solved).\n" );
			break;

		case HST_IDENTITY:
			myPrintf( "Hessian is identity matrix.\n" );
			break;

		case HST_POSDEF:
			myPrintf( "Hessian matrix is (strictly) positive definite.\n" );
			break;

		case HST_POSDEF_NULLSPACE:
			myPrintf( "Hessian matrix is positive definite on null space of active constraints.\n" );
			break;

		case HST_SEMIDEF:
			myPrintf( "Hessian matrix is positive semi-definite.\n" );
			break;

		case HST_INDEF:
			myPrintf( "Hessian matrix is indefinite.\n" );
			break;

		default:
			myPrintf( "Hessian matrix has unknown type.\n" );
			break;
	}

	if ( infeasible == BT_TRUE )
		myPrintf( "QP was found to be infeasible.\n" );
	else
		myPrintf( "QP seems to be feasible.\n" );

	if ( unbounded == BT_TRUE )
		myPrintf( "QP was found to be unbounded from below.\n" );
	else
		myPrintf( "QP seems to be bounded from below.\n" );

	myPrintf( "\n" );


	/* 4) QP object properties. */
	switch ( status )
	{
		case QPS_NOTINITIALISED:
			myPrintf( "Status of QP object: freshly instantiated or reset.\n" );
			break;

		case QPS_PREPARINGAUXILIARYQP:
			myPrintf( "Status of QP object: an auxiliary QP is currently setup.\n" );
			break;

		case QPS_AUXILIARYQPSOLVED:
			myPrintf( "Status of QP object: an auxilary QP was solved.\n" );
			break;

		case QPS_PERFORMINGHOMOTOPY:
			myPrintf( "Status of QP object: a homotopy step is performed.\n" );
			break;

		case QPS_HOMOTOPYQPSOLVED:
			myPrintf( "Status of QP object: an intermediate QP along the homotopy path was solved.\n" );
			break;

		case QPS_SOLVED:
			myPrintf( "Status of QP object: solution of the actual QP was found.\n" );
			break;
	}

	switch ( options.printLevel )
	{
		case PL_DEBUG_ITER:
			myPrintf( "Print level of QP object is set to display a tabular output for debugging.\n" );
			break;

		case PL_TABULAR:
			myPrintf( "Print level of QP object is set to display a tabular output.\n" );
			break;

		case PL_LOW:
			myPrintf( "Print level of QP object is low, i.e. only error are printed.\n" );
			break;

		case PL_MEDIUM:
			myPrintf( "Print level of QP object is medium, i.e. error and warnings are printed.\n" );
			break;

		case PL_HIGH:
			myPrintf( "Print level of QP object is high, i.e. all available output is printed.\n" );
			break;

		default:
			break;
	}

	myPrintf( "\n" );

	#endif /* __SUPPRESSANYOUTPUT__ */

	return SUCCESSFUL_RETURN;
}

returnValue QProblem::getFreeVariablesFlags( BooleanType* varIsFree )
{
	int_t nV  = getNV( );
	for ( int_t i=0; i<nV; i++ )
		varIsFree[i] = BT_FALSE;

	int_t nFR  = getNFR( );
	int_t* FR_idx;
	bounds.getFree( )->getNumberArray( &FR_idx );

	for ( int_t i=0; i<nFR; i++ )
		varIsFree[FR_idx[i]] = BT_TRUE;

	return SUCCESSFUL_RETURN;
}


/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/

/*
 *	c l e a r
 */
returnValue QProblem::clear( )
{
	if ( ( freeConstraintMatrix == BT_TRUE ) && ( A != 0 ) )
	{
		delete A;
		A = 0;
	}

	if ( lbA != 0 )
	{
		delete[] lbA;
		lbA = 0;
	}

	if ( ubA != 0 )
	{
		delete[] ubA;
		ubA = 0;
	}

	if ( T != 0 )
	{
		delete[] T;
		T = 0;
	}

	if ( Q != 0 )
	{
		delete[] Q;
		Q = 0;
	}

	if ( Ax != 0 )
	{
		delete[] Ax;
		Ax = 0;
	}

	if ( Ax_l != 0 )
	{
		delete[] Ax_l;
		Ax_l = 0;
	}

	if ( Ax_u != 0 )
	{
		delete[] Ax_u;
		Ax_u = 0;
	}

	if ( tempA != 0 )
	{
		delete[] tempA;
		tempA = 0;
	}

	if ( ZFR_delta_xFRz != 0 )
	{
		delete[] ZFR_delta_xFRz;
		ZFR_delta_xFRz = 0;
	}

	if ( delta_xFRy != 0 )
	{
		delete[] delta_xFRy;
		delta_xFRy = 0;
	}

	if ( delta_xFRz != 0 )
	{
		delete[] delta_xFRz;
		delta_xFRz = 0;
	}

	if ( tempB != 0 )
	{
		delete[] tempB;
		tempB = 0;
	}

	if ( delta_yAC_TMP != 0 )
	{
		delete[] delta_yAC_TMP;
		delta_yAC_TMP = 0;
	}

	if ( tempC != 0 )
	{
		delete[] tempC;
		tempC = 0;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	c o p y
 */
returnValue QProblem::copy(	const QProblem& rhs
							)
{
	uint_t _nV = (uint_t)rhs.getNV( );
	uint_t _nC = (uint_t)rhs.getNC( );

	constraints = rhs.constraints;

	if ( ( freeConstraintMatrix == BT_TRUE ) && ( A != 0 ) )
	{
		delete A;
		A = 0;
	}

	freeConstraintMatrix = rhs.freeConstraintMatrix;

	if ( freeConstraintMatrix == BT_TRUE )
		A = rhs.A->duplicate();
	else
		A = rhs.A;

	if ( rhs.lbA != 0 )
	{
		lbA = new real_t[_nC];
		setLBA( rhs.lbA );
	}
	else
		lbA = 0;

	if ( rhs.ubA != 0 )
	{
		ubA = new real_t[_nC];
		setUBA( rhs.ubA );
	}
	else
		ubA = 0;

	if ( rhs.y != 0 )
	{
		delete[] y; /* y of no constraints version too short! */
		y = new real_t[_nV+_nC];
		memcpy( y,rhs.y,(_nV+_nC)*sizeof(real_t) );
	}
	else
		y = 0;

	sizeT = rhs.sizeT;

	if ( rhs.T != 0 )
	{
		T = new real_t[sizeT*sizeT];
		memcpy( T,rhs.T,((uint_t)(sizeT*sizeT))*sizeof(real_t) );
	}
	else
		T = 0;

	if ( rhs.Q != 0 )
	{
		Q = new real_t[_nV*_nV];
		memcpy( Q,rhs.Q,_nV*_nV*sizeof(real_t) );
	}
	else
		Q = 0;

	if ( rhs.Ax != 0 )
	{
		Ax = new real_t[_nC];
		memcpy( Ax,rhs.Ax,_nC*sizeof(real_t) );
	}
	else
		Ax = 0;

	if ( rhs.Ax_l != 0 )
	{
		Ax_l = new real_t[_nC];
		memcpy( Ax_l,rhs.Ax_l,_nC*sizeof(real_t) );
	}
	else
		Ax_l = 0;

	if ( rhs.Ax_u != 0 )
	{
		Ax_u = new real_t[_nC];
		memcpy( Ax_u,rhs.Ax_u,_nC*sizeof(real_t) );
	}
	else
		Ax_u = 0;

	if ( rhs.constraintProduct != 0 )
		constraintProduct = rhs.constraintProduct;
	else
		constraintProduct = 0;

	tempA = new real_t[_nV];			/* nFR */
	ZFR_delta_xFRz = new real_t[_nV];	/* nFR */
	delta_xFRz = new real_t[_nV];		/* nZ */

	if ( _nC > 0 )
	{
		delta_xFRy = new real_t[_nC];		/* nAC */
		tempB = new real_t[_nC];			/* nAC */
		delta_yAC_TMP = new real_t[_nC];    /* nAC */
		tempC = new real_t[_nC];            /* nAC */
	}
	else
	{
		delta_xFRy = 0;
		tempB = 0;
		delta_yAC_TMP = 0;
		tempC = 0;
	}

	return SUCCESSFUL_RETURN;
}



/*
 *	s o l v e I n i t i a l Q P
 */
returnValue QProblem::solveInitialQP(	const real_t* const xOpt, const real_t* const yOpt,
										const Bounds* const guessedBounds, const Constraints* const guessedConstraints,
										const real_t* const _R,
										int_t& nWSR, real_t* const cputime
										)
{
	int_t i,j;

	/* some definitions */
	int_t nV = getNV( );
	int_t nC = getNC( );

	//writeQpDataIntoMatFile( "qpData.mat" );

	/* start runtime measurement */
	real_t starttime = 0.0;
	if ( cputime != 0 )
		starttime = getCPUtime( );

	status = QPS_NOTINITIALISED;

	/* I) ANALYSE QP DATA: */
	/* 1) Check if Hessian happens to be the identity matrix. */
	if ( determineHessianType( ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* 2) Setup type of bounds and constraints (i.e. unbounded, implicitly fixed etc.). */
	if ( setupSubjectToType( ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	status = QPS_PREPARINGAUXILIARYQP;


	/* II) SETUP AUXILIARY QP WITH GIVEN OPTIMAL SOLUTION: */
	/* 1) Setup bounds and constraints data structure. */
	if ( bounds.setupAllFree( ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	if ( constraints.setupAllInactive( ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* 2) Setup optimal primal/dual solution for auxiliary QP. */
	if ( setupAuxiliaryQPsolution( xOpt,yOpt ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* 3) Obtain linear independent working set for auxiliary QP. */
	Bounds auxiliaryBounds( nV );
	Constraints auxiliaryConstraints( nC );

	if ( obtainAuxiliaryWorkingSet(	xOpt,yOpt,guessedBounds,guessedConstraints,
									&auxiliaryBounds,&auxiliaryConstraints ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* 4) Setup working set of auxiliary QP and setup matrix factorisations. */
	/* a) Regularise Hessian if necessary. */
	if ( ( hessianType == HST_ZERO ) || ( hessianType == HST_SEMIDEF ) )
	{
		if ( regulariseHessian( ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_INIT_FAILED_REGULARISATION );
	}

	/* b) TQ factorisation. */
	if ( setupTQfactorisation( ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED_TQ );

	/* c) Working set of auxiliary QP. */
	if ( setupAuxiliaryWorkingSet( &auxiliaryBounds,&auxiliaryConstraints,BT_TRUE ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* d) Copy external Cholesky factor if provided */
	haveCholesky = BT_FALSE;

	if ( _R != 0 )
	{
		if ( options.initialStatusBounds != ST_INACTIVE )
		{
			THROWWARNING( RET_NO_CHOLESKY_WITH_INITIAL_GUESS );
		}
		else
		{
			if ( _R == R )
			{
				/* Cholesky factor read from file and already loaded into R. */
				haveCholesky = BT_TRUE;
			}
			else if ( ( xOpt == 0 ) && ( yOpt == 0 ) && ( guessedBounds == 0 ) && ( guessedConstraints == 0 ) )
			{
				for( i=0; i<nV; ++i )
					for( j=i; j<nV; ++j )
						RR(i,j) = _R[i*nV+j];
				haveCholesky = BT_TRUE;
			}
			else
			{
				THROWWARNING( RET_NO_CHOLESKY_WITH_INITIAL_GUESS );
			}
		}
	}

	/* 5) Store original QP formulation... */
	real_t* g_original = new real_t[nV];
	real_t* lb_original = new real_t[nV];
	real_t* ub_original = new real_t[nV];
	real_t* lbA_original = new real_t[nC];
	real_t* ubA_original = new real_t[nC];

	for( i=0; i<nV; ++i )
	{
		g_original[i] = g[i];
		lb_original[i] = lb[i];
		ub_original[i] = ub[i];
	}

	for( i=0; i<nC; ++i )
	{
		lbA_original[i] = lbA[i];
		ubA_original[i] = ubA[i];
	}

	/* ... and setup QP data of an auxiliary QP having an optimal solution
	 * as specified by the user (or xOpt = yOpt = 0, by default). */
	if ( setupAuxiliaryQPgradient( ) != SUCCESSFUL_RETURN )
	{
		delete[] ubA_original; delete[] lbA_original; delete[] ub_original; delete[] lb_original; delete[] g_original;
		return THROWERROR( RET_INIT_FAILED );
	}

	if ( setupAuxiliaryQPbounds( &auxiliaryBounds,&auxiliaryConstraints,BT_TRUE ) != SUCCESSFUL_RETURN )
	{
		delete[] ubA_original; delete[] lbA_original; delete[] ub_original; delete[] lb_original; delete[] g_original;
		return THROWERROR( RET_INIT_FAILED );
	}

	status = QPS_AUXILIARYQPSOLVED;


	if ( options.enableRamping == BT_TRUE )
		performRamping( );


	/* III) SOLVE ACTUAL INITIAL QP: */
	/* Allow only remaining CPU time for usual hotstart. */
	if ( cputime != 0 )
		*cputime -= getCPUtime( ) - starttime;

	/* Use hotstart method to find the solution of the original initial QP,... */
	returnValue returnvalue = hotstart( g_original,lb_original,ub_original,lbA_original,ubA_original, nWSR,cputime );

	/* ... deallocate memory,... */
	delete[] ubA_original; delete[] lbA_original; delete[] ub_original; delete[] lb_original; delete[] g_original;

	/* ... check for infeasibility and unboundedness... */
	if ( isInfeasible( ) == BT_TRUE )
		return THROWERROR( RET_INIT_FAILED_INFEASIBILITY );

	if ( isUnbounded( ) == BT_TRUE )
		return THROWERROR( RET_INIT_FAILED_UNBOUNDEDNESS );

	/* ... and internal errors. */
	if ( ( returnvalue != SUCCESSFUL_RETURN ) && ( returnvalue != RET_MAX_NWSR_REACHED ) )
		return THROWERROR( RET_INIT_FAILED_HOTSTART );


	/* stop runtime measurement */
	if ( cputime != 0 )
		*cputime = getCPUtime( ) - starttime;

	THROWINFO( RET_INIT_SUCCESSFUL );

	return returnvalue;
}


/*
 *	s o l v e Q P
 */
returnValue QProblem::solveQP(	const real_t* const g_new,
								const real_t* const lb_new, const real_t* const ub_new,
								const real_t* const lbA_new, const real_t* const ubA_new,
								int_t& nWSR, real_t* const cputime, int_t nWSRperformed,
								BooleanType isFirstCall
								)
{
	int_t iter;
	int_t nV  = getNV( );
	int_t nC  = getNC( );

	returnValue returnvalue;

	/* consistency check */
	if ( ( getStatus( ) == QPS_NOTINITIALISED )       ||
		 ( getStatus( ) == QPS_PREPARINGAUXILIARYQP ) ||
		 ( getStatus( ) == QPS_PERFORMINGHOMOTOPY )   )
	{
		return THROWERROR( RET_HOTSTART_FAILED_AS_QP_NOT_INITIALISED );
	}

	/* start runtime measurement */
	real_t starttime = 0.0;
	if ( cputime != 0 )
		starttime = getCPUtime( );

	/* AW: Remove bounds if they were active before but are now infinity */
	status = QPS_PERFORMINGHOMOTOPY; // AW TODO: Not sure if this is too early, but otherwise removeBounds will complain
	returnvalue = updateActivitiesForHotstart( lb_new, ub_new, lbA_new, ubA_new );
	if ( returnvalue != SUCCESSFUL_RETURN )
	{
		THROWERROR( RET_HOTSTART_FAILED );
		return returnvalue;
	}

	/* I) PREPARATIONS */
	/* 1) Allocate delta vectors of gradient and (constraints') bounds,
	 *    index arrays and step direction arrays. */
	real_t* delta_xFR = new real_t[nV];
	real_t* delta_xFX = new real_t[nV];
	real_t* delta_yAC = new real_t[nC];
	real_t* delta_yFX = new real_t[nV];

	real_t* delta_g   = new real_t[nV];
	real_t* delta_lb  = new real_t[nV];
	real_t* delta_ub  = new real_t[nV];
	real_t* delta_lbA = new real_t[nC];
	real_t* delta_ubA = new real_t[nC];

	BooleanType Delta_bC_isZero, Delta_bB_isZero;

	int_t BC_idx;
	SubjectToStatus BC_status;
	BooleanType BC_isBound;

	real_t homotopyLength;

	#ifndef __SUPPRESSANYOUTPUT__
	char messageString[MAX_STRING_LENGTH];
	#endif


	/* 2) Update type of bounds and constraints, e.g.
	 *    a former equality constraint might have become a normal one etc. */

  // (ckirches) disabled this, as inactive but tight bounds may become inactive equalities
    //            which would then never become active again!
/*
	if ( setupSubjectToType( lb_new,ub_new,lbA_new,ubA_new ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_HOTSTART_FAILED );
*/

	/* 3) Reset status flags. */
	infeasible = BT_FALSE;
	unbounded  = BT_FALSE;


	/* II) MAIN HOMOTOPY LOOP */
	for( iter=nWSRperformed; iter<nWSR; ++iter )
	{
		tabularOutput.idxAddB = tabularOutput.idxRemB = tabularOutput.idxAddC = tabularOutput.idxRemC = -1;
		tabularOutput.excAddB = tabularOutput.excRemB = tabularOutput.excAddC = tabularOutput.excRemC = 0;

		if ( isCPUtimeLimitExceeded( cputime,starttime,iter-nWSRperformed ) == BT_TRUE )
		{
			/* If CPU time limit is exceeded, stop homotopy loop immediately!
			* Assign number of working set recalculations (runtime measurement is stopped later). */
			nWSR = iter;
			break;
		}

		status = QPS_PERFORMINGHOMOTOPY;

		#ifndef __SUPPRESSANYOUTPUT__
		if ( isFirstCall == BT_TRUE )
			snprintf( messageString,MAX_STRING_LENGTH,"%d ...",(int)iter );
		else
			snprintf( messageString,MAX_STRING_LENGTH,"%d* ...",(int)iter );
		getGlobalMessageHandler( )->throwInfo( RET_ITERATION_STARTED,messageString,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
		#endif

		/* 2) Determination of shift direction of the gradient and the (constraints') bounds. */
		returnvalue = determineDataShift(	g_new,lbA_new,ubA_new,lb_new,ub_new,
											delta_g,delta_lbA,delta_ubA,delta_lb,delta_ub,
											Delta_bC_isZero, Delta_bB_isZero
											);
		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			delete[] delta_yAC; delete[] delta_yFX; delete[] delta_xFX; delete[] delta_xFR;
			delete[] delta_ub; delete[] delta_lb; delete[] delta_ubA; delete[] delta_lbA; delete[] delta_g;

			/* Assign number of working set recalculations and stop runtime measurement. */
			nWSR = iter;
			if ( cputime != 0 )
				*cputime = getCPUtime( ) - starttime;

			THROWERROR( RET_SHIFT_DETERMINATION_FAILED );
			return returnvalue;
		}

		/* 3) Determination of step direction of X and Y. */
		returnvalue = determineStepDirection(	delta_g,delta_lbA,delta_ubA,delta_lb,delta_ub,
												Delta_bC_isZero, Delta_bB_isZero,
												delta_xFX,delta_xFR,delta_yAC,delta_yFX
												);
		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			delete[] delta_yAC; delete[] delta_yFX; delete[] delta_xFX; delete[] delta_xFR;
			delete[] delta_ub; delete[] delta_lb; delete[] delta_ubA; delete[] delta_lbA; delete[] delta_g;

			/* Assign number of working set recalculations and stop runtime measurement. */
			nWSR = iter;
			if ( cputime != 0 )
				*cputime = getCPUtime( ) - starttime;

			THROWERROR( RET_STEPDIRECTION_DETERMINATION_FAILED );
			return returnvalue;
		}

		/* 4) Determination of step length TAU.
		 *    This step along the homotopy path is also taken (without changing working set). */
		returnvalue = performStep(	delta_g, delta_lbA,delta_ubA,delta_lb,delta_ub,
									delta_xFX,delta_xFR,delta_yAC,delta_yFX,
									BC_idx,BC_status,BC_isBound
									);
		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			delete[] delta_yAC; delete[] delta_yFX; delete[] delta_xFX; delete[] delta_xFR;
			delete[] delta_ub; delete[] delta_lb; delete[] delta_ubA; delete[] delta_lbA; delete[] delta_g;

			/* Assign number of working set recalculations and stop runtime measurement. */
			nWSR = iter;
			if ( cputime != 0 )
				*cputime = getCPUtime( ) - starttime;

			THROWERROR( RET_STEPLENGTH_DETERMINATION_FAILED );
			return returnvalue;
		}

		/* 5) Termination criterion. */
		nV = getNV( );
		nC = getNC( );

		homotopyLength = getRelativeHomotopyLength( g_new,lb_new,ub_new,lbA_new,ubA_new );
		if ( homotopyLength <= options.terminationTolerance )
		{
			status = QPS_SOLVED;

			THROWINFO( RET_OPTIMAL_SOLUTION_FOUND );

			if ( printIteration( iter,BC_idx,BC_status,BC_isBound,homotopyLength,isFirstCall ) != SUCCESSFUL_RETURN )
				THROWERROR( RET_PRINT_ITERATION_FAILED ); /* do not pass this as return value! */

			nWSR = iter;
			if ( cputime != 0 )
				*cputime = getCPUtime( ) - starttime;

			delete[] delta_yAC; delete[] delta_yFX; delete[] delta_xFX; delete[] delta_xFR;
			delete[] delta_ub; delete[] delta_lb; delete[] delta_ubA; delete[] delta_lbA; delete[] delta_g;

			return SUCCESSFUL_RETURN;
		}

		/* 6) Change active set. */
		returnvalue = changeActiveSet( BC_idx,BC_status,BC_isBound );
		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			delete[] delta_yAC; delete[] delta_yFX; delete[] delta_xFX; delete[] delta_xFR;
			delete[] delta_ub; delete[] delta_lb; delete[] delta_ubA; delete[] delta_lbA; delete[] delta_g;

			/* Assign number of working set recalculations and stop runtime measurement. */
			nWSR = iter;
			if ( cputime != 0 )
				*cputime = getCPUtime( ) - starttime;

			/* Checks for infeasibility... */
			if ( isInfeasible( ) == BT_TRUE )
			{
				status = QPS_HOMOTOPYQPSOLVED;
				return setInfeasibilityFlag( RET_HOTSTART_STOPPED_INFEASIBILITY );
			}

			/* ...unboundedness... */
			if ( unbounded == BT_TRUE ) /* not necessary since objective function convex! */
				return THROWERROR( RET_HOTSTART_STOPPED_UNBOUNDEDNESS );

			/* ... and throw unspecific error otherwise */
			THROWERROR( RET_HOMOTOPY_STEP_FAILED );
			return returnvalue;
		}

		/* 6a) Possibly refactorise projected Hessian from scratch. */
		if ( ( options.enableCholeskyRefactorisation > 0 ) && ( (iter % options.enableCholeskyRefactorisation) == 0 ) )
		{
			returnvalue = computeProjectedCholesky( );
			if (returnvalue != SUCCESSFUL_RETURN)
			{
				delete[] delta_yAC; delete[] delta_yFX; delete[] delta_xFX; delete[] delta_xFR;
				delete[] delta_ub; delete[] delta_lb; delete[] delta_ubA; delete[] delta_lbA; delete[] delta_g;
				return returnvalue;
			}
		}

		/* 7) Output information of successful QP iteration. */
		status = QPS_HOMOTOPYQPSOLVED;

		if ( printIteration( iter,BC_idx,BC_status,BC_isBound,homotopyLength,isFirstCall ) != SUCCESSFUL_RETURN )
			THROWERROR( RET_PRINT_ITERATION_FAILED ); /* do not pass this as return value! */

		/* 8) Perform Ramping Strategy on zero homotopy step or drift correction (if desired). */
		if (BC_status != ST_UNDEFINED)
		{
			if ( ( tau <= EPS ) && ( options.enableRamping == BT_TRUE ) )
				performRamping( );
			else
			if ( (options.enableDriftCorrection > 0)
			  && ((iter+1) % options.enableDriftCorrection == 0) )
				performDriftCorrection( );  /* always returns SUCCESSFUL_RETURN */
		}
		else // AW: Added this.  Otherwise, I observed that the gradient might become incorrect
		{
			if ( (options.enableDriftCorrection > 0)
			  && ((iter+1) % options.enableDriftCorrection == 0) )
				performDriftCorrection( );  /* always returns SUCCESSFUL_RETURN */
		}
	}

	delete[] delta_yAC; delete[] delta_yFX; delete[] delta_xFX; delete[] delta_xFR;
	delete[] delta_ub; delete[] delta_lb; delete[] delta_ubA; delete[] delta_lbA; delete[] delta_g;

	/* stop runtime measurement */
	if ( cputime != 0 )
		*cputime = getCPUtime( ) - starttime;


	/* if program gets to here, output information that QP could not be solved
	 * within the given maximum numbers of working set changes */
	if ( options.printLevel == PL_HIGH )
	{
		#ifndef __SUPPRESSANYOUTPUT__
		snprintf( messageString,MAX_STRING_LENGTH,"(nWSR = %d)",(int)iter );
		return getGlobalMessageHandler( )->throwWarning( RET_MAX_NWSR_REACHED,messageString,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
		#else
		return RET_MAX_NWSR_REACHED;
		#endif
	}
	else
	{
		return RET_MAX_NWSR_REACHED;
	}
}


/*
 *	s o l v e R e g u l a r i s e d Q P
 */
returnValue QProblem::solveRegularisedQP(	const real_t* const g_new,
											const real_t* const lb_new, const real_t* const ub_new,
											const real_t* const lbA_new, const real_t* const ubA_new,
											int_t& nWSR, real_t* const cputime, int_t nWSRperformed,
											BooleanType isFirstCall
											)
{
	int_t i, step;
	int_t nV = getNV( );


	/* Perform normal QP solution if QP has not been regularised. */
	if ( usingRegularisation( ) == BT_FALSE )
		return solveQP( g_new,lb_new,ub_new,lbA_new,ubA_new, nWSR,cputime,nWSRperformed,isFirstCall );


	/* I) SOLVE USUAL REGULARISED QP */
	returnValue returnvalue;

	int_t nWSR_max   = nWSR;
	int_t nWSR_total = nWSRperformed;

	real_t cputime_total = 0.0;
	real_t cputime_cur   = 0.0;

	if ( cputime == 0 )
	{
		returnvalue = solveQP( g_new,lb_new,ub_new,lbA_new,ubA_new, nWSR,0,nWSRperformed,isFirstCall );
	}
	else
	{
		cputime_cur = *cputime;
		returnvalue = solveQP( g_new,lb_new,ub_new,lbA_new,ubA_new, nWSR,&cputime_cur,nWSRperformed,isFirstCall );
	}
	nWSR_total     = nWSR;
	cputime_total += cputime_cur;
	isFirstCall    = BT_FALSE;

	/* Only continue if QP solution has been successful. */
	if ( returnvalue != SUCCESSFUL_RETURN )
	{
		if ( cputime != 0 )
			*cputime = cputime_total;

		if ( returnvalue == RET_MAX_NWSR_REACHED )
			THROWWARNING( RET_NO_REGSTEP_NWSR );

		return returnvalue;
	}


	/* II) PERFORM SUCCESSIVE REGULARISATION STEPS */
	real_t* gMod = new real_t[nV];

	for( step=0; step<options.numRegularisationSteps; ++step )
	{
		/* 1) Modify gradient: gMod = g - eps*xOpt
		 *    (assuming regularisation matrix to be regVal*Id). */
		for( i=0; i<nV; ++i )
			gMod[i] = g_new[i] - regVal*x[i];

		/* 2) Solve regularised QP with modified gradient allowing
		 *    only as many working set recalculations and CPU time
		 *    as have been left from previous QP solutions. */
		nWSR = nWSR_max;

		if ( cputime == 0 )
		{
			returnvalue = solveQP( gMod,lb_new,ub_new,lbA_new,ubA_new, nWSR,0,nWSR_total,isFirstCall );
		}
		else
		{
			cputime_cur = *cputime - cputime_total;
			returnvalue = solveQP( gMod,lb_new,ub_new,lbA_new,ubA_new, nWSR,&cputime_cur,nWSR_total,isFirstCall );
		}

		nWSR_total     = nWSR;
		cputime_total += cputime_cur;

		/* Only continue if QP solution has been successful. */
		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			delete[] gMod;

			if ( cputime != 0 )
				*cputime = cputime_total;

			if ( returnvalue == RET_MAX_NWSR_REACHED )
				THROWWARNING( RET_FEWER_REGSTEPS_NWSR );

			return returnvalue;
		}
	}

	for( i=0; i<nV; ++i )
		g[i] = g_new[i];

	delete[] gMod;

	if ( cputime != 0 )
		*cputime = cputime_total;

	return SUCCESSFUL_RETURN;
}

/*
 *	s e t u p S u b j e c t T o T y p e
 */
returnValue QProblem::updateActivitiesForHotstart( const real_t* const lb_new, const real_t* const ub_new,
												   const real_t* const lbA_new, const real_t* const ubA_new
												   )
{
	int_t i;
	int_t nV = getNV( );

	returnValue returnvalue;

	if ( QProblemB::setupSubjectToType( lb_new,ub_new ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_SETUPSUBJECTTOTYPE_FAILED );

	for ( i=0; i<nV; i++ )
	{
		if ( lb_new[i] <= -INFTY && bounds.getStatus(i) == ST_LOWER )
		{
			returnvalue = removeBound( i, BT_TRUE, BT_FALSE, options.enableNZCTests );
			if ( returnvalue != SUCCESSFUL_RETURN )
				return returnvalue;
			g[i] -= y[i];
			y[i] = 0.0;
		}
		if ( ub_new[i] >= INFTY && bounds.getStatus(i) == ST_UPPER )
		{
			returnvalue = removeBound( i, BT_TRUE, BT_FALSE, options.enableNZCTests );
			if ( returnvalue != SUCCESSFUL_RETURN )
				return returnvalue;
			g[i] -= y[i];
			y[i] = 0.0;
		}
		if ( lb_new[i] > -INFTY && lb[i] <= -INFTY )
		{
			/* Now a lower bound has become finite.  To avoid numerical issues, adjust lb */
			lb[i] = x[i] - options.boundRelaxation;
		}
		if ( ub_new[i] < INFTY && ub[i] >= INFTY )
		{
			/* Now a lower bound has become finite.  To avoid numerical issues, adjust lb */
			ub[i] = x[i] + options.boundRelaxation;
		}
	}

	for ( i=0; i<nV; i++ )
	{
	  if ( bounds.getType(i) == ST_EQUALITY ) // ?? && lb[i] != ub[i]
		{
			/* AW: Are the following two lines OK? */
			lb[i] = x[i];
			ub[i] = x[i];
			if (  bounds.getStatus(i) == ST_INACTIVE )
			{
				returnvalue = addBound_checkLI(i);
				if ( returnvalue == RET_LINEARLY_INDEPENDENT )
				{
					returnvalue = addBound( i,ST_LOWER, BT_TRUE );
					if ( returnvalue != SUCCESSFUL_RETURN )
						return returnvalue;
				}
				/* AW: Check: This allows to have variables that are
				   equalities to be in the set of free variables. */
			}
		}
	}


	// AW TODO: We could also implement something here for the constraints

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p S u b j e c t T o T y p e
 */
returnValue QProblem::setupSubjectToType( )
{
	return setupSubjectToType( lb,ub,lbA,ubA );
}


/*
 *	s e t u p S u b j e c t T o T y p e
 */
returnValue QProblem::setupSubjectToType(	const real_t* const lb_new, const real_t* const ub_new,
											const real_t* const lbA_new, const real_t* const ubA_new
											)
{
	int_t i;
	int_t nC = getNC( );


	/* I) SETUP SUBJECTTOTYPE FOR BOUNDS */
	if ( QProblemB::setupSubjectToType( lb_new,ub_new ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_SETUPSUBJECTTOTYPE_FAILED );


	/* II) SETUP SUBJECTTOTYPE FOR CONSTRAINTS */
	/* 1) Check if lower constraints' bounds are present. */
	constraints.setNoLower( BT_TRUE );
	if ( lbA_new != 0 )
	{
		for( i=0; i<nC; ++i )
		{
			if ( lbA_new[i] > -INFTY )
			{
				constraints.setNoLower( BT_FALSE );
				break;
			}
		}
	}

	/* 2) Check if upper constraints' bounds are present. */
	constraints.setNoUpper( BT_TRUE );
	if ( ubA_new != 0 )
	{
		for( i=0; i<nC; ++i )
		{
			if ( ubA_new[i] < INFTY )
			{
				constraints.setNoUpper( BT_FALSE );
				break;
			}
		}
	}

	/* 3) Determine implicit equality constraints and unbounded constraints. */
	if ( ( lbA_new != 0 ) && ( ubA_new != 0 ) )
	{
		for( i=0; i<nC; ++i )
		{
			if ( constraints.getType(i) == ST_DISABLED )
				continue;

			if ( ( lbA_new[i] < -INFTY+options.boundTolerance ) && ( ubA_new[i] > INFTY-options.boundTolerance )
					&& (options.enableFarBounds == BT_FALSE))
			{
				constraints.setType( i,ST_UNBOUNDED );
			}
			else
			{
				if ( options.enableEqualities && lbA[i] > ubA[i] - options.boundTolerance
				                              && lbA_new[i] > ubA_new[i] - options.boundTolerance)
					constraints.setType( i,ST_EQUALITY );
				else
					constraints.setType( i,ST_BOUNDED );
			}
		}
	}
	else
	{
		if ( ( lbA_new == 0 ) && ( ubA_new == 0 ) )
		{
			for( i=0; i<nC; ++i )
			{
				if ( constraints.getType(i) != ST_DISABLED )
					constraints.setType( i,ST_UNBOUNDED );
			}
		}
		else
		{
			for( i=0; i<nC; ++i )
			{
				if ( constraints.getType(i) != ST_DISABLED )
					constraints.setType( i,ST_BOUNDED );
			}
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	c o m p u t e P r o j e c t e d C h o l e s k y
 */
returnValue QProblem::computeProjectedCholesky( )
{
	int_t i, j;
	int_t nV  = getNV( );
	int_t nZ  = getNZ( );

	SymSparseMat* Id;

	/* Revert to unprotected Cholesky decomposition */
	if ( getNFX() + getNAC() == 0 )
		return QProblemB::computeCholesky( );

	/* 1) Initialises R with all zeros. */
	for( i=0; i<nV*nV; ++i )
		R[i] = 0.0;

	/* Do not do anything for empty null spaces (important for LP case, HST_ZERO !)*/
	if ( nZ == 0 ) // nZ == nV - getNFX() - getNAC()
		return SUCCESSFUL_RETURN;

	/* 2) Calculate Cholesky decomposition of projected Hessian Z'*H*Z. */
	int_t* FR_idx;
	bounds.getFree( )->getNumberArray( &FR_idx );

	int_t* AC_idx;
	constraints.getActive( )->getNumberArray( &AC_idx );

	/* calculate Z'*H*Z */
	switch ( hessianType )
	{
		case HST_ZERO:
			if ( usingRegularisation() == BT_TRUE )
			{
				Id = createDiagSparseMat( nV, regVal );
				Id->bilinear(bounds.getFree(), nZ, Q, nV, R, nV);
				delete Id;
			}
			else
			{
				/* Code should not get here, as  nZ == 0  always holds for an LP (without regularisation)! */
				if ( nZ > 0 )
					return THROWERROR( RET_UNKNOWN_BUG );
			}
			break;

		case HST_IDENTITY:
			Id = createDiagSparseMat( nV, 1.0 );
			Id->bilinear(bounds.getFree(), nZ, Q, nV, R, nV);
			delete Id;
			break;

		default:
			if ( getNAC() == 0 ) {
				/* make Z trivial */
				for ( j=0; j < nZ; ++j ) {
					for ( i=0; i < nV; ++i )
						QQ(i,j) = 0.0;
					QQ(FR_idx[j],j) = 1.0;
				}
				/* now Z is trivial, and so is Z'HZ */
				int_t nFR = getNFR ();
				for ( j=0; j < nFR; ++j )
					H->getCol (FR_idx[j], bounds.getFree (), 1.0, &R[j*nV]);
			} else {
				/* this is expensive if Z is large! */
				H->bilinear(bounds.getFree(), nZ, Q, nV, R, nV);
			}
	}

	/* R'*R = Z'*H*Z */
	la_int_t info = 0;
	la_uint_t _nZ = (la_uint_t)nZ, _nV = (la_uint_t)nV;

	POTRF( "U", &_nZ, R, &_nV, &info );

	/* <0 = invalid call, =0 ok, >0 not spd */
	if (info > 0) {
		if ( R[0] < 0.0 )
		{
			/* Cholesky decomposition has tunneled a negative
			 * diagonal element. */
			options.epsRegularisation = getMin( -R[0]+options.epsRegularisation,getSqrt(getAbs(options.epsRegularisation)) );
		}

		hessianType = HST_SEMIDEF;
		return RET_HESSIAN_NOT_SPD;
	}

	/* zero first subdiagonal to make givens updates work */
	for (i=0;i<nZ-1;++i)
		RR(i+1,i) = 0.0;

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p I n i t i a l C h o l e s k y
 */
returnValue QProblem::setupInitialCholesky( )
{
	returnValue returnvalueCholesky;

	/* If regularisation shall be used, always regularise at beginning
	 * if initial working set is not empty. */
	if ( ( getNV() != getNFR()-getNFV() ) && ( options.enableRegularisation == BT_TRUE ) )
		if ( regulariseHessian( ) != SUCCESSFUL_RETURN )
			return RET_INIT_FAILED_REGULARISATION;

	/* Factorise projected Hessian
	 * now handles all special cases (no active bounds/constraints, no nullspace) */
	returnvalueCholesky = computeProjectedCholesky( );

	/* If Hessian is not positive definite, regularise and try again. */
	if ( returnvalueCholesky == RET_HESSIAN_NOT_SPD )
	{
		if ( regulariseHessian( ) != SUCCESSFUL_RETURN )
			return RET_INIT_FAILED_REGULARISATION;

		returnvalueCholesky = computeProjectedCholesky( );
	}

	if ( returnvalueCholesky != SUCCESSFUL_RETURN )
		return RET_INIT_FAILED_CHOLESKY;

	haveCholesky = BT_TRUE;
	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p T Q f a c t o r i s a t i o n
 */
returnValue QProblem::setupTQfactorisation( )
{
	int_t i, ii;
	int_t nV  = getNV( );
	int_t nFR = getNFR( );

	int_t* FR_idx;
	bounds.getFree( )->getNumberArray( &FR_idx );

	/* 1) Set Q to unity matrix. */
	for( i=0; i<nV*nV; ++i )
		Q[i] = 0.0;

	for( i=0; i<nFR; ++i )
	{
		ii = FR_idx[i];
		QQ(ii,i) = 1.0;
	}

 	/* 2) Set T to zero matrix. */
	for( i=0; i<sizeT*sizeT; ++i )
		T[i] = 0.0;

	return SUCCESSFUL_RETURN;
}


/*
 *	o b t a i n A u x i l i a r y W o r k i n g S e t
 */
returnValue QProblem::obtainAuxiliaryWorkingSet(	const real_t* const xOpt, const real_t* const yOpt,
													const Bounds* const guessedBounds, const Constraints* const guessedConstraints,
													Bounds* auxiliaryBounds, Constraints* auxiliaryConstraints
													) const
{
	int_t i = 0;
	int_t nV = getNV( );
	int_t nC = getNC( );


	/* 1) Ensure that desiredBounds is allocated (and different from guessedBounds). */
	if ( ( auxiliaryBounds == 0 ) || ( auxiliaryBounds == guessedBounds ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	if ( ( auxiliaryConstraints == 0 ) || ( auxiliaryConstraints == guessedConstraints ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );


	SubjectToStatus guessedStatus;

	/* 2) Setup working set of bounds for auxiliary initial QP. */
	if ( QProblemB::obtainAuxiliaryWorkingSet( xOpt,yOpt,guessedBounds, auxiliaryBounds ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );

	/* 3) Setup working set of constraints for auxiliary initial QP. */
	if ( guessedConstraints != 0 )
	{
		/* If an initial working set is specific, use it!
		 * Moreover, add all equality constraints if specified. */
		for( i=0; i<nC; ++i )
		{
			/* Add constraint only if it is not (going to be) disabled! */
			guessedStatus = guessedConstraints->getStatus( i );

			#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
			if ( constraints.getType( i ) == ST_EQUALITY )
			{
				if ( auxiliaryConstraints->setupConstraint( i,ST_LOWER ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
			}
			else
			#endif
			{
				if ( auxiliaryConstraints->setupConstraint( i,guessedStatus ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
			}
		}
	}
	else	/* No initial working set specified. */
	{
		/* Obtain initial working set by "clipping". */
		if ( ( xOpt != 0 ) && ( yOpt == 0 ) )
		{
			for( i=0; i<nC; ++i )
			{
				if ( Ax[i] - lbA[i] <= options.boundTolerance )
				{
					if ( auxiliaryConstraints->setupConstraint( i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
					continue;
				}

				if ( ubA[i] - Ax_u[i] <= options.boundTolerance )
				{
					if ( auxiliaryConstraints->setupConstraint( i,ST_UPPER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
					continue;
				}

				/* Moreover, add all equality constraints if specified. */
				#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
				if ( constraints.getType( i ) == ST_EQUALITY )
				{
					if ( auxiliaryConstraints->setupConstraint( i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
				else
				#endif
				{
					if ( auxiliaryConstraints->setupConstraint( i,ST_INACTIVE ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
			}
		}

		/* Obtain initial working set in accordance to sign of dual solution vector. */
		if ( yOpt != 0 )
		{
			for( i=0; i<nC; ++i )
			{
				if ( yOpt[nV+i] > EPS )
				{
					if ( auxiliaryConstraints->setupConstraint( i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
					continue;
				}

				if ( yOpt[nV+i] < -EPS )
				{
					if ( auxiliaryConstraints->setupConstraint( i,ST_UPPER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
					continue;
				}

				/* Moreover, add all equality constraints if specified. */
				#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
				if ( constraints.getType( i ) == ST_EQUALITY )
				{
					if ( auxiliaryConstraints->setupConstraint( i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
				else
				#endif
				{
					if ( auxiliaryConstraints->setupConstraint( i,ST_INACTIVE ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
			}
		}

		/* If xOpt and yOpt are null pointer and no initial working is specified,
		 * start with empty working set (or implicitly fixed bounds and equality constraints only)
		 * for auxiliary QP. */
		if ( ( xOpt == 0 ) && ( yOpt == 0 ) )
		{
			for( i=0; i<nC; ++i )
			{
				/* Only add all equality constraints if specified. */
				#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
				if ( constraints.getType( i ) == ST_EQUALITY )
				{
					if ( auxiliaryConstraints->setupConstraint( i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
				else
				#endif
				{
					if ( auxiliaryConstraints->setupConstraint( i,ST_INACTIVE ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
			}
		}
	}

	return SUCCESSFUL_RETURN;
}



/*
 *	s e t u p A u x i l i a r y W o r k i n g S e t
 */
returnValue QProblem::setupAuxiliaryWorkingSet(	const Bounds* const auxiliaryBounds,
												const Constraints* const auxiliaryConstraints,
												BooleanType setupAfresh
												)
{
	int_t i;
	int_t nV = getNV( );
	int_t nC = getNC( );
	BooleanType WSisTrivial = BT_TRUE;

	/* consistency checks */
	if ( auxiliaryBounds != 0 )
	{
		for( i=0; i<nV; ++i )
			if ( ( bounds.getStatus( i ) == ST_UNDEFINED ) || ( auxiliaryBounds->getStatus( i ) == ST_UNDEFINED ) )
				return THROWERROR( RET_UNKNOWN_BUG );
	}
	else
	{
		return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	if ( auxiliaryConstraints != 0 )
	{
		for( i=0; i<nC; ++i )
			if ( ( constraints.getStatus( i ) == ST_UNDEFINED ) || ( auxiliaryConstraints->getStatus( i ) == ST_UNDEFINED ) )
				return THROWERROR( RET_UNKNOWN_BUG );
	}
	else
	{
		return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	/* Check for trivial working set (all and only bounds active) */
	for (i = 0; i < nV; i++)
		if (auxiliaryBounds->getStatus(i) == ST_INACTIVE)
		{
			WSisTrivial = BT_FALSE;
			break;
		}
	for (i = 0; i < nC; i++)
		// (ckirches) here we chose to ignore an invalid ST_INACTIVE on
		//            constraints that are ST_EQUALITies or may just have become equalities
		if ( (constraints.getType(i) == ST_EQUALITY) // NOT auxiliaryConstraints here
			|| (auxiliaryConstraints->getStatus(i) != ST_INACTIVE) )
		{
			WSisTrivial = BT_FALSE;
			break;
		}

	if (WSisTrivial == BT_TRUE)
	{
		for (i = 0; i < nV; i++)
			if (bounds.getStatus(i) == ST_INACTIVE)
				bounds.moveFreeToFixed(i, auxiliaryBounds->getStatus(i));

		return SUCCESSFUL_RETURN;
	}


	/* I) SETUP CHOLESKY FLAG:
	 *    Cholesky decomposition shall only be updated if working set
	 *    shall be updated (i.e. NOT setup afresh!) */
	BooleanType updateCholesky;
	if ( setupAfresh == BT_TRUE )
		updateCholesky = BT_FALSE;
	else
		updateCholesky = BT_TRUE;


	BooleanType was_fulli = options.enableFullLITests;
	real_t backupEpsLITests = options.epsLITests;

	options.enableFullLITests = BT_FALSE;
	/* options.epsLITests = 1e-1; */

	/* II) REMOVE FORMERLY ACTIVE (CONSTRAINTS') BOUNDS (IF NECESSARY): */
	if ( setupAfresh == BT_FALSE )
	{
		/* 1) Remove all active constraints that shall be inactive or disabled AND
		*    all active constraints that are active at the wrong bound. */
		for( i=0; i<nC; ++i )
		{
			if ( ( constraints.getStatus( i ) == ST_LOWER ) && ( auxiliaryConstraints->getStatus( i ) != ST_LOWER ) )
				if ( removeConstraint( i,updateCholesky,BT_FALSE,options.enableNZCTests ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SETUP_WORKINGSET_FAILED );

			if ( ( constraints.getStatus( i ) == ST_UPPER ) && ( auxiliaryConstraints->getStatus( i ) != ST_UPPER ) )
				if ( removeConstraint( i,updateCholesky,BT_FALSE,options.enableNZCTests ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
		}

		/* 2) Remove all active bounds that shall be inactive AND
		*    all active bounds that are active at the wrong bound. */
		for( i=0; i<nV; ++i )
		{
			if ( ( bounds.getStatus( i ) == ST_LOWER ) && ( auxiliaryBounds->getStatus( i ) != ST_LOWER ) )
				if ( removeBound( i,updateCholesky,BT_FALSE,options.enableNZCTests ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SETUP_WORKINGSET_FAILED );

			if ( ( bounds.getStatus( i ) == ST_UPPER ) && ( auxiliaryBounds->getStatus( i ) != ST_UPPER ) )
				if ( removeBound( i,updateCholesky,BT_FALSE,options.enableNZCTests ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
		}
	}


	/* III) ADD NEWLY ACTIVE (CONSTRAINTS') BOUNDS: */

	/* 1) Add all equality bounds. */
	for( i=0; i<nV; ++i )
	{
		//if ( ( bounds.getType( i ) == ST_EQUALITY ) && ( ( bounds.getStatus( i ) == ST_INACTIVE ) && ( auxiliaryBounds->getStatus( i ) != ST_INACTIVE ) ) )

		// (ckirches) force equalities active

		if ( ( bounds.getType( i ) == ST_EQUALITY ) && ( bounds.getStatus( i ) == ST_INACTIVE ) )
		{
            // assert ( auxiliaryBounds->getStatus( i ) != ST_INACTIVE );
			/* No check for linear independence necessary. */
			if ( addBound( i,ST_LOWER,updateCholesky ) != SUCCESSFUL_RETURN ) // was auxiliaryBounds->getStatus( i )
				return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
		}
	}

	/* 2) Add all equality constraints. */
	for( i=0; i<nC; ++i )
	{
        //if ( ( constraints.getType( i ) == ST_EQUALITY ) && ( ( constraints.getStatus( i ) == ST_INACTIVE ) && ( auxiliaryConstraints->getStatus( i ) != ST_INACTIVE ) ) )

		// (ckirches) force equalities active

		if ( ( constraints.getType( i ) == ST_EQUALITY ) && ( constraints.getStatus( i ) == ST_INACTIVE ) )
		{
            // assert ( auxiliaryConstraints->getStatus( i ) != ST_INACTIVE );
			/* Add constraint only if it is linearly independent from the current working set. */
			if ( addConstraint_checkLI( i ) == RET_LINEARLY_INDEPENDENT )
			{
				if ( addConstraint( i,ST_LOWER,updateCholesky ) != SUCCESSFUL_RETURN )  // was auxiliaryConstraints->getStatus( i )
					return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
			}
			else
			{
				/* Equalities are not linearly independent! */
				constraints.setType(i, ST_BOUNDED);
			}
		}
	}


	/* 3) Add all inactive bounds that shall be active AND
	 *    all formerly active bounds that have been active at the wrong bound. */
	for( i=0; i<nV; ++i )
	{
		if ( ( bounds.getType( i ) != ST_EQUALITY ) && ( ( bounds.getStatus( i ) == ST_INACTIVE ) && ( auxiliaryBounds->getStatus( i ) != ST_INACTIVE ) ) )
		{
			/* Add bound only if it is linearly independent from the current working set. */
			if ( addBound_checkLI( i ) == RET_LINEARLY_INDEPENDENT )
			{
				if ( addBound( i,auxiliaryBounds->getStatus( i ),updateCholesky ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
			}
		}
	}

	/* 4) Add all inactive constraints that shall be active AND
	 *    all formerly active constraints that have been active at the wrong bound. */
	for( i=0; i<nC; ++i )
	{
		if ( ( constraints.getType( i ) != ST_EQUALITY ) && ( auxiliaryConstraints->getStatus( i ) != ST_INACTIVE ) )
		{
			/* formerly inactive */
			if ( constraints.getStatus( i ) == ST_INACTIVE )
			{
				/* Add constraint only if it is linearly independent from the current working set. */
				if ( addConstraint_checkLI( i ) == RET_LINEARLY_INDEPENDENT )
				{
					if ( addConstraint( i,auxiliaryConstraints->getStatus( i ),updateCholesky ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
				}
			}
		}
	}

	options.enableFullLITests = was_fulli;
	options.epsLITests = backupEpsLITests;

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y Q P s o l u t i o n
 */
returnValue QProblem::setupAuxiliaryQPsolution(	const real_t* const xOpt, const real_t* const yOpt
												)
{
	int_t i, j;
	int_t nV = getNV( );
	int_t nC = getNC( );


	/* Setup primal/dual solution vector for auxiliary initial QP:
	 * if a null pointer is passed, a zero vector is assigned;
	 *  old solution vector is kept if pointer to internal solution vector is passed. */
	if ( xOpt != 0 )
	{
		if ( xOpt != x )
			for( i=0; i<nV; ++i )
				x[i] = xOpt[i];

		A->times(1, 1.0, x, nV, 0.0, Ax, nC);

		for ( j=0; j<nC; ++j )
		{
			Ax_l[j] = Ax[j];
			Ax_u[j] = Ax[j];
		}
	}
	else
	{
		for( i=0; i<nV; ++i )
			x[i] = 0.0;

		for ( j=0; j<nC; ++j )
		{
			Ax[j] = 0.0;
			Ax_l[j] = 0.0;
			Ax_u[j] = 0.0;
		}
	}

	if ( yOpt != 0 )
	{
		if ( yOpt != y )
			for( i=0; i<nV+nC; ++i )
				y[i] = yOpt[i];
	}
	else
	{
		for( i=0; i<nV+nC; ++i )
			y[i] = 0.0;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y Q P g r a d i e n t
 */
returnValue QProblem::setupAuxiliaryQPgradient( )
{
	int_t i;
	int_t nV = getNV( );
	int_t nC = getNC( );


	/* Setup gradient vector: g = -H*x + [Id A]'*[yB yC]
	 *                          = yB - H*x + A'*yC. */
	switch ( hessianType )
	{
		case HST_ZERO:
			if ( usingRegularisation( ) == BT_FALSE )
				for ( i=0; i<nV; ++i )
					g[i] = y[i];
			else
				for ( i=0; i<nV; ++i )
					g[i] = y[i] - regVal*x[i];
			break;

		case HST_IDENTITY:
			for ( i=0; i<nV; ++i )
				g[i] = y[i] - x[i];
			break;

		default:
			/* y'*Id */
			for ( i=0; i<nV; ++i )
				g[i] = y[i];

			/* - H*x */
			H->times(1, -1.0, x, nV, 1.0, g, nV);
			break;
	}

	/* + A'*yC */
	A->transTimes(1, 1.0, y + nV, nC, 1.0, g, nV);

	return SUCCESSFUL_RETURN;
}


/*
 *	a r e B o u n d s C o n s i s t e n t
 */
returnValue QProblem::areBoundsConsistent(	const real_t* const lb_new, const real_t* const ub_new,
											const real_t* const lbA_new, const real_t* const ubA_new) const
{
	if (QProblemB::areBoundsConsistent(lb_new, ub_new) == RET_QP_INFEASIBLE)
		return RET_QP_INFEASIBLE;

	if (lbA_new && ubA_new) {
		for (int_t i = 0; i < getNC(); ++i) {
			if (lbA_new[i] > ubA_new[i]+EPS) {
				return RET_QP_INFEASIBLE;
			}
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y Q P b o u n d s
 */
returnValue QProblem::setupAuxiliaryQPbounds(	const Bounds* const auxiliaryBounds,
												const Constraints* const auxiliaryConstraints,
												BooleanType useRelaxation
												)
{
	int_t i;
	int_t nV = getNV( );
	int_t nC = getNC( );


	/* 1) Setup bound vectors. */
	for ( i=0; i<nV; ++i )
	{
		switch ( bounds.getStatus( i ) )
		{
			case ST_INACTIVE:
				if ( useRelaxation == BT_TRUE )
				{
					if ( bounds.getType( i ) == ST_EQUALITY )
					{
						lb[i] = x[i];
						ub[i] = x[i];
					}
					else
					{
						/* If a bound is inactive although it was supposed to be
						* active by the auxiliaryBounds, it could not be added
						* due to linear dependence. Thus set it "strongly inactive". */
						if ( auxiliaryBounds->getStatus( i ) == ST_LOWER )
							lb[i] = x[i];
						else
							lb[i] = x[i] - options.boundRelaxation;

						if ( auxiliaryBounds->getStatus( i ) == ST_UPPER )
							ub[i] = x[i];
						else
							ub[i] = x[i] + options.boundRelaxation;
					}
				}
				break;

			case ST_LOWER:
				lb[i] = x[i];
				if ( bounds.getType( i ) == ST_EQUALITY )
				{
					ub[i] = x[i];
				}
				else
				{
					if ( useRelaxation == BT_TRUE )
						ub[i] = x[i] + options.boundRelaxation;
				}
				break;

			case ST_UPPER:
				ub[i] = x[i];
				if ( bounds.getType( i ) == ST_EQUALITY )
				{
					lb[i] = x[i];
				}
				else
				{
					if ( useRelaxation == BT_TRUE )
						lb[i] = x[i] - options.boundRelaxation;
				}
				break;

            case ST_INFEASIBLE_LOWER:
			case ST_INFEASIBLE_UPPER:
                break;

			default:
				return THROWERROR( RET_UNKNOWN_BUG );
		}
	}

	/* 2) Setup constraints vectors. */
	for ( i=0; i<nC; ++i )
	{
		switch ( constraints.getStatus( i ) )
		{
			case ST_INACTIVE:
				if ( useRelaxation == BT_TRUE )
				{
					if ( constraints.getType( i ) == ST_EQUALITY )
					{
						lbA[i] = Ax_l[i];
						ubA[i] = Ax_u[i];
					}
					else
					{
						/* If a constraint is inactive although it was supposed to be
						* active by the auxiliaryConstraints, it could not be added
						* due to linear dependence. Thus set it "strongly inactive". */
						if ( auxiliaryConstraints->getStatus( i ) == ST_LOWER )
							lbA[i] = Ax_l[i];
						else
							lbA[i] = Ax_l[i] - options.boundRelaxation;

						if ( auxiliaryConstraints->getStatus( i ) == ST_UPPER )
							ubA[i] = Ax_u[i];
						else
							ubA[i] = Ax_u[i] + options.boundRelaxation;
					}
				}
				break;

			case ST_LOWER:
				lbA[i] = Ax_l[i];
				if ( constraints.getType( i ) == ST_EQUALITY )
				{
					ubA[i] = Ax_l[i];
				}
				else
				{
					if ( useRelaxation == BT_TRUE )
						ubA[i] = Ax_l[i] + options.boundRelaxation;
				}
				break;

			case ST_UPPER:
				ubA[i] = Ax_u[i];
				if ( constraints.getType( i ) == ST_EQUALITY )
				{
					lbA[i] = Ax_u[i];
				}
				else
				{
					if ( useRelaxation == BT_TRUE )
						lbA[i] = Ax_u[i] - options.boundRelaxation;
				}
				break;

            case ST_INFEASIBLE_LOWER:
			case ST_INFEASIBLE_UPPER:
                break;

			default:
				return THROWERROR( RET_UNKNOWN_BUG );
		}
		Ax_l[i] = Ax_l[i] - lbA[i];
		Ax_u[i] = ubA[i] - Ax_u[i];
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	a d d C o n s t r a i n t
 */
returnValue QProblem::addConstraint(	int_t number, SubjectToStatus C_status,
										BooleanType updateCholesky,
										BooleanType ensureLI
										)
{
	int_t i, j, ii;

	/* consistency checks */
	if ( constraints.getStatus( number ) != ST_INACTIVE )
		return THROWERROR( RET_CONSTRAINT_ALREADY_ACTIVE );

	if ( ( constraints.getNC( ) - getNAC( ) ) == constraints.getNUC( ) )
		return THROWERROR( RET_ALL_CONSTRAINTS_ACTIVE );

	if ( ( getStatus( ) == QPS_NOTINITIALISED )    ||
		 ( getStatus( ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( getStatus( ) == QPS_HOMOTOPYQPSOLVED )  ||
		 ( getStatus( ) == QPS_SOLVED )            )
	{
		return THROWERROR( RET_UNKNOWN_BUG );
	}


	/* I) ENSURE LINEAR INDEPENDENCE OF THE WORKING SET,
	 *    i.e. remove a constraint or bound if linear dependence occurs. */
	/* check for LI only if Cholesky decomposition shall be updated! */
	if ( updateCholesky == BT_TRUE && ensureLI == BT_TRUE )
	{
		returnValue ensureLIreturnvalue = addConstraint_ensureLI( number,C_status );

		switch ( ensureLIreturnvalue )
		{
			case SUCCESSFUL_RETURN:
				break;

			case RET_LI_RESOLVED:
				break;

			case RET_ENSURELI_FAILED_NOINDEX:
				return RET_ADDCONSTRAINT_FAILED_INFEASIBILITY;

			case RET_ENSURELI_FAILED_CYCLING:
				return RET_ADDCONSTRAINT_FAILED_INFEASIBILITY;

			case RET_ENSURELI_DROPPED:
				return SUCCESSFUL_RETURN;

			default:
				return THROWERROR( RET_ENSURELI_FAILED );
		}
	}

	/* some definitions */
	int_t nV  = getNV( );
	int_t nFR = getNFR( );
	int_t nAC = getNAC( );
	int_t nZ  = getNZ( );

	int_t tcol = sizeT - nAC;


	int_t* FR_idx;
	bounds.getFree( )->getNumberArray( &FR_idx );

	real_t* aFR = new real_t[nFR];
	real_t* wZ = new real_t[nZ];
	for( i=0; i<nZ; ++i )
		wZ[i] = 0.0;


	/* II) ADD NEW ACTIVE CONSTRAINT TO MATRIX T: */
	/* 1) Add row [wZ wY] = aFR'*[Z Y] to the end of T: assign aFR. */
	A->getRow(number, bounds.getFree(), 1.0, aFR);

	/* calculate wZ */
	for( i=0; i<nFR; ++i )
	{
		ii = FR_idx[i];
		for( j=0; j<nZ; ++j )
			wZ[j] += aFR[i] * QQ(ii,j);
	}

	/* 2) Calculate wY and store it directly into T. */
	if ( nAC > 0 )
	{
		for( j=0; j<nAC; ++j )
			TT(nAC,tcol+j) = 0.0;
		for( i=0; i<nFR; ++i )
		{
			ii = FR_idx[i];
			for( j=0; j<nAC; ++j )
				TT(nAC,tcol+j) += aFR[i] * QQ(ii,nZ+j);
		}
	}

	delete[] aFR;


	real_t c, s, nu;

	if ( nZ > 0 )
	{
		/* II) RESTORE TRIANGULAR FORM OF T: */
		/*     Use column-wise Givens rotations to restore reverse triangular form
		*      of T, simultanenous change of Q (i.e. Z) and R. */
		for( j=0; j<nZ-1; ++j )
		{
			computeGivens( wZ[j+1],wZ[j], wZ[j+1],wZ[j],c,s );
			nu = s/(1.0+c);

			for( i=0; i<nFR; ++i )
			{
				ii = FR_idx[i];
				applyGivens( c,s,nu,QQ(ii,1+j),QQ(ii,j), QQ(ii,1+j),QQ(ii,j) );
			}

			if ( ( updateCholesky == BT_TRUE ) &&
				 ( hessianType != HST_ZERO )   && ( hessianType != HST_IDENTITY ) )
			{
				for( i=0; i<=j+1; ++i )
					applyGivens( c,s,nu,RR(i,1+j),RR(i,j), RR(i,1+j),RR(i,j) );
			}
		}

		TT(nAC,tcol-1) = wZ[nZ-1];


		if ( ( updateCholesky == BT_TRUE ) &&
			 ( hessianType != HST_ZERO )   && ( hessianType != HST_IDENTITY ) )
		{
			/* III) RESTORE TRIANGULAR FORM OF R:
			 *      Use row-wise Givens rotations to restore upper triangular form of R. */
			for( i=0; i<nZ-1; ++i )
			{
				computeGivens( RR(i,i),RR(1+i,i), RR(i,i),RR(1+i,i),c,s );
				nu = s/(1.0+c);

				for( j=(1+i); j<(nZ-1); ++j ) /* last column of R is thrown away */
					applyGivens( c,s,nu,RR(i,j),RR(1+i,j), RR(i,j),RR(1+i,j) );
			}
			/* last column of R is thrown away */
			for( i=0; i<nZ; ++i )
				RR(i,nZ-1) = 0.0;
		}
	}

	delete[] wZ;


	/* IV) UPDATE INDICES */
	tabularOutput.idxAddC = number;
	if ( constraints.moveInactiveToActive( number,C_status ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_ADDCONSTRAINT_FAILED );


	return SUCCESSFUL_RETURN;
}



/*
 *	a d d C o n s t r a i n t _ c h e c k L I
 */
returnValue QProblem::addConstraint_checkLI( int_t number )
{
	returnValue returnvalue = RET_LINEARLY_DEPENDENT;

	int_t i, j, ii;
	int_t nV  = getNV( );
	int_t nFR = getNFR( );
	int_t nZ  = getNZ( );
	int_t nC  = getNC( );
	int_t nAC = getNAC();
	int_t nFX = getNFX();
	int_t *FR_idx;

	bounds.getFree( )->getNumberArray( &FR_idx );


	if (options.enableFullLITests)
	{
		/*
		 * expensive LI test. Backsolve with refinement using special right
		 * hand side. This gives an estimate for what should be considered
		 * "zero". We then check linear independence relative to this estimate.
		 */

		int_t *FX_idx, *AC_idx, *IAC_idx;

		real_t *delta_g   = new real_t[nV];
		real_t *delta_xFX = new real_t[nFX];
		real_t *delta_xFR = new real_t[nFR];
		real_t *delta_yAC = new real_t[nAC];
		real_t *delta_yFX = new real_t[nFX];

		bounds.getFixed( )->getNumberArray( &FX_idx );
		constraints.getActive( )->getNumberArray( &AC_idx );
		constraints.getInactive( )->getNumberArray( &IAC_idx );

		int_t dim = (nC>nV)?nC:nV;
		real_t *nul = new real_t[dim];
		for (ii = 0; ii < dim; ++ii)
			nul[ii]=0.0;

		A->getRow (number, 0, 1.0, delta_g);

		// AW: I think original line overwrote correct return value
		// original: returnvalue = determineStepDirection ( delta_g,
		returnValue dsdreturnvalue = determineStepDirection ( delta_g,
											  nul, nul, nul, nul,
											  BT_FALSE, BT_FALSE,
											  delta_xFX, delta_xFR, delta_yAC, delta_yFX);
		if (dsdreturnvalue!=SUCCESSFUL_RETURN)
			returnvalue = dsdreturnvalue;

		delete[] nul;

		/* compute the weight in inf-norm */
		real_t weight = 0.0;
		for (ii = 0; ii < nAC; ++ii)
		{
			real_t a = getAbs (delta_yAC[ii]);
			if (weight < a) weight = a;
		}
		for (ii = 0; ii < nFX; ++ii)
		{
			real_t a = getAbs (delta_yFX[ii]);
			if (weight < a) weight = a;
		}

		/* look at the "zero" in a relative inf-norm */
		real_t zero = 0.0;
		for (ii = 0; ii < nFX; ++ii)
		{
			real_t a = getAbs (delta_xFX[ii]);
			if (zero < a) zero = a;
		}
		for (ii = 0; ii < nFR; ++ii)
		{
			real_t a = getAbs (delta_xFR[ii]);
			if (zero < a) zero = a;
		}

		/* relative test against zero in inf-norm */
		if (zero > options.epsLITests * weight)
			returnvalue = RET_LINEARLY_INDEPENDENT;

		delete[] delta_yFX;
		delete[] delta_yAC;
		delete[] delta_xFR;
		delete[] delta_xFX;
		delete[] delta_g;

	}
	else
	{
		/*
		 * cheap LI test for constraint. Check if constraint <number> is
		 * linearly independent from the the active ones (<=> is element of null
		 * space of Afr).
		 */

		real_t *Arow = new real_t[nFR];
		A->getRow(number, bounds.getFree(), 1.0, Arow);

		real_t sum, l2;

		l2  = 0.0;
		for (i = 0; i < nFR; i++)
			l2  += Arow[i]*Arow[i];

		for( j=0; j<nZ; ++j )
		{
			sum = 0.0;
			for( i=0; i<nFR; ++i )
			{
				ii = FR_idx[i];
				sum += Arow[i] * QQ(ii,j);
			}

			if ( getAbs( sum ) > options.epsLITests*l2 )
			{
				/*fprintf(stdFile, "LI test: |sum| = %9.2e, l2 = %9.2e, var = %d\n", getAbs(sum), l2, jj+1); */
				returnvalue = RET_LINEARLY_INDEPENDENT;
				break;
			}
		}

		delete[] Arow;
	}

	return THROWINFO( returnvalue );
}


/*
 *	a d d C o n s t r a i n t _ e n s u r e L I
 */
returnValue QProblem::addConstraint_ensureLI( int_t number, SubjectToStatus C_status )
{
	int_t i, j, ii, jj;
	int_t nV  = getNV( );
	int_t nFR = getNFR( );
	int_t nFX = getNFX( );
	int_t nAC = getNAC( );
	int_t nZ  = getNZ( );


	/* I) Check if new constraint is linearly independent from the active ones. */
	returnValue returnvalueCheckLI = addConstraint_checkLI( number );

	if ( returnvalueCheckLI == RET_INDEXLIST_CORRUPTED )
		return THROWERROR( RET_ENSURELI_FAILED );

	if ( returnvalueCheckLI == RET_LINEARLY_INDEPENDENT )
		return SUCCESSFUL_RETURN;


 	/* II) NEW CONSTRAINT IS LINEARLY DEPENDENT: */
	/* 1) Determine coefficients of linear combination,
	 *    cf. M.J. Best. Applied Mathematics and Parallel Computing, chapter:
	 *    An Algorithm for the Solution of the Parametric Quadratic Programming
	 *    Problem, pages 57-76. Physica-Verlag, Heidelberg, 1996. */
	int_t* FR_idx;
	bounds.getFree( )->getNumberArray( &FR_idx );

	int_t* FX_idx;
	bounds.getFixed( )->getNumberArray( &FX_idx );

	real_t* xiC = new real_t[nAC];
	real_t* xiC_TMP = new real_t[nAC];
	real_t* xiB = new real_t[nFX];
	real_t* Arow = new real_t[nFR];
	real_t* num = new real_t[nV];

	returnValue returnvalue = SUCCESSFUL_RETURN;

	real_t y_min = options.maxDualJump;
	int_t y_min_number = -1;
	int_t y_min_number_bound = -1;
	BooleanType y_min_isBound = BT_FALSE;

	A->getRow(number, bounds.getFree(), C_status == ST_LOWER ? 1.0 : -1.0, Arow);

	/* 2) Calculate xiC */
	if ( nAC > 0 )
	{
		for( i=0; i<nAC; ++i )
		{
			xiC_TMP[i] = 0.0;
			for( j=0; j<nFR; ++j )
			{
				jj = FR_idx[j];
				xiC_TMP[i] += QQ(jj,nZ+i) * Arow[j];
			}
		}

		if ( backsolveT( xiC_TMP, BT_TRUE, xiC ) != SUCCESSFUL_RETURN )
		{
			returnvalue = RET_ENSURELI_FAILED_TQ;
			goto farewell;
		}
	}

	/* 3) Calculate xiB. */
	int_t* AC_idx;
	constraints.getActive( )->getNumberArray( &AC_idx );

	A->getRow(number, bounds.getFixed(), C_status == ST_LOWER ? 1.0 : -1.0, xiB);
	A->transTimes(constraints.getActive(), bounds.getFixed(), 1, -1.0, xiC, nAC, 1.0, xiB, nFX);

	/* III) DETERMINE CONSTRAINT/BOUND TO BE REMOVED. */

	/* 1) Constraints. */
	for( i=0; i<nAC; ++i )
	{
		ii = AC_idx[i];
		num[i] = y[nV+ii];
	}

	performRatioTest (nAC, AC_idx, &constraints, num, xiC, options.epsNum, options.epsDen, y_min, y_min_number);

	/* 2) Bounds. */
	for( i=0; i<nFX; ++i )
	{
		ii = FX_idx[i];
		num[i] = y[ii];
	}

	performRatioTest (nFX, FX_idx, &bounds, num, xiB, options.epsNum, options.epsDen, y_min, y_min_number_bound);

	if ( y_min_number_bound >= 0 )
	{
		y_min_number = y_min_number_bound;
		y_min_isBound = BT_TRUE;
	}

	#ifndef __SUPPRESSANYOUTPUT__
	/* setup output preferences */
	char messageString[MAX_STRING_LENGTH];
	#endif

	/* IV) REMOVE CONSTRAINT/BOUND FOR RESOLVING LINEAR DEPENDENCE: */
	if ( y_min_number >= 0 )
	{
		/* Update Lagrange multiplier... */
		for( i=0; i<nAC; ++i )
		{
			ii = AC_idx[i];
			y[nV+ii] -= y_min * xiC[i];
		}
		for( i=0; i<nFX; ++i )
		{
			ii = FX_idx[i];
			y[ii] -= y_min * xiB[i];
		}

		/* ... also for newly active constraint... */
		if ( C_status == ST_LOWER )
			y[nV+number] = y_min;
		else
			y[nV+number] = -y_min;

		/* ... and for constraint to be removed. */
		if ( y_min_isBound == BT_TRUE )
		{
			#ifndef __SUPPRESSANYOUTPUT__
			snprintf( messageString,MAX_STRING_LENGTH,"bound no. %d.",(int)y_min_number );
			getGlobalMessageHandler( )->throwInfo( RET_REMOVE_FROM_ACTIVESET,messageString,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
			#endif

			if ( removeBound( y_min_number,BT_TRUE,BT_FALSE,BT_FALSE ) != SUCCESSFUL_RETURN )
			{
				returnvalue = RET_REMOVE_FROM_ACTIVESET_FAILED;
				goto farewell;
			}
			tabularOutput.excRemB = 1;

			y[y_min_number] = 0.0;
		}
		else
		{
			#ifndef __SUPPRESSANYOUTPUT__
			snprintf( messageString,MAX_STRING_LENGTH,"constraint no. %d.",(int)y_min_number );
			getGlobalMessageHandler( )->throwInfo( RET_REMOVE_FROM_ACTIVESET,messageString,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
			#endif

			if ( removeConstraint( y_min_number,BT_TRUE,BT_FALSE,BT_FALSE ) != SUCCESSFUL_RETURN )
			{
				returnvalue = RET_REMOVE_FROM_ACTIVESET_FAILED;
				goto farewell;
			}
			tabularOutput.excRemC = 1;

			y[nV+y_min_number] = 0.0;
		}
	}
	else
	{
		if (options.enableDropInfeasibles == BT_TRUE) {
			/* dropping of infeasible constraints according to drop priorities */
			returnvalue = dropInfeasibles (number, C_status, BT_FALSE, xiB, xiC);
		}
		else
		{
			/* no constraint/bound can be removed => QP is infeasible! */
			returnvalue = RET_ENSURELI_FAILED_NOINDEX;
			setInfeasibilityFlag( returnvalue );
		}
	}

farewell:
	delete[] num;
	delete[] Arow;
	delete[] xiB;
	delete[] xiC_TMP;
	delete[] xiC;

	getGlobalMessageHandler( )->throwInfo( RET_LI_RESOLVED,0,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );

	return ( (returnvalue != SUCCESSFUL_RETURN) && (returnvalue != RET_ENSURELI_FAILED_NOINDEX ) ) ? THROWERROR (returnvalue) : returnvalue;
}



/*
 *	a d d B o u n d
 */
returnValue QProblem::addBound(	int_t number, SubjectToStatus B_status,
								BooleanType updateCholesky,
								BooleanType ensureLI
								)
{
	int_t i, j, ii;

	/* consistency checks */
	if ( bounds.getStatus( number ) != ST_INACTIVE )
		return THROWERROR( RET_BOUND_ALREADY_ACTIVE );

	if ( getNFR( ) == bounds.getNUV( ) )
		return THROWERROR( RET_ALL_BOUNDS_ACTIVE );

	if ( ( getStatus( ) == QPS_NOTINITIALISED )    ||
		 ( getStatus( ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( getStatus( ) == QPS_HOMOTOPYQPSOLVED )  ||
 		 ( getStatus( ) == QPS_SOLVED )            )
	{
		return THROWERROR( RET_UNKNOWN_BUG );
	}


	/* I) ENSURE LINEAR INDEPENDENCE OF THE WORKING SET,
	 *    i.e. remove a constraint or bound if linear dependence occurs. */
	/* check for LI only if Cholesky decomposition shall be updated! */
	if ( ( updateCholesky == BT_TRUE ) && ( ensureLI == BT_TRUE ) )
	{
		returnValue ensureLIreturnvalue = addBound_ensureLI( number,B_status );

		switch ( ensureLIreturnvalue )
		{
			case SUCCESSFUL_RETURN:
				break;

			case RET_LI_RESOLVED:
				break;

			case RET_ENSURELI_FAILED_NOINDEX:
				return RET_ADDBOUND_FAILED_INFEASIBILITY;

			case RET_ENSURELI_FAILED_CYCLING:
				return RET_ADDBOUND_FAILED_INFEASIBILITY;

			case RET_ENSURELI_DROPPED:
				return SUCCESSFUL_RETURN;

			default:
				return THROWERROR( RET_ENSURELI_FAILED );
		}
	}


	/* some definitions */
	int_t nV  = getNV( );
	int_t nFR = getNFR( );
	int_t nAC = getNAC( );
	int_t nZ  = getNZ( );

	int_t tcol = sizeT - nAC;


	/* II) SWAP INDEXLIST OF FREE VARIABLES:
	 *     move the variable to be fixed to the end of the list of free variables. */
	int_t lastfreenumber = bounds.getFree( )->getLastNumber( );
	if ( lastfreenumber != number )
		if ( bounds.swapFree( number,lastfreenumber ) != SUCCESSFUL_RETURN )
			THROWERROR( RET_ADDBOUND_FAILED );


	int_t* FR_idx;
	bounds.getFree( )->getNumberArray( &FR_idx );

	real_t* w = new real_t[nFR];


	/* III) ADD NEW ACTIVE BOUND TO TOP OF MATRIX T: */
	/* 1) add row [wZ wY] = [Z Y](number) at the top of T: assign w */
	for( i=0; i<nFR; ++i )
		w[i] = QQ(FR_idx[nFR-1],i);


	/* 2) Use column-wise Givens rotations to restore reverse triangular form
	 *    of the first row of T, simultanenous change of Q (i.e. Z) and R. */
	real_t c, s, nu;

	for( j=0; j<nZ-1; ++j )
	{
		computeGivens( w[j+1],w[j], w[j+1],w[j],c,s );
		nu = s/(1.0+c);

		for( i=0; i<nFR; ++i )
		{
			ii = FR_idx[i];
			applyGivens( c,s,nu,QQ(ii,1+j),QQ(ii,j), QQ(ii,1+j),QQ(ii,j) );
		}

		if ( ( updateCholesky == BT_TRUE ) &&
			 ( hessianType != HST_ZERO )   && ( hessianType != HST_IDENTITY ) )
		{
			for( i=0; i<=j+1; ++i )
				applyGivens( c,s,nu,RR(i,1+j),RR(i,j), RR(i,1+j),RR(i,j) );
		}
	}


	if ( nAC > 0 )	  /* ( nAC == 0 ) <=> ( nZ == nFR ) <=> Y and T are empty => nothing to do */
	{
		/* store new column a in a temporary vector instead of shifting T one column to the left */
		real_t* tmp = new real_t[nAC];
		for( i=0; i<nAC; ++i )
			tmp[i] = 0.0;

		{
			j = nZ-1;

			computeGivens( w[j+1],w[j], w[j+1],w[j],c,s );
			nu = s/(1.0+c);

			for( i=0; i<nFR; ++i )
			{
				ii = FR_idx[i];
				applyGivens( c,s,nu,QQ(ii,1+j),QQ(ii,j), QQ(ii,1+j),QQ(ii,j) );
			}

			applyGivens( c,s,nu,TT(nAC-1,tcol),tmp[nAC-1], tmp[nAC-1],TT(nAC-1,tcol) );
		}

		for( j=nZ; j<nFR-1; ++j )
		{
			computeGivens( w[j+1],w[j], w[j+1],w[j],c,s );
			nu = s/(1.0+c);

			for( i=0; i<nFR; ++i )
			{
				ii = FR_idx[i];
				applyGivens( c,s,nu,QQ(ii,1+j),QQ(ii,j), QQ(ii,1+j),QQ(ii,j) );
			}

			for( i=(nFR-2-j); i<nAC; ++i )
				applyGivens( c,s,nu,TT(i,1+tcol-nZ+j),tmp[i], tmp[i],TT(i,1+tcol-nZ+j) );
		}

		delete[] tmp;
	}

	delete[] w;


	if ( ( updateCholesky == BT_TRUE ) &&
		 ( hessianType != HST_ZERO )   && ( hessianType != HST_IDENTITY ) )
	{
		/* IV) RESTORE TRIANGULAR FORM OF R:
		 *     use row-wise Givens rotations to restore upper triangular form of R */
		for( i=0; i<nZ-1; ++i )
		{
			computeGivens( RR(i,i),RR(1+i,i), RR(i,i),RR(1+i,i),c,s );
			nu = s/(1.0+c);

			for( j=(1+i); j<nZ-1; ++j ) /* last column of R is thrown away */
				applyGivens( c,s,nu,RR(i,j),RR(1+i,j), RR(i,j),RR(1+i,j) );
		}
		/* last column of R is thrown away */
		for( i=0; i<nZ; ++i )
			RR(i,nZ-1) = 0.0;
	}


	/* V) UPDATE INDICES */
	tabularOutput.idxAddB = number;
	if ( bounds.moveFreeToFixed( number,B_status ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_ADDBOUND_FAILED );

	return SUCCESSFUL_RETURN;
}


/*
 *	a d d B o u n d _ c h e c k L I
 */
returnValue QProblem::addBound_checkLI( int_t number )
{
	int_t i, ii;
	int_t nV  = getNV( );  /* for QQ() macro */
	int_t nFR = getNFR( );
	int_t nAC = getNAC();
	int_t nFX = getNFX();
	int_t nC  = getNC( );
	returnValue returnvalue = RET_LINEARLY_DEPENDENT;

	if (options.enableFullLITests)
	{
		/*
		 * expensive LI test. Backsolve with refinement using special right
		 * hand side. This gives an estimate for what should be considered
		 * "zero". We then check linear independence relative to this estimate.
		 */

		/*
		 * expensive LI test. Backsolve with refinement using special right
		 * hand side. This gives an estimate for what should be considered
		 * "zero". We then check linear independence relative to this estimate.
		 */

		real_t *delta_g   = new real_t[nV];
		real_t *delta_xFX = new real_t[nFX];
		real_t *delta_xFR = new real_t[nFR];
		real_t *delta_yAC = new real_t[nAC];
		real_t *delta_yFX = new real_t[nFX];

		for (ii = 0; ii < nV; ++ii)
			delta_g[ii] = 0.0;
		delta_g[number] = 1.0;	/* sign doesn't matter here */

		int_t dim = (nC>nV)?nC:nV;
		real_t *nul = new real_t[dim];
		for (ii = 0; ii < dim; ++ii)
			nul[ii]=0.0;

		returnValue dsdReturnValue = determineStepDirection (
				delta_g, nul, nul, nul, nul, BT_FALSE, BT_FALSE,
				delta_xFX, delta_xFR, delta_yAC, delta_yFX);
		if (dsdReturnValue != SUCCESSFUL_RETURN)
			returnvalue = dsdReturnValue;

		/* compute the weight in inf-norm */
		real_t weight = 0.0;
		for (ii = 0; ii < nAC; ++ii)
		{
			real_t a = getAbs (delta_yAC[ii]);
			if (weight < a) weight = a;
		}
		for (ii = 0; ii < nFX; ++ii)
		{
			real_t a = getAbs (delta_yFX[ii]);
			if (weight < a) weight = a;
		}

		/* look at the "zero" in a relative inf-norm */
		real_t zero = 0.0;
		for (ii = 0; ii < nFX; ++ii)
		{
			real_t a = getAbs (delta_xFX[ii]);
			if (zero < a) zero = a;
		}
		for (ii = 0; ii < nFR; ++ii)
		{
			real_t a = getAbs (delta_xFR[ii]);
			if (zero < a) zero = a;
		}

		/* relative test against zero in inf-norm */
		if (zero > options.epsLITests * weight)
			returnvalue = RET_LINEARLY_INDEPENDENT;

		delete[] nul;
		delete[] delta_yFX;
		delete[] delta_yAC;
		delete[] delta_xFR;
		delete[] delta_xFX;
		delete[] delta_g;

	}
	else
	{
		/*
		 * cheap LI test for simple bound. Check if constraint <number> is
		 * linearly independent from the the active ones (<=> is element of null
		 * space of Afr).
		 */

		/* some definitions */
		int_t nZ  = getNZ( );

		for( i=0; i<nZ; ++i )
			if ( getAbs( QQ(number,i) ) > options.epsLITests )
			{
				returnvalue = RET_LINEARLY_INDEPENDENT;
				break;
			}
	}

	return THROWINFO( returnvalue );
}


/*
 *	a d d B o u n d _ e n s u r e L I
 */
returnValue QProblem::addBound_ensureLI( int_t number, SubjectToStatus B_status )
{
	int_t i, ii;
	int_t nV  = getNV( );
	int_t nFX = getNFX( );
	int_t nAC = getNAC( );
	int_t nZ  = getNZ( );


	/* I) Check if new constraint is linearly independent from the active ones. */
	returnValue returnvalueCheckLI = addBound_checkLI( number );

	if ( returnvalueCheckLI == RET_INDEXLIST_CORRUPTED )
		return THROWERROR( RET_ENSURELI_FAILED );

	if ( returnvalueCheckLI == RET_LINEARLY_INDEPENDENT )
		return SUCCESSFUL_RETURN;


 	/* II) NEW BOUND IS LINEARLY DEPENDENT: */
	/* 1) Determine coefficients of linear combination,
	 *    cf. M.J. Best. Applied Mathematics and Parallel Computing, chapter:
	 *    An Algorithm for the Solution of the Parametric Quadratic Programming
	 *    Problem, pages 57-76. Physica-Verlag, Heidelberg, 1996. */
	int_t* FR_idx;
	bounds.getFree( )->getNumberArray( &FR_idx );

	int_t* FX_idx;
	bounds.getFixed( )->getNumberArray( &FX_idx );

	int_t* AC_idx;
	constraints.getActive( )->getNumberArray( &AC_idx );

	real_t* xiC = new real_t[nAC];
	real_t* xiC_TMP = new real_t[nAC];
	real_t* xiB = new real_t[nFX];
	real_t* num = new real_t[nV];

	real_t y_min = options.maxDualJump;
	int_t y_min_number = -1;
	int_t y_min_number_bound = -1;
	BooleanType y_min_isBound = BT_FALSE;

	returnValue returnvalue = SUCCESSFUL_RETURN;


	/* 2) Calculate xiC. */
	if ( nAC > 0 )
	{
		if ( B_status == ST_LOWER )
		{
			for( i=0; i<nAC; ++i )
				xiC_TMP[i] = QQ(number,nZ+i);
		}
		else
		{
			for( i=0; i<nAC; ++i )
				xiC_TMP[i] = -QQ(number,nZ+i);
		}

		if ( backsolveT( xiC_TMP, BT_TRUE, xiC ) != SUCCESSFUL_RETURN )
		{
			returnvalue = RET_ENSURELI_FAILED_TQ;
			goto farewell;
		}
	}

	/* 3) Calculate xiB. */
	A->transTimes(constraints.getActive(), bounds.getFixed(), 1, -1.0, xiC, nAC, 0.0, xiB, nFX);


	/* III) DETERMINE CONSTRAINT/BOUND TO BE REMOVED. */

	/* 1) Constraints. */
	for( i=0; i<nAC; ++i )
	{
		ii = AC_idx[i];
		num[i] = y[nV+ii];
	}

	performRatioTest( nAC,AC_idx,&constraints, num,xiC, options.epsNum,options.epsDen, y_min,y_min_number );

	/* 2) Bounds. */
	for( i=0; i<nFX; ++i )
	{
		ii = FX_idx[i];
		num[i] = y[ii];
	}

	performRatioTest( nFX,FX_idx,&bounds, num,xiB, options.epsNum,options.epsDen, y_min,y_min_number_bound );

	if ( y_min_number_bound >= 0 )
	{
		y_min_number = y_min_number_bound;
		y_min_isBound = BT_TRUE;
	}

	/* IV) REMOVE CONSTRAINT/BOUND FOR RESOLVING LINEAR DEPENDENCE: */
	#ifndef __SUPPRESSANYOUTPUT__
	char messageString[MAX_STRING_LENGTH];
	#endif

	if ( y_min_number >= 0 )
	{
		/* Update Lagrange multiplier... */
		for( i=0; i<nAC; ++i )
		{
			ii = AC_idx[i];
			y[nV+ii] -= y_min * xiC[i];
		}
		for( i=0; i<nFX; ++i )
		{
			ii = FX_idx[i];
			y[ii] -= y_min * xiB[i];
		}

		/* ... also for newly active bound ... */
		if ( B_status == ST_LOWER )
			y[number] = y_min;
		else
			y[number] = -y_min;

		/* ... and for bound to be removed. */
		if ( y_min_isBound == BT_TRUE )
		{
			#ifndef __SUPPRESSANYOUTPUT__
			snprintf( messageString,MAX_STRING_LENGTH,"bound no. %d.",(int)y_min_number );
			getGlobalMessageHandler( )->throwInfo( RET_REMOVE_FROM_ACTIVESET,messageString,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
			#endif

			if ( removeBound( y_min_number,BT_TRUE,BT_FALSE,BT_FALSE ) != SUCCESSFUL_RETURN )
			{
				returnvalue = RET_REMOVE_FROM_ACTIVESET_FAILED;
				goto farewell;
			}
			tabularOutput.excRemB = 1;

			y[y_min_number] = 0.0;
		}
		else
		{
			#ifndef __SUPPRESSANYOUTPUT__
			snprintf( messageString,MAX_STRING_LENGTH,"constraint no. %d.",(int)y_min_number );
			getGlobalMessageHandler( )->throwInfo( RET_REMOVE_FROM_ACTIVESET,messageString,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
			#endif

			if ( removeConstraint( y_min_number,BT_TRUE,BT_FALSE,BT_FALSE ) != SUCCESSFUL_RETURN )
			{
				returnvalue = RET_REMOVE_FROM_ACTIVESET_FAILED;
				goto farewell;
			}
			tabularOutput.excRemC = 1;

			y[nV+y_min_number] = 0.0;
		}
	}
	else
	{
		if (options.enableDropInfeasibles == BT_TRUE) {
			/* dropping of infeasible constraints according to drop priorities */
			returnvalue = dropInfeasibles (number, B_status, BT_TRUE, xiB, xiC);
		}
		else
		{
			/* no constraint/bound can be removed => QP is infeasible! */
			returnvalue = RET_ENSURELI_FAILED_NOINDEX;
			setInfeasibilityFlag( returnvalue );
		}
	}

farewell:
	delete[] num;
	delete[] xiB;
	delete[] xiC_TMP;
	delete[] xiC;

	getGlobalMessageHandler( )->throwInfo( RET_LI_RESOLVED,0,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );

	return ( (returnvalue != SUCCESSFUL_RETURN) && (returnvalue != RET_ENSURELI_FAILED_NOINDEX ) ) ? THROWERROR (returnvalue) : returnvalue;
}



/*
 *	r e m o v e C o n s t r a i n t
 */
returnValue QProblem::removeConstraint(	int_t number,
										BooleanType updateCholesky,
										BooleanType allowFlipping,
										BooleanType ensureNZC
										)
{
	int_t i, j, ii, jj;
	returnValue returnvalue = SUCCESSFUL_RETURN;
	BooleanType hasFlipped = BT_FALSE;

	/* consistency check */
	if ( ( getStatus( ) == QPS_NOTINITIALISED )    ||
		 ( getStatus( ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( getStatus( ) == QPS_HOMOTOPYQPSOLVED )  ||
 		 ( getStatus( ) == QPS_SOLVED )            )
	{
		return THROWERROR( RET_UNKNOWN_BUG );
	}

	/* some definitions */
	int_t nV  = getNV( );
	int_t nFR = getNFR( );
	int_t nAC = getNAC( );
	int_t nZ  = getNZ( );

	int_t tcol = sizeT - nAC;
	int_t number_idx = constraints.getActive( )->getIndex( number );

	int_t addIdx;
	BooleanType addBoundNotConstraint;
	SubjectToStatus addStatus;
	BooleanType exchangeHappened = BT_FALSE;


	/* consistency checks */
	if ( constraints.getStatus( number ) == ST_INACTIVE )
		return THROWERROR( RET_CONSTRAINT_NOT_ACTIVE );

	if ( ( number_idx < 0 ) || ( number_idx >= nAC ) )
		return THROWERROR( RET_CONSTRAINT_NOT_ACTIVE );


	int_t* FR_idx;
	bounds.getFree( )->getNumberArray( &FR_idx );

	/* N) PERFORM ZERO CURVATURE TEST. */
	if (ensureNZC == BT_TRUE)
	{
		returnvalue = ensureNonzeroCurvature(BT_FALSE, number, exchangeHappened, addBoundNotConstraint, addIdx, addStatus);

		if (returnvalue != SUCCESSFUL_RETURN)
			return returnvalue;
	}

	/* save index sets and decompositions for flipping bounds strategy */
	if ( ( exchangeHappened == BT_FALSE ) && ( options.enableFlippingBounds == BT_TRUE ) && ( allowFlipping == BT_TRUE ) )
		flipper.set( &bounds,R,&constraints,Q,T );

	/* I) REMOVE <number>th ROW FROM T,
	 *    i.e. shift rows number+1 through nAC  upwards (instead of the actual
	 *    constraint number its corresponding index within matrix A is used). */
	if ( number_idx < nAC-1 )
	{
		for( i=(number_idx+1); i<nAC; ++i )
			for( j=(nAC-i-1); j<nAC; ++j )
				TT(i-1,tcol+j) = TT(i,tcol+j);
		/* gimmick: write zeros into the last row of T */
		for( j=0; j<nAC; ++j )
			TT(nAC-1,tcol+j) = 0.0;


		/* II) RESTORE TRIANGULAR FORM OF T,
		 *     use column-wise Givens rotations to restore reverse triangular form
		 *     of T simultanenous change of Q (i.e. Y). */
		real_t c, s, nu;

		for( j=(nAC-2-number_idx); j>=0; --j )
		{
			computeGivens( TT(nAC-2-j,tcol+1+j),TT(nAC-2-j,tcol+j), TT(nAC-2-j,tcol+1+j),TT(nAC-2-j,tcol+j),c,s );
			nu = s/(1.0+c);

			for( i=(nAC-j-1); i<(nAC-1); ++i )
				applyGivens( c,s,nu,TT(i,tcol+1+j),TT(i,tcol+j), TT(i,tcol+1+j),TT(i,tcol+j) );

			for( i=0; i<nFR; ++i )
			{
				ii = FR_idx[i];
				applyGivens( c,s,nu,QQ(ii,nZ+1+j),QQ(ii,nZ+j), QQ(ii,nZ+1+j),QQ(ii,nZ+j) );
			}
		}
	}
	else
	{
		/* gimmick: write zeros into the last row of T */
		for( j=0; j<nAC; ++j )
			TT(nAC-1,tcol+j) = 0.0;
	}


	if ( ( updateCholesky == BT_TRUE ) &&
		 ( hessianType != HST_ZERO )   && ( hessianType != HST_IDENTITY ) )
	{
		/* III) UPDATE CHOLESKY DECOMPOSITION,
		 *      calculate new additional column (i.e. [r sqrt(rho2)]')
		 *      of the Cholesky factor R. */
		real_t* Hz = new real_t[nFR];
		real_t* z = new real_t[nFR];
		real_t rho2 = 0.0;

		/* 1) Calculate Hz = H*z, where z is the new rightmost column of Z
		 *    (i.e. the old leftmost column of Y).  */
		for( j=0; j<nFR; ++j )
			z[j] = QQ(FR_idx[j],nZ);
		H->times(bounds.getFree(), bounds.getFree(), 1, 1.0, z, nFR, 0.0, Hz, nFR);
		delete[] z;

		if ( nZ > 0 )
		{
			real_t* ZHz = new real_t[nZ];
			for ( i=0; i<nZ; ++i )
				ZHz[i] = 0.0;
			real_t* r = new real_t[nZ];

			/* 2) Calculate ZHz = Z'*Hz (old Z). */
			for( j=0; j<nFR; ++j )
			{
				jj = FR_idx[j];

				for( i=0; i<nZ; ++i )
					ZHz[i] += QQ(jj,i) * Hz[j];
			}

			/* 3) Calculate r = R^-T * ZHz. */
			if ( backsolveR( ZHz,BT_TRUE,r ) != SUCCESSFUL_RETURN )
			{
				delete[] Hz; delete[] r; delete[] ZHz;
				return THROWERROR( RET_REMOVECONSTRAINT_FAILED );
			}

			/* 4) Calculate rho2 = rho^2 = z'*Hz - r'*r
			 *    and store r into R. */
			for( i=0; i<nZ; ++i )
			{
				rho2 -= r[i]*r[i];
				RR(i,nZ) = r[i];
			}

			delete[] r; delete[] ZHz;
		}

		/* 5) Store rho into R. */
		for( j=0; j<nFR; ++j )
			rho2 += QQ(FR_idx[j],nZ) * Hz[j];

		delete[] Hz;

		if ( ( options.enableFlippingBounds == BT_TRUE ) && ( allowFlipping == BT_TRUE ) && ( exchangeHappened == BT_FALSE ) )
		{
			if ( rho2 > options.epsFlipping )
				RR(nZ,nZ) = getSqrt( rho2 );
			else
			{
				hessianType = HST_SEMIDEF;

				flipper.get( &bounds,R,&constraints,Q,T );
				constraints.flipFixed(number);
				tabularOutput.idxAddC = number;
				tabularOutput.excAddC = 2;

				switch (constraints.getStatus(number))
				{
					case ST_LOWER:
						lbA[number] = ubA[number]; Ax_l[number] = -Ax_u[number]; break;
					case ST_UPPER:
						ubA[number] = lbA[number]; Ax_u[number] = -Ax_l[number]; break;
					default:
						return THROWERROR( RET_MOVING_BOUND_FAILED );
				}

				hasFlipped = BT_TRUE;
			}
		}
		else if ( exchangeHappened == BT_FALSE )
		{
			if ( rho2 > ZERO )
				RR(nZ,nZ) = getSqrt( rho2 );
			else
			{
				if ( allowFlipping == BT_FALSE )
				{
					RR(nZ,nZ) = 100.0*EPS;
				}
				else
				{
					hessianType = HST_SEMIDEF;
					return THROWERROR( RET_HESSIAN_NOT_SPD );
				}
			}
		}
	}


	/* IV) UPDATE INDICES */
	tabularOutput.idxRemC = number;
	if ( hasFlipped == BT_FALSE )
	{
		if ( constraints.moveActiveToInactive( number ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_REMOVECONSTRAINT_FAILED );
	}

	if (exchangeHappened == BT_TRUE)
	{
		/* add bound or constraint */

		/* hessianType = HST_SEMIDEF; */
		RR(nZ,nZ) = 0.0;

		if ( addBoundNotConstraint )
		{
			addBound(addIdx, addStatus, BT_TRUE, BT_FALSE);
			tabularOutput.excAddB = 1;
		}
		else
		{
			addConstraint(addIdx, addStatus, BT_TRUE, BT_FALSE);
			tabularOutput.excAddC = 1;
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	r e m o v e B o u n d
 */
returnValue QProblem::removeBound(	int_t number,
									BooleanType updateCholesky,
									BooleanType allowFlipping,
									BooleanType ensureNZC
									)
{
	int_t i, j, ii, jj;
	returnValue returnvalue = SUCCESSFUL_RETURN;
	int_t addIdx;
	BooleanType addBoundNotConstraint;
	SubjectToStatus addStatus;
	BooleanType exchangeHappened = BT_FALSE;


	/* consistency checks */
	if ( bounds.getStatus( number ) == ST_INACTIVE )
		return THROWERROR( RET_BOUND_NOT_ACTIVE );

	if ( ( getStatus( ) == QPS_NOTINITIALISED )    ||
		 ( getStatus( ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( getStatus( ) == QPS_HOMOTOPYQPSOLVED )  ||
 		 ( getStatus( ) == QPS_SOLVED )            )
	{
		return THROWERROR( RET_UNKNOWN_BUG );
	}

	/* some definitions */
	int_t nV  = getNV( );
	int_t nFR = getNFR( );
	int_t nAC = getNAC( );
	int_t nZ  = getNZ( );

	int_t tcol = sizeT - nAC;

	/* 0) PERFORM ZERO CURVATURE TEST. */
	if (ensureNZC == BT_TRUE)
	{
		returnvalue = ensureNonzeroCurvature(BT_TRUE, number, exchangeHappened, addBoundNotConstraint, addIdx, addStatus);

		if (returnvalue != SUCCESSFUL_RETURN)
			return returnvalue;
	}

	/* save index sets and decompositions for flipping bounds strategy */
	if ( ( options.enableFlippingBounds == BT_TRUE ) && ( allowFlipping == BT_TRUE ) && ( exchangeHappened == BT_FALSE ) )
		flipper.set( &bounds,R,&constraints,Q,T );

	/* I) UPDATE INDICES */
	tabularOutput.idxRemB = number;
	if ( bounds.moveFixedToFree( number ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_REMOVEBOUND_FAILED );

	int_t* FR_idx;
	bounds.getFree( )->getNumberArray( &FR_idx );

	/* I) APPEND <nFR+1>th UNITY VECTOR TO Q. */
	int_t nnFRp1 = FR_idx[nFR];
	for( i=0; i<nFR; ++i )
	{
		ii = FR_idx[i];
		QQ(ii,nFR) = 0.0;
		QQ(nnFRp1,i) = 0.0;
	}
	QQ(nnFRp1,nFR) = 1.0;

	if ( nAC > 0 )
	{
		/* store new column a in a temporary vector instead of shifting T one column to the left and appending a */
		int_t* AC_idx;
		constraints.getActive( )->getNumberArray( &AC_idx );

		real_t* tmp = new real_t[nAC];
		A->getCol(number, constraints.getActive(), 1.0, tmp);


		/* II) RESTORE TRIANGULAR FORM OF T,
		 *     use column-wise Givens rotations to restore reverse triangular form
		 *     of T = [T A(:,number)], simultanenous change of Q (i.e. Y and Z). */
		real_t c, s, nu;

		for( j=(nAC-1); j>=0; --j )
		{
			computeGivens( tmp[nAC-1-j],TT(nAC-1-j,tcol+j),TT(nAC-1-j,tcol+j),tmp[nAC-1-j],c,s );
			nu = s/(1.0+c);

			for( i=(nAC-j); i<nAC; ++i )
				applyGivens( c,s,nu,tmp[i],TT(i,tcol+j),TT(i,tcol+j),tmp[i] );

			for( i=0; i<=nFR; ++i )
			{
				ii = FR_idx[i];
				/* nZ+1+nAC = nFR+1  /  nZ+(1) = nZ+1 */
				applyGivens( c,s,nu,QQ(ii,nZ+1+j),QQ(ii,nZ+j),QQ(ii,nZ+1+j),QQ(ii,nZ+j) );
			}
		}

		delete[] tmp;
	}


	if ( ( updateCholesky == BT_TRUE ) &&
	     ( hessianType != HST_ZERO ) && ( hessianType != HST_IDENTITY ) )
	{
		/* III) UPDATE CHOLESKY DECOMPOSITION,
		 *      calculate new additional column (i.e. [r sqrt(rho2)]')
		 *      of the Cholesky factor R: */
		real_t z2 = QQ(nnFRp1,nZ);
		real_t rho2 = H->diag(nnFRp1) * z2*z2;

		if ( nFR > 0 )
		{
			/* Attention: Index list of free variables has already grown by one! */
			real_t* Hz = new real_t[nFR+1];
			real_t* z = new real_t[nFR+1];
			/* 1) Calculate R'*r = Zfr'*Hfr*z1 + z2*Zfr'*h1 =: Zfr'*Hz + z2*Zfr'*h1 =: rhs and
			 *    rho2 = z1'*Hfr*z1 + 2*z2*h1'*z1 + h2*z2^2 - r'*r =: z1'*Hz + 2*z2*h1'*z1 + h2*z2^2 - r'r */
			for( j=0; j<nFR; ++j )
				z[j] = QQ(FR_idx[j],nZ);
			z[nFR] = 0.0;

			H->times(bounds.getFree(), bounds.getFree(), 1, 1.0, z, nFR+1, 0.0, Hz, nFR+1);
			H->getCol(nnFRp1, bounds.getFree(), 1.0, z);

			if ( nZ > 0 )
			{
				real_t* r = new real_t[nZ];
				real_t* rhs = new real_t[nZ];
				for( i=0; i<nZ; ++i )
					rhs[i] = 0.0;

				/* 2) Calculate rhs. */
				for( j=0; j<nFR; ++j )
				{
					jj = FR_idx[j];
					for( i=0; i<nZ; ++i )
										/* Zfr' * ( Hz + z2*h1 ) */
						rhs[i] += QQ(jj,i) * ( Hz[j] + z2 * z[j] );
				}

				/* 3) Calculate r = R^-T * rhs. */
				if ( backsolveR( rhs,BT_TRUE,BT_TRUE,r ) != SUCCESSFUL_RETURN )
				{
					delete[] z;
					delete[] Hz; delete[] r; delete[] rhs;
					return THROWERROR( RET_REMOVEBOUND_FAILED );
				}


				/* 4) Calculate rho2 = rho^2 = z'*Hz - r'*r
				 *    and store r into R. */
				for( i=0; i<nZ; ++i )
				{
					rho2 -= r[i]*r[i];
					RR(i,nZ) = r[i];
				}

				delete[] rhs; delete[] r;
			}

			for( j=0; j<nFR; ++j )
			{
				jj = FR_idx[j];
							/* z1' * ( Hz + 2*z2*h1 ) */
				rho2 += QQ(jj,nZ) * ( Hz[j] + 2.0*z2*z[j] );
			}

			delete[] z;
			delete[] Hz;
		}

		/* 5) Store rho into R. */
		if ( ( options.enableFlippingBounds == BT_TRUE ) && ( allowFlipping == BT_TRUE ) && ( exchangeHappened == BT_FALSE ) )
		{
			if ( rho2 > options.epsFlipping )
				RR(nZ,nZ) = getSqrt( rho2 );
			else
			{
				if ( hessianType != HST_ZERO )
					hessianType = HST_SEMIDEF;

				flipper.get( &bounds,R,&constraints,Q,T );
				bounds.flipFixed(number);
				tabularOutput.idxAddB = number;
				tabularOutput.excAddB = 2;

				switch (bounds.getStatus(number))
				{
					case ST_LOWER:
						lb[number] = ub[number];
						break;
					case ST_UPPER:
						ub[number] = lb[number];
						break;
					default: return THROWERROR( RET_MOVING_BOUND_FAILED );
				}

			}
		}
		else if ( exchangeHappened == BT_FALSE )
		{
			if ( rho2 > ZERO )
				RR(nZ,nZ) = getSqrt( rho2 );
			else
			{
				if ( allowFlipping == BT_FALSE )
					RR(nZ,nZ) = 100.0*EPS;
				else
				{
					hessianType = HST_SEMIDEF;
					return THROWERROR( RET_HESSIAN_NOT_SPD );
				}
			}
		}
		else
		{
			/* add bound or constraint */

			/* hessianType = HST_SEMIDEF; */
			RR(nZ,nZ) = 0.0;

			if ( addBoundNotConstraint )
			{
				addBound(addIdx, addStatus, BT_TRUE, BT_FALSE);
				tabularOutput.excAddB = 1;
			}
			else
			{
				addConstraint(addIdx, addStatus, BT_TRUE, BT_FALSE);
				tabularOutput.excAddC = 1;
			}
		}
	}

	return SUCCESSFUL_RETURN;
}



returnValue QProblem::performPlainRatioTest(	int_t nIdx,
												const int_t* const idxList,
												const real_t* const num,
												const real_t* const den,
												real_t epsNum,
												real_t epsDen,
												real_t& t,
												int_t& BC_idx
												) const
{
	int_t i;
	for (i = 0; i < nIdx; i++)
		if ( (num[i] > epsNum) && (den[i] > epsDen) && (t * den[i] > num[i]) )
		{
			t = num[i] / den[i];
			BC_idx = idxList[i];
		}

	return SUCCESSFUL_RETURN;
}


returnValue QProblem::ensureNonzeroCurvature(	BooleanType removeBoundNotConstraint,
												int_t remIdx,
												BooleanType &exchangeHappened,
												BooleanType &addBoundNotConstraint,
												int_t &addIdx,
												SubjectToStatus &addStatus
												)
{
	int_t i, ii;
	int_t addLBndIdx = -1, addLCnstrIdx = -1, addUBndIdx = -1, addUCnstrIdx = -1; /* exchange indices */
	int_t *FX_idx, *AC_idx, *IAC_idx;
	returnValue returnvalue = SUCCESSFUL_RETURN;

	int_t nV  = getNV( );
	int_t nFR = getNFR( );
	int_t nAC = getNAC( );
	int_t nC  = getNC( );
	int_t nFX = getNFX( );
	int_t nIAC = getNIAC( );

	int_t* FR_idx;
	bounds.getFree( )->getNumberArray( &FR_idx );

// 	real_t *delta_g   = new real_t[nV];
	real_t *delta_xFX = new real_t[nFX];
	real_t *delta_xFR = new real_t[nFR];
	real_t *delta_yAC = new real_t[nAC];
	real_t *delta_yFX = new real_t[nFX];

	bounds.getFixed( )->getNumberArray( &FX_idx );
	constraints.getActive( )->getNumberArray( &AC_idx );
	constraints.getInactive( )->getNumberArray( &IAC_idx );

	addBoundNotConstraint = BT_TRUE;
	addStatus = ST_INACTIVE;
	exchangeHappened = BT_FALSE;

	if (removeBoundNotConstraint)
	{
		int_t dim = nV < nC ? nC : nV;
		real_t *nul = new real_t[dim];
		real_t *ek = new real_t[nV]; /* minus e_k (bound k is removed) */
		for (ii = 0; ii < dim; ++ii)
			nul[ii]=0.0;
		for (ii = 0; ii < nV; ++ii)
			ek[ii]=0.0;
		ek[remIdx] = bounds.getStatus(remIdx) == ST_LOWER ? 1.0 : -1.0;

		returnvalue = determineStepDirection (nul, nul, nul, ek, ek,
											  BT_FALSE, BT_FALSE,
											  delta_xFX, delta_xFR, delta_yAC, delta_yFX);
		delete[] ek;
		delete[] nul;
	}
	else
	{
		real_t *nul = new real_t[nV];
		real_t *ek = new real_t[nC]; /* minus e_k (constraint k is removed) */
		for (ii = 0; ii < nV; ++ii)
			nul[ii]=0.0;
		for (ii = 0; ii < nC; ++ii)
			ek[ii]=0.0;
		ek[remIdx] = constraints.getStatus(remIdx) == ST_LOWER ? 1.0 : -1.0;

		returnvalue = determineStepDirection (nul,
											  ek, ek, nul, nul,
											  BT_FALSE, BT_TRUE,
											  delta_xFX, delta_xFR, delta_yAC, delta_yFX);
		delete[] ek;
		delete[] nul;
	}

	/* compute the weight in inf-norm */
	real_t normXi = 0.0;
	for (ii = 0; ii < nAC; ++ii)
	{
		real_t a = getAbs (delta_yAC[ii]);
		if (normXi < a) normXi = a;
	}
	for (ii = 0; ii < nFX; ++ii)
	{
		real_t a = getAbs (delta_yFX[ii]);
		if (normXi < a) normXi = a;
	}

	/* look at the "zero" in a relative inf-norm */
	real_t normS = 0.0;
	for (ii = 0; ii < nFX; ++ii)
	{
		real_t a = getAbs (delta_xFX[ii]);
		if (normS < a) normS = a;
	}
	for (ii = 0; ii < nFR; ++ii)
	{
		real_t a = getAbs (delta_xFR[ii]);
		if (normS < a) normS = a;
	}

	/* relative test against zero in inf-norm */
	if (normXi < options.epsNZCTests * normS)
	{
		/* determine jump in x via ratio tests */
		real_t sigmaLBnd, sigmaLCnstr, sigmaUBnd, sigmaUCnstr, sigma;

		/* bounds */

		/* compress x-u */
		real_t *x_W = new real_t[getMax(1,nFR)];
		for (i = 0; i < nFR; i++)
		{
			ii = FR_idx[i];
			x_W[i] = ub[ii] - x[ii];
		}
		/* performRatioTest( nFR,FR_idx,&bounds, x_W,delta_xFR, options.epsNum,options.epsDen, sigmaUBnd,addUBndIdx ); */
		sigmaUBnd = options.maxPrimalJump;
		addUBndIdx = -1;
		performPlainRatioTest(nFR, FR_idx, x_W, delta_xFR, options.epsNum, options.epsDen, sigmaUBnd, addUBndIdx);
		if (removeBoundNotConstraint == BT_TRUE && bounds.getStatus(remIdx) == ST_LOWER)
		{
			/* also consider bound which is to be removed */
			real_t one = 1.0;
			x_W[0] = ub[remIdx] - x[remIdx];
			performPlainRatioTest(1, &remIdx, x_W, &one, options.epsNum, options.epsDen, sigmaUBnd, addUBndIdx);
		}

		/* compress x-l */
		for (i = 0; i < nFR; i++)
		{
			ii = FR_idx[i];
			x_W[i] = x[ii] - lb[ii];
		}
		for (i = 0; i < nFR; i++)
			delta_xFR[i] = -delta_xFR[i];
		/* performRatioTest( nFR,FR_idx,&bounds, x_W,delta_xFR, options.epsNum,options.epsDen, sigmaLBnd,addLBndIdx ); */
		sigmaLBnd = options.maxPrimalJump;
		addLBndIdx = -1;
		performPlainRatioTest(nFR, FR_idx, x_W, delta_xFR, options.epsNum, options.epsDen, sigmaLBnd, addLBndIdx);
		if (removeBoundNotConstraint == BT_TRUE && bounds.getStatus(remIdx) == ST_UPPER)
		{
			/* also consider bound which is to be removed */
			real_t one = 1.0;
			x_W[0] = x[remIdx] - lb[remIdx];
			performPlainRatioTest(1, &remIdx, x_W, &one, options.epsNum, options.epsDen, sigmaLBnd, addLBndIdx);
		}
		for (i = 0; i < nFR; i++)
			delta_xFR[i] = -delta_xFR[i];

		delete[] x_W;

		/* constraints */

		/* compute As (compressed to inactive constraints) */
		real_t *As = new real_t[nIAC];
		A->times(constraints.getInactive(), bounds.getFixed(), 1, 1.0, delta_xFX, nFX, 0.0, As, nIAC);
		A->times(constraints.getInactive(), bounds.getFree(), 1, 1.0, delta_xFR, nFR, 1.0, As, nIAC);

		/* compress Ax_u */
		real_t *Ax_W = new real_t[nIAC];
		for (i = 0; i < nIAC; i++)
		{
			ii = IAC_idx[i];
			Ax_W[i] = Ax_u[ii];
		}
		/* performRatioTest( nIAC,IAC_idx,&constraints, Ax_W,As, options.epsNum,options.epsDen, sigmaUCnstr,addUCnstrIdx ); */
		sigmaUCnstr = options.maxPrimalJump;
		addUCnstrIdx = -1;
		performPlainRatioTest(nIAC, IAC_idx, Ax_W, As, options.epsNum, options.epsDen, sigmaUCnstr, addUCnstrIdx);
		if (removeBoundNotConstraint == BT_FALSE && constraints.getStatus(remIdx) == ST_LOWER)
		{
			/* also consider constraint which is to be removed */
			real_t one = 1.0;
			performPlainRatioTest(1, &remIdx, &Ax_u[remIdx], &one, options.epsNum, options.epsDen, sigmaUCnstr, addUCnstrIdx);
		}

		/* compress Ax_l */
		for (i = 0; i < nIAC; i++)
		{
			ii = IAC_idx[i];
			Ax_W[i] = Ax_l[ii];
		}
		for (i = 0; i < nIAC; i++)
			As[i] = -As[i];
		/* performRatioTest( nIAC,IAC_idx,&constraints, Ax_W,As, options.epsNum,options.epsDen, sigmaLCnstr,addLCnstrIdx ); */
		sigmaLCnstr = options.maxPrimalJump;
		addLCnstrIdx = -1;
		performPlainRatioTest(nIAC, IAC_idx, Ax_W, As, options.epsNum, options.epsDen, sigmaLCnstr, addLCnstrIdx);
		if (removeBoundNotConstraint == BT_FALSE && constraints.getStatus(remIdx) == ST_UPPER)
		{
			/* also consider constraint which is to be removed */
			real_t one = 1.0;
			performPlainRatioTest(1, &remIdx, &Ax_l[remIdx], &one, options.epsNum, options.epsDen, sigmaLCnstr, addLCnstrIdx);
		}

		/* perform primal jump */
		sigma = options.maxPrimalJump;
		if (sigmaUCnstr < sigma) { sigma = sigmaUCnstr; addStatus = ST_UPPER; addBoundNotConstraint = BT_FALSE; addIdx = addUCnstrIdx; }
		if (sigmaLCnstr < sigma) { sigma = sigmaLCnstr; addStatus = ST_LOWER; addBoundNotConstraint = BT_FALSE; addIdx = addLCnstrIdx; }
		if (sigmaUBnd < sigma) { sigma = sigmaUBnd; addStatus = ST_UPPER; addBoundNotConstraint = BT_TRUE; addIdx = addUBndIdx; }
		if (sigmaLBnd < sigma) { sigma = sigmaLBnd; addStatus = ST_LOWER; addBoundNotConstraint = BT_TRUE; addIdx = addLBndIdx; }

		if (sigma >= options.maxPrimalJump)
		{
			unbounded = BT_TRUE;
			returnvalue = RET_HOTSTART_STOPPED_UNBOUNDEDNESS;
		}
		else
		{
			for (i = 0; i < nFR; i++)
				x[FR_idx[i]] += sigma * delta_xFR[i];

			for (i = 0; i < nFX; i++)
				x[FX_idx[i]] += sigma * delta_xFX[i];

			/* update Ax, Ax_u, and Ax_l */
			A->times(1, 1.0, x, nV, 0.0, Ax, nC);
			for (i = 0; i < nC; i++) Ax_u[i] = ubA[i] - Ax[i];
			for (i = 0; i < nC; i++) Ax_l[i] = Ax[i] - lbA[i];

			/* change working set later */
			exchangeHappened = BT_TRUE;
		}

		delete[] Ax_W;
		delete[] As;
	}

	delete[] delta_yFX;
	delete[] delta_yAC;
	delete[] delta_xFR;
	delete[] delta_xFX;
// 	delete[] delta_g;

	return returnvalue;
}



/*
 *	b a c k s o l v e T
 */
returnValue QProblem::backsolveT( const real_t* const b, BooleanType transposed, real_t* const a ) const
{
	int_t i, j;
	int_t nT = getNAC( );
	int_t tcol = sizeT - nT;

	real_t sum;

	/* nothing to do */
	if ( nT <= 0 )
		return SUCCESSFUL_RETURN;


	/* Solve Ta = b, where T might be transposed. */
	if ( transposed == BT_FALSE )
	{
		/* solve Ta = b */
		for( i=0; i<nT; ++i )
		{
			sum = b[i];
			for( j=0; j<i; ++j )
				sum -= TT(i,sizeT-1-j) * a[nT-1-j];

			if ( getAbs( TT(i,sizeT-1-i) ) > EPS )
				a[nT-1-i] = sum / TT(i,sizeT-1-i);
			else
				return THROWERROR( RET_DIV_BY_ZERO );
		}
	}
	else
	{
		/* solve T^T*a = b */
		for( i=0; i<nT; ++i )
		{
			sum = b[i];
			for( j=0; j<i; ++j )
				sum -= TT(nT-1-j,tcol+i) * a[nT-1-j];

			if ( getAbs( TT(nT-1-i,tcol+i) ) > EPS )
				a[nT-1-i] = sum / TT(nT-1-i,tcol+i);
			else
				return THROWERROR( RET_DIV_BY_ZERO );
		}
	}


	return SUCCESSFUL_RETURN;
}


/*
 *	d e t e r m i n e D a t a S h i f t
 */
returnValue QProblem::determineDataShift(	const real_t* const g_new, const real_t* const lbA_new, const real_t* const ubA_new,
											const real_t* const lb_new, const real_t* const ub_new,
											real_t* const delta_g, real_t* const delta_lbA, real_t* const delta_ubA,
											real_t* const delta_lb, real_t* const delta_ub,
											BooleanType& Delta_bC_isZero, BooleanType& Delta_bB_isZero
											)
{
	int_t i, ii;
	int_t nC  = getNC( );
	int_t nAC = getNAC( );

	int_t* FX_idx;
	int_t* AC_idx;

	bounds.getFixed( )->getNumberArray( &FX_idx );
	constraints.getActive( )->getNumberArray( &AC_idx );



	/* I) DETERMINE DATA SHIFT FOR BOUNDS */
	QProblemB::determineDataShift(	g_new,lb_new,ub_new,
									delta_g,delta_lb,delta_ub,
									Delta_bB_isZero );


	/* II) DETERMINE DATA SHIFT FOR CONSTRAINTS */
	/* 1) Calculate shift directions. */
	for( i=0; i<nC; ++i )
	{
		/* if lower constraints' bounds are to be disabled or do not exist, shift them to -infinity */
		if ( lbA_new != 0 )
			delta_lbA[i] = lbA_new[i] - lbA[i];
		else
			delta_lbA[i] = -INFTY - lbA[i];
	}

	for( i=0; i<nC; ++i )
	{
		/* if upper constraints' bounds are to be disabled or do not exist, shift them to infinity */
		if ( ubA_new != 0 )
			delta_ubA[i] = ubA_new[i] - ubA[i];
		else
			delta_ubA[i] = INFTY - ubA[i];
	}

	/* 2) Determine if active constraints' bounds are to be shifted. */
	Delta_bC_isZero = BT_TRUE;

	for ( i=0; i<nAC; ++i )
	{
		ii = AC_idx[i];

		if ( ( getAbs( delta_lbA[ii] ) > EPS ) || ( getAbs( delta_ubA[ii] ) > EPS ) )
		{
			Delta_bC_isZero = BT_FALSE;
			break;
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	d e t e r m i n e S t e p D i r e c t i o n
 */
returnValue QProblem::determineStepDirection(	const real_t* const delta_g, const real_t* const delta_lbA, const real_t* const delta_ubA,
												const real_t* const delta_lb, const real_t* const delta_ub,
												BooleanType Delta_bC_isZero, BooleanType Delta_bB_isZero,
												real_t* const delta_xFX, real_t* const delta_xFR,
												real_t* const delta_yAC, real_t* const delta_yFX
												)
{
	int_t i, j, ii, jj, r;
	int_t nV  = getNV( );
	int_t nFR = getNFR( );
	int_t nFX = getNFX( );
	int_t nAC = getNAC( );
	int_t nZ  = getNZ( );

	int_t* FR_idx;
	int_t* FX_idx;
	int_t* AC_idx;

	bounds.getFree( )->getNumberArray( &FR_idx );
	bounds.getFixed( )->getNumberArray( &FX_idx );
	constraints.getActive( )->getNumberArray( &AC_idx );


	/* I) DETERMINE delta_xFX (this is exact, does not need refinement) */
	if ( Delta_bB_isZero == BT_FALSE )
	{
		for( i=0; i<nFX; ++i )
		{
			ii = FX_idx[i];

			if ( bounds.getStatus( ii ) == ST_LOWER )
				delta_xFX[i] = delta_lb[ii];
			else
				delta_xFX[i] = delta_ub[ii];
		}
	}
	else
	{
		for( i=0; i<nFX; ++i )
			delta_xFX[i] = 0.0;
	}


	/* tempA and tempB hold the residuals in gFR and bA (= lbA or ubA)
	 * delta_xFR, delta_yAC hold the steps that get refined */
	for ( i=0; i<nFR; ++i )
	{
		ii = FR_idx[i];
		tempA[i] = delta_g[ii];
		delta_xFR[i] = 0.0;
	}
	for ( i=0; i<nAC; ++i )
		delta_yAC[i] = 0.0;
	if ( Delta_bC_isZero == BT_FALSE )
	{
		for ( i=0; i<nAC; ++i )
		{
			ii = AC_idx[i];
			if ( constraints.getStatus( ii ) == ST_LOWER )
				tempB[i] = delta_lbA[ii];
			else
				tempB[i] = delta_ubA[ii];
		}
	}
	else
	{
		for ( i=0; i<nAC; ++i )
			tempB[i] = 0.0;
	}

	/* Iterative refinement loop for delta_xFRz, delta_xFRy, delta_yAC */
	for ( r=0; r<=options.numRefinementSteps; ++r )
	{
		/* II) DETERMINE delta_xFR */
		if ( nFR > 0 )
		{
			for( i=0; i<nFR; ++i )
				delta_xFR_TMP[i] = 0.0;

			/* 1) Determine delta_xFRy. */
			if ( nAC > 0 )
			{
				if ( ( Delta_bC_isZero == BT_TRUE ) && ( Delta_bB_isZero == BT_TRUE ) )
				{
					for( i=0; i<nAC; ++i )
						delta_xFRy[i] = 0.0;
				}
				else
				{
					/* compute bA - A * delta_xFX. tempB already holds bA->
					 * in refinements r>=1, delta_xFX is exactly zero */
					if ( ( Delta_bB_isZero == BT_FALSE ) && ( r == 0 ) )
						A->times(constraints.getActive(), bounds.getFixed(), 1, -1.0, delta_xFX, nFX, 1.0, tempB, nAC);

					if ( backsolveT( tempB, BT_FALSE, delta_xFRy ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_STEPDIRECTION_FAILED_TQ );

					for( i=0; i<nFR; ++i )
					{
						ii = FR_idx[i];
						for( j=0; j<nAC; ++j )
							delta_xFR_TMP[i] += QQ(ii,nZ+j) * delta_xFRy[j];
					}
				}
			}


			/* 2) Determine delta_xFRz. */
			for( i=0; i<nZ; ++i )
				delta_xFRz[i] = 0.0;

			if ( ( hessianType == HST_ZERO ) || ( hessianType == HST_IDENTITY ) )
			{
				/* compute Z*delta_gFR [/eps] (delta_gFR is stored in tempA) */
				for( j=0; j<nFR; ++j )
				{
					jj = FR_idx[j];
					for( i=0; i<nZ; ++i )
						delta_xFRz[i] -= QQ(jj,i) * tempA[j];
				}

				if ( hessianType == HST_ZERO )
				{
					if ( usingRegularisation( ) == BT_TRUE )
					{
						for( i=0; i<nZ; ++i )
							delta_xFRz[i] /= regVal;
					}
					else
					{
						/* When solving LPs without regularisation, iterates must always be at a vertex. */
						if ( nZ > 0 )
							return THROWERROR( RET_UNKNOWN_BUG );
					}
				}
			}
			else
			{
				/* compute HMX*delta_xFX. DESTROY delta_gFR that was in tempA */
				if ( ( Delta_bB_isZero == BT_FALSE ) && ( r == 0 ) )
					H->times(bounds.getFree(), bounds.getFixed(), 1, 1.0, delta_xFX, nFX, 1.0, tempA, nFR);

				/* compute HFR*delta_xFRy */
				if ( ( nAC > 0 ) && ( ( Delta_bC_isZero == BT_FALSE ) || ( Delta_bB_isZero == BT_FALSE ) ) )
					H->times(bounds.getFree(), bounds.getFree(), 1, 1.0, delta_xFR_TMP, nFR, 1.0, tempA, nFR);

				/* compute ZFR_delta_xFRz = (Z'*HFR*Z) \ Z * (HFR*delta_xFR + HMX*delta_xFX + delta_gFR) */
				if ( nZ > 0 )
				{
					for( j=0; j<nFR; ++j )
					{
						jj = FR_idx[j];
						for( i=0; i<nZ; ++i )
							delta_xFRz[i] -= QQ(jj,i) * tempA[j];
					}

					if ( backsolveR( delta_xFRz,BT_TRUE,delta_xFRz ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_STEPDIRECTION_FAILED_CHOLESKY );

					if ( backsolveR( delta_xFRz,BT_FALSE,delta_xFRz ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_STEPDIRECTION_FAILED_CHOLESKY );
				}
			}

			/* compute Z * ZFR_delta_xFRz */
			if ( nZ > 0 )
			{
				for( i=0; i<nFR; ++i )
				{
					ZFR_delta_xFRz[i] = 0.0;

					ii = FR_idx[i];
					for( j=0; j<nZ; ++j )
						ZFR_delta_xFRz[i] += QQ(ii,j) * delta_xFRz[j];

					delta_xFR_TMP[i] += ZFR_delta_xFRz[i];
				}
			}
		}

		/* III) DETERMINE delta_yAC */
		if ( nAC > 0 ) /* => ( nFR = nZ + nAC > 0 ) */
		{
			if ( ( hessianType == HST_ZERO ) || ( hessianType == HST_IDENTITY ) )
			{
				/* if zero:     delta_yAC = (T')^-1 * ( Yfr*delta_gFR + eps*delta_xFRy ),
				 * if identity: delta_yAC = (T')^-1 * ( Yfr*delta_gFR +     delta_xFRy )
				 *
				 * DESTROY residual_bA that was stored in tempB
				 * If we come here, residual_gFR in tempA is STILL VALID
				 */
				if ( hessianType == HST_IDENTITY )
				{
					for( j=0; j<nAC; ++j )
						tempB[j] = delta_xFRy[j];
				}
				else /* hessianType == HST_ZERO */
				{
					if ( usingRegularisation( ) == BT_TRUE )
					{
						for( j=0; j<nAC; ++j )
							tempB[j] = regVal*delta_xFRy[j];
					}
					else
					{
						for( j=0; j<nAC; ++j )
							tempB[j] = 0.0;
					}
				}

				for( j=0; j<nAC; ++j )
				{
					for( i=0; i<nFR; ++i )
					{
						ii = FR_idx[i];
						tempB[j] += QQ(ii,nZ+j) * tempA[i];
					}
				}
			}
			else
			{
				/* Compute HFR * delta_xFR + HMX*delta_xFX
				 * Here, tempA holds (HFR*delta_xFRy + HMX*delta_xFX) */
				if ( nZ > 0 )
					H->times(bounds.getFree(), bounds.getFree(), 1, 1.0, ZFR_delta_xFRz, nFR, 1.0, tempA, nFR);

				for( i=0; i<nAC; ++i)
				{
					tempB[i] = 0.0;
					for( j=0; j<nFR; ++j )
					{
						jj = FR_idx[j];
						tempB[i] += QQ(jj,nZ+i) * tempA[j];
					}
				}
			}

			if ( backsolveT( tempB,BT_TRUE,delta_yAC_TMP ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_STEPDIRECTION_FAILED_TQ );
		}

		/* refine the solution found so far */
		for ( i=0; i<nFR; ++i )
			delta_xFR[i] += delta_xFR_TMP[i];
		for ( i=0; i<nAC; ++i )
			delta_yAC[i] += delta_yAC_TMP[i];

		if ( options.numRefinementSteps > 0 )
		{
			/* compute residuals in tempA and tempB, and max-norm */
			for ( i=0; i<nFR; ++i )
			{
				ii = FR_idx[i];
				tempA[i] = delta_g[ii];
			}

			switch ( hessianType )
			{
				case HST_ZERO:
					if ( usingRegularisation( ) == BT_TRUE )
						for ( i=0; i<nFR; ++i )
							tempA[i] += regVal*delta_xFR[i];
					break;

				case HST_IDENTITY:
					for ( i=0; i<nFR; ++i )
						tempA[i] += delta_xFR[i];
					break;

				default:
					H->times(bounds.getFree(), bounds.getFree(),  1, 1.0, delta_xFR, nFR, 1.0, tempA, nFR);
					H->times(bounds.getFree(), bounds.getFixed(), 1, 1.0, delta_xFX, nFX, 1.0, tempA, nFR);
					break;
			}

			A->transTimes(constraints.getActive(), bounds.getFree(), 1, -1.0, delta_yAC, nAC, 1.0, tempA, nFR);
			real_t rnrm = 0.0;
			for ( i=0; i<nFR; ++i )
				if (rnrm < getAbs (tempA[i]))
					rnrm = getAbs (tempA[i]);

			if (!Delta_bC_isZero)
			{
				for ( i=0; i<nAC; ++i )
				{
					ii = AC_idx[i];
					if ( constraints.getStatus( ii ) == ST_LOWER )
						tempB[i] = delta_lbA[ii];
					else
						tempB[i] = delta_ubA[ii];
				}
			}
			else
			{
				for ( i=0; i<nAC; ++i )
					tempB[i] = 0.0;
			}
			A->times(constraints.getActive(), bounds.getFree(), 1, -1.0, delta_xFR, nFR, 1.0, tempB, nAC);
			A->times(constraints.getActive(), bounds.getFixed(), 1, -1.0, delta_xFX, nFX, 1.0, tempB, nAC);
			for ( i=0; i<nAC; ++i )
				if (rnrm < getAbs (tempB[i]))
					rnrm = getAbs (tempB[i]);

			/* early termination of residual norm small enough */
			if ( rnrm < options.epsIterRef )
				break;
		}
	} /* end of refinement loop for delta_xFRz, delta_xFRy, delta_yAC */


	/* IV) DETERMINE delta_yFX */
	if ( nFX > 0 )
	{
		for( i=0; i<nFX; ++i )
			delta_yFX[i] = delta_g[FX_idx[i]];

		A->transTimes(constraints.getActive(), bounds.getFixed(), 1, -1.0, delta_yAC, nAC, 1.0, delta_yFX, nFX);

		switch( hessianType )
		{
			case HST_ZERO:
				if ( usingRegularisation( ) == BT_TRUE )
					for( i=0; i<nFX; ++i )
						delta_yFX[i] += regVal*delta_xFX[i];
				break;

			case HST_IDENTITY:
				for( i=0; i<nFX; ++i )
					delta_yFX[i] += 1.0 * delta_xFX[i];
				break;

			default:
				H->times(bounds.getFixed(), bounds.getFree(), 1, 1.0, delta_xFR, nFR, 1.0, delta_yFX, nFX);
				H->times(bounds.getFixed(), bounds.getFixed(), 1, 1.0, delta_xFX, nFX, 1.0, delta_yFX, nFX);
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	p e r f o r m S t e p
 */
returnValue QProblem::performStep(	const real_t* const delta_g,
									const real_t* const delta_lbA, const real_t* const delta_ubA,
									const real_t* const delta_lb, const real_t* const delta_ub,
									const real_t* const delta_xFX, const real_t* const delta_xFR,
									const real_t* const delta_yAC, const real_t* const delta_yFX,
									int_t& BC_idx, SubjectToStatus& BC_status, BooleanType& BC_isBound
									)
{
	int_t i, j, ii, jj;
	int_t nV  = getNV( );
	int_t nC  = getNC( );
	int_t nFR = getNFR( );
	int_t nFX = getNFX( );
	int_t nAC = getNAC( );
	int_t nIAC = getNIAC( );

	int_t* FR_idx;
	int_t* FX_idx;
	int_t* AC_idx;
	int_t* IAC_idx;

	bounds.getFree( )->getNumberArray( &FR_idx );
	bounds.getFixed( )->getNumberArray( &FX_idx );
	constraints.getActive( )->getNumberArray( &AC_idx );
	constraints.getInactive( )->getNumberArray( &IAC_idx );

	/* initialise maximum steplength array */
	tau = 1.0;
	BC_idx = -1;
	BC_status = ST_UNDEFINED;

	int_t BC_idx_tmp = -1;

	real_t* num = new real_t[ getMax( nV,nC ) ];
	real_t* den = new real_t[ getMax( nV,nC ) ];

	real_t* delta_Ax_l = new real_t[nC];
	real_t* delta_Ax_u = new real_t[nC];
	real_t* delta_Ax   = new real_t[nC];

	real_t* delta_x = new real_t[nV];
	for( j=0; j<nFR; ++j )
	{
		jj = FR_idx[j];
		delta_x[jj] = delta_xFR[j];
	}
	for( j=0; j<nFX; ++j )
	{
		jj = FX_idx[j];
		delta_x[jj] = delta_xFX[j];
	}


	/* I) DETERMINE MAXIMUM DUAL STEPLENGTH: */
	/* 1) Ensure that active dual constraints' bounds remain valid
	 *    (ignoring inequality constraints).  */
	for( i=0; i<nAC; ++i )
	{
		ii = AC_idx[i];

		num[i] = y[nV+ii];
		den[i] = -delta_yAC[i];
	}

	performRatioTest( nAC,AC_idx,&constraints, num,den, options.epsNum,options.epsDen, tau,BC_idx_tmp );

	if ( BC_idx_tmp >= 0 )
	{
		BC_idx = BC_idx_tmp;
		BC_status = ST_INACTIVE;
		BC_isBound = BT_FALSE;
	}


	/* 2) Ensure that active dual bounds remain valid
	 *    (ignoring implicitly fixed variables). */
	for( i=0; i<nFX; ++i )
	{
		ii = FX_idx[i];
		num[i] = y[ii];
		den[i] = -delta_yFX[i];
	}

	performRatioTest( nFX,FX_idx,&bounds, num,den, options.epsNum,options.epsDen, tau,BC_idx_tmp );

	if ( BC_idx_tmp >= 0 )
	{
		BC_idx = BC_idx_tmp;
		BC_status = ST_INACTIVE;
		BC_isBound = BT_TRUE;
	}


 	/* II) DETERMINE MAXIMUM PRIMAL STEPLENGTH */
	/* 1) Ensure that inactive constraints' bounds remain valid
	 *    (ignoring unbounded constraints). */

	/* calculate product A*x */
	if ( constraintProduct == 0 )
	{
		A->times(constraints.getInactive(), 0, 1, 1.0, delta_x, nV, 0.0, delta_Ax, nC, BT_FALSE);
	}
	else
	{
		for( i=0; i<nIAC; ++i )
		{
			ii = IAC_idx[i];

			if ( constraints.getType( ii ) != ST_UNBOUNDED )
			{
				if ( (*constraintProduct)( ii,delta_x, &(delta_Ax[ii]) ) != 0 )
				{
					delete[] den; delete[] num;
					delete[] delta_Ax; delete[] delta_Ax_u; delete[] delta_Ax_l; delete[] delta_x;
					return THROWERROR( RET_ERROR_IN_CONSTRAINTPRODUCT );
				}
			}
		}
	}

	if ( constraints.hasNoLower( ) == BT_FALSE )
	{
		for( i=0; i<nIAC; ++i )
		{
			ii = IAC_idx[i];
			num[i] = getMax( Ax_l[ii],0.0 );
			den[i] = delta_lbA[ii] - delta_Ax[ii];
		}

		performRatioTest( nIAC,IAC_idx,&constraints, num,den, options.epsNum,options.epsDen, tau,BC_idx_tmp );

		if ( BC_idx_tmp >= 0 )
		{
			BC_idx = BC_idx_tmp;
			BC_status = ST_LOWER;
			BC_isBound = BT_FALSE;
		}
	}

	if ( constraints.hasNoUpper( ) == BT_FALSE )
	{
		for( i=0; i<nIAC; ++i )
		{
			ii = IAC_idx[i];
			num[i] = getMax( Ax_u[ii],0.0 );
			den[i] = delta_Ax[ii] - delta_ubA[ii];
		}

		performRatioTest( nIAC,IAC_idx,&constraints, num,den, options.epsNum,options.epsDen, tau,BC_idx_tmp );

		if ( BC_idx_tmp >= 0 )
		{
			BC_idx = BC_idx_tmp;
			BC_status = ST_UPPER;
			BC_isBound = BT_FALSE;
		}
	}


	for( i=0; i<nIAC; ++i )
	{
		ii = IAC_idx[i];

		if ( constraints.getType( ii ) != ST_UNBOUNDED )
		{
			delta_Ax_l[ii] = delta_Ax[ii] - delta_lbA[ii];
			delta_Ax_u[ii] = delta_ubA[ii] - delta_Ax[ii];
		}
	}


	/* 2) Ensure that inactive bounds remain valid
	 *    (ignoring unbounded variables). */
	/* inactive lower bounds */
	if ( bounds.hasNoLower( ) == BT_FALSE )
	{
		for( i=0; i<nFR; ++i )
		{
			ii = FR_idx[i];
			num[i] = getMax( x[ii] - lb[ii],0.0 );
			den[i] = delta_lb[ii] - delta_xFR[i];
		}

		performRatioTest( nFR,FR_idx,&bounds, num,den, options.epsNum,options.epsDen, tau,BC_idx_tmp );

		if ( BC_idx_tmp >= 0 )
		{
			BC_idx = BC_idx_tmp;
			BC_status = ST_LOWER;
			BC_isBound = BT_TRUE;
		}
	}

	/* inactive upper bounds */
	if ( bounds.hasNoUpper( ) == BT_FALSE )
	{
		for( i=0; i<nFR; ++i )
		{
			ii = FR_idx[i];
			num[i] = getMax( ub[ii] - x[ii],0.0 );
			den[i] = delta_xFR[i] - delta_ub[ii];
		}

		performRatioTest( nFR,FR_idx,&bounds, num,den, options.epsNum,options.epsDen, tau,BC_idx_tmp );

		if ( BC_idx_tmp >= 0 )
		{
			BC_idx = BC_idx_tmp;
			BC_status = ST_UPPER;
			BC_isBound = BT_TRUE;
		}
	}

	delete[] den;
	delete[] num;
	delete[] delta_x;


	#ifndef __SUPPRESSANYOUTPUT__
	char messageString[MAX_STRING_LENGTH];

	if ( BC_status == ST_UNDEFINED )
		snprintf( messageString,MAX_STRING_LENGTH,"Stepsize is %.15e!",tau );
	else
		snprintf( messageString,MAX_STRING_LENGTH,"Stepsize is %.15e! (idx = %d, isBound = %d, status = %d)",tau,(int)BC_idx,(int)BC_isBound,(int)BC_status );

	getGlobalMessageHandler( )->throwInfo( RET_STEPSIZE_NONPOSITIVE,messageString,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
	#endif


	/* III) PERFORM STEP ALONG HOMOTOPY PATH */
	if ( tau > ZERO )
	{
		/* 1) Perform step in primal and dual space... */
		for( i=0; i<nFR; ++i )
		{
			ii = FR_idx[i];
			x[ii] += tau*delta_xFR[i];
		}

		for( i=0; i<nFX; ++i )
		{
			ii = FX_idx[i];
			x[ii] += tau*delta_xFX[i];
			y[ii] += tau*delta_yFX[i];
		}

		for( i=0; i<nAC; ++i )
		{
			ii = AC_idx[i];
			y[nV+ii] += tau*delta_yAC[i];
		}

		/* 2) Shift QP data. */
		for( i=0; i<nV; ++i )
		{
			g[i]  += tau*delta_g[i];
			lb[i] += tau*delta_lb[i];
			ub[i] += tau*delta_ub[i];
		}

		for( i=0; i<nC; ++i )
		{
			lbA[i] += tau*delta_lbA[i];
			ubA[i] += tau*delta_ubA[i];
		}

		A->times( constraints.getActive(),0, 1, 1.0, x, nV, 0.0, Ax, nC, BT_FALSE );
		for( i=0; i<nAC; ++i )
		{
			ii = AC_idx[i];
			Ax_u[ii] = ubA[ii] - Ax[ii];
			Ax_l[ii] = Ax[ii] - lbA[ii];
		}
		for( i=0; i<nIAC; ++i )
		{
			ii = IAC_idx[i];
			if ( constraints.getType( ii ) != ST_UNBOUNDED )
			{
				Ax[ii]   += tau*delta_Ax[ii];
				Ax_l[ii] += tau*delta_Ax_l[ii];
				Ax_u[ii] += tau*delta_Ax_u[ii];
			}
		}
	}
	else
	{
		/* print a stepsize warning if stepsize is zero */
		#ifndef __SUPPRESSANYOUTPUT__
		snprintf( messageString,MAX_STRING_LENGTH,"Stepsize is %.15e",tau );
		getGlobalMessageHandler( )->throwWarning( RET_STEPSIZE,messageString,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
		#endif
	}

	delete[] delta_Ax; delete[] delta_Ax_u; delete[] delta_Ax_l;

	return SUCCESSFUL_RETURN;
}


/*
 *	c h a n g e A c t i v e S e t
 */
returnValue QProblem::changeActiveSet( int_t BC_idx, SubjectToStatus BC_status, BooleanType BC_isBound )
{
	int_t nV = getNV( );

	#ifndef __SUPPRESSANYOUTPUT__
	char messageString[MAX_STRING_LENGTH];
	#endif

	switch ( BC_status )
	{
		/* Optimal solution found as no working set change detected. */
		case ST_UNDEFINED:
			return SUCCESSFUL_RETURN;

		/* Remove one variable from active set. */
		case ST_INACTIVE:
			if ( BC_isBound == BT_TRUE )
			{
				#ifndef __SUPPRESSANYOUTPUT__
				snprintf( messageString,MAX_STRING_LENGTH,"bound no. %d.",(int)BC_idx );
				getGlobalMessageHandler( )->throwInfo( RET_REMOVE_FROM_ACTIVESET,messageString,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
				#endif

				if ( removeBound( BC_idx,BT_TRUE,BT_TRUE,options.enableNZCTests ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_REMOVE_FROM_ACTIVESET_FAILED );

				y[BC_idx] = 0.0;
			}
			else
			{
				#ifndef __SUPPRESSANYOUTPUT__
				snprintf( messageString,MAX_STRING_LENGTH,"constraint no. %d.",(int)BC_idx );
				getGlobalMessageHandler( )->throwInfo( RET_REMOVE_FROM_ACTIVESET,messageString,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
				#endif

				if ( removeConstraint( BC_idx,BT_TRUE,BT_TRUE,options.enableNZCTests ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_REMOVE_FROM_ACTIVESET_FAILED );

				y[nV+BC_idx] = 0.0;
			}
			break;


		/* Add one variable to active set. */
		default:
			returnValue returnvalue;
			if ( BC_isBound == BT_TRUE )
			{
				#ifndef __SUPPRESSANYOUTPUT__
				if ( BC_status == ST_LOWER )
					snprintf( messageString,MAX_STRING_LENGTH,"lower bound no. %d.",(int)BC_idx );
				else
					snprintf( messageString,MAX_STRING_LENGTH,"upper bound no. %d.",(int)BC_idx );
				getGlobalMessageHandler( )->throwInfo( RET_ADD_TO_ACTIVESET,messageString,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
				#endif

				returnvalue = addBound( BC_idx,BC_status,BT_TRUE );
				if ( returnvalue == RET_ADDBOUND_FAILED_INFEASIBILITY )
					return returnvalue;
				if ( returnvalue != SUCCESSFUL_RETURN )
					return THROWERROR( RET_ADD_TO_ACTIVESET_FAILED );
			}
			else
			{
				#ifndef __SUPPRESSANYOUTPUT__
				if ( BC_status == ST_LOWER )
					snprintf( messageString,MAX_STRING_LENGTH,"lower constraint's bound no. %d.",(int)BC_idx );
				else
					snprintf( messageString,MAX_STRING_LENGTH,"upper constraint's bound no. %d.",(int)BC_idx );
				getGlobalMessageHandler( )->throwInfo( RET_ADD_TO_ACTIVESET,messageString,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
				#endif

				returnvalue = addConstraint( BC_idx,BC_status,BT_TRUE );
				if ( returnvalue == RET_ADDCONSTRAINT_FAILED_INFEASIBILITY )
					return returnvalue;
				if ( returnvalue != SUCCESSFUL_RETURN )
					return THROWERROR( RET_ADD_TO_ACTIVESET_FAILED );
			}
	}

	return SUCCESSFUL_RETURN;
}



/*
 * g e t R e l a t i v e H o m o t o p y L e n g t h
 */
real_t QProblem::getRelativeHomotopyLength(	const real_t* const g_new, const real_t* const lb_new, const real_t* const ub_new,
											const real_t* const lbA_new, const real_t* const ubA_new
											)
{
	int_t i;
	int_t nC = getNC( );
	real_t len = QProblemB::getRelativeHomotopyLength( g_new,lb_new,ub_new );
	real_t d, s;

	/*fprintf( stdFile, "len in homotopyLength = %.3e\n",len ); */

	/* lower constraint bounds */
	if ( lbA_new != 0 )
	{
		for (i = 0; i < nC; i++)
		{
			s = getAbs(lbA_new[i]);
			if (s < 1.0) s = 1.0;
			d = getAbs(lbA_new[i] - lbA[i]) / s;
			if (d > len) len = d;
		}
	}
	/*fprintf( stdFile, "len in homotopyLength = %.3e\n",len ); */

	/* upper constraint bounds */
	if ( ubA_new != 0 )
	{
		for (i = 0; i < nC; i++)
		{
			s = getAbs(ubA_new[i]);
			if (s < 1.0) s = 1.0;
			d = getAbs(ubA_new[i] - ubA[i]) / s;
			if (d > len) len = d;
		}
	}
	/*fprintf( stdFile, "len in homotopyLength = %.3e\n",len ); */

	return len;
}


/*
 * p e r f o r m R a m p i n g
 */
returnValue QProblem::performRamping( )
{
	int_t nV = getNV( ), nC = getNC( ), bstat, cstat, i, nRamp;
	real_t tP, rampValP, tD, rampValD, sca;

	/* compute number of values in ramp */
	nRamp = nV + nC + nC + nV;

	/* ramp inactive variable bounds and active dual bound variables */
	for (i = 0; i < nV; i++)
	{
		switch (bounds.getType(i))
		{
			case ST_EQUALITY:
				lb[i] = x[i]; ub[i] = x[i];  /* reestablish exact feasibility */
				continue;

			case ST_BOUNDED:
				tP = static_cast<real_t>((i+rampOffset) % nRamp) / static_cast<real_t>(nRamp-1);
				rampValP = (1.0-tP) * ramp0 + tP * ramp1;
				tD = static_cast<real_t>((nV+nC+nC+i+rampOffset) % nRamp) / static_cast<real_t>(nRamp-1);
				rampValD = (1.0-tD) * ramp0 + tD * ramp1;
				bstat = bounds.getStatus(i);
				if (bstat != ST_LOWER) { sca = getMax(getAbs(x[i]), 1.0); lb[i] = x[i] - sca * rampValP; }
				if (bstat != ST_UPPER) { sca = getMax(getAbs(x[i]), 1.0); ub[i] = x[i] + sca * rampValP; }
				if (bstat == ST_LOWER) { lb[i] = x[i]; y[i] = +rampValD; }
				if (bstat == ST_UPPER) { ub[i] = x[i]; y[i] = -rampValD; }
				if (bstat == ST_INACTIVE) y[i] = 0.0; /* reestablish exact complementarity */
				break;

			case ST_UNBOUNDED:
			case ST_DISABLED:
			default:
				 continue;
		}
	}

	/* ramp inactive constraints and active dual constraint variables */
	for (i = 0; i < nC; i++)
	{
		switch (constraints.getType(i))
		{
			case ST_EQUALITY:
				lbA[i] = Ax[i]; ubA[i] = Ax[i];  /* reestablish exact feasibility */
				continue;

			case ST_BOUNDED:
				tP = static_cast<real_t>((nV+i+rampOffset) % nRamp) / static_cast<real_t>(nRamp-1);
				rampValP = (1.0-tP) * ramp0 + tP * ramp1;
				tD = static_cast<real_t>((nV+nC+i+rampOffset) % nRamp) / static_cast<real_t>(nRamp-1);
				rampValD = (1.0-tD) * ramp0 + tD * ramp1;
				cstat = constraints.getStatus(i);
				if (cstat != ST_LOWER) { sca = getMax(getAbs(Ax[i]), 1.0); lbA[i] = Ax[i] - sca * rampValP; }
				if (cstat != ST_UPPER) { sca = getMax(getAbs(Ax[i]), 1.0); ubA[i] = Ax[i] + sca * rampValP; }
				if (cstat == ST_LOWER) { lbA[i] = Ax[i]; y[nV+i] = +rampValD; }
				if (cstat == ST_UPPER) { ubA[i] = Ax[i]; y[nV+i] = -rampValD; }
				if (cstat == ST_INACTIVE) y[nV+i] = 0.0; /* reestablish exact complementarity */

				Ax_l[i] = Ax[i] - lbA[i];
				Ax_u[i] = ubA[i] - Ax[i];
				break;

			case ST_UNBOUNDED:
			case ST_DISABLED:
			default:
				continue;
		}
	}

	/* reestablish exact stationarity */
	setupAuxiliaryQPgradient( );

	/* advance ramp offset to avoid Ramping cycles */
	rampOffset++;

	return SUCCESSFUL_RETURN;
}


/*
 * u p d a t e F a r B o u n d s
 */
returnValue QProblem::updateFarBounds(	real_t curFarBound, int_t nRamp,
										const real_t* const lb_new, real_t* const lb_new_far,
										const real_t* const ub_new, real_t* const ub_new_far,
										const real_t* const lbA_new, real_t* const lbA_new_far,
										const real_t* const ubA_new, real_t* const ubA_new_far
										) const
{
	int_t i;
	real_t rampVal, t;
	int_t nV = getNV( );
	int_t nC = getNC( );

    returnValue returnvalue = QProblemB::updateFarBounds(	curFarBound,nRamp,
															lb_new,lb_new_far, ub_new,ub_new_far
															);
    if ( returnvalue != SUCCESSFUL_RETURN )
        return returnvalue;

	if ( options.enableRamping == BT_TRUE )
	{
		for ( i=0; i<nC; ++i )
		{
			t = static_cast<real_t>((nV+i + rampOffset) % nRamp) / static_cast<real_t>(nRamp-1);
			rampVal = curFarBound * (1.0 + (1.0-t)*ramp0 + t*ramp1);

			if ( lbA_new == 0 )
				lbA_new_far[i] = -rampVal;
			else
				lbA_new_far[i] = getMax( -rampVal,lbA_new[i] );

			if ( ubA_new == 0 )
				ubA_new_far[i] = rampVal;
			else
				ubA_new_far[i] = getMin( rampVal,ubA_new[i] );
		}
	}
	else
	{
		for ( i=0; i<nC; ++i )
		{
			if ( lbA_new == 0 )
				lbA_new_far[i] = -curFarBound;
			else
				lbA_new_far[i] = getMax( -curFarBound,lbA_new[i] );

			if ( ubA_new == 0 )
				ubA_new_far[i] = curFarBound;
			else
				ubA_new_far[i] = getMin( curFarBound,ubA_new[i] );
		}
	}

	return SUCCESSFUL_RETURN;
}



/*
 * p e r f o r m D r i f t C o r r e c t i o n
 */
returnValue QProblem::performDriftCorrection( )
{
	int_t i;
	int_t nV = getNV ();
	int_t nC = getNC ();

	for ( i=0; i<nV; ++i )
	{
		switch ( bounds.getType ( i ) )
		{
			case ST_BOUNDED:
				switch ( bounds.getStatus ( i ) )
				{
					case ST_LOWER:
						lb[i] = x[i];
						ub[i] = getMax (ub[i], x[i]);
						y[i] = getMax (y[i], 0.0);
						break;
					case ST_UPPER:
						lb[i] = getMin (lb[i], x[i]);
						ub[i] = x[i];
						y[i] = getMin (y[i], 0.0);
						break;
					case ST_INACTIVE:
						lb[i] = getMin (lb[i], x[i]);
						ub[i] = getMax (ub[i], x[i]);
						y[i] = 0.0;
						break;
					case ST_UNDEFINED:
					case ST_INFEASIBLE_LOWER:
					case ST_INFEASIBLE_UPPER:
						break;
				}
				break;
			case ST_EQUALITY:
				lb[i] = x[i];
				ub[i] = x[i];
				break;
			case ST_UNBOUNDED:
			case ST_UNKNOWN:
            case ST_DISABLED:
				break;
		}
	}

	for ( i=0; i<nC; ++i )
	{
		switch ( constraints.getType ( i ) )
		{
			case ST_BOUNDED:
				switch ( constraints.getStatus ( i ) )
				{
					case ST_LOWER:
						lbA[i] = Ax[i];
						Ax_l[i] = 0.0;
						ubA[i] = getMax (ubA[i], Ax[i]);
						Ax_u[i] = ubA[i] - Ax[i];
						y[i+nV] = getMax (y[i+nV], 0.0);
						break;
					case ST_UPPER:
						lbA[i] = getMin (lbA[i], Ax[i]);
						Ax_l[i] = Ax[i] - lbA[i];
						ubA[i] = Ax[i];
						Ax_u[i] = 0.0;
						y[i+nV] = getMin (y[i+nV], 0.0);
						break;
					case ST_INACTIVE:
						lbA[i] = getMin (lbA[i], Ax[i]);
						Ax_l[i] = Ax[i] - lbA[i];
						ubA[i] = getMax (ubA[i], Ax[i]);
						Ax_u[i] = ubA[i] - Ax[i];
						y[i+nV] = 0.0;
						break;
					case ST_UNDEFINED:
					case ST_INFEASIBLE_LOWER:
					case ST_INFEASIBLE_UPPER:
						break;
				}
				break;
			case ST_EQUALITY:
				lbA[i] = Ax[i];
				Ax_l[i] = 0.0;
				ubA[i] = Ax[i];
				Ax_u[i] = 0.0;
				break;
			case ST_UNBOUNDED:
			case ST_UNKNOWN:
            case ST_DISABLED:
				break;
		}
	}

	return setupAuxiliaryQPgradient( );
}


/*
 *	s e t u p A u x i l i a r y Q P
 */
returnValue QProblem::setupAuxiliaryQP( const Bounds* const guessedBounds, const Constraints* const guessedConstraints )
{
	int_t i, j;
	int_t nV = getNV( );
	int_t nC = getNC( );

	/* consistency check */
	if ( ( guessedBounds == 0 ) || ( guessedConstraints == 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* nothing to do */
	if ( ( guessedBounds == &bounds ) && ( guessedConstraints == &constraints ) )
		return SUCCESSFUL_RETURN;

	status = QPS_PREPARINGAUXILIARYQP;


	/* I) SETUP WORKING SET ... */
	if ( shallRefactorise( guessedBounds,guessedConstraints ) == BT_TRUE )
	{
		/* ... WITH REFACTORISATION: */
		/* 1) Reset bounds/constraints ... */
		bounds.init( nV );
		constraints.init( nC );

		/*    ... and set them up afresh. */
		if ( setupSubjectToType( ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		if ( bounds.setupAllFree( ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		if ( constraints.setupAllInactive( ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		/* 2) Setup TQ factorisation. */
		if ( setupTQfactorisation( ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		/* 3) Setup guessed working sets afresh. */
		if ( setupAuxiliaryWorkingSet( guessedBounds,guessedConstraints,BT_TRUE ) != SUCCESSFUL_RETURN )
			THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		/* 4) Computes Cholesky decomposition of projected Hessian
		 *    This now handles all special cases (no active bounds/constraints, no nullspace) */
		if ( computeProjectedCholesky( ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
	}
	else
	{
		/* ... WITHOUT REFACTORISATION: */
		if ( setupAuxiliaryWorkingSet( guessedBounds,guessedConstraints,BT_FALSE ) != SUCCESSFUL_RETURN )
			THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
	}


	/* II) SETUP AUXILIARY QP DATA: */
	/* 1) Ensure that dual variable is zero for free bounds and inactive constraints. */
	for ( i=0; i<nV; ++i )
		if ( bounds.getStatus( i ) == ST_INACTIVE )
			y[i] = 0.0;

	for ( i=0; i<nC; ++i )
		if ( constraints.getStatus( i ) == ST_INACTIVE )
			y[nV+i] = 0.0;

	/* 2) Setup gradient and (constraints') bound vectors. */
	if ( setupAuxiliaryQPgradient( ) != SUCCESSFUL_RETURN )
		THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

	A->times(1, 1.0, x, nV, 0.0, Ax, nC);
	for ( j=0; j<nC; ++j )
	{
		Ax_l[j] = Ax[j];
		Ax_u[j] = Ax[j];
	}

	/* (also sets Ax_l and Ax_u) */
	if ( setupAuxiliaryQPbounds( 0,0,BT_FALSE ) != SUCCESSFUL_RETURN )
		THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

	return SUCCESSFUL_RETURN;
}


/*
 *	s h a l l R e f a c t o r i s e
 */

BooleanType QProblem::shallRefactorise(	const Bounds* const guessedBounds,
										const Constraints* const guessedConstraints
										) const
{
	int_t i;
	int_t nV = getNV( );
	int_t nC = getNC( );

	/* always refactorise if Hessian is not known to be positive definite */
	if ( ( hessianType == HST_SEMIDEF ) || ( hessianType == HST_INDEF ) )
		return BT_TRUE;

	/* 1) Determine number of bounds that have same status
	 *    in guessed AND current bounds.*/
	int_t differenceNumberBounds = 0;

	for( i=0; i<nV; ++i )
		if ( guessedBounds->getStatus( i ) != bounds.getStatus( i ) )
			++differenceNumberBounds;

	/* 2) Determine number of constraints that have same status
	 *    in guessed AND current constraints.*/
	int_t differenceNumberConstraints = 0;

	for( i=0; i<nC; ++i )
		if ( guessedConstraints->getStatus( i ) != constraints.getStatus( i ) )
			++differenceNumberConstraints;

	/* 3) Decide wheter to refactorise or not. */
	if ( 2*(differenceNumberBounds+differenceNumberConstraints) > guessedConstraints->getNAC( )+guessedBounds->getNFX( ) )
		return BT_TRUE;
	else
		return BT_FALSE;
}


/*
 *	s e t u p Q P d a t a
 */
returnValue QProblem::setupQPdata(	SymmetricMatrix *_H, const real_t* const _g, Matrix *_A,
									const real_t* const _lb, const real_t* const _ub,
									const real_t* const _lbA, const real_t* const _ubA
									)
{
	int_t nC = getNC( );

#ifdef __WRITE_DATA_FILES__
	{
		int_t i;
		const double Infinity = 1e20;
		int_t nV = getNV( );
		GlobalOutputFileCounter++;
		char buf[256];
		snprintf(buf,256,"QP%d_setupQPdata.dat",GlobalOutputFileCounter);
		MyPrintf("+++ Writing output file %s\n", buf);

		FILE* output_file = fopen(buf,"w");

		fprintf(output_file,"nVar = %d\n", nV);
		fprintf(output_file,"nCon = %d\n", nC);

		_H->writeToFile(output_file,"H_");
		for (i=0; i<nV; i++) {
			fprintf(output_file,"g[%d] = %23.16e\n",i,_g[i]);
		}
		_A->writeToFile(output_file,"A_");
		for (i=0; i<nV; i++) {
			fprintf(output_file,"lb[%d] = %23.16e\n",i,getMax(-Infinity,_lb[i]));
		}
		for (i=0; i<nV; i++) {
			fprintf(output_file,"ub[%d] = %23.16e\n",i,getMin(Infinity,_ub[i]));
		}

		for (i=0; i<nC; i++) {
			fprintf(output_file,"lbA[%d] = %23.16e\n",i,getMax(-Infinity,_lbA[i]));
		}
		for (i=0; i<nC; i++) {
			fprintf(output_file,"ubA[%d] = %23.16e\n",i,getMin(Infinity,_ubA[i]));
		}
		fclose(output_file);
  }
#endif

	/* 1) Load Hessian matrix as well as lower and upper bounds vectors. */
	if ( QProblemB::setupQPdata( _H,_g,_lb,_ub ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	if ( ( nC > 0 ) && ( _A == 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	if ( nC > 0 )
	{
		/* 2) Setup lower/upper constraints' bounds vector. */
		setLBA( _lbA );
		setUBA( _ubA );

		/* 3) Only load constraint matrix after setting up vectors! */
		setA( _A );
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p Q P d a t a
 */
returnValue QProblem::setupQPdata(	const real_t* const _H, const real_t* const _g, const real_t* const _A,
									const real_t* const _lb, const real_t* const _ub,
									const real_t* const _lbA, const real_t* const _ubA
									)
{
	int_t nC = getNC( );

	/* 1) Load Hessian matrix as well as lower and upper bounds vectors. */
	if ( QProblemB::setupQPdata( _H,_g,_lb,_ub ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	if ( ( nC > 0 ) && ( _A == 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	if ( nC > 0 )
	{
		/* 2) Setup lower/upper constraints' bounds vector. */
		setLBA( _lbA );
		setUBA( _ubA );

		/* 3) Only load constraint matrix after setting up vectors! */
		setA( _A );
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p Q P d a t a F r o m F i l e
 */
returnValue QProblem::setupQPdataFromFile(	const char* const H_file, const char* const g_file, const char* const A_file,
											const char* const lb_file, const char* const ub_file,
											const char* const lbA_file, const char* const ubA_file
											)
{
	int_t i;
	int_t nV = getNV( );
	int_t nC = getNC( );

	returnValue returnvalue;


	/* 1) Load Hessian matrix as well as lower and upper bounds vectors from files. */
	returnvalue = QProblemB::setupQPdataFromFile( H_file,g_file,lb_file,ub_file );
	if ( returnvalue != SUCCESSFUL_RETURN )
		return THROWERROR( returnvalue );

	if ( ( nC > 0 ) && ( A_file == 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	if ( nC > 0 )
	{
		/* 2) Load lower constraints' bounds vector from file. */
		if ( lbA_file != 0 )
		{
			returnvalue = readFromFile( lbA, nC, lbA_file );
			if ( returnvalue != SUCCESSFUL_RETURN )
				return THROWERROR( returnvalue );
		}
		else
		{
			/* if no lower constraints' bounds are specified, set them to -infinity */
			for( i=0; i<nC; ++i )
				lbA[i] = -INFTY;
		}

		/* 3) Load upper constraints' bounds vector from file. */
		if ( ubA_file != 0 )
		{
			returnvalue = readFromFile( ubA, nC, ubA_file );
			if ( returnvalue != SUCCESSFUL_RETURN )
				return THROWERROR( returnvalue );
		}
		else
		{
			/* if no upper constraints' bounds are specified, set them to infinity */
			for( i=0; i<nC; ++i )
				ubA[i] = INFTY;
		}

		/* 4) Only load constraint matrix from file after setting up vectors! */
		real_t* _A = new real_t[nC * nV];
		returnvalue = readFromFile( _A, nC,nV, A_file );
		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			delete[] _A;
			return THROWERROR( returnvalue );
		}
		setA( _A );
		A->doFreeMemory( );
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	l o a d Q P v e c t o r s F r o m F i l e
 */
returnValue QProblem::loadQPvectorsFromFile(	const char* const g_file, const char* const lb_file, const char* const ub_file,
												const char* const lbA_file, const char* const ubA_file,
												real_t* const g_new, real_t* const lb_new, real_t* const ub_new,
												real_t* const lbA_new, real_t* const ubA_new
												) const
{
	int_t nC = getNC( );

	returnValue returnvalue;


	/* 1) Load gradient vector as well as lower and upper bounds vectors from files. */
	returnvalue = QProblemB::loadQPvectorsFromFile( g_file,lb_file,ub_file, g_new,lb_new,ub_new );
	if ( returnvalue != SUCCESSFUL_RETURN )
		return THROWERROR( returnvalue );

	if ( nC > 0 )
	{
		/* 2) Load lower constraints' bounds vector from file. */
		if ( lbA_file != 0 )
		{
			if ( lbA_new != 0 )
			{
				returnvalue = readFromFile( lbA_new, nC, lbA_file );
				if ( returnvalue != SUCCESSFUL_RETURN )
					return THROWERROR( returnvalue );
			}
			else
			{
				/* If filename is given, storage must be provided! */
				return THROWERROR( RET_INVALID_ARGUMENTS );
			}
		}

		/* 3) Load upper constraints' bounds vector from file. */
		if ( ubA_file != 0 )
		{
			if ( ubA_new != 0 )
			{
				returnvalue = readFromFile( ubA_new, nC, ubA_file );
				if ( returnvalue != SUCCESSFUL_RETURN )
					return THROWERROR( returnvalue );
			}
			else
			{
				/* If filename is given, storage must be provided! */
				return THROWERROR( RET_INVALID_ARGUMENTS );
			}
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t I t e r a t i o n
 */
returnValue QProblem::printIteration( 	int_t iter,
										int_t BC_idx,	SubjectToStatus BC_status, BooleanType BC_isBound, real_t homotopyLength,
										BooleanType isFirstCall
		  								)
{
	#ifndef __SUPPRESSANYOUTPUT__

	/* consistency check */
	if ( iter < 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	int_t i;
	int_t nV = getNV();
	int_t nC = getNC();
	int_t nAC = getNAC();

	real_t stat, bfeas, cfeas, bcmpl, ccmpl, Tmaxomin;
	real_t *grad = 0;
	real_t *AX = 0;
	real_t Tmin, Tmax;

	char myPrintfString[MAX_STRING_LENGTH];
	char info[MAX_STRING_LENGTH];
	const char excStr[] = " ef";

	switch ( options.printLevel )
	{
		case PL_DEBUG_ITER:
			grad = new real_t[nV];
			AX = new real_t[nC];
			stat = bfeas = cfeas = bcmpl = ccmpl = Tmaxomin = 0.0;

			/* stationarity */
			for (i = 0; i < nV; i++) grad[i] = g[i] - y[i];

			switch ( hessianType )
			{
				case HST_ZERO:
					for( i=0; i<nV; ++i )
						grad[i] += regVal * x[i];
					break;

				case HST_IDENTITY:
					for( i=0; i<nV; ++i )
						grad[i] += 1.0 * x[i];
					break;

				default:
					H->times(1, 1.0, x, nV, 1.0, grad, nV);
					break;
			}

			A->transTimes(1, -1.0, y+nV, nC, 1.0, grad, nV);
			for (i = 0; i < nV; i++) if (getAbs(grad[i]) > stat) stat = getAbs(grad[i]);

			/* feasibility */
			for (i = 0; i < nV; i++) if (lb[i] - x[i] > bfeas) bfeas = lb[i] - x[i];
			for (i = 0; i < nV; i++) if (x[i] - ub[i] > bfeas) bfeas = x[i] - ub[i];
			A->times(1, 1.0, x, nV, 0.0, AX, nC);
			for (i = 0; i < nC; i++) if (lbA[i] - AX[i] > cfeas) cfeas = lbA[i] - AX[i];
			for (i = 0; i < nC; i++) if (AX[i] - ubA[i] > cfeas) cfeas = AX[i] - ubA[i];

			/* complementarity */
			for (i = 0; i < nV; i++) if (y[i] > +EPS && getAbs((lb[i] - x[i])*y[i]) > bcmpl) bcmpl = getAbs((lb[i] - x[i])*y[i]);
			for (i = 0; i < nV; i++) if (y[i] < -EPS && getAbs((ub[i] - x[i])*y[i]) > bcmpl) bcmpl = getAbs((ub[i] - x[i])*y[i]);
			for (i = 0; i < nC; i++) if (y[nV+i] > +EPS && getAbs((lbA[i]-AX[i])*y[nV+i]) > ccmpl) ccmpl = getAbs((lbA[i]-AX[i])*y[nV+i]);
			for (i = 0; i < nC; i++) if (y[nV+i] < -EPS && getAbs((ubA[i]-AX[i])*y[nV+i]) > ccmpl) ccmpl = getAbs((ubA[i]-AX[i])*y[nV+i]);

			Tmin = 1.0e16; Tmax = 0.0;
			for (i = 0; i < nAC; i++)
				if (getAbs(TT(i,sizeT-i-1)) < Tmin)
					Tmin = getAbs(TT(i,sizeT-i-1));
				else if (getAbs(TT(i,sizeT-i-1)) > Tmax)
					Tmax = getAbs(TT(i,sizeT-i-1));
			Tmaxomin = Tmax/Tmin;

			if ( (iter % 10 == 0) && ( isFirstCall == BT_TRUE ) )
			{
				snprintf( myPrintfString,MAX_STRING_LENGTH, "\n%5s %4s %4s %4s %4s %9s %9s %9s %9s %9s %9s %9s %9s\n",
						"iter", "addB", "remB", "addC", "remC", "hom len", "tau", "stat",
						"bfeas", "cfeas", "bcmpl", "ccmpl", "Tmin");
				myPrintf( myPrintfString );
			}

			if ( isFirstCall == BT_TRUE )
				snprintf( myPrintfString,MAX_STRING_LENGTH, "%5d ",(int)iter );
			else
				snprintf( myPrintfString,MAX_STRING_LENGTH, "%5d*",(int)iter );
			myPrintf( myPrintfString );

			if (tabularOutput.idxAddB >= 0)
			{
				snprintf( myPrintfString,MAX_STRING_LENGTH, "%4d ",(int)(tabularOutput.idxAddB) );
				myPrintf( myPrintfString );
			}
			else
			{
				myPrintf( "     " );
			}

			if (tabularOutput.idxRemB >= 0)
			{
				snprintf( myPrintfString,MAX_STRING_LENGTH, "%4d ",(int)(tabularOutput.idxRemB) );
				myPrintf( myPrintfString );
			}
			else
			{
				myPrintf( "     " );
			}

			if (tabularOutput.idxAddC >= 0)
			{
				snprintf( myPrintfString,MAX_STRING_LENGTH, "%4d ",(int)(tabularOutput.idxAddC) );
				myPrintf( myPrintfString );
			}
			else
			{
				myPrintf( "     " );
			}

			if (tabularOutput.idxRemC >= 0)
			{
				snprintf( myPrintfString,MAX_STRING_LENGTH, "%4d ",(int)(tabularOutput.idxRemC) );
				myPrintf( myPrintfString );
			}
			else
			{
				myPrintf( "     " );
			}

			snprintf( myPrintfString,MAX_STRING_LENGTH, "%9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n",
					homotopyLength, tau, stat, bfeas, cfeas, bcmpl, ccmpl, Tmin);
			myPrintf( myPrintfString );

			delete[] AX;
			delete[] grad;
			break;

		case PL_TABULAR:
			if ( (iter % 10 == 0) && ( isFirstCall == BT_TRUE ) )
			{
				snprintf( myPrintfString,MAX_STRING_LENGTH, "\n%5s %6s %6s %6s %6s %9s %9s\n",
						"iter", "addB", "remB", "addC", "remC", "hom len", "tau" );
				myPrintf( myPrintfString );
			}
			if ( isFirstCall == BT_TRUE )
				snprintf( myPrintfString,MAX_STRING_LENGTH, "%5d ",(int)iter);
			else
				snprintf( myPrintfString,MAX_STRING_LENGTH, "%5d*",(int)iter);
			myPrintf( myPrintfString );

			if (tabularOutput.idxAddB >= 0)
			{
				snprintf( myPrintfString,MAX_STRING_LENGTH, "%5d%c ",(int)(tabularOutput.idxAddB), excStr[tabularOutput.excAddB]);
				myPrintf( myPrintfString );
			}
			else
			{
				myPrintf( "       " );
			}

			if (tabularOutput.idxRemB >= 0)
			{
				snprintf( myPrintfString,MAX_STRING_LENGTH, "%5d%c ",(int)(tabularOutput.idxRemB), excStr[tabularOutput.excRemB]);
				myPrintf( myPrintfString );
			}
			else
			{
				myPrintf( "       " );
			}

			if (tabularOutput.idxAddC >= 0)
			{
				snprintf( myPrintfString,MAX_STRING_LENGTH, "%5d%c ",(int)(tabularOutput.idxAddC), excStr[tabularOutput.excAddC]);
				myPrintf( myPrintfString );
			}
			else
			{
				myPrintf( "       " );
			}

			if (tabularOutput.idxRemC >= 0)
			{
				snprintf( myPrintfString,MAX_STRING_LENGTH, "%5d%c ",(int)(tabularOutput.idxRemC), excStr[tabularOutput.excRemC]);
				myPrintf( myPrintfString );
			}
			else
			{
				myPrintf( "       " );
			}

			snprintf( myPrintfString,MAX_STRING_LENGTH, "%9.2e %9.2e\n", homotopyLength, tau);
			myPrintf( myPrintfString );
			break;

		case PL_MEDIUM:
			/* 1) Print header at first iteration. */
 			if ( ( iter == 0 ) && ( isFirstCall == BT_TRUE ) )
			{
				snprintf( myPrintfString,MAX_STRING_LENGTH,"\n\n####################   qpOASES  --  QP NO. %3.0d   #####################\n\n",(int)count );
				myPrintf( myPrintfString );

				myPrintf( "    Iter   |    StepLength    |       Info       |   nFX   |   nAC    \n" );
				myPrintf( " ----------+------------------+------------------+---------+--------- \n" );
			}

			/* 2) Print iteration line. */
			if ( BC_status == ST_UNDEFINED )
			{
				if ( hessianType == HST_ZERO )
					snprintf( info,3,"LP" );
				else
					snprintf( info,3,"QP" );

				if ( isFirstCall == BT_TRUE )
					snprintf( myPrintfString,MAX_STRING_LENGTH,"   %5.1d   |   %1.6e   |    %s SOLVED     |  %4.1d   |  %4.1d   \n", (int)iter,tau,info,(int)getNFX( ),(int)getNAC( ) );
				else
					snprintf( myPrintfString,MAX_STRING_LENGTH,"   %5.1d*  |   %1.6e   |    %s SOLVED     |  %4.1d   |  %4.1d   \n", (int)iter,tau,info,(int)getNFX( ),(int)getNAC( ) );
				myPrintf( myPrintfString );
			}
			else
			{
				if ( BC_status == ST_INACTIVE )
					snprintf( info,5,"REM " );
				else
					snprintf( info,5,"ADD " );

				if ( BC_isBound == BT_TRUE )
					snprintf( &(info[4]),4,"BND" );
				else
					snprintf( &(info[4]),4,"CON" );

				snprintf( myPrintfString,MAX_STRING_LENGTH,"   %5.1d   |   %1.6e   |   %s %4.1d   |  %4.1d   |  %4.1d   \n", (int)iter,tau,info,(int)BC_idx,(int)getNFX( ),(int)getNAC( ) );
				myPrintf( myPrintfString );
			}
			break;

		default:
			/* nothing to display */
			break;
	}

	#endif /* __SUPPRESSANYOUTPUT__ */

	return SUCCESSFUL_RETURN;
}


inline real_t abs (real_t x) { return (x>0)?x:-x; }

/*
 * d r o p I n f e a s i b l e s
 */
returnValue QProblem::dropInfeasibles( int_t BC_number, SubjectToStatus BC_status, BooleanType BC_isBound,
										real_t *xiB, real_t *xiC )
{
	int_t i;

	int_t nAC                   = getNAC ();
	int_t nFX                   = getNFX ();
	int_t blockingPriority      = (BC_isBound) ? options.dropBoundPriority : options.dropIneqConPriority;
	int_t y_min_number          = -1;
	BooleanType y_min_isBound = BC_isBound;
	int_t y_min_priority        = blockingPriority;

	int_t* AC_idx;
	constraints.getActive( )->getNumberArray( &AC_idx );

	int_t* FX_idx;
	bounds.getFixed( )->getNumberArray( &FX_idx );

	if (options.dropEqConPriority <= y_min_priority)
	{
		// look for an equality constraint we can drop according to priorities
		for ( i = 0; i < nAC; ++i )
			if ( (constraints.getType (i) == ST_EQUALITY)
				&& (getAbs (xiC[i]) > options.epsDen) )
			{
				y_min_number = AC_idx[i];
				y_min_isBound = BT_FALSE;
				y_min_priority = options.dropEqConPriority;
				break;
			}
	}

	if (options.dropIneqConPriority <= y_min_priority)
	{
		// look for an inequality constraint we can drop according to priorities
		for ( i = 0; i < nAC; ++i )
			if ( (constraints.getType (i) == ST_BOUNDED)
				&& (getAbs (xiC[i]) > options.epsDen) )
			{
				y_min_number = AC_idx[i];
				y_min_isBound = BT_FALSE;
				y_min_priority = options.dropIneqConPriority;
				break;
			}
	}

	if (options.dropBoundPriority <= y_min_priority)
	{
		// look for a simple bound we can drop according to priorities
		for ( i = 0; i < nFX; ++i )
			if (getAbs (xiB[i]) > options.epsDen)
			{
				y_min_number = FX_idx[i];
				y_min_isBound = BT_TRUE;
				y_min_priority = options.dropBoundPriority;
				break;
			}
	}

	if (y_min_number >= 0) {

		// drop active equality or active bound we have found
		if (y_min_isBound) {
			SubjectToStatus status_ = bounds.getStatus (y_min_number);
			removeBound (y_min_number, BT_TRUE, BT_FALSE, BT_FALSE);
			bounds.setStatus (y_min_number, (status_ == ST_LOWER) ? ST_INFEASIBLE_LOWER : ST_INFEASIBLE_UPPER);
			// TODO: fix duals y[]
			/* fprintf (stdFile, "Dropping bounds %d for %s %d\n", y_min_number, BC_isBound?"bound":"constraint", BC_number); */
		} else {
			SubjectToStatus status_ = constraints.getStatus (y_min_number);
			removeConstraint (y_min_number, BT_TRUE, BT_FALSE, BT_FALSE);
			constraints.setStatus (y_min_number, (status_ == ST_LOWER) ? ST_INFEASIBLE_LOWER : ST_INFEASIBLE_UPPER);
			// TODO: fix duals y[]
			/* fprintf (stdFile, "Dropping constraint %d for %s %d\n", y_min_number, BC_isBound?"bound":"constraint", BC_number); */
		}

		// ... now return, add the blocking constraint, and continue solving QP with dropped bound/constraint
		return SUCCESSFUL_RETURN;

	} else {

		// nothing found, then drop the blocking (still inactive) constraint
		if (BC_isBound)
			bounds.setStatus (BC_number, (BC_status == ST_LOWER) ? ST_INFEASIBLE_LOWER : ST_INFEASIBLE_UPPER);
		else
			constraints.setStatus (BC_number, (BC_status == ST_LOWER) ? ST_INFEASIBLE_LOWER : ST_INFEASIBLE_UPPER);

		/* fprintf (stdFile, "Dropping %s %d itself\n", BC_isBound?"bound":"constraint", BC_number); */

		// ... now return, and continue solving QP with dropped bound/constraint
		return RET_ENSURELI_DROPPED;
	}
}



/*
 *  w r i t e Q p D a t a I n t o M a t F i l e
 */
returnValue QProblem::writeQpDataIntoMatFile(	const char* const filename
												) const
{
	#ifndef __SUPPRESSANYOUTPUT__

	FILE* matFile;
	matFile = fopen( filename,"w+" );

	if ( matFile == 0 )
		return RET_UNABLE_TO_OPEN_FILE;

	int_t nV = getNV();
	int_t nC = getNC();

	real_t* Hfull = H->full();
	writeIntoMatFile( matFile, Hfull, nV,nV, "H"   );
	delete[] Hfull;

	writeIntoMatFile( matFile, g,     nV,1,  "g"   );

	real_t* Afull = A->full();
	writeIntoMatFile( matFile, Afull, nC,nV, "A"   );
	delete[] Afull;

	writeIntoMatFile( matFile, lb,    nV,1,  "lb"  );
	writeIntoMatFile( matFile, ub,    nV,1,  "ub"  );
	writeIntoMatFile( matFile, lbA,   nC,1,  "lbA" );
	writeIntoMatFile( matFile, ubA,   nC,1,  "ubA" );

	fclose( matFile );

	return SUCCESSFUL_RETURN;
    
	#else /* __SUPPRESSANYOUTPUT__ */

	return RET_NOT_YET_IMPLEMENTED;

	#endif /* __SUPPRESSANYOUTPUT__ */
}


/*
 *  w r i t e Q p W o r k s p a c e I n t o M a t F i l e
 */
returnValue QProblem::writeQpWorkspaceIntoMatFile(	const char* const filename
													)
{
	#ifndef __SUPPRESSANYOUTPUT__

	FILE* matFile;
	matFile = fopen( filename,"w+" );

	if ( matFile == 0 )
		return RET_UNABLE_TO_OPEN_FILE;

	int_t nV = getNV();
	int_t nC = getNC();
	int_t nFR  = getNFR();
	int_t nFX  = getNFX();
	int_t nAC  = getNAC();
	int_t nIAC = getNIAC();


	writeIntoMatFile( matFile, T, sizeT,sizeT, "T" );
	writeIntoMatFile( matFile, Q, nV,nV, "Q" );

	writeIntoMatFile( matFile, Ax, nC,1, "Ax" );
	writeIntoMatFile( matFile, Ax_l, nC,1, "Ax_l" );
	writeIntoMatFile( matFile, Ax_u, nC,1, "Ax_u" );


	int_t *FR_idx, *FX_idx, *AC_idx, *IAC_idx;
	bounds.getFree( )->getNumberArray( &FR_idx );
	bounds.getFixed( )->getNumberArray( &FX_idx );
	constraints.getActive( )->getNumberArray( &AC_idx );
	constraints.getInactive( )->getNumberArray( &IAC_idx );

	writeIntoMatFile( matFile, FR_idx,  nFR, 1, "FR_idx"  );
	writeIntoMatFile( matFile, FX_idx,  nFX, 1, "FX_idx"  );
	writeIntoMatFile( matFile, AC_idx,  nAC, 1, "AC_idx"  );
	writeIntoMatFile( matFile, IAC_idx, nIAC,1, "IAC_idx" );

	fclose( matFile );

	return SUCCESSFUL_RETURN;
    
	#else /* __SUPPRESSANYOUTPUT__ */

	return RET_NOT_YET_IMPLEMENTED;

	#endif /* __SUPPRESSANYOUTPUT__ */
}


END_NAMESPACE_QPOASES


/*
 *	end of file
 */
