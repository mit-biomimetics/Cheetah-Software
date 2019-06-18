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
 *	\file src/QProblemB.cpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.2
 *	\date 2007-2017
 *
 *	Implementation of the QProblemB class which is able to use the newly
 *	developed online active set strategy for parametric quadratic programming.
 */


#include <qpOASES/QProblemB.hpp>
#include <qpOASES/LapackBlasReplacement.hpp>


BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/


/*
 *	Q P r o b l e m B
 */
QProblemB::QProblemB( )
{
	/* print copyright notice */
	if (options.printLevel != PL_NONE)
		printCopyrightNotice( );

	/* reset global message handler */
	getGlobalMessageHandler( )->reset( );

	freeHessian = BT_FALSE;
	H = 0;

	g = 0;
	lb = 0;
	ub = 0;

	R = 0;
	haveCholesky = BT_FALSE;

	x = 0;
	y = 0;

	tau = 0.0;

	hessianType = HST_UNKNOWN;
	regVal = 0.0;

	infeasible  = BT_FALSE;
	unbounded   = BT_FALSE;

	status = QPS_NOTINITIALISED;

	count = 0;

	ramp0 = options.initialRamping;
	ramp1 = options.finalRamping;
	rampOffset = 0;

	delta_xFR_TMP = 0;

	setPrintLevel( options.printLevel );
}


/*
 *	Q P r o b l e m B
 */
QProblemB::QProblemB( int_t _nV, HessianType _hessianType, BooleanType allocDenseMats )
{
	int_t i;

	/* print copyright notice */
	if (options.printLevel != PL_NONE)
		printCopyrightNotice( );

	/* consistency check */
	if ( _nV <= 0 )
	{
		_nV = 1;
		THROWERROR( RET_INVALID_ARGUMENTS );
	}

	/* reset global message handler */
	getGlobalMessageHandler( )->reset( );

	freeHessian = BT_FALSE;
	H = 0;

	g = new real_t[_nV];
	for( i=0; i<_nV; ++i ) g[i] = 0.0;

	lb = new real_t[_nV];
	for( i=0; i<_nV; ++i ) lb[i] = 0.0;

	ub = new real_t[_nV];
	for( i=0; i<_nV; ++i ) ub[i] = 0.0;

	bounds.init( _nV );

	if ( allocDenseMats == BT_TRUE )
	{
		R = new real_t[_nV*_nV];
		for( i=0; i<_nV*_nV; ++i ) R[i] = 0.0;
	}
	else R = 0;
	haveCholesky = BT_FALSE;

	x = new real_t[_nV];
	for( i=0; i<_nV; ++i ) x[i] = 0.0;

	y = new real_t[_nV];
	for( i=0; i<_nV; ++i ) y[i] = 0.0;

	tau = 0.0;

	hessianType = _hessianType;
	regVal = 0.0;

	infeasible  = BT_FALSE;
	unbounded   = BT_FALSE;

	status = QPS_NOTINITIALISED;

	count = 0;

	ramp0 = options.initialRamping;
	ramp1 = options.finalRamping;
	rampOffset = 0;

	delta_xFR_TMP = new real_t[_nV];

	setPrintLevel( options.printLevel );

	flipper.init( (uint_t)_nV );
}


/*
 *	Q P r o b l e m B
 */
QProblemB::QProblemB( const QProblemB& rhs )
{
	freeHessian = BT_FALSE;
	H = 0;

	copy( rhs );
}


/*
 *	~ Q P r o b l e m B
 */
QProblemB::~QProblemB( )
{
	clear( );

	/* reset global message handler */
	getGlobalMessageHandler( )->reset( );
}


/*
 *	o p e r a t o r =
 */
QProblemB& QProblemB::operator=( const QProblemB& rhs )
{
	if ( this != &rhs )
	{
		clear( );
		copy( rhs );
	}

	return *this;
}


/*
 *	r e s e t
 */
returnValue QProblemB::reset( )
{
	int_t i;
	int_t nV = getNV( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Reset bounds. */
	bounds.init( nV );

	/* 2) Reset Cholesky decomposition. */
	if ( R!=0 )
		for( i=0; i<nV*nV; ++i )
			R[i] = 0.0;

	haveCholesky = BT_FALSE;

	/* 3) Reset steplength and status flags. */
	tau = 0.0;

	hessianType = HST_UNKNOWN;
	regVal = 0.0;

	infeasible  = BT_FALSE;
	unbounded   = BT_FALSE;

	status = QPS_NOTINITIALISED;

	ramp0 = options.initialRamping;
	ramp1 = options.finalRamping;
	rampOffset = 0;

	/* 4) Reset flipper object */
	flipper.init( (uint_t)nV );

	return SUCCESSFUL_RETURN;
}


/*
 *	i n i t
 */
returnValue QProblemB::init( 	SymmetricMatrix *_H, const real_t* const _g,
								const real_t* const _lb, const real_t* const _ub,
								int_t& nWSR, real_t* const cputime,
								const real_t* const xOpt, const real_t* const yOpt,
								const Bounds* const guessedBounds,
								const real_t* const _R
								)
{
	int_t i;
	int_t nV = getNV( );

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

	/* exclude this possibility in order to avoid inconsistencies */
	if ( ( xOpt == 0 ) && ( yOpt != 0 ) && ( guessedBounds != 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	if ( ( _R != 0 ) && ( ( xOpt != 0 ) || ( yOpt != 0 ) || ( guessedBounds != 0 ) ) )
		return THROWERROR( RET_NO_CHOLESKY_WITH_INITIAL_GUESS );

	/* 2) Setup QP data. */
	if ( setupQPdata( _H,_g,_lb,_ub ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 3) Call to main initialisation routine. */
	return solveInitialQP( xOpt,yOpt,guessedBounds,_R, nWSR,cputime );
}


/*
 *	i n i t
 */
returnValue QProblemB::init( 	const real_t* const _H, const real_t* const _g,
								const real_t* const _lb, const real_t* const _ub,
								int_t& nWSR, real_t* const cputime,
								const real_t* const xOpt, const real_t* const yOpt,
								const Bounds* const guessedBounds,
								const real_t* const _R
								)
{
	int_t i;
	int_t nV = getNV( );

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

	/* exclude this possibility in order to avoid inconsistencies */
	if ( ( xOpt == 0 ) && ( yOpt != 0 ) && ( guessedBounds != 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	if ( ( _R != 0 ) && ( ( xOpt != 0 ) || ( yOpt != 0 ) || ( guessedBounds != 0 ) ) )
		return THROWERROR( RET_NO_CHOLESKY_WITH_INITIAL_GUESS );

	/* 2) Setup QP data. */
	if ( setupQPdata( _H,_g,_lb,_ub ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 3) Call to main initialisation routine. */
	return solveInitialQP( xOpt,yOpt,guessedBounds,_R, nWSR,cputime );
}


/*
 *	i n i t
 */
returnValue QProblemB::init( 	const char* const H_file, const char* const g_file,
								const char* const lb_file, const char* const ub_file,
								int_t& nWSR, real_t* const cputime,
								const real_t* const xOpt, const real_t* const yOpt,
								const Bounds* const guessedBounds,
								const char* const R_file
								)
{
	int_t i;
	int_t nV = getNV( );

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

	/* exclude this possibility in order to avoid inconsistencies */
	if ( ( xOpt == 0 ) && ( yOpt != 0 ) && ( guessedBounds != 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	if ( ( R_file != 0 ) && ( ( xOpt != 0 ) || ( yOpt != 0 ) || ( guessedBounds != 0 ) ) )
		return THROWERROR( RET_NO_CHOLESKY_WITH_INITIAL_GUESS );

	/* 2) Setup QP data from files. */
	if ( setupQPdataFromFile( H_file,g_file,lb_file,ub_file ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_UNABLE_TO_READ_FILE );

	if ( R_file == 0 )
	{
		/* 3) Call to main initialisation routine. */
		return solveInitialQP( xOpt,yOpt,guessedBounds,0, nWSR,cputime );
	}
	else
	{
		/* Also read Cholesky factor from file and store it directly into R [thus... */
		returnValue returnvalue = readFromFile( R, nV,nV, R_file );
		if ( returnvalue != SUCCESSFUL_RETURN )
			return THROWWARNING( returnvalue );

		/* 3) Call to main initialisation routine. ...passing R here!] */
		return solveInitialQP( xOpt,yOpt,guessedBounds,R, nWSR,cputime );
	}
}



/*
 *	h o t s t a r t
 */
returnValue QProblemB::hotstart(	const real_t* const g_new,
									const real_t* const lb_new, const real_t* const ub_new,
									int_t& nWSR, real_t* const cputime,
									const Bounds* const guessedBounds
									)
{
	int_t i, nActiveFar;
	int_t nV = getNV ();
	real_t starttime = 0.0;
	real_t auxTime = 0.0;

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );


	/* Possibly update working set according to guess for working set of bounds. */
	if ( guessedBounds != 0 )
	{
		if ( cputime != 0 )
			starttime = getCPUtime( );

		if ( setupAuxiliaryQP( guessedBounds ) != SUCCESSFUL_RETURN )
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

	/* Simple check for consistency of bounds */
	if ( areBoundsConsistent( lb_new,ub_new ) != SUCCESSFUL_RETURN )
		return setInfeasibilityFlag(returnvalue,BT_TRUE);

	++count;


	int_t nWSR_max = nWSR;
	int_t nWSR_performed = 0;

	real_t cputime_remaining = INFTY;
	real_t cputime_needed = 0.0;

	real_t farbound = options.initialFarBounds;

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
		returnvalue = solveRegularisedQP(	g_new,lb_new,ub_new,
											nWSR,cputime,0,
											isFirstCall
											);
	}
	else
	{
		real_t *ub_new_far = new real_t[nV];
		real_t *lb_new_far = new real_t[nV];

		/* possibly extend initial far bounds to largest bound/constraint data */
		if (ub_new)
			for (i = 0; i < nV; i++)
				if ((ub_new[i] < INFTY) && (ub_new[i] > farbound)) farbound = ub_new[i];
		if (lb_new)
			for (i = 0; i < nV; i++)
				if ((lb_new[i] > -INFTY) && (lb_new[i] < -farbound)) farbound = -lb_new[i];

		updateFarBounds(	farbound,nV,
							lb_new,lb_new_far, ub_new,ub_new_far
							);

		for ( ;; )
		{
			nWSR = nWSR_max;
			if ( cputime != 0 )
				cputime_remaining = *cputime - cputime_needed;

			/* Automatically call standard solveQP if regularisation is not active. */
			returnvalue = solveRegularisedQP(	g_new,lb_new_far,ub_new_far,
												nWSR,&cputime_remaining,nWSR_performed,
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
					goto farewell;
				}

				updateFarBounds(	farbound,nV,
									lb_new,lb_new_far, ub_new,ub_new_far
									);
			}
			else if ( status == QPS_SOLVED )
			{
				real_t tol = farbound/options.growFarBounds * options.boundTolerance;
				nActiveFar = 0;
				for ( i=0; i<nV; ++i )
				{
					if ( ( ( lb_new == 0 ) || ( lb_new_far[i] > lb_new[i] ) ) && ( getAbs ( lb_new_far[i] - x[i] ) < tol ) )
						++nActiveFar;
					if ( ( ( ub_new == 0 ) || ( ub_new_far[i] < ub_new[i] ) ) && ( getAbs ( ub_new_far[i] - x[i] ) < tol ) )
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

				updateFarBounds(	farbound,nV,
									lb_new,lb_new_far, ub_new,ub_new_far
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
			delete[] lb_new_far; delete[] ub_new_far;
	}

	return ( returnvalue != SUCCESSFUL_RETURN ) ? THROWERROR( returnvalue ) : returnvalue;
}


/*
 *	h o t s t a r t
 */
returnValue QProblemB::hotstart(	const char* const g_file,
									const char* const lb_file, const char* const ub_file,
									int_t& nWSR, real_t* const cputime,
									const Bounds* const guessedBounds
									)
{
	int_t nV  = getNV( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* consistency check */
	if ( g_file == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );


	/* 1) Allocate memory (if bounds exist). */
	real_t* g_new  = new real_t[nV];
	real_t* lb_new = ( lb_file != 0 ) ? new real_t[nV] : 0;
	real_t* ub_new = ( ub_file != 0 ) ? new real_t[nV] : 0;


	/* 2) Load new QP vectors from file. */
	returnValue returnvalue;
	returnvalue = loadQPvectorsFromFile(	g_file,lb_file,ub_file,
											g_new,lb_new,ub_new
											);
	if ( returnvalue != SUCCESSFUL_RETURN )
	{
		if ( ub_file != 0 )
			delete[] ub_new;
		if ( lb_file != 0 )
			delete[] lb_new;
		delete[] g_new;

		return THROWERROR( RET_UNABLE_TO_READ_FILE );
	}


	/* 3) Actually perform hotstart. */
	returnvalue = hotstart(	g_new,lb_new,ub_new,
							nWSR,cputime,
							guessedBounds
							);


	/* 4) Free memory. */
	if ( ub_file != 0 )
		delete[] ub_new;
	if ( lb_file != 0 )
		delete[] lb_new;
	delete[] g_new;

	return returnvalue;
}


/*
 *	g e t W o r k i n g S e t
 */
returnValue QProblemB::getWorkingSet( real_t* workingSet )
{
	return getWorkingSetBounds( workingSet );
}


/*
 *	g e t W o r k i n g S e t B o u n d s
 */
returnValue QProblemB::getWorkingSetBounds( real_t* workingSetB )
{
	int_t i;
	int_t nV = this->getNV();

	if ( workingSetB == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	for ( i=0; i<nV; ++i )
	{
		switch ( bounds.getStatus(i) )
		{
			case ST_LOWER: workingSetB[i] = -1.0; break;
			case ST_UPPER: workingSetB[i] = +1.0; break;
			default:       workingSetB[i] =  0.0; break;
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	g e t W o r k i n g S e t C o n s t r a i n t s
 */
returnValue QProblemB::getWorkingSetConstraints( real_t* workingSetC )
{
	if ( workingSetC == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );
	else
		return SUCCESSFUL_RETURN;
}


/*
 *	g e t N Z
 */
int_t QProblemB::getNZ( ) const
{
	/* if no constraints are present: nZ=nFR */
	return getNFR( );
}


/*
 *	g e t O b j V a l
 */
real_t QProblemB::getObjVal( ) const
{
	real_t objVal;

	/* calculated optimal objective function value
	 * only if current QP has been solved */
	if ( ( getStatus( ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( getStatus( ) == QPS_HOMOTOPYQPSOLVED )  ||
		 ( getStatus( ) == QPS_SOLVED ) )
	{
		objVal = getObjVal( x );
	}
	else
	{
		objVal = INFTY;
	}

	return objVal;
}


/*
 *	g e t O b j V a l
 */
real_t QProblemB::getObjVal( const real_t* const _x ) const
{
	int_t i;
	int_t nV = getNV( );

	if ( nV == 0 )
		return 0.0;

	real_t objVal = 0.0;

	for( i=0; i<nV; ++i )
		objVal += _x[i]*g[i];

	switch ( hessianType )
	{
		case HST_ZERO:
			break;

		case HST_IDENTITY:
			for( i=0; i<nV; ++i )
				objVal += 0.5*_x[i]*_x[i];
			break;

		default:
			real_t *Hx = new real_t[nV];
			H->times(1, 1.0, _x, nV, 0.0, Hx, nV);
			for( i=0; i<nV; ++i )
				objVal += 0.5*_x[i]*Hx[i];
			delete[] Hx;
			break;
	}

	/* When using regularisation, the objective function value
	 * needs to be modified as follows:
	 * objVal = objVal - 0.5*_x*(Hmod-H)*_x - _x'*(gMod-g)
	 *        = objVal - 0.5*_x*eps*_x * - _x'*(-eps*_x)
	 *        = objVal + 0.5*_x*eps*_x */
	if ( usingRegularisation( ) == BT_TRUE )
	{
		for( i=0; i<nV; ++i )
			objVal += 0.5*_x[i]*regVal*_x[i];
	}

	return objVal;
}


/*
 *	g e t P r i m a l S o l u t i o n
 */
returnValue QProblemB::getPrimalSolution( real_t* const xOpt ) const
{
	int_t i;

	/* return optimal primal solution vector
	 * only if current QP has been solved */
	if ( ( getStatus( ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( getStatus( ) == QPS_HOMOTOPYQPSOLVED )  ||
		 ( getStatus( ) == QPS_SOLVED ) )
	{
		for( i=0; i<getNV( ); ++i )
			xOpt[i] = x[i];

		return SUCCESSFUL_RETURN;
	}
	else
	{
		return RET_QP_NOT_SOLVED;
	}
}


/*
 *	g e t D u a l S o l u t i o n
 */
returnValue QProblemB::getDualSolution( real_t* const yOpt ) const
{
	int_t i;

	for( i=0; i<getNV( ); ++i )
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
 *	s e t P r i n t L e v e l
 */
returnValue QProblemB::setPrintLevel( PrintLevel _printLevel )
{
	#ifndef __SUPPRESSANYOUTPUT__
		#ifndef __MATLAB__
			if ( ( options.printLevel == PL_HIGH ) && ( options.printLevel != _printLevel ) )
				THROWINFO( RET_PRINTLEVEL_CHANGED );
		#endif /* __MATLAB__ */
		options.printLevel = _printLevel;
	#else
	options.printLevel = PL_NONE;
	#endif /* __SUPPRESSANYOUTPUT__ */

	/* update message handler preferences */
 	switch ( options.printLevel )
 	{
 		case PL_NONE:
 			getGlobalMessageHandler( )->setErrorVisibilityStatus( VS_HIDDEN );
			getGlobalMessageHandler( )->setWarningVisibilityStatus( VS_HIDDEN );
			getGlobalMessageHandler( )->setInfoVisibilityStatus( VS_HIDDEN );
			break;

		case PL_TABULAR:
		case PL_LOW:
			getGlobalMessageHandler( )->setErrorVisibilityStatus( VS_VISIBLE );
			getGlobalMessageHandler( )->setWarningVisibilityStatus( VS_HIDDEN );
			getGlobalMessageHandler( )->setInfoVisibilityStatus( VS_HIDDEN );
			break;

		case PL_DEBUG_ITER:
		case PL_MEDIUM:
			getGlobalMessageHandler( )->setErrorVisibilityStatus( VS_VISIBLE );
			getGlobalMessageHandler( )->setWarningVisibilityStatus( VS_VISIBLE );
			getGlobalMessageHandler( )->setInfoVisibilityStatus( VS_HIDDEN );
			break;

		default: /* PL_HIGH */
			getGlobalMessageHandler( )->setErrorVisibilityStatus( VS_VISIBLE );
			getGlobalMessageHandler( )->setWarningVisibilityStatus( VS_VISIBLE );
			getGlobalMessageHandler( )->setInfoVisibilityStatus( VS_VISIBLE );
			break;
 	}

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t P r o p e r t i e s
 */
returnValue QProblemB::printProperties( )
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


	/* 2) Further properties. */
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


	/* 3) QP object properties. */
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


returnValue QProblemB::printOptions( ) const
{
	return options.print( );
}



/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/

/*
 *	c l e a r
 */
returnValue QProblemB::clear( )
{
	if ( ( freeHessian == BT_TRUE ) && ( H != 0 ) )
	{
		delete H;
		H = 0;
	}

	if ( g != 0 )
	{
		delete[] g;
		g = 0;
	}

	if ( lb != 0 )
	{
		delete[] lb;
		lb = 0;
	}

	if ( ub != 0 )
	{
		delete[] ub;
		ub = 0;
	}

	if ( R != 0 )
	{
		delete[] R;
		R = 0;
	}

	if ( x != 0 )
	{
		delete[] x;
		x = 0;
	}

	if ( y != 0 )
	{
		delete[] y;
		y = 0;
	}

	if ( delta_xFR_TMP != 0 )
	{
		delete[] delta_xFR_TMP;
		delta_xFR_TMP = 0;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	c o p y
 */
returnValue QProblemB::copy(	const QProblemB& rhs
								)
{
	uint_t _nV = (uint_t)rhs.getNV( );

	bounds = rhs.bounds;

	freeHessian = rhs.freeHessian;

	if ( freeHessian == BT_TRUE )
		H = (SymmetricMatrix *)(rhs.H->duplicateSym());
	else
		H = rhs.H;

	if ( rhs.g != 0 )
	{
		g = new real_t[_nV];
		setG( rhs.g );
	}
	else
		g = 0;

	if ( rhs.lb != 0 )
	{
		lb = new real_t[_nV];
		setLB( rhs.lb );
	}
	else
		lb = 0;

	if ( rhs.ub != 0 )
	{
		ub = new real_t[_nV];
		setUB( rhs.ub );
	}
	else
		ub = 0;

	if ( rhs.R != 0 )
	{
		R = new real_t[_nV*_nV];
		memcpy( R,rhs.R,_nV*_nV*sizeof(real_t) );
	}
	else
		R = 0;

	haveCholesky = rhs.haveCholesky;

	if ( rhs.x != 0 )
	{
		x = new real_t[_nV];
		memcpy( x,rhs.x,_nV*sizeof(real_t) );
	}
	else
		x = 0;

	if ( rhs.y != 0 )
	{
		y = new real_t[_nV];
		memcpy( y,rhs.y,_nV*sizeof(real_t) );
	}
	else
		y = 0;

	tau = rhs.tau;

	hessianType = rhs.hessianType;
	regVal = rhs.regVal;

	infeasible = rhs.infeasible;
	unbounded = rhs.unbounded;

	status = rhs.status;

	count = rhs.count;

	ramp0 = rhs.ramp0;
	ramp1 = rhs.ramp1;
	// AW: Following line seemed to be missing
	rampOffset = rhs.rampOffset;

	delta_xFR_TMP = new real_t[_nV];	/* nFR */

	options = rhs.options;
	setPrintLevel( options.printLevel );

	flipper = rhs.flipper;

	return SUCCESSFUL_RETURN;
}


/*
 *	d e t e r m i n e H e s s i a n T y p e
 */
returnValue QProblemB::determineHessianType( )
{
	int_t i;
	int_t nV = getNV( );
	real_t curDiag;

	/* if Hessian type has been set by user, do NOT change it! */
	switch ( hessianType )
	{
		case HST_ZERO:
			/* ensure regularisation as default options do not always solve LPs */
			if ( options.enableRegularisation == BT_FALSE )
			{
				options.enableRegularisation = BT_TRUE;
				options.numRegularisationSteps = 1;
			}
			return SUCCESSFUL_RETURN;

		case HST_IDENTITY:
			return SUCCESSFUL_RETURN;

		case HST_POSDEF:
        case HST_POSDEF_NULLSPACE:
        case HST_SEMIDEF:
		case HST_INDEF:
			/* if H == 0, continue to reset hessianType to HST_ZERO
			 *  to avoid segmentation faults! */
			if ( H != 0 )
				return SUCCESSFUL_RETURN;

		default:
			/* HST_UNKNOWN, continue */
			break;
	}

	/* if Hessian has not been allocated, assume it to be all zeros! */
	if ( H == 0 )
	{
		hessianType = HST_ZERO;
		THROWINFO( RET_ZERO_HESSIAN_ASSUMED );

		/* ensure regularisation as default options do not always solve LPs */
		if ( options.enableRegularisation == BT_FALSE )
		{
			options.enableRegularisation = BT_TRUE;
			options.numRegularisationSteps = 1;
		}

		return SUCCESSFUL_RETURN;
	}

	/* 1) If Hessian has outer-diagonal elements,
	 *    Hessian is assumed to be positive definite. */
	hessianType = HST_POSDEF;
	if ( H->isDiag() == BT_FALSE )
		return SUCCESSFUL_RETURN;

	/* 2) Otherwise it is diagonal and test for identity or zero matrix is performed. */
	BooleanType isIdentity = BT_TRUE;
	BooleanType isZero = BT_TRUE;

	for ( i=0; i<nV; ++i )
	{
        curDiag = H->diag(i);
        if ( curDiag >= INFTY )
            return RET_DIAGONAL_NOT_INITIALISED;

		if ( curDiag < -ZERO )
		{
			hessianType = HST_INDEF;
			if ( options.enableFlippingBounds == BT_FALSE )
				return THROWERROR( RET_HESSIAN_INDEFINITE );
			else
				return SUCCESSFUL_RETURN;
		}

		if ( getAbs( curDiag - 1.0 ) > EPS )
			isIdentity = BT_FALSE;

		if ( getAbs( curDiag ) > EPS )
			isZero = BT_FALSE;
	}

	if ( isIdentity == BT_TRUE )
		hessianType = HST_IDENTITY;

	if ( isZero == BT_TRUE )
	{
		hessianType = HST_ZERO;

		/* ensure regularisation as default options do not always solve LPs */
		if ( options.enableRegularisation == BT_FALSE )
		{
			options.enableRegularisation = BT_TRUE;
			options.numRegularisationSteps = 1;
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p S u b j e c t T o T y p e
 */
returnValue QProblemB::setupSubjectToType( )
{
	return setupSubjectToType( lb,ub );
}


/*
 *	s e t u p S u b j e c t T o T y p e
 */
returnValue QProblemB::setupSubjectToType( const real_t* const lb_new, const real_t* const ub_new )
{
	int_t i;
	int_t nV = getNV( );


	/* 1) Check if lower bounds are present. */
	bounds.setNoLower( BT_TRUE );
	if ( lb_new != 0 )
	{
		for( i=0; i<nV; ++i )
		{
			if ( lb_new[i] > -INFTY )
			{
				bounds.setNoLower( BT_FALSE );
				break;
			}
		}
	}

	/* 2) Check if upper bounds are present. */
	bounds.setNoUpper( BT_TRUE );
	if ( ub_new != 0 )
	{
		for( i=0; i<nV; ++i )
		{
			if ( ub_new[i] < INFTY )
			{
				bounds.setNoUpper( BT_FALSE );
				break;
			}
		}
	}

	/* 3) Determine implicitly fixed and unbounded variables. */
	if ( ( lb_new != 0 ) && ( ub_new != 0 ) )
	{
		for( i=0; i<nV; ++i )
		{
			if ( ( lb_new[i] < -INFTY+options.boundTolerance ) && ( ub_new[i] > INFTY-options.boundTolerance )
					&& (options.enableFarBounds == BT_FALSE))
			{
				bounds.setType( i,ST_UNBOUNDED );
			}
			else
			{
				if (options.enableEqualities
						&& lb[i] > ub[i] - options.boundTolerance
						&& lb_new[i] > ub_new[i] - options.boundTolerance)
					bounds.setType( i,ST_EQUALITY );
				else
					bounds.setType( i,ST_BOUNDED );
			}
		}
	}
	else
	{
		if ( ( lb_new == 0 ) && ( ub_new == 0 ) )
		{
			for( i=0; i<nV; ++i )
				bounds.setType( i,ST_UNBOUNDED );
		}
		else
		{
			for( i=0; i<nV; ++i )
				bounds.setType( i,ST_BOUNDED );
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	c o m p u t e C h o l e s k y
 */
returnValue QProblemB::computeCholesky( )
{
	int_t i, j;
	int_t nV  = getNV( );
	int_t nFR = getNFR( );

	/* 1) Initialises R with all zeros. */
	for( i=0; i<nV*nV; ++i )
		R[i] = 0.0;

	/* 2) Calculate Cholesky decomposition of H (projected to free variables). */
	switch ( hessianType )
	{
		case HST_ZERO:

			/* if Hessian is zero matrix and it has been regularised,
			 * its Cholesky factor is the identity matrix scaled by sqrt(eps). */
			if ( usingRegularisation( ) == BT_TRUE )
			{
				for( i=0; i<nV; ++i )
					RR(i,i) = getSqrt( regVal );
			}
			else
			{
				return THROWERROR( RET_CHOLESKY_OF_ZERO_HESSIAN );
			}
			break;


		case HST_IDENTITY:

			/* if Hessian is identity, so is its Cholesky factor. */
			for( i=0; i<nV; ++i )
				RR(i,i) = 1.0;
			break;


		default:

			if ( nFR > 0 )
			{
				int_t* FR_idx;
				bounds.getFree( )->getNumberArray( &FR_idx );

				/* get H */
				for ( j=0; j < nFR; ++j )
					H->getCol (FR_idx[j], bounds.getFree (), 1.0, &(R[j*nV]) );

				/* R'*R = H */
				la_int_t info = 0;
				la_uint_t _nFR = (la_uint_t)nFR, _nV = (la_uint_t)nV;

				POTRF( "U", &_nFR, R, &_nV, &info );

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
				for ( i=0; i<nFR-1; ++i )
					RR(i+1,i) = 0.0;
			}
			break;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p I n i t i a l C h o l e s k y
 */
returnValue QProblemB::setupInitialCholesky( )
{
	returnValue returnvalueCholesky;

	/* If regularisation shall be used, always regularise at beginning
	 * if initial working set is not empty. */
	if ( ( getNV() != getNFR()-getNFV() ) && ( options.enableRegularisation == BT_TRUE ) )
		if ( regulariseHessian( ) != SUCCESSFUL_RETURN )
			return RET_INIT_FAILED_REGULARISATION;

	returnvalueCholesky = computeCholesky( );

	/* If Hessian is not positive definite, regularise and try again. */
	if ( returnvalueCholesky == RET_HESSIAN_NOT_SPD )
	{
		if ( regulariseHessian( ) != SUCCESSFUL_RETURN )
			return RET_INIT_FAILED_REGULARISATION;

		returnvalueCholesky = computeCholesky( );
	}

	if ( returnvalueCholesky != SUCCESSFUL_RETURN )
		return RET_INIT_FAILED_CHOLESKY;

	haveCholesky = BT_TRUE;
	return SUCCESSFUL_RETURN;
}


/*
 *	o b t a i n A u x i l i a r y W o r k i n g S e t
 */
returnValue QProblemB::obtainAuxiliaryWorkingSet(	const real_t* const xOpt, const real_t* const yOpt,
													const Bounds* const guessedBounds, Bounds* auxiliaryBounds
													) const
{
	int_t i = 0;
	int_t nV = getNV( );


	/* 1) Ensure that desiredBounds is allocated (and different from guessedBounds). */
	if ( ( auxiliaryBounds == 0 ) || ( auxiliaryBounds == guessedBounds ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );


	/* 2) Setup working set for auxiliary initial QP. */
	if ( guessedBounds != 0 )
	{
		/* If an initial working set is specific, use it!
		 * Moreover, add all implictly fixed variables if specified. */
		for( i=0; i<nV; ++i )
		{
			#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
			if ( bounds.getType( i ) == ST_EQUALITY )
			{
				if ( auxiliaryBounds->setupBound( i,ST_LOWER ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
			}
			else
			#endif
			{
				if ( auxiliaryBounds->setupBound( i,guessedBounds->getStatus( i ) ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
			}
		}
	}
	else	/* No initial working set specified. */
	{
		if ( ( xOpt != 0 ) && ( yOpt == 0 ) )
		{
			/* Obtain initial working set by "clipping". */
			for( i=0; i<nV; ++i )
			{
				if ( xOpt[i] <= lb[i] + options.boundTolerance )
				{
					if ( auxiliaryBounds->setupBound( i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
					continue;
				}

				if ( xOpt[i] >= ub[i] - options.boundTolerance )
				{
					if ( auxiliaryBounds->setupBound( i,ST_UPPER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
					continue;
				}

				/* Moreover, add all implictly fixed variables if specified. */
				#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
				if ( bounds.getType( i ) == ST_EQUALITY )
				{
					if ( auxiliaryBounds->setupBound( i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
				else
				#endif
				{
					if ( auxiliaryBounds->setupBound( i,ST_INACTIVE ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
			}
		}

		if ( yOpt != 0 )
		{
			/* Obtain initial working set in accordance to sign of dual solution vector. */
			for( i=0; i<nV; ++i )
			{
				if ( yOpt[i] > EPS )
				{
					if ( auxiliaryBounds->setupBound( i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
					continue;
				}

				if ( yOpt[i] < -EPS )
				{
					if ( auxiliaryBounds->setupBound( i,ST_UPPER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
					continue;
				}

				/* Moreover, add all implictly fixed variables if specified. */
				#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
				if ( bounds.getType( i ) == ST_EQUALITY )
				{
					if ( auxiliaryBounds->setupBound( i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
				else
				#endif
				{
					if ( auxiliaryBounds->setupBound( i,ST_INACTIVE ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
			}
		}

		/* If xOpt and yOpt are null pointer and no initial working is specified,
		 * start with empty working set (or implicitly fixed bounds only)
		 * for auxiliary QP. */
		if ( ( xOpt == 0 ) && ( yOpt == 0 ) )
		{
			for( i=0; i<nV; ++i )
			{
				switch( bounds.getType( i ) )
				{
					case ST_UNBOUNDED:
						if ( auxiliaryBounds->setupBound( i,ST_INACTIVE ) != SUCCESSFUL_RETURN )
							return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
						break;

					/* Only add all implictly fixed variables if specified. */
					#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
					case ST_EQUALITY:
						if ( auxiliaryBounds->setupBound( i,ST_LOWER ) != SUCCESSFUL_RETURN )
							return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
						break;
					#endif

					default:
						if ( auxiliaryBounds->setupBound( i,options.initialStatusBounds ) != SUCCESSFUL_RETURN )
							return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
						break;
				}
			}
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	b a c k s o l v e R
 */
returnValue QProblemB::backsolveR(	const real_t* const b, BooleanType transposed,
									real_t* const a
									) const
{
	/* Call standard backsolve procedure (i.e. removingBound == BT_FALSE). */
	return backsolveR( b,transposed,BT_FALSE,a );
}


/*
 *	b a c k s o l v e R
 */
returnValue QProblemB::backsolveR(	const real_t* const b, BooleanType transposed,
									BooleanType removingBound,
									real_t* const a
									) const
{
	int_t i, j;
	int_t nV = getNV( );
	int_t nR = getNZ( );

	real_t sum;

	/* if backsolve is called while removing a bound, reduce nZ by one. */
	if ( removingBound == BT_TRUE )
		--nR;

	/* nothing to do */
	if ( nR <= 0 )
		return SUCCESSFUL_RETURN;


	/* Solve Ra = b, where R might be transposed. */
	if ( transposed == BT_FALSE )
	{
		/* solve Ra = b */
		for( i=(nR-1); i>=0; --i )
		{
			sum = b[i];
			for( j=(i+1); j<nR; ++j )
				sum -= RR(i,j) * a[j];

			if ( getAbs( RR(i,i) ) >= ZERO*getAbs( sum ) )
				a[i] = sum / RR(i,i);
			else
				return THROWERROR( RET_DIV_BY_ZERO );
		}
	}
	else
	{
		/* solve R^T*a = b */
		for( i=0; i<nR; ++i )
		{
			sum = b[i];
			for( j=0; j<i; ++j )
				sum -= RR(j,i) * a[j];

			if ( getAbs( RR(i,i) ) >= ZERO*getAbs( sum ) )
				a[i] = sum / RR(i,i);
			else
				return THROWERROR( RET_DIV_BY_ZERO );
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	d e t e r m i n e D a t a S h i f t
 */
returnValue QProblemB::determineDataShift(	const real_t* const g_new, const real_t* const lb_new, const real_t* const ub_new,
											real_t* const delta_g, real_t* const delta_lb, real_t* const delta_ub,
											BooleanType& Delta_bB_isZero
											)
{
	int_t i, ii;
	int_t nV  = getNV( );
	int_t nFX = getNFX( );

	int_t* FX_idx;
	bounds.getFixed( )->getNumberArray( &FX_idx );


	/* 1) Calculate shift directions. */
	for( i=0; i<nV; ++i )
		delta_g[i]  = g_new[i]  - g[i];

	if ( lb_new != 0 )
	{
		for( i=0; i<nV; ++i )
			delta_lb[i] = lb_new[i] - lb[i];
	}
	else
	{
		/* if no lower bounds exist, assume the new lower bounds to be -infinity */
		for( i=0; i<nV; ++i )
			delta_lb[i] = -INFTY - lb[i];
	}

	if ( ub_new != 0 )
	{
		for( i=0; i<nV; ++i )
			delta_ub[i] = ub_new[i] - ub[i];
	}
	else
	{
		/* if no upper bounds exist, assume the new upper bounds to be infinity */
		for( i=0; i<nV; ++i )
			delta_ub[i] = INFTY - ub[i];
	}

	/* 2) Determine if active bounds are to be shifted. */
	Delta_bB_isZero = BT_TRUE;

	for ( i=0; i<nFX; ++i )
	{
		ii = FX_idx[i];

		if ( ( getAbs( delta_lb[ii] ) > EPS ) || ( getAbs( delta_ub[ii] ) > EPS ) )
		{
			Delta_bB_isZero = BT_FALSE;
			break;
		}
	}

	return SUCCESSFUL_RETURN;
}



/*
 *	s e t u p Q P d a t a
 */
returnValue QProblemB::setupQPdata(	SymmetricMatrix *_H, const real_t* const _g,
									const real_t* const _lb, const real_t* const _ub
									)
{
	/* 1) Setup Hessian matrix. */
	setH( _H );

	/* 2) Setup gradient vector. */
	if ( _g == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );
	else
		setG( _g );

	/* 3) Setup lower/upper bounds vector. */
	setLB( _lb );
	setUB( _ub );

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p Q P d a t a
 */
returnValue QProblemB::setupQPdata(	const real_t* const _H, const real_t* const _g,
									const real_t* const _lb, const real_t* const _ub
									)
{
	/* 1) Setup Hessian matrix. */
	setH( _H );

	/* 2) Setup gradient vector. */
	if ( _g == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );
	else
		setG( _g );

	/* 3) Setup lower/upper bounds vector. */
	setLB( _lb );
	setUB( _ub );

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p Q P d a t a F r o m F i l e
 */
returnValue QProblemB::setupQPdataFromFile(	const char* const H_file, const char* const g_file,
											const char* const lb_file, const char* const ub_file
											)
{
	int_t i;
	int_t nV = getNV( );

	returnValue returnvalue;


	/* 1) Load Hessian matrix from file. */
	if ( H_file != 0 )
	{
		real_t* _H = new real_t[nV * nV];
		returnvalue = readFromFile( _H, nV,nV, H_file );
		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			delete[] _H;
			return THROWERROR( returnvalue );
		}
		setH( _H );
		H->doFreeMemory( );
	}
	else
	{
		real_t* _H = 0;
		setH( _H );
	}

	/* 2) Load gradient vector from file. */
	if ( g_file == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	returnvalue = readFromFile( g, nV, g_file );
	if ( returnvalue != SUCCESSFUL_RETURN )
		return THROWERROR( returnvalue );

	/* 3) Load lower bounds vector from file. */
	if ( lb_file != 0 )
	{
		returnvalue = readFromFile( lb, nV, lb_file );
		if ( returnvalue != SUCCESSFUL_RETURN )
			return THROWERROR( returnvalue );
	}
	else
	{
		/* if no lower bounds are specified, set them to -infinity */
		for( i=0; i<nV; ++i )
			lb[i] = -INFTY;
	}

	/* 4) Load upper bounds vector from file. */
	if ( ub_file != 0 )
	{
		returnvalue = readFromFile( ub, nV, ub_file );
		if ( returnvalue != SUCCESSFUL_RETURN )
			return THROWERROR( returnvalue );
	}
	else
	{
		/* if no upper bounds are specified, set them to infinity */
		for( i=0; i<nV; ++i )
			ub[i] = INFTY;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	l o a d Q P v e c t o r s F r o m F i l e
 */
returnValue QProblemB::loadQPvectorsFromFile(	const char* const g_file, const char* const lb_file, const char* const ub_file,
												real_t* const g_new, real_t* const lb_new, real_t* const ub_new
												) const
{
	int_t nV = getNV( );

	returnValue returnvalue;


	/* 1) Load gradient vector from file. */
	if ( ( g_file != 0 ) && ( g_new != 0 ) )
	{
		returnvalue = readFromFile( g_new, nV, g_file );
		if ( returnvalue != SUCCESSFUL_RETURN )
			return THROWERROR( returnvalue );
	}
	else
	{
		/* At least gradient vector needs to be specified! */
		return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	/* 2) Load lower bounds vector from file. */
	if ( lb_file != 0 )
	{
		if ( lb_new != 0 )
		{
			returnvalue = readFromFile( lb_new, nV, lb_file );
			if ( returnvalue != SUCCESSFUL_RETURN )
				return THROWERROR( returnvalue );
		}
		else
		{
			/* If filename is given, storage must be provided! */
			return THROWERROR( RET_INVALID_ARGUMENTS );
		}
	}

	/* 3) Load upper bounds vector from file. */
	if ( ub_file != 0 )
	{
		if ( ub_new != 0 )
		{
			returnvalue = readFromFile( ub_new, nV, ub_file );
			if ( returnvalue != SUCCESSFUL_RETURN )
				return THROWERROR( returnvalue );
		}
		else
		{
			/* If filename is given, storage must be provided! */
			return THROWERROR( RET_INVALID_ARGUMENTS );
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t I n f e a s i b i l i t y F l a g
 */
returnValue QProblemB::setInfeasibilityFlag(	returnValue returnvalue,
												BooleanType doThrowError
												)
{
	infeasible = BT_TRUE;

	if ( ( doThrowError == BT_TRUE ) || ( options.enableFarBounds == BT_FALSE ) )
		THROWERROR( returnvalue );

	return returnvalue;
}


/*
 *	a r e B o u n d s C o n s i s t e n t
 */
returnValue QProblemB::areBoundsConsistent(	const real_t* const lb_new, const real_t* const ub_new ) const
{
	if (lb_new && ub_new) {
		for (int_t i = 0; i < getNV(); ++i) {
			if (lb_new[i] > ub_new[i]+EPS) {
				return RET_QP_INFEASIBLE;
			}
		}
	}
	return SUCCESSFUL_RETURN;
}



/*
 *	i s C P U t i m e L i m i t E x c e e d e d
 */
BooleanType QProblemB::isCPUtimeLimitExceeded(	const real_t* const cputime,
												real_t starttime,
												int_t nWSR
												) const
{
	/* Always perform next QP iteration if no CPU time limit is given. */
	if ( cputime == 0 )
		return BT_FALSE;

	/* Always perform first QP iteration. */
	if ( nWSR <= 0 )
		return BT_FALSE;

	real_t elapsedTime = getCPUtime( ) - starttime;
	real_t timePerIteration = elapsedTime / ((real_t) nWSR);

	/* Determine if next QP iteration exceed CPU time limit
	 * considering the (current) average CPU time per iteration. */
	if ( ( elapsedTime + timePerIteration*1.25 ) <= ( *cputime ) )
		return BT_FALSE;
	else
		return BT_TRUE;
}


/*
 *	r e g u l a r i s e H e s s i a n
 */
returnValue QProblemB::regulariseHessian( )
{
	/* Do nothing if Hessian regularisation is disbaled! */
	if ( options.enableRegularisation == BT_FALSE )
		return SUCCESSFUL_RETURN;

	/* Regularisation of identity Hessian not possible. */
	if ( hessianType == HST_IDENTITY )
		return THROWERROR( RET_CANNOT_REGULARISE_IDENTITY );

	/* Determine regularisation parameter. */
	if ( usingRegularisation( ) == BT_TRUE )
		return SUCCESSFUL_RETURN; /*THROWERROR( RET_HESSIAN_ALREADY_REGULARISED );*/
	else
	{
		/* Regularisation of zero Hessian is done implicitly. */
		if ( hessianType == HST_ZERO )
		{
			regVal = getNorm( g,getNV() ) * options.epsRegularisation;
		}
		else
		{
			regVal = H->getNorm() * options.epsRegularisation;

			if ( H->addToDiag( regVal ) == RET_NO_DIAGONAL_AVAILABLE )
				return THROWERROR( RET_CANNOT_REGULARISE_SPARSE );
		}

		THROWINFO( RET_USING_REGULARISATION );
	}

	return SUCCESSFUL_RETURN;
}



/*
 *  c r e a t e D i a g S p a r s e M a t
 */
SymSparseMat* QProblemB::createDiagSparseMat( int_t n, real_t diagVal )
{
	real_t* M_val = new real_t[n];
	sparse_int_t* M_jc = new sparse_int_t[n+1];
	sparse_int_t* M_ir = new sparse_int_t[n+1];

	for( int_t ii=0; ii<n; ++ii )
	{
		M_val[ii] = diagVal;
		M_jc[ii] = (sparse_int_t)ii;
		M_ir[ii] = (sparse_int_t)ii;
	}
	M_jc[n] = (sparse_int_t)n;
	M_ir[n] = (sparse_int_t)n;

	SymSparseMat* M = new SymSparseMat( n,n, M_ir,M_jc,M_val );
	M->createDiagInfo( );
	M->doFreeMemory( );

	return M;
}



/*
 *	p e r f o r m R a t i o T e s t
 */
returnValue QProblemB::performRatioTest(	int_t nIdx,
											const int_t* const idxList,
											const SubjectTo* const subjectTo,
											const real_t* const num,
											const real_t* const den,
											real_t epsNum,
											real_t epsDen,
											real_t& t,
											int_t& BC_idx
											) const
{
	int_t i, ii;

	BC_idx = -1;

	for( i=0; i<nIdx; ++i )
	{
		ii = idxList[i];

		if ( subjectTo->getType( ii ) != ST_EQUALITY )
		{
			if ( ( subjectTo->getStatus( ii ) == ST_LOWER ) || ( subjectTo->getStatus( ii ) == ST_INACTIVE ) )
			{
				if ( isBlocking( num[i],den[i],epsNum,epsDen,t ) == BT_TRUE )
				{
					t = num[i] / den[i];
					BC_idx = ii;
				}
			}
			else
			if ( subjectTo->getStatus( ii ) == ST_UPPER )
			{
				if ( isBlocking( -num[i],-den[i],epsNum,epsDen,t ) == BT_TRUE )
				{
					t = num[i] / den[i];
					BC_idx = ii;
				}
			}
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 * g e t R e l a t i v e H o m o t o p y L e n g t h
 */
real_t QProblemB::getRelativeHomotopyLength(	const real_t* const g_new, const real_t* const lb_new, const real_t* const ub_new
												)
{
	int_t i;
	int_t nV = getNV( );
	real_t d, s, len = 0.0;

	/* gradient */
	for (i = 0; i < nV; i++)
	{
		s = getAbs(g_new[i]);
		if (s < 1.0) s = 1.0;
		d = getAbs(g_new[i] - g[i]) / s;
		if (d > len) len = d;
	}
	/*fprintf( stderr, "homLen = %e\n", len );*/

	/* lower bounds */
	if ( lb_new != 0 )
	{
		for (i = 0; i < nV; i++)
		{
			s = getAbs(lb_new[i]);
			if (s < 1.0) s = 1.0;
			d = getAbs(lb_new[i] - lb[i]) / s;
			if (d > len) len = d;
		}
	}
	/*fprintf( stderr, "homLen = %e\n", len );*/

	/* upper bounds */
	if ( ub_new != 0 )
	{
		for (i = 0; i < nV; i++)
		{
			s = getAbs(ub_new[i]);
			if (s < 1.0) s = 1.0;
			d = getAbs(ub_new[i] - ub[i]) / s;
			if (d > len) len = d;
		}
	}
	/*fprintf( stderr, "homLen = %e\n", len );*/

	return len;
}


/*
 * p e r f o r m R a m p i n g
 */
returnValue QProblemB::performRamping( )
{
	int_t nV = getNV( ), bstat, i;
	real_t t, rampVal;

	/* ramp inactive bounds and active dual variables */
	for (i = 0; i < nV; i++)
	{
		switch (bounds.getType(i))
		{
			case ST_EQUALITY: lb[i] = x[i]; ub[i] = x[i]; continue; /* reestablish exact feasibility */
			case ST_UNBOUNDED: continue;
			case ST_DISABLED: continue;
			default: break;
		}

		t = static_cast<real_t>((i + rampOffset) % nV) / static_cast<real_t>(nV-1);
		rampVal = (1.0-t) * ramp0 + t * ramp1;
		bstat = bounds.getStatus(i);
		if (bstat != ST_LOWER) { lb[i] = x[i] - rampVal; }
		if (bstat != ST_UPPER) { ub[i] = x[i] + rampVal; }
		if (bstat == ST_LOWER) { lb[i] = x[i]; y[i] = +rampVal; }
		if (bstat == ST_UPPER) { ub[i] = x[i]; y[i] = -rampVal; }
		if (bstat == ST_INACTIVE) y[i] = 0.0; /* reestablish exact complementarity */
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
returnValue QProblemB::updateFarBounds(	real_t curFarBound, int_t nRamp,
											const real_t* const lb_new, real_t* const lb_new_far,
											const real_t* const ub_new, real_t* const ub_new_far
											) const
{
	int_t i;
	real_t rampVal, t;
	int_t nV = getNV( );

	if ( options.enableRamping == BT_TRUE )
	{
		for ( i=0; i<nV; ++i )
		{
			t = static_cast<real_t>((i + rampOffset) % nRamp) / static_cast<real_t>(nRamp-1);
			rampVal = curFarBound * (1.0 + (1.0-t)*ramp0 + t*ramp1);

			if ( lb_new == 0 )
				lb_new_far[i] = -rampVal;
			else
				lb_new_far[i] = getMax( -rampVal,lb_new[i] );

			if ( ub_new == 0 )
				ub_new_far[i] = rampVal;
			else
				ub_new_far[i] = getMin( rampVal,ub_new[i] );
		}
	}
	else
	{
		for ( i=0; i<nV; ++i )
		{
			if ( lb_new == 0 )
				lb_new_far[i] = -curFarBound;
			else
				lb_new_far[i] = getMax( -curFarBound,lb_new[i] );

			if ( ub_new == 0 )
				ub_new_far[i] = curFarBound;
			else
				ub_new_far[i] = getMin( curFarBound,ub_new[i] );
		}
	}

	return SUCCESSFUL_RETURN;
}





/*****************************************************************************
 *  P R I V A T E                                                            *
 *****************************************************************************/

/*
 *	s o l v e I n i t i a l Q P
 */
returnValue QProblemB::solveInitialQP(	const real_t* const xOpt, const real_t* const yOpt,
										const Bounds* const guessedBounds,
										const real_t* const _R,
										int_t& nWSR, real_t* const cputime
										)
{
	int_t i,j;
	int_t nV = getNV( );


	/* start runtime measurement */
	real_t starttime = 0.0;
	if ( cputime != 0 )
		starttime = getCPUtime( );


	status = QPS_NOTINITIALISED;

	/* I) ANALYSE QP DATA: */
	/* 1) Check if Hessian happens to be the identity matrix. */
	if ( determineHessianType( ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* 2) Setup type of bounds (i.e. unbounded, implicitly fixed etc.). */
	if ( setupSubjectToType( ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	status = QPS_PREPARINGAUXILIARYQP;


	/* II) SETUP AUXILIARY QP WITH GIVEN OPTIMAL SOLUTION: */
	/* 1) Setup bounds data structure. */
	if ( bounds.setupAllFree( ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* 2) Setup optimal primal/dual solution for auxiliary QP. */
	if ( setupAuxiliaryQPsolution( xOpt,yOpt ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* 3) Obtain linear independent working set for auxiliary QP. */
	Bounds auxiliaryBounds( nV );
	if ( obtainAuxiliaryWorkingSet( xOpt,yOpt,guessedBounds, &auxiliaryBounds ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* 4) Setup working set of auxiliary QP and possibly cholesky decomposition. */
	/* a) Working set of auxiliary QP. */
	if ( setupAuxiliaryWorkingSet( &auxiliaryBounds,BT_TRUE ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* b) Regularise Hessian if necessary. */
	if ( ( hessianType == HST_ZERO ) || ( hessianType == HST_SEMIDEF ) )
	{
		if ( regulariseHessian( ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_INIT_FAILED_REGULARISATION );
	}

	/* c) Copy external Cholesky factor if provided */
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
			else if ( ( xOpt == 0 ) && ( yOpt == 0 ) && ( guessedBounds == 0 ) )
			{
				for( i=0; i<nV; ++i )
					for( j=i; j<nV; ++j )
						RR(i,j) = _R[i*nV+j];
				haveCholesky = BT_TRUE;
			}
		}
	}

	/* 5) Store original QP formulation... */
	real_t* g_original = new real_t[nV];
	real_t* lb_original = new real_t[nV];
	real_t* ub_original = new real_t[nV];

	for( i=0; i<nV; ++i )
	{
		g_original[i]  = g[i];
		lb_original[i] = lb[i];
		ub_original[i] = ub[i];
	}

	/* ... and setup QP data of an auxiliary QP having an optimal solution
	 * as specified by the user (or xOpt = yOpt = 0, by default). */
	if ( setupAuxiliaryQPgradient( ) != SUCCESSFUL_RETURN )
	{
		delete[] ub_original; delete[] lb_original; delete[] g_original;
		return THROWERROR( RET_INIT_FAILED );
	}

	if ( setupAuxiliaryQPbounds( BT_TRUE ) != SUCCESSFUL_RETURN )
	{
 		delete[] ub_original; delete[] lb_original; delete[] g_original;
		return THROWERROR( RET_INIT_FAILED );
	}

	status = QPS_AUXILIARYQPSOLVED;


	/* III) SOLVE ACTUAL INITIAL QP: */

	/* Allow only remaining CPU time for usual hotstart. */
	if ( cputime != 0 )
		*cputime -= getCPUtime( ) - starttime;

	/* Use hotstart method to find the solution of the original initial QP,... */
	returnValue returnvalue = hotstart( g_original,lb_original,ub_original, nWSR,cputime );

	/* ... deallocate memory,... */
	delete[] ub_original; delete[] lb_original; delete[] g_original;

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
returnValue QProblemB::solveQP(	const real_t* const g_new,
								const real_t* const lb_new, const real_t* const ub_new,
								int_t& nWSR, real_t* const cputime, int_t nWSRperformed,
								BooleanType isFirstCall
								)
{
	int_t iter;
	int_t nV  = getNV( );

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


	/* I) PREPARATIONS */
	/* 1) Allocate delta vectors of gradient and bounds,
	 *    index arrays and step direction arrays. */
	real_t* delta_xFR = new real_t[nV];
	real_t* delta_xFX = new real_t[nV];
	real_t* delta_yFX = new real_t[nV];

	real_t* delta_g  = new real_t[nV];
	real_t* delta_lb = new real_t[nV];
	real_t* delta_ub = new real_t[nV];

	returnValue returnvalue;
	BooleanType Delta_bB_isZero;

	int_t BC_idx;
	SubjectToStatus BC_status;

	real_t homotopyLength;

	#ifndef __SUPPRESSANYOUTPUT__
	char messageString[MAX_STRING_LENGTH];
	#endif

	/* 2) Update type of bounds, e.g. a formerly implicitly fixed
	 *    variable might have become a normal one etc. */
	if ( setupSubjectToType( lb_new,ub_new ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_HOTSTART_FAILED );

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
			/* Assign number of working set recalculations and stop runtime measurement. */
			nWSR = iter;
			if ( cputime != 0 )
				*cputime = getCPUtime( ) - starttime;

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

		/* 2) Initialise shift direction of the gradient and the bounds. */
		returnvalue = determineDataShift(	g_new,lb_new,ub_new,
											delta_g,delta_lb,delta_ub,
											Delta_bB_isZero
											);
		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			delete[] delta_yFX; delete[] delta_xFX; delete[] delta_xFR;
			delete[] delta_ub; delete[] delta_lb; delete[] delta_g;

			/* Assign number of working set recalculations and stop runtime measurement. */
			nWSR = iter;
			if ( cputime != 0 )
				*cputime = getCPUtime( ) - starttime;

			THROWERROR( RET_SHIFT_DETERMINATION_FAILED );
			return returnvalue;
		}

		/* 3) Determination of step direction of X and Y. */
		returnvalue = determineStepDirection(	delta_g,delta_lb,delta_ub,
												Delta_bB_isZero,
												delta_xFX,delta_xFR,delta_yFX
												);
		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			delete[] delta_yFX; delete[] delta_xFX; delete[] delta_xFR;
			delete[] delta_ub; delete[] delta_lb; delete[] delta_g;

			/* Assign number of working set recalculations and stop runtime measurement. */
			nWSR = iter;
			if ( cputime != 0 )
				*cputime = getCPUtime( ) - starttime;

			THROWERROR( RET_STEPDIRECTION_DETERMINATION_FAILED );
			return returnvalue;
		}


		/* 4) Determination of step length TAU.
		 *    This step along the homotopy path is also taken (without changing working set). */
		returnvalue = performStep(	delta_g,delta_lb,delta_ub,
									delta_xFX,delta_xFR,delta_yFX,
									BC_idx,BC_status
									);
		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			delete[] delta_yFX; delete[] delta_xFX; delete[] delta_xFR;
			delete[] delta_ub; delete[] delta_lb; delete[] delta_g;

			/* Assign number of working set recalculations and stop runtime measurement. */
			nWSR = iter;
			if ( cputime != 0 )
				*cputime = getCPUtime( ) - starttime;

			THROWERROR( RET_STEPLENGTH_DETERMINATION_FAILED );
			return returnvalue;
		}

		/* 5) Termination criterion. */
		homotopyLength = getRelativeHomotopyLength(g_new, lb_new, ub_new);
		if ( homotopyLength <= options.terminationTolerance )
		{
			status = QPS_SOLVED;

			THROWINFO( RET_OPTIMAL_SOLUTION_FOUND );

			if ( printIteration( iter,BC_idx,BC_status,homotopyLength,isFirstCall ) != SUCCESSFUL_RETURN )
				THROWERROR( RET_PRINT_ITERATION_FAILED ); /* do not pass this as return value! */

			nWSR = iter;
			if ( cputime != 0 )
				*cputime = getCPUtime( ) - starttime;

			delete[] delta_yFX; delete[] delta_xFX; delete[] delta_xFR;
			delete[] delta_ub; delete[] delta_lb; delete[] delta_g;

			return SUCCESSFUL_RETURN;
		}


		/* 6) Change active set. */
		returnvalue = changeActiveSet( BC_idx,BC_status );
		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			delete[] delta_yFX; delete[] delta_xFX; delete[] delta_xFR;
			delete[] delta_ub; delete[] delta_lb; delete[] delta_g;

			/* Assign number of working set recalculations and stop runtime measurement. */
			nWSR = iter;
			if ( cputime != 0 )
				*cputime = getCPUtime( ) - starttime;

			/* checks for infeasibility... */
			if ( infeasible == BT_TRUE )
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
			returnvalue = computeCholesky( );
			if (returnvalue != SUCCESSFUL_RETURN)
			{
				delete[] delta_yFX; delete[] delta_xFX; delete[] delta_xFR;
				delete[] delta_ub; delete[] delta_lb; delete[] delta_g;
				return returnvalue;
			}
		}


		/* 7) Perform Ramping Strategy on zero homotopy step or drift correction (if desired). */
		 if ( ( tau <= EPS ) && ( options.enableRamping == BT_TRUE ) )
			performRamping( );
		else
		if ( (options.enableDriftCorrection > 0) && ((iter+1) % options.enableDriftCorrection == 0) )
			performDriftCorrection( );  /* always returns SUCCESSFUL_RETURN */

		/* 8) Output information of successful QP iteration. */
		status = QPS_HOMOTOPYQPSOLVED;

		if ( printIteration( iter,BC_idx,BC_status,homotopyLength,isFirstCall ) != SUCCESSFUL_RETURN )
			THROWERROR( RET_PRINT_ITERATION_FAILED ); /* do not pass this as return value! */
	}

	delete[] delta_yFX; delete[] delta_xFX; delete[] delta_xFR;
	delete[] delta_ub; delete[] delta_lb; delete[] delta_g;

	/* stop runtime measurement */
	if ( cputime != 0 )
		*cputime = getCPUtime( ) - starttime;


	/* if programm gets to here, output information that QP could not be solved
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
returnValue QProblemB::solveRegularisedQP(	const real_t* const g_new,
											const real_t* const lb_new, const real_t* const ub_new,
											int_t& nWSR, real_t* const cputime, int_t nWSRperformed,
											BooleanType isFirstCall
											)
{
	int_t i, step;
	int_t nV = getNV( );


	/* Perform normal QP solution if QP has not been regularised. */
	if ( usingRegularisation( ) == BT_FALSE )
		return solveQP( g_new,lb_new,ub_new, nWSR,cputime,nWSRperformed,isFirstCall );


	/* I) SOLVE USUAL REGULARISED QP */
	returnValue returnvalue;

	int_t nWSR_max   = nWSR;
	int_t nWSR_total = nWSRperformed;

	real_t cputime_total = 0.0;
	real_t cputime_cur   = 0.0;

	if ( cputime == 0 )
	{
		returnvalue = solveQP( g_new,lb_new,ub_new, nWSR,0,nWSRperformed,isFirstCall );
	}
	else
	{
		cputime_cur = *cputime;
		returnvalue = solveQP( g_new,lb_new,ub_new, nWSR,&cputime_cur,nWSRperformed,isFirstCall );
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
		if ( cputime == 0 )
		{
			nWSR = nWSR_max;
			returnvalue = solveQP( gMod,lb_new,ub_new, nWSR,0,nWSR_total,isFirstCall );
		}
		else
		{
			nWSR = nWSR_max;
			cputime_cur = *cputime - cputime_total;
			returnvalue = solveQP( gMod,lb_new,ub_new, nWSR,&cputime_cur,nWSR_total,isFirstCall );
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
 *	s e t u p A u x i l i a r y W o r k i n g S e t
 */
returnValue QProblemB::setupAuxiliaryWorkingSet( 	const Bounds* const auxiliaryBounds,
													BooleanType setupAfresh
													)
{
	int_t i;
	int_t nV = getNV( );

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


	/* I) SETUP CHOLESKY FLAG:
	 *    Cholesky decomposition shall only be updated if working set
	 *    shall be updated (i.e. NOT setup afresh!) */
	BooleanType updateCholesky;
	if ( setupAfresh == BT_TRUE )
		updateCholesky = BT_FALSE;
	else
		updateCholesky = BT_TRUE;


	/* II) REMOVE FORMERLY ACTIVE BOUNDS (IF NECESSARY): */
	if ( setupAfresh == BT_FALSE )
	{
		/* Remove all active bounds that shall be inactive AND
		*  all active bounds that are active at the wrong bound. */
		for( i=0; i<nV; ++i )
		{
			if ( ( bounds.getStatus( i ) == ST_LOWER ) && ( auxiliaryBounds->getStatus( i ) != ST_LOWER ) )
				if ( removeBound( i,updateCholesky ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SETUP_WORKINGSET_FAILED );

			if ( ( bounds.getStatus( i ) == ST_UPPER ) && ( auxiliaryBounds->getStatus( i ) != ST_UPPER ) )
				if ( removeBound( i,updateCholesky ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
		}
	}


	/* III) ADD NEWLY ACTIVE BOUNDS: */
	/*      Add all inactive bounds that shall be active AND
	 *      all formerly active bounds that have been active at the wrong bound. */
	for( i=0; i<nV; ++i )
	{
		if ( ( bounds.getStatus( i ) == ST_INACTIVE ) && ( auxiliaryBounds->getStatus( i ) != ST_INACTIVE ) )
		{
			if ( addBound( i,auxiliaryBounds->getStatus( i ),updateCholesky ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y Q P s o l u t i o n
 */
returnValue QProblemB::setupAuxiliaryQPsolution(	const real_t* const xOpt, const real_t* const yOpt
													)
{
	int_t i;
	int_t nV = getNV( );


	/* Setup primal/dual solution vectors for auxiliary initial QP:
	 * if a null pointer is passed, a zero vector is assigned;
	 * old solution vector is kept if pointer to internal solution vector is passed. */
	if ( xOpt != 0 )
	{
		if ( xOpt != x )
			for( i=0; i<nV; ++i )
				x[i] = xOpt[i];
	}
	else
	{
		for( i=0; i<nV; ++i )
			x[i] = 0.0;
	}

	if ( yOpt != 0 )
	{
		if ( yOpt != y )
			for( i=0; i<nV; ++i )
				y[i] = yOpt[i];
	}
	else
	{
		for( i=0; i<nV; ++i )
			y[i] = 0.0;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y Q P g r a d i e n t
 */
returnValue QProblemB::setupAuxiliaryQPgradient( )
{
	int_t i;
	int_t nV = getNV( );

	/* Setup gradient vector: g = -H*x + y'*Id. */
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

			/* -H*x */
			H->times(1, -1.0, x, nV, 1.0, g, nV);

			break;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y Q P b o u n d s
 */
returnValue QProblemB::setupAuxiliaryQPbounds( BooleanType useRelaxation )
{
	int_t i;
	int_t nV = getNV( );


	/* Setup bound vectors. */
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
						lb[i] = x[i] - options.boundRelaxation;
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

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y Q P
 */
returnValue QProblemB::setupAuxiliaryQP( const Bounds* const guessedBounds )
{
	int_t i;
	int_t nV = getNV( );

	/* nothing to do */
	if ( guessedBounds == &bounds )
		return SUCCESSFUL_RETURN;

	status = QPS_PREPARINGAUXILIARYQP;


	/* I) SETUP WORKING SET ... */
	if ( shallRefactorise( guessedBounds ) == BT_TRUE )
	{
		/* ... WITH REFACTORISATION: */
		/* 1) Reset bounds ... */
		bounds.init( nV );

		/*    ... and set them up afresh. */
		if ( setupSubjectToType( ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		if ( bounds.setupAllFree( ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		/* 2) Setup guessed working set afresh. */
		if ( setupAuxiliaryWorkingSet( guessedBounds,BT_TRUE ) != SUCCESSFUL_RETURN )
			THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		/* 3) Calculate Cholesky decomposition. */
		if ( computeCholesky( ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
	}
	else
	{
		/* ... WITHOUT REFACTORISATION: */
		if ( setupAuxiliaryWorkingSet( guessedBounds,BT_FALSE ) != SUCCESSFUL_RETURN )
			THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
	}


	/* II) SETUP AUXILIARY QP DATA: */
	/* 1) Ensure that dual variable is zero for free bounds. */
	for ( i=0; i<nV; ++i )
		if ( bounds.getStatus( i ) == ST_INACTIVE )
			y[i] = 0.0;

	/* 2) Setup gradient and bound vectors. */
	if ( setupAuxiliaryQPgradient( ) != SUCCESSFUL_RETURN )
		THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

	if ( setupAuxiliaryQPbounds( BT_FALSE ) != SUCCESSFUL_RETURN )
		THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

	return SUCCESSFUL_RETURN;
}


/*
 *	d e t e r m i n e S t e p D i r e c t i o n
 */
returnValue QProblemB::determineStepDirection(	const real_t* const delta_g, const real_t* const delta_lb, const real_t* const delta_ub,
												BooleanType Delta_bB_isZero,
												real_t* const delta_xFX, real_t* const delta_xFR,
												real_t* const delta_yFX
												)
{
	int_t i, ii;
	int_t r;
	int_t nFR = getNFR( );
	int_t nFX = getNFX( );

	int_t* FR_idx;
	int_t* FX_idx;

	bounds.getFree( )->getNumberArray( &FR_idx );
	bounds.getFixed( )->getNumberArray( &FX_idx );


	/* This routine computes
	 * delta_xFX := delta_b
	 * delta_xFR := R \ R' \ -( delta_g + HMX*delta_xFX )
	 * delta_yFX := HMX'*delta_xFR + HFX*delta_xFX  { + eps*delta_xFX }
	 */

	/* I) DETERMINE delta_xFX := delta_{l|u}b */
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


	/* delta_xFR_TMP holds the residual, initialized with right hand side
	 * delta_xFR holds the step that gets refined incrementally */
	for ( i=0; i<nFR; ++i )
	{
		ii = FR_idx[i];
		delta_xFR_TMP[i] = - delta_g[ii];
		delta_xFR[i] = 0.0;
	}


	/* Iterative refinement loop for delta_xFR */
	for ( r=0; r<=options.numRefinementSteps; ++r )
	{
		/* II) DETERMINE delta_xFR */
		if ( nFR > 0 )
		{
			/* Add - HMX*delta_xFX
			 * This is skipped if delta_b=0 or mixed part HM=0 (H=0 or H=Id) */
			if ( ( hessianType != HST_ZERO ) && ( hessianType != HST_IDENTITY ) && ( Delta_bB_isZero == BT_FALSE ) && ( r == 0 ) )
				H->times(bounds.getFree(), bounds.getFixed(), 1, -1.0, delta_xFX, nFX, 1.0, delta_xFR_TMP, nFR);

			/* Determine R' \ ( - HMX*delta_xFX - delta_gFR ) where R'R = HFR */
			if ( backsolveR( delta_xFR_TMP,BT_TRUE,delta_xFR_TMP ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_STEPDIRECTION_FAILED_CHOLESKY );

			/* Determine HFR \ ( - HMX*delta_xFX - delta_gFR ) */
			if ( backsolveR( delta_xFR_TMP,BT_FALSE,delta_xFR_TMP ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_STEPDIRECTION_FAILED_CHOLESKY );
		}

		/* refine solution found for delta_xFR so far */
		for ( i=0; i<nFR; ++i )
			delta_xFR[i] += delta_xFR_TMP[i];

		if ( options.numRefinementSteps > 0 )
		{
			real_t rnrm = 0.0;
			/* compute new residual in delta_xFR_TMP:
			 * residual := - HFR*delta_xFR - HMX*delta_xFX - delta_gFR
			 * set to -delta_gFR */
			for ( i=0; i<nFR; ++i )
			{
				ii = FR_idx[i];
				delta_xFR_TMP[i] = -delta_g[ii];
			}
			/* add - HFR*delta_xFR */
			switch ( hessianType )
			{
				case HST_ZERO:
					break;

				case HST_IDENTITY:
					for ( i=0; i<nFR; ++i )
					{
						delta_xFR_TMP[i] -= delta_xFR[i];

						/* compute max norm */
						if (rnrm < getAbs (delta_xFR_TMP[i]))
							rnrm = getAbs (delta_xFR_TMP[i]);
					}
					break;

				default:
					H->times(bounds.getFree(), bounds.getFree(),  1, -1.0, delta_xFR, nFR, 1.0, delta_xFR_TMP, nFR);
					H->times(bounds.getFree(), bounds.getFixed(), 1, -1.0, delta_xFX, nFX, 1.0, delta_xFR_TMP, nFR);

					/* compute max norm */
					for ( i=0; i<nFR; ++i )
						if (rnrm < getAbs (delta_xFR_TMP[i]))
							rnrm = getAbs (delta_xFR_TMP[i]);

					break;
			}

			/* early termination of residual norm small enough */
			if ( rnrm < options.epsIterRef )
				break;
		}

	} /* end of refinement loop for delta_xFR */

	/* III) DETERMINE delta_yFX */
	if ( nFX > 0 )
	{
		if ( ( hessianType == HST_ZERO ) || ( hessianType == HST_IDENTITY ) )
		{
			for( i=0; i<nFX; ++i )
			{
				/* set to delta_g */
				ii = FX_idx[i];
				delta_yFX[i] = delta_g[ii];

				/* add HFX*delta_xFX = {0|I}*delta_xFX */
				if ( hessianType == HST_ZERO )
				{
					if ( usingRegularisation( ) == BT_TRUE )
						delta_yFX[i] += regVal*delta_xFX[i];
				}
				else
					delta_yFX[i] += delta_xFX[i];
			}
		}
		else
		{
			for( i=0; i<nFX; ++i )
			{
				/* set to delta_g */
				ii = FX_idx[i];
				delta_yFX[i] = delta_g[ii];
			}
			H->times(bounds.getFixed(), bounds.getFree(), 1, 1.0, delta_xFR, nFR, 1.0, delta_yFX, nFX);
			if (Delta_bB_isZero == BT_FALSE)
				H->times(bounds.getFixed(), bounds.getFixed(), 1, 1.0, delta_xFX, nFX, 1.0, delta_yFX, nFX);
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	p e r f o r m S t e p
 */
returnValue QProblemB::performStep(	const real_t* const delta_g,
									const real_t* const delta_lb, const real_t* const delta_ub,
									const real_t* const delta_xFX,
									const real_t* const delta_xFR,
									const real_t* const delta_yFX,
									int_t& BC_idx, SubjectToStatus& BC_status
									)
{
	int_t i, ii;
	int_t nV = getNV( );
	int_t nFR = getNFR( );
	int_t nFX = getNFX( );

	int_t* FR_idx;
	int_t* FX_idx;

	bounds.getFree( )->getNumberArray( &FR_idx );
	bounds.getFixed( )->getNumberArray( &FX_idx );

	tau = 1.0;
	BC_idx = -1;
	BC_status = ST_UNDEFINED;

	int_t BC_idx_tmp = -1;

	real_t* num = new real_t[nV];
	real_t* den = new real_t[nV];


	/* I) DETERMINE MAXIMUM DUAL STEPLENGTH, i.e. ensure that
	 *    active dual bounds remain valid (ignoring implicitly fixed variables): */
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
	}


	/* II) DETERMINE MAXIMUM PRIMAL STEPLENGTH, i.e. ensure that
	 *     inactive bounds remain valid (ignoring unbounded variables). */
	/* 1) Inactive lower bounds. */
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
		}
	}

	/* 2) Inactive upper bounds. */
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
		}
	}

	delete[] den;
	delete[] num;


	#ifndef __SUPPRESSANYOUTPUT__
	char messageString[MAX_STRING_LENGTH];

	if ( BC_status == ST_UNDEFINED )
		snprintf( messageString,MAX_STRING_LENGTH,"Stepsize is %.15e!",tau );
	else
		snprintf( messageString,MAX_STRING_LENGTH,"Stepsize is %.15e! (idx = %d, status = %d)",tau,(int)BC_idx,(int)BC_status );

	getGlobalMessageHandler( )->throwInfo( RET_STEPSIZE_NONPOSITIVE,messageString,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
	#endif


	/* III) PERFORM STEP ALONG HOMOTOPY PATH */
	if ( tau > ZERO )
	{
		/* 1) Perform step in primal und dual space. */
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

		/* 2) Shift QP data. */
		for( i=0; i<nV; ++i )
		{
			g[i]  += tau*delta_g[i];
			lb[i] += tau*delta_lb[i];
			ub[i] += tau*delta_ub[i];
		}
	}
	else
	{
		/* print a warning if stepsize is zero */
		#ifndef __SUPPRESSANYOUTPUT__
		snprintf( messageString,MAX_STRING_LENGTH,"Stepsize is %.15e",tau );
		getGlobalMessageHandler( )->throwWarning( RET_STEPSIZE,messageString,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
		#endif
	}


	return SUCCESSFUL_RETURN;
}


/*
 *	c h a n g e A c t i v e S e t
 */
returnValue QProblemB::changeActiveSet( int_t BC_idx, SubjectToStatus BC_status )
{
	#ifndef __SUPPRESSANYOUTPUT__
	char messageString[MAX_STRING_LENGTH];
	#endif

	/* IV) UPDATE ACTIVE SET */
	switch ( BC_status )
	{
		/* Optimal solution found as no working set change detected. */
		case ST_UNDEFINED:
			return RET_OPTIMAL_SOLUTION_FOUND;


		/* Remove one variable from active set. */
		case ST_INACTIVE:
			#ifndef __SUPPRESSANYOUTPUT__
			snprintf( messageString,MAX_STRING_LENGTH,"bound no. %d.",(int)BC_idx );
			getGlobalMessageHandler( )->throwInfo( RET_REMOVE_FROM_ACTIVESET,messageString,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
			#endif

			if ( removeBound( BC_idx,BT_TRUE ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_REMOVE_FROM_ACTIVESET_FAILED );

			y[BC_idx] = 0.0;
			break;


		/* Add one variable to active set. */
		default:
			#ifndef __SUPPRESSANYOUTPUT__
			if ( BC_status == ST_LOWER )
				snprintf( messageString,MAX_STRING_LENGTH,"lower bound no. %d.",(int)BC_idx );
			else
				snprintf( messageString,MAX_STRING_LENGTH,"upper bound no. %d.",(int)BC_idx );
				getGlobalMessageHandler( )->throwInfo( RET_ADD_TO_ACTIVESET,messageString,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
			#endif

			if ( addBound( BC_idx,BC_status,BT_TRUE ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_ADD_TO_ACTIVESET_FAILED );
			break;
	}

	return SUCCESSFUL_RETURN;
}



/*
 * p e r f o r m D r i f t C o r r e c t i o n
 */
returnValue QProblemB::performDriftCorrection( )
{
	int_t i;
	int_t nV = getNV ();

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

	return setupAuxiliaryQPgradient( );
}


/*
 *	s h a l l R e f a c t o r i s e
 */
BooleanType QProblemB::shallRefactorise( const Bounds* const guessedBounds ) const
{
	int_t i;
	int_t nV = getNV( );

	/* always refactorise if Hessian is not known to be positive definite */
	if ( ( hessianType == HST_SEMIDEF ) || ( hessianType == HST_INDEF ) )
		return BT_TRUE;

	/* 1) Determine number of bounds that have same status
	 *    in guessed AND current bounds.*/
	int_t differenceNumber = 0;

	for( i=0; i<nV; ++i )
		if ( guessedBounds->getStatus( i ) != bounds.getStatus( i ) )
			++differenceNumber;

	/* 2) Decide wheter to refactorise or not. */
	if ( 2*differenceNumber > guessedBounds->getNFX( ) )
		return BT_TRUE;
	else
		return BT_FALSE;
}


/*
 *	a d d B o u n d
 */
returnValue QProblemB::addBound(	int_t number, SubjectToStatus B_status,
									BooleanType updateCholesky
									)
{
	int_t i, j;
	int_t nV  = getNV( );
	int_t nFR = getNFR( );


	/* consistency check */
	if ( ( getStatus( ) == QPS_NOTINITIALISED )    ||
		 ( getStatus( ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( getStatus( ) == QPS_HOMOTOPYQPSOLVED )  ||
		 ( getStatus( ) == QPS_SOLVED )            )
	{
		return THROWERROR( RET_UNKNOWN_BUG );
	}

	/* Perform cholesky updates only if QProblemB has been initialised! */
	if ( getStatus( ) == QPS_PREPARINGAUXILIARYQP )
	{
		/* UPDATE INDICES */
		if ( bounds.moveFreeToFixed( number,B_status ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_ADDBOUND_FAILED );

		return SUCCESSFUL_RETURN;
	}


	/* I) PERFORM CHOLESKY UPDATE: */
	if ( ( updateCholesky == BT_TRUE ) &&
		 ( hessianType != HST_ZERO )   && ( hessianType != HST_IDENTITY ) )
	{
		/* 1) Index of variable to be added within the list of free variables. */
		int_t number_idx = bounds.getFree( )->getIndex( number );

		real_t c, s, nu;

		/* 2) Use row-wise Givens rotations to restore upper triangular form of R. */
		for( i=number_idx+1; i<nFR; ++i )
		{
			computeGivens( RR(i-1,i),RR(i,i), RR(i-1,i),RR(i,i),c,s );
			nu = s/(1.0+c);

			for( j=(1+i); j<nFR; ++j ) /* last column of R is thrown away */
				applyGivens( c,s,nu,RR(i-1,j),RR(i,j), RR(i-1,j),RR(i,j) );
		}

		/* 3) Delete <number_idx>th column and ... */
		for( i=0; i<nFR-1; ++i )
			for( j=number_idx+1; j<nFR; ++j )
				RR(i,j-1) = RR(i,j);
		/* ... last column of R. */
		for( i=0; i<nFR; ++i )
			RR(i,nFR-1) = 0.0;
	}

	/* II) UPDATE INDICES */
	tabularOutput.idxAddB = number;
	if ( bounds.moveFreeToFixed( number,B_status ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_ADDBOUND_FAILED );


	return SUCCESSFUL_RETURN;
}


/*
 *	r e m o v e B o u n d
 */
returnValue QProblemB::removeBound(	int_t number,
									BooleanType updateCholesky
									)
{
	int_t i;
	int_t nV  = getNV( );
	int_t nFR = getNFR( );


	/* consistency check */
	if ( ( getStatus( ) == QPS_NOTINITIALISED )    ||
		 ( getStatus( ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( getStatus( ) == QPS_HOMOTOPYQPSOLVED )  ||
		 ( getStatus( ) == QPS_SOLVED )            )
	{
		return THROWERROR( RET_UNKNOWN_BUG );
	}

	/* save index sets and decompositions for flipping bounds strategy */
	if ( options.enableFlippingBounds == BT_TRUE )
		flipper.set( &bounds,R );

	/* I) UPDATE INDICES */
	tabularOutput.idxRemB = number;
	if ( bounds.moveFixedToFree( number ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_REMOVEBOUND_FAILED );

	/* Perform cholesky updates only if QProblemB has been initialised! */
	if ( getStatus( ) == QPS_PREPARINGAUXILIARYQP )
		return SUCCESSFUL_RETURN;


	/* II) PERFORM CHOLESKY UPDATE */
	if ( ( updateCholesky == BT_TRUE ) &&
		 ( hessianType != HST_ZERO )   && ( hessianType != HST_IDENTITY ) )
	{
		int_t* FR_idx;
		bounds.getFree( )->getNumberArray( &FR_idx );

		/* 1) Calculate new column of cholesky decomposition. */
		real_t* rhs = new real_t[nFR+1];
		real_t* r   = new real_t[nFR];

		real_t r0;
		switch ( hessianType )
		{
			case HST_ZERO: /* TODO: Code can/should? never get here!! */
				if ( usingRegularisation( ) == BT_FALSE )
					r0 = 0.0;
				else
					r0 = regVal;
				for( i=0; i<nFR; ++i )
					rhs[i] = 0.0;
				break;

			case HST_IDENTITY:
				r0 = 1.0;
				for( i=0; i<nFR; ++i )
					rhs[i] = 0.0;
				break;

			default:
				H->getRow(number, bounds.getFree(), 1.0, rhs);
				r0 = H->diag(number);
				break;
		}

		if ( backsolveR( rhs,BT_TRUE,BT_TRUE,r ) != SUCCESSFUL_RETURN )
		{
			delete[] rhs; delete[] r;
			return THROWERROR( RET_REMOVEBOUND_FAILED );
		}

		for( i=0; i<nFR; ++i )
			r0 -= r[i]*r[i];

		/* 2) Store new column into R. */
		for( i=0; i<nFR; ++i )
			RR(i,nFR) = r[i];

		if ( options.enableFlippingBounds == BT_TRUE )
		{
			if ( r0 > options.epsFlipping )
				RR(nFR,nFR) = getSqrt( r0 );
			else
			{
				hessianType = HST_SEMIDEF;

				flipper.get( &bounds,R );
				bounds.flipFixed(number);

				switch (bounds.getStatus(number))
				{
					case ST_LOWER: lb[number] = ub[number]; break;
					case ST_UPPER: ub[number] = lb[number]; break;
					default: delete[] rhs; delete[] r; return THROWERROR( RET_MOVING_BOUND_FAILED );
				}

			}
		}
		else
		{
			if ( r0 > ZERO )
				RR(nFR,nFR) = getSqrt( r0 );
			else
			{
				delete[] rhs; delete[] r;

				hessianType = HST_SEMIDEF;
				return THROWERROR( RET_HESSIAN_NOT_SPD );
			}
		}

		delete[] rhs; delete[] r;
	}

	if ( ( hessianType == HST_ZERO ) && ( options.enableFlippingBounds == BT_TRUE ) )
	{
		flipper.get( &bounds,R );
		bounds.flipFixed(number);

		switch (bounds.getStatus(number))
		{
			case ST_LOWER: lb[number] = ub[number]; break;
			case ST_UPPER: ub[number] = lb[number]; break;
			default: return THROWERROR( RET_MOVING_BOUND_FAILED );
		}

	}

	return SUCCESSFUL_RETURN;
}



/*
 *	p r i n t I t e r a t i o n
 */
returnValue QProblemB::printIteration( 	int_t iter,
										int_t BC_idx,	SubjectToStatus BC_status, real_t homotopyLength,
										BooleanType isFirstCall
										)
{
	#ifndef __SUPPRESSANYOUTPUT__

	/* consistency check */
	if ( iter < 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	int_t i;
	int_t nV = getNV();
	real_t stat, bfeas, bcmpl;
	real_t *grad = 0;

	char myPrintfString[MAX_STRING_LENGTH];
	char info[MAX_STRING_LENGTH];
	const char excStr[] = " ef";

	switch ( options.printLevel )
	{
		case PL_DEBUG_ITER:
			grad = new real_t[nV];
			stat = bfeas = bcmpl = 0.0;

			/* stationarity */
			for (i = 0; i < nV; i++) grad[i] = g[i] - y[i];
			H->times(1, 1.0, x, nV, 1.0, grad, nV);
			for (i = 0; i < nV; i++) if (getAbs(grad[i]) > stat) stat = getAbs(grad[i]);

			/* feasibility */
			for (i = 0; i < nV; i++) if (lb[i] - x[i] > bfeas) bfeas = lb[i] - x[i];
			for (i = 0; i < nV; i++) if (x[i] - ub[i] > bfeas) bfeas = x[i] - ub[i];

			/* complementarity */
			for (i = 0; i < nV; i++) if (y[i] > +EPS && getAbs((lb[i] - x[i])*y[i]) > bcmpl) bcmpl = getAbs((lb[i] - x[i])*y[i]);
			for (i = 0; i < nV; i++) if (y[i] < -EPS && getAbs((ub[i] - x[i])*y[i]) > bcmpl) bcmpl = getAbs((ub[i] - x[i])*y[i]);

			if ( (iter % 10 == 0) && ( isFirstCall == BT_TRUE ) )
			{
				snprintf( myPrintfString,MAX_STRING_LENGTH, "\n%5s %4s %4s %9s %9s %9s %9s %9s\n",
						"iter", "addB", "remB", "hom len", "tau", "stat", "bfeas", "bcmpl");
			}
			myPrintf( myPrintfString );

			snprintf( myPrintfString,MAX_STRING_LENGTH, "%5d ",(int)iter );
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

			snprintf( myPrintfString,MAX_STRING_LENGTH, "%9.2e %9.2e %9.2e %9.2e %9.2e\n",
					homotopyLength, tau, stat, bfeas, bcmpl);
			myPrintf( myPrintfString );

			delete[] grad;
			break;

		case PL_TABULAR:
			if ( (iter % 10 == 0) && ( isFirstCall == BT_TRUE ) )
			{
				snprintf( myPrintfString,MAX_STRING_LENGTH, "\n%5s %6s %6s %9s %9s\n",
						"iter", "addB", "remB", "hom len", "tau");
				myPrintf( myPrintfString );
			}

			snprintf( myPrintfString,MAX_STRING_LENGTH, "%5d ",(int)iter );
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

			snprintf( myPrintfString,MAX_STRING_LENGTH, "%9.2e %9.2e\n", homotopyLength, tau);
			myPrintf( myPrintfString );
			break;

		case PL_MEDIUM:
			/* 1) Print header at first iteration. */
 			if ( ( iter == 0 ) && ( isFirstCall == BT_TRUE ) )
			{
				snprintf( myPrintfString,MAX_STRING_LENGTH,"\n\n#################   qpOASES  --  QP NO. %3.0d   ##################\n\n",(int)count );
				myPrintf( myPrintfString );

				myPrintf( "    Iter   |    StepLength    |       Info       |   nFX    \n" );
				myPrintf( " ----------+------------------+------------------+--------- \n" );
			}

			/* 2) Print iteration line. */
			if ( BC_status == ST_UNDEFINED )
			{
				if ( hessianType == HST_ZERO )
					snprintf( info,3,"LP" );
				else
					snprintf( info,3,"QP" );

				if ( isFirstCall == BT_TRUE )
					snprintf( myPrintfString,MAX_STRING_LENGTH,"   %5.1d   |   %1.6e   |    %s SOLVED     |  %4.1d   \n", (int)iter,tau,info,(int)getNFX( ) );
				else
					snprintf( myPrintfString,MAX_STRING_LENGTH,"   %5.1d*  |   %1.6e   |    %s SOLVED     |  %4.1d   \n", (int)iter,tau,info,(int)getNFX( ) );
				myPrintf( myPrintfString );
			}
			else
			{
				if ( BC_status == ST_INACTIVE )
					snprintf( info,8,"REM BND" );
				else
					snprintf( info,8,"ADD BND" );

				snprintf( myPrintfString,MAX_STRING_LENGTH,"   %5.1d   |   %1.6e   |   %s %4.1d   |  %4.1d   \n", (int)iter,tau,info,(int)BC_idx,(int)getNFX( ) );
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



END_NAMESPACE_QPOASES


/*
 *	end of file
 */
