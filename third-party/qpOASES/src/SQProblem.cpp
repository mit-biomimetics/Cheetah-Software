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
 *	\file src/SQProblem.cpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.2
 *	\date 2007-2017
 *
 *	Implementation of the SQProblem class which is able to use the newly
 *	developed online active set strategy for parametric quadratic programming
 *	with varying matrices.
 */


#include <qpOASES/SQProblem.hpp>


BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/


/*
 *	S Q P r o b l e m
 */
SQProblem::SQProblem( ) : QProblem( )
{
}


/*
 *	S Q P r o b l e m
 */
SQProblem::SQProblem( int_t _nV, int_t _nC, HessianType _hessianType, BooleanType allocDenseMats ) 
	: QProblem( _nV,_nC,_hessianType,allocDenseMats )
{
}


/*
 *	S Q P r o b l e m
 */
SQProblem::SQProblem( const SQProblem& rhs ) : QProblem( rhs )
{
}


/*
 *	~ S Q P r o b l e m
 */
SQProblem::~SQProblem( )
{
}


/*
 *	o p e r a t o r =
 */
SQProblem& SQProblem::operator=( const SQProblem& rhs )
{
	if ( this != &rhs )
	{
		QProblem::operator=( rhs );
	}

	return *this;
}



/*
 * h o t s t a r t
 */
returnValue SQProblem::hotstart(	SymmetricMatrix *H_new, const real_t* const g_new, Matrix *A_new,
									const real_t* const lb_new, const real_t* const ub_new,
									const real_t* const lbA_new, const real_t* const ubA_new,
									int_t& nWSR, real_t* const cputime,
									const Bounds* const guessedBounds, const Constraints* const guessedConstraints
									)
{
	if ( ( getStatus( ) == QPS_NOTINITIALISED )       ||
		 ( getStatus( ) == QPS_PREPARINGAUXILIARYQP ) ||
		 ( getStatus( ) == QPS_PERFORMINGHOMOTOPY )   )
	{
		return THROWERROR( RET_HOTSTART_FAILED_AS_QP_NOT_INITIALISED );
	}


	real_t starttime = 0.0;
	real_t auxTime = 0.0;

    if ( cputime != 0 )
        starttime = getCPUtime( );


	/* I) UPDATE QP MATRICES AND VECTORS */
	if ( setupNewAuxiliaryQP( H_new,A_new,lb_new,ub_new,lbA_new,ubA_new ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );


	/* II) PERFORM USUAL HOMOTOPY */

	/* Allow only remaining CPU time for usual hotstart. */
	if ( cputime != 0 )
	{
		auxTime = getCPUtime( ) - starttime;
		*cputime -= auxTime;
	}

	returnValue returnvalue = QProblem::hotstart(	g_new,lb_new,ub_new,lbA_new,ubA_new,
													nWSR,cputime,
													guessedBounds,guessedConstraints
													);
	
	if ( cputime != 0 )
		*cputime += auxTime;

	return returnvalue;
}


/*
 *	h o t s t a r t
 */
returnValue SQProblem::hotstart(	const real_t* const H_new, const real_t* const g_new, const real_t* const A_new,
									const real_t* const lb_new, const real_t* const ub_new,
									const real_t* const lbA_new, const real_t* const ubA_new,
									int_t& nWSR, real_t* const cputime,
									const Bounds* const guessedBounds, const Constraints* const guessedConstraints
									)
{
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


	/* I) UPDATE QP MATRICES AND VECTORS */
	if ( setupNewAuxiliaryQP( H_new,A_new,lb_new,ub_new,lbA_new,ubA_new ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );


	/* II) PERFORM USUAL HOMOTOPY */

	/* Allow only remaining CPU time for usual hotstart. */
	if ( cputime != 0 )
		*cputime -= getCPUtime( ) - starttime;

	returnValue returnvalue = QProblem::hotstart(	g_new,lb_new,ub_new,lbA_new,ubA_new,
													nWSR,cputime,
													guessedBounds,guessedConstraints
													);


	/* stop runtime measurement */
	if ( cputime != 0 )
		*cputime = getCPUtime( ) - starttime;

	return returnvalue;
}


/*
 *	h o t s t a r t
 */
returnValue SQProblem::hotstart(	const char* const H_file, const char* const g_file, const char* const A_file,
									const char* const lb_file, const char* const ub_file,
									const char* const lbA_file, const char* const ubA_file,
									int_t& nWSR, real_t* const cputime,
									const Bounds* const guessedBounds, const Constraints* const guessedConstraints
									)
{
	int_t nV = getNV( );
	int_t nC = getNC( );

	returnValue returnvalue;

	/* consistency checks */
	if ( ( H_file == 0 ) || ( g_file == 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	if ( ( nC > 0 ) && ( A_file == 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );


	/* 1) Load new QP matrices from files. */
	real_t* H_new  = new real_t[nV*nV];
	real_t* A_new  = new real_t[nC*nV];

	if ( readFromFile( H_new, nV,nV, H_file ) != SUCCESSFUL_RETURN )
	{
		delete[] A_new;
		delete[] H_new;
		return THROWERROR( RET_UNABLE_TO_READ_FILE );
	}

	if ( readFromFile( A_new, nC,nV, A_file ) != SUCCESSFUL_RETURN )
	{
		delete[] A_new;
		delete[] H_new;
		return THROWERROR( RET_UNABLE_TO_READ_FILE );
	}

	/* 2) Load new QP vectors from files. */
	real_t* g_new   = new real_t[nV];
	real_t* lb_new  = ( lb_file != 0 )  ? new real_t[nV] : 0;
	real_t* ub_new  = ( ub_file != 0 )  ? new real_t[nV] : 0;
	real_t* lbA_new = ( lbA_file != 0 ) ? new real_t[nC] : 0;
	real_t* ubA_new = ( ubA_file != 0 ) ? new real_t[nC] : 0;

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
		delete[] A_new;
		delete[] H_new;

		return THROWERROR( RET_UNABLE_TO_READ_FILE );
	}

	/* 3) Actually perform hotstart. */
	returnvalue = hotstart(	H_new,g_new,A_new,lb_new,ub_new,lbA_new,ubA_new,
							nWSR,cputime,
							guessedBounds,guessedConstraints
							);

	if ( ubA_file != 0 )
		delete[] ubA_new;
	if ( lbA_file != 0 )
		delete[] lbA_new;
	if ( ub_file != 0 )
		delete[] ub_new;
	if ( lb_file != 0 )
		delete[] lb_new;
	delete[] g_new;
	delete[] A_new;
	delete[] H_new;

	return returnvalue;
}


/*
 *	h o t s t a r t
 */
returnValue SQProblem::hotstart(	const real_t* const g_new,
									const real_t* const lb_new, const real_t* const ub_new,
									const real_t* const lbA_new, const real_t* const ubA_new,
									int_t& nWSR, real_t* const cputime,
									const Bounds* const guessedBounds, const Constraints* const guessedConstraints
									)
{
	/* Call to hotstart function for fixed QP matrices. */
	return QProblem::hotstart( g_new,lb_new,ub_new,lbA_new,ubA_new, nWSR,cputime, guessedBounds,guessedConstraints );
}


/*
 *	h o t s t a r t
 */
returnValue SQProblem::hotstart(	const char* const g_file,
									const char* const lb_file, const char* const ub_file,
									const char* const lbA_file, const char* const ubA_file,
									int_t& nWSR, real_t* const cputime,
									const Bounds* const guessedBounds, const Constraints* const guessedConstraints
									)
{
	/* Call to hotstart function for fixed QP matrices. */
	return QProblem::hotstart( g_file,lb_file,ub_file,lbA_file,ubA_file, nWSR,cputime, guessedBounds,guessedConstraints );
}



#ifdef __MATLAB__
returnValue SQProblem::resetMatrixPointers( )
{
	H = 0;
	A = 0;

	return SUCCESSFUL_RETURN;
}
#endif



/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/

/*
 *	s e t u p N e w A u x i l i a r y Q P
 */
returnValue SQProblem::setupNewAuxiliaryQP(	SymmetricMatrix *H_new, Matrix *A_new,
											const real_t *lb_new, const real_t *ub_new, const real_t *lbA_new, const real_t *ubA_new
											)
{
	int_t i;
	int_t nV = getNV( );
	int_t nC = getNC( );
	returnValue returnvalue;

	if ( ( getStatus( ) == QPS_NOTINITIALISED )       ||
		 ( getStatus( ) == QPS_PREPARINGAUXILIARYQP ) ||
		 ( getStatus( ) == QPS_PERFORMINGHOMOTOPY )   )
	{
		return THROWERROR( RET_UPDATEMATRICES_FAILED_AS_QP_NOT_SOLVED );
	}

	status = QPS_PREPARINGAUXILIARYQP;


	/* I) SETUP NEW QP MATRICES AND VECTORS: */
	/* 1) Shift constraints' bounds vectors by (A_new - A)'*x_opt to ensure
	 *    that old optimal solution remains feasible for new QP data. */
	/*    Firstly, shift by -A'*x_opt and ... */
	if ( nC > 0 )
	{
		if ( A_new == 0 )
			return THROWERROR( RET_INVALID_ARGUMENTS );

		for ( i=0; i<nC; ++i )
		{
			lbA[i] = -Ax_l[i];
			ubA[i] =  Ax_u[i];
		}

		/* Set constraint matrix as well as ... */
		setA( A_new );

		/* ... secondly, shift by +A_new'*x_opt. */
		for ( i=0; i<nC; ++i )
		{
			lbA[i] += Ax[i];
			ubA[i] += Ax[i];
		}

		/* update constraint products. */
		for ( i=0; i<nC; ++i )
		{
			Ax_u[i] = ubA[i] - Ax[i];
			Ax_l[i] = Ax[i] - lbA[i];
		}
	}

	/* 2) Set new Hessian matrix, determine Hessian type and
	 *    regularise new Hessian matrix if necessary. */
	/* a) Setup new Hessian matrix and determine its type. */
	if ( H_new != 0 )
	{
		setH( H_new );

		hessianType = HST_UNKNOWN;
		if ( determineHessianType( ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		/* b) Regularise new Hessian if necessary. */
		if ( ( hessianType == HST_ZERO ) ||
			 ( hessianType == HST_SEMIDEF ) ||
			 ( usingRegularisation( ) == BT_TRUE ) )
		{
			regVal = 0.0; /* reset previous regularisation */

			if ( regulariseHessian( ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
		}
	}
	else
	{
		/* if no Hessian is specified, keep previous trivial Hessian (HST_ZERO or HST_IDENTITY),
		   otherwise abort */
		if ( H != 0 )
			return THROWERROR( RET_NO_HESSIAN_SPECIFIED );
	}

	/* 3) Setup QP gradient. */
	if ( setupAuxiliaryQPgradient( ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );


	/* II) SETUP WORKING SETS AND MATRIX FACTORISATIONS: */
	/* 1) Make a copy of current bounds/constraints ... */
	Bounds      oldBounds      = bounds;
	Constraints oldConstraints = constraints;

    /* we're trying to find an active set with positive definite null
     * space Hessian twice:
     * - first for the current active set including all equalities
     * - second after moving all inactive variables to a bound
     *   (depending on Options). This creates an empty null space and
     *   is guaranteed to succeed. Thus this loop will exit after n_try=1.
     */
    int_t n_try;
    for (n_try = 0; n_try < 2; ++n_try) {

        if (n_try > 0) {
            // the current active set leaves an indefinite null space Hessian
            // move all inactive variables to a bound, creating an empty null space
            for (int_t ii = 0; ii < nV; ++ii)
                if (oldBounds.getStatus (ii) == ST_INACTIVE)
                    oldBounds.setStatus (ii, options.initialStatusBounds);
        }

        /*    ... reset them ... */
        bounds.init( nV );
        constraints.init( nC );

        /*    ... and set them up afresh. */
        if ( setupSubjectToType(lb_new,ub_new,lbA_new,ubA_new ) != SUCCESSFUL_RETURN )
            return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

        if ( bounds.setupAllFree( ) != SUCCESSFUL_RETURN )
            return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

        if ( constraints.setupAllInactive( ) != SUCCESSFUL_RETURN )
            return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

        /* 2) Setup TQ factorisation. */
        if ( setupTQfactorisation( ) != SUCCESSFUL_RETURN )
            return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		// check for equalities that have become bounds ...
		for (int_t ii = 0; ii < nC; ++ii) {
			if (oldConstraints.getType (ii) == ST_EQUALITY && constraints.getType (ii) == ST_BOUNDED) {
				if (oldConstraints.getStatus (ii) == ST_LOWER && y[nV+ii] < 0.0)
					oldConstraints.setStatus (ii, ST_UPPER);
				else if (oldConstraints.getStatus (ii) == ST_UPPER && y[nV+ii] > 0.0)
					oldConstraints.setStatus (ii, ST_LOWER);
			}
		}

		// ... and do the same also for the bounds!
		for (int_t ii = 0; ii < nV; ++ii) {
			if (oldBounds.getType(ii) == ST_EQUALITY
					&& bounds.getType(ii) == ST_BOUNDED) {
				if (oldBounds.getStatus(ii) == ST_LOWER && y[ii] < 0.0)
					oldBounds.setStatus(ii, ST_UPPER);
				else if (oldBounds.getStatus(ii) == ST_UPPER && y[ii] > 0.0)
					oldBounds.setStatus(ii, ST_LOWER);
			}
		}

        /* 3) Setup old working sets afresh (updating TQ factorisation). */
        if ( setupAuxiliaryWorkingSet( &oldBounds,&oldConstraints,BT_TRUE ) != SUCCESSFUL_RETURN )
            return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

        /* Factorise projected Hessian
        * this now handles all special cases (no active bounds/constraints, no nullspace) */
        returnvalue = computeProjectedCholesky( );

        /* leave the loop if decomposition was successful, i.e. we have
         * found an active set with positive definite null space Hessian */
        if ( returnvalue == SUCCESSFUL_RETURN )
            break;
    }

    /* adjust lb/ub if we changed the old active set in the second try
     */
    if (n_try > 0) {
		// as per setupAuxiliaryQPbounds assumptions ... oh the troubles
		for (int_t ii = 0; ii < nC; ++ii)
			Ax_l[ii] = Ax_u[ii] = Ax[ii];
        setupAuxiliaryQPbounds( &bounds, &constraints, BT_FALSE );
	}

	status = QPS_AUXILIARYQPSOLVED;

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p N e w A u x i l i a r y Q P
 */
returnValue SQProblem::setupNewAuxiliaryQP(	const real_t* const H_new, const real_t* const A_new,
											const real_t *lb_new, const real_t *ub_new, const real_t *lbA_new, const real_t *ubA_new
											)
{
	int_t nV = getNV( );
	int_t nC = getNC( );

	DenseMatrix *dA = 0;
	SymDenseMat *sH = 0;

	if ( A_new != 0 )
	{
		dA = new DenseMatrix(nC, nV, nV, (real_t*) A_new);
	}
	else
	{
		if ( nC > 0 )
			return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	if ( H_new != 0 )
		sH = new SymDenseMat(nV, nV, nV, (real_t*) H_new);

	returnValue returnvalue = setupNewAuxiliaryQP( sH,dA, lb_new,ub_new,lbA_new,ubA_new );

	if ( H_new != 0 )
		freeHessian = BT_TRUE;
	freeConstraintMatrix = BT_TRUE;

	return returnvalue;
}


END_NAMESPACE_QPOASES


/*
 *	end of file
 */
