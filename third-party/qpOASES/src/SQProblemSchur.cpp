/*
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2014 by Hans Joachim Ferreau, Andreas Potschka,
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
 *	\file src/SQProblemSchur.cpp
 *	\author Andreas Waechter and Dennis Janka, based on QProblem.cpp by Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.2
 *	\date 2012-2017
 *
 *	Implementation of the SQProblemSchur class which is able to use the newly
 *	developed online active set strategy for parametric quadratic programming.
 *	This implementation uses a Schur complement approach to solve the linear
 *	systems.
 */

#include <qpOASES/SQProblemSchur.hpp>


#ifndef __MATLAB__
# include <cstdarg>
void MyPrintf(const char* pformat, ... )
{
  va_list ap;
  va_start(ap, pformat);

  vfprintf(stdout, pformat, ap);

  va_end(ap);
}
#else
# include <mex.h>
# define MyPrintf mexPrintf
#endif


BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/


/*
 *	Q P r o b l e m
 */
SQProblemSchur::SQProblemSchur( ) : SQProblem( )
{
#ifdef SOLVER_MA57
	sparseSolver = new Ma57SparseSolver();
#elif defined SOLVER_MA27
	sparseSolver = new Ma27SparseSolver();
#elif defined SOLVER_NONE
	sparseSolver = new DummySparseSolver();
#endif

	nSmax = 0;
	nS = -1;
	S = 0;
	Q_ = 0;
	R_ = 0;
	detS = 0.0;
	rcondS = 0.0;
	schurUpdateIndex = 0;
	schurUpdate = 0;
	numFactorizations = 0;

	M_physicallength = 0;
	M_vals = 0;
	M_ir = 0;
	M_jc = 0;
}


/*
 *	Q P r o b l e m
 */
SQProblemSchur::SQProblemSchur( int_t _nV, int_t _nC, HessianType _hessianType, int_t maxSchurUpdates ) 
	: SQProblem( _nV,_nC,_hessianType, BT_FALSE )
{
	/* The interface to the sparse linear solver.  In the long run,
	   different linear solvers might be optionally chosen. */
#ifdef SOLVER_MA57
	sparseSolver = new Ma57SparseSolver();
#elif defined SOLVER_MA27
	sparseSolver = new Ma27SparseSolver();
#elif defined SOLVER_NONE
	sparseSolver = new DummySparseSolver();
#endif

	nSmax = maxSchurUpdates;
	nS = -1;
	if ( nSmax > 0 )
	{
		S = new real_t[nSmax*nSmax];
		schurUpdateIndex = new int_t[nSmax];
		schurUpdate = new SchurUpdateType[nSmax];
		Q_ = new real_t[nSmax*nSmax];
		R_ = new real_t[nSmax*nSmax];
		M_physicallength = 10*nSmax;  /* TODO: Decide good default. */
		M_vals = new real_t[M_physicallength];
		M_ir = new sparse_int_t[M_physicallength];
		M_jc = new sparse_int_t[nSmax+1];
		detS = 1.0;
		rcondS = 1.0;
	}
	else
	{
		S = 0;
		Q_ = 0;
		R_ = 0;
		detS = 0.0;
		rcondS = 0.0;
		schurUpdateIndex = 0;
		schurUpdate = 0;
		M_physicallength = 0;
		M_vals = 0;
		M_ir = 0;
		M_jc = 0;
	}
	numFactorizations = 0;
}


/*
 *	Q P r o b l e m
 */
SQProblemSchur::SQProblemSchur( const SQProblemSchur& rhs ) : SQProblem( rhs )
{
#ifdef SOLVER_MA57
	sparseSolver = new Ma57SparseSolver();
#elif defined SOLVER_MA27
	sparseSolver = new Ma27SparseSolver();
#elif defined SOLVER_NONE
	sparseSolver = new DummySparseSolver();
#endif
	copy( rhs );
}


/*
 *	~ Q P r o b l e m
 */
SQProblemSchur::~SQProblemSchur( )
{
	delete sparseSolver;

	clear( );
}


/*
 *	o p e r a t o r =
 */
SQProblemSchur& SQProblemSchur::operator=( const SQProblemSchur& rhs )
{
	if ( this != &rhs )
	{
		clear( );
		SQProblem::operator=( rhs );
		copy( rhs );
	}
	return *this;
}


/*
 *	r e s e t
 */
returnValue SQProblemSchur::reset( )
{
	/* AW: We probably want to avoid resetting factorization in QProblem */
	if ( SQProblem::reset( ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_RESET_FAILED );

	sparseSolver->reset();
	nS = -1;

	return SUCCESSFUL_RETURN;
}


/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/

/*
 *	c l e a r
 */
returnValue SQProblemSchur::clear( )
{
	nSmax = 0;
	nS = -1;
	detS = 0.0;
	rcondS = 0.0;
	numFactorizations = 0;
	delete [] S; S=0;
	delete [] Q_; Q_=0;
	delete [] R_; R_=0;
	delete [] schurUpdateIndex; schurUpdateIndex=0;
	delete [] schurUpdate; schurUpdate=0;
	M_physicallength = 0;
	delete [] M_vals; M_vals=0;
	delete [] M_ir; M_ir=0;
	delete [] M_jc; M_jc=0;

	return SUCCESSFUL_RETURN;
}


/*
 *	c o p y
 */
returnValue SQProblemSchur::copy(	const SQProblemSchur& rhs
									)
{
	int_t i, j, length;

	*sparseSolver = *(rhs.sparseSolver);

	nS = rhs.nS;
	nSmax = rhs.nSmax;
	if ( nSmax > 0 )
	{
		detS = rhs.detS;
		rcondS = rhs.rcondS;
		S = new real_t[nSmax*nSmax];
		Q_ = new real_t[nSmax*nSmax];
		R_ = new real_t[nSmax*nSmax];
		schurUpdateIndex = new int_t[nSmax];
		schurUpdate = new SchurUpdateType[nSmax];

		if ( nS>0 )
		{
			for ( i=0; i<nS; i++)
				for ( j=0; j<nS; j++)
				{
					S[i*nSmax + j] = rhs.S[i*nSmax + j];
					Q_[i*nSmax + j] = rhs.Q_[i*nSmax + j];
					R_[i*nSmax + j] = rhs.R_[i*nSmax + j];
				}

			memcpy( schurUpdateIndex, rhs.schurUpdateIndex, ((unsigned int)nS)*sizeof(int_t));
			memcpy( schurUpdate, rhs.schurUpdate, ((unsigned int)nS)*sizeof(SchurUpdateType));
		}

		M_physicallength = rhs.M_physicallength;
		if ( M_physicallength>0 )
		{
			M_vals = new real_t[M_physicallength];
			M_ir = new sparse_int_t[M_physicallength];
			M_jc = new sparse_int_t[nSmax+1];

			if ( nS>0 )
			{
				memcpy(M_jc, rhs.M_jc, ((unsigned int)(nS+1))*sizeof(sparse_int_t));
				length = M_jc[nS];
				memcpy(M_vals, rhs.M_vals, ((unsigned int)length)*sizeof(real_t));
				memcpy(M_ir, rhs.M_ir, ((unsigned int)length)*sizeof(sparse_int_t));
			}
			else if ( nS==0 )
				M_jc[0] = rhs.M_jc[0];
		}
	}
	else
	{
		S = 0;
		Q_ = 0;
		R_ = 0;
		detS = 0.0;
		rcondS = 0.0;
		schurUpdateIndex = 0;
		schurUpdate = 0;
		M_physicallength = 0;
		M_vals = 0;
		M_ir = 0;
		M_jc = 0;
	}
	numFactorizations = rhs.numFactorizations;

	boundsFreeStart = rhs.boundsFreeStart;
	constraintsActiveStart = rhs.constraintsActiveStart;

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y Q P
 */
returnValue SQProblemSchur::setupAuxiliaryQP(	SymmetricMatrix *H_new,
												Matrix *A_new,
												const real_t *lb_new,
												const real_t *ub_new,
												const real_t *lbA_new,
												const real_t *ubA_new
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

	/* 2) Set new Hessian matrix,determine Hessian type and
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
		if ( H != 0 )
			return THROWERROR( RET_NO_HESSIAN_SPECIFIED );
	}

	/* 3) Setup QP gradient. */
	if ( setupAuxiliaryQPgradient( ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

	/* II) SETUP WORKING SET AND MATRIX FACTORISATION: */

	/* 1) Check if current active set is linearly independent and has the correct inertia */
	returnvalue = resetSchurComplement( BT_FALSE );
	int_t neig = sparseSolver->getNegativeEigenvalues( );

	if ( returnvalue == SUCCESSFUL_RETURN && neig == getNAC( ) )
	{
		/* a) This means the proposed working set is linearly independent and
		 *    leaves no zero curvature exposed in the nullspace and can be used to start QP solve. */
		if ( options.printLevel == PL_HIGH )
			MyPrintf( "In hotstart for new matrices, old working set is linearly independent and has correct inertia.\n");

		status = QPS_AUXILIARYQPSOLVED;
		return SUCCESSFUL_RETURN;
	}
	else if ( returnvalue == SUCCESSFUL_RETURN && neig > getNAC( ) )
	{
		/* b) KKT matrix has too many negative eigenvalues. Try to correct the inertia by adding bounds (reduce nullspace dimension). */
		if ( options.printLevel == PL_HIGH )
			MyPrintf( "WARNING: In hotstart for new matrices, reduced Hessian for initial working set has %i negative eigenvalues, should be %i.\n", neig, getNAC( ) );

		/* If enabling inertia correction is disabled, exit here */
		if ( options.enableInertiaCorrection )
		{
			returnvalue = correctInertia();
			if ( returnvalue == SUCCESSFUL_RETURN )
			{
				status = QPS_AUXILIARYQPSOLVED;
				return SUCCESSFUL_RETURN;
			}
		}
		else
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
	}

	/* 2) If inertia correction has failed or factorization yielded some other error,
	 *    try to rebuild the active set with all simple bounds set according to initialStatusBounds
	 *    (Note: in exact arithmetic, this cannot happen) */
	if ( options.printLevel == PL_HIGH )
		MyPrintf( "WARNING: hotstart for old active set failed. Trying to rebuild a working set.\n");

	Bounds      oldBounds      = bounds;
	Constraints oldConstraints = constraints;

	/* Move all inactive variables to a bound */
	for ( i=0; i<nV; i++ )
	{
		#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
		if ( bounds.getType( i ) == ST_EQUALITY )
		{
			if ( oldBounds.setStatus( i,ST_LOWER ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
		}
		else
		#endif
		{
			if ( oldBounds.getStatus( i ) == ST_INACTIVE )
				if ( oldBounds.setStatus( i, options.initialStatusBounds ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
		}
	}

	/* Set all equalities active */
	#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
	for( i=0; i<nC; ++i )
	{
		if ( constraints.getType( i ) == ST_EQUALITY )
			if ( oldConstraints.setStatus( i,ST_LOWER ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
	}
	#endif

	/* Set all inequalities inactive */
	for( i=0; i<nC; ++i )
	{
		if ( constraints.getType( i ) != ST_EQUALITY )
			if ( oldConstraints.setStatus( i,ST_INACTIVE ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
	}

	/* Reset bounds and constraints */
	bounds.init( nV );
	constraints.init( nC );

	if ( setupSubjectToType(lb_new,ub_new,lbA_new,ubA_new ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

	if ( bounds.setupAllFree( ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

	if ( constraints.setupAllInactive( ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

	/* Setup working sets afresh. */
	if ( setupAuxiliaryWorkingSet( &oldBounds,&oldConstraints,BT_TRUE ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

	/* adjust lb/ub */
	for (int_t ii = 0; ii < nC; ++ii)
		Ax_l[ii] = Ax_u[ii] = Ax[ii];
	setupAuxiliaryQPbounds (&bounds, &constraints, BT_FALSE);

	status = QPS_AUXILIARYQPSOLVED;

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y W o r k i n g S e t
 */
returnValue SQProblemSchur::setupAuxiliaryWorkingSet(	const Bounds* const auxiliaryBounds,
														const Constraints* const auxiliaryConstraints,
														BooleanType setupAfresh
														)
{
	int_t i;
	int_t nV = getNV( );
	int_t nC = getNC( );

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

	/* I.) REMOVE INEQUALITY BOUNDS/CONSTRAINTS */

	/* I.1) Remove inequality bounds that are active now but shall be
	 *      inactive or active at the other bound according to auxiliaryBounds */
	for( i=0; i<nV; ++i )
	{
		if ( ( bounds.getStatus( i ) == ST_LOWER ) && ( auxiliaryBounds->getStatus( i ) != ST_LOWER ) )
			if ( bounds.moveFixedToFree( i ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_WORKINGSET_FAILED );

		if ( ( bounds.getStatus( i ) == ST_UPPER ) && ( auxiliaryBounds->getStatus( i ) != ST_UPPER ) )
			if ( bounds.moveFixedToFree( i ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
	}

	/* I.2.) Remove inequality constraints that are active now but shall be
	 *       inactive or active at the other bound according to auxiliaryConstraints */
	for( i=0; i<nC; ++i )
	{
		if ( ( constraints.getStatus( i ) == ST_LOWER ) && ( auxiliaryConstraints->getStatus( i ) != ST_LOWER ) )
			if ( constraints.moveActiveToInactive( i ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_WORKINGSET_FAILED );

		if ( ( constraints.getStatus( i ) == ST_UPPER ) && ( auxiliaryConstraints->getStatus( i ) != ST_UPPER ) )
			if ( constraints.moveActiveToInactive( i ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
	}

	/* II.) ADD BOUNDS/CONSTRAINTS */

	/* II.1.) Add bounds according to auxiliaryBounds */
	for( i=0; i<nV; ++i )
	{
		if ( ( bounds.getStatus( i ) == ST_INACTIVE ) && ( auxiliaryBounds->getStatus( i ) != ST_INACTIVE ) )
			if ( bounds.moveFreeToFixed( i, auxiliaryBounds->getStatus( i ) ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
	}

	/* II.2.) Add constraints according to auxiliaryConstraints */
	for( i=0; i<nC; ++i )
	{
		if ( ( constraints.getStatus( i ) == ST_INACTIVE ) && ( auxiliaryConstraints->getStatus( i ) != ST_INACTIVE ) )
			if ( constraints.moveInactiveToActive( i,auxiliaryConstraints->getStatus( i ) ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
	}

	/* III) FACTORIZATION */

	/* III.1.) Factorize (resolves linear dependency) */
	if( resetSchurComplement( BT_FALSE ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_SETUP_WORKINGSET_FAILED );

	/* III.2.) Check if inertia is correct. If so, we now have a linearly independent working set with a pos def reduced Hessian */
	int_t neig = sparseSolver->getNegativeEigenvalues( );
	if ( neig == getNAC( ) )
	{
		/* We now have a linearly independent working set with a pos def reduced Hessian.
		 * We need to correct the QP bounds and gradient after this. */
		return SUCCESSFUL_RETURN;
	}

	/* IV.) INERTIA CORRECTION IF NECESSARY */

	/* We now have a fresh factorization and can start the usual inertia correction routine */
	if ( options.printLevel == PL_HIGH )
		MyPrintf( "WARNING: In setupAuxiliaryWorkingSet: Initial working set reduced Hessian has %i negative eigenvalues, should be %i.\n", neig, getNAC( ) );

	if ( options.enableInertiaCorrection == BT_TRUE )
		return correctInertia( );
	else
		return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
}


/*
 *	c h o l e s k y D e c o m p o s i t i o n P r o j e c t e d
 */
returnValue SQProblemSchur::computeProjectedCholesky( )
{
	return SUCCESSFUL_RETURN;
}


/*
 *	c o m p u t e I n i t i a l C h o l e s k y
 */
returnValue SQProblemSchur::computeInitialCholesky( )
{
	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p T Q f a c t o r i s a t i o n
 */
returnValue SQProblemSchur::setupTQfactorisation( )
{
	return SUCCESSFUL_RETURN;
}


/*
 *	a d d C o n s t r a i n t
 */
returnValue SQProblemSchur::addConstraint(	int_t number, 
											SubjectToStatus C_status,
											BooleanType updateCholesky,
											BooleanType ensureLI
											)
{
	int_t idxDeleted = -1;

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
	if ( ensureLI == BT_TRUE )
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

	/* IV) UPDATE INDICES */
	tabularOutput.idxAddC = number;
	if ( constraints.moveInactiveToActive( number,C_status ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_ADDCONSTRAINT_FAILED );

	/* Also update the Schur complement. */

	/* First check if this constraint had been removed before. In that
	   case delete this constraint from the Schur complement. */
	bool found = false;
	for ( int_t i=0; i<nS; i++ )
	{
		if ( schurUpdate[i] == SUT_ConRemoved && number == schurUpdateIndex[i] )
		{
			if ( deleteFromSchurComplement( i ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_ADDCONSTRAINT_FAILED );
			found = true;
			idxDeleted = i;
			break;
		}
	}
	if ( !found )
	{
		if ( nS < 0 || nS==nSmax )
		{
			/* The schur complement has become too large, reset. */
			/* Correct inertia if necessary. */
			returnValue retval = resetSchurComplement( BT_TRUE );
			if ( retval != SUCCESSFUL_RETURN )
			{
				if ( retval == RET_KKT_MATRIX_SINGULAR && options.printLevel == PL_HIGH )
					MyPrintf( "In addConstraint: KKT matrix singular when resetting Schur complement\n" );
				else if ( options.printLevel == PL_HIGH )
					MyPrintf( "In addConstraint, resetSchurComplement failed with retval = %d\n", retval);
				return THROWERROR( RET_ADDCONSTRAINT_FAILED );
			}
			found = true;
		}
		else
		{
			/* If the constraint was not yet in Schur complement, add it now. */
			int_t nFRStart = boundsFreeStart.getLength();
			int_t* FR_idxStart;
			boundsFreeStart.getNumberArray( &FR_idxStart );

			sparse_int_t* MNpos = new sparse_int_t[nFRStart+nS]; // This is an overestimate
			real_t* MNvals = new real_t[nFRStart+nS];

			int_t* irn = new int_t[nFRStart+nS];
			int_t* jcn = new int_t[nFRStart+nS];
			real_t* vals = new real_t[nFRStart+nS];
			int_t* icolsNumber = new int_t[nFRStart+nS];
			int_t* icolsSIdx = new int_t[nS];

			for ( int_t i=0; i<nFRStart; i++)
				icolsNumber[i] = FR_idxStart[i];

			int_t icolsLength = nFRStart;
			for ( int_t i=0; i<nS; i++)
				if ( schurUpdate[i] == SUT_VarFreed )
				{
					icolsNumber[icolsLength] = schurUpdateIndex[i];
					icolsSIdx[icolsLength-nFRStart] = i;
					icolsLength++;
				}

			if ( constraintProduct != 0 )
			{
				MyPrintf( "In SQProblemSchur::addConstraint, constraintProduct not yet implemented.\n");
				return THROWERROR(RET_NOT_YET_IMPLEMENTED);
			}
			int_t numNonzerosA;
			A->getSparseSubmatrix( 1, &number, icolsLength, icolsNumber, 0, 0, numNonzerosA, irn, jcn, vals );
			delete [] irn;

			int_t numNonzerosM = 0;
			int_t numNonzerosN = 0;
			for ( int_t i=0; i<numNonzerosA; i++ )
				if ( jcn[i] < nFRStart )
				{
					MNpos[numNonzerosM] = jcn[i];
					MNvals[numNonzerosM] = vals[i];
					numNonzerosM++;
				}
				else
				{
					MNpos[nFRStart+numNonzerosN] = icolsSIdx[jcn[i]-nFRStart];
					MNvals[nFRStart+numNonzerosN] = vals[i];
					numNonzerosN++;
				}

			returnValue returnvalue = addToSchurComplement( number, SUT_ConAdded, numNonzerosM, MNpos, MNvals, numNonzerosN, MNpos+nFRStart, MNvals+nFRStart, 0.0 );

			delete [] icolsSIdx;
			delete [] icolsNumber;
			delete [] vals;
			delete [] jcn;
			delete [] MNvals;
			delete [] MNpos;

			if ( returnvalue != SUCCESSFUL_RETURN )
				return THROWERROR( RET_ADDCONSTRAINT_FAILED );

			found = true;
		}
	}

	if ( !found )
		return THROWERROR( RET_ADDCONSTRAINT_FAILED );

	updateSchurQR( idxDeleted );

	/* If reciprocal of condition number becomes to small, refactorize KKT matrix */
	if( rcondS < options.rcondSMin )
	{
		returnValue retval = resetSchurComplement( BT_TRUE );
		if ( retval != SUCCESSFUL_RETURN )
		{
			if ( retval == RET_KKT_MATRIX_SINGULAR  && options.printLevel == PL_HIGH )
				MyPrintf( "In addConstraint: KKT matrix singular when resetting Schur complement\n" );
			else if ( options.printLevel == PL_HIGH )
				MyPrintf( "In addConstraint, resetSchurComplement failed with retval = %d\n", retval);
			return THROWERROR( RET_ADDCONSTRAINT_FAILED );
		}
	}

	return SUCCESSFUL_RETURN;
}



/*
 *	a d d C o n s t r a i n t _ c h e c k L I
 */
returnValue SQProblemSchur::addConstraint_checkLI( int_t number )
{
	/* Get space for the multipliers xi in linear independence test */
	int_t nAC = getNAC();
	int_t nFX = getNFX();
	real_t *xiC = new real_t[nAC];
	real_t *xiB = new real_t[nFX];

	/* I) Check if new constraint is linearly independent from the active ones. */
	returnValue returnvalueCheckLI = addConstraint_checkLISchur( number, xiC, xiB );

	delete [] xiB;
	delete [] xiC;

	return returnvalueCheckLI;
}


/*
 *	a d d C o n s t r a i n t _ c h e c k L I S c h u r
 */
returnValue SQProblemSchur::addConstraint_checkLISchur( int_t number, real_t* xiC, real_t* xiB )
{
	returnValue returnvalue = RET_LINEARLY_DEPENDENT;

	int_t ii;
	int_t nV  = getNV( );
	int_t nFR = getNFR( );
	int_t nC  = getNC( );
	int_t nAC = getNAC();
	int_t nFX = getNFX();
	int_t *FR_idx;

	bounds.getFree( )->getNumberArray( &FR_idx );

	/* For the Schur complement version we only use options.enableFullLITests = TRUE */
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
		real_t *delta_yAC = xiC;
		real_t *delta_yFX = xiB;

		bounds.getFixed( )->getNumberArray( &FX_idx );
		constraints.getActive( )->getNumberArray( &AC_idx );
		constraints.getInactive( )->getNumberArray( &IAC_idx );

		int_t dim = (nC>nV)?nC:nV;
		real_t *nul = new real_t[dim];
		for (ii = 0; ii < dim; ++ii)
			nul[ii]=0.0;

		A->getRow (number, 0, 1.0, delta_g);

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

		delete[] delta_xFR;
		delete[] delta_xFX;
		delete[] delta_g;

	}
	return THROWINFO( returnvalue );
}


/*
 *	a d d C o n s t r a i n t _ e n s u r e L I
 */
returnValue SQProblemSchur::addConstraint_ensureLI( int_t number, SubjectToStatus C_status )
{
	/* Get space for the multipliers xi in linear independence test */
	int_t nAC = getNAC();
	int_t nFX = getNFX();
	real_t *xiC = new real_t[nAC];
	real_t *xiB = new real_t[nFX];

	/* I) Check if new constraint is linearly independent from the active ones. */
	returnValue returnvalueCheckLI = addConstraint_checkLISchur( number, xiC, xiB );

	if ( returnvalueCheckLI == RET_INDEXLIST_CORRUPTED )
	{
		delete [] xiB;
		delete [] xiC;
		return THROWERROR( RET_ENSURELI_FAILED );
	}

	if ( returnvalueCheckLI == RET_LINEARLY_INDEPENDENT )
	{
		delete [] xiB;
		delete [] xiC;
		return SUCCESSFUL_RETURN;
	}

 	/* II) NEW BOUND IS LINEARLY DEPENDENT: */
	/* 1) Coefficients of linear combination, have already been computed, but we need to correct the sign.  */
	int_t i, ii;

	if ( C_status != ST_LOWER )
	{
		for( i=0; i<nAC; ++i )
			xiC[i] = -xiC[i];
		for( i=0; i<nFX; ++i )
			xiB[i] = -xiB[i];
	}

	int_t nV  = getNV( );

	int_t* FX_idx;
	bounds.getFixed( )->getNumberArray( &FX_idx );

	int_t* AC_idx;
	constraints.getActive( )->getNumberArray( &AC_idx );

	real_t* num = new real_t[nV];

	real_t y_min = options.maxDualJump;
	int_t y_min_number = -1;
	int_t y_min_number_bound = -1;
	BooleanType y_min_isBound = BT_FALSE;

	returnValue returnvalue = SUCCESSFUL_RETURN;

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

	#ifndef __XPCTARGET__
	/* setup output preferences */
	char messageString[80];
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
			#ifndef __XPCTARGET__
			snprintf( messageString,80,"bound no. %d.",(int)y_min_number );
			getGlobalMessageHandler( )->throwInfo( RET_REMOVE_FROM_ACTIVESET,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
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
			#ifndef __XPCTARGET__
			snprintf( messageString,80,"constraint no. %d.",(int)y_min_number );
			getGlobalMessageHandler( )->throwInfo( RET_REMOVE_FROM_ACTIVESET,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
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
	delete [] xiB;
	delete [] xiC;

	getGlobalMessageHandler( )->throwInfo( RET_LI_RESOLVED,0,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );

	return (returnvalue != SUCCESSFUL_RETURN) ? THROWERROR (returnvalue) : returnvalue;
}


/*
 *	a d d B o u n d
 */
returnValue SQProblemSchur::addBound(	int_t number,
										SubjectToStatus B_status,
										BooleanType updateCholesky,
										BooleanType ensureLI
										)
{
	int_t idxDeleted = -1;

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
	if ( ensureLI == BT_TRUE )
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

	/* II) UPDATE INDICES */
	tabularOutput.idxAddB = number;
	if ( bounds.moveFreeToFixed( number,B_status ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_ADDBOUND_FAILED );

	/* Also update the Schur complement. */

	/* First check if this variable had been freed before. In that
	   case delete this variable from the Schur complement. */
	bool found = false;
	for ( int_t i=0; i<nS; i++ )
	{
		if ( schurUpdate[i] == SUT_VarFreed && number == schurUpdateIndex[i] )
		{
			if ( deleteFromSchurComplement( i ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_ADDBOUND_FAILED );
			found = true;
			idxDeleted = i;
			break;
		}
	}
	if ( !found )
	{
		if ( nS < 0 || nS==nSmax )
		{
			/* The schur complement has become too large, reset. */
			/* Correct inertia if necessary. */
			returnValue retval = resetSchurComplement( BT_TRUE );
			if ( retval != SUCCESSFUL_RETURN )
			{
				if ( retval == RET_KKT_MATRIX_SINGULAR && options.printLevel == PL_HIGH )
					MyPrintf( "In addBound: KKT matrix singular when resetting Schur complement\n" );
				else if ( options.printLevel == PL_HIGH )
					MyPrintf( "In addBound, resetSchurComplement failed with retval = %d\n", retval);
				return THROWERROR( RET_ADDBOUND_FAILED );
			}
			found = true;
		}
		else
		{
			/* If the variable was not yet in Schur complement, add it now. */
			int_t nFRStart = boundsFreeStart.getLength();
			int_t* FR_idxStart;
			boundsFreeStart.getNumberArray( &FR_idxStart );
			for ( int_t i=0; i<nFRStart; i++ )
				if ( FR_idxStart[i] == number )
				{
					real_t one = 1.0;
					sparse_int_t pos = i;
					if ( addToSchurComplement( number, SUT_VarFixed, 1, &pos, &one, 0, 0, 0, 0.0 ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_ADDBOUND_FAILED );
					found = true;
					break;
				}
		}
	}

	if ( !found )
		return THROWERROR( RET_ADDBOUND_FAILED );

	updateSchurQR( idxDeleted );

	/* If reciprocal of condition number becomes to small, refactorize KKT matrix */
	if( rcondS < options.rcondSMin )
	{
		returnValue retval = resetSchurComplement( BT_TRUE );
		if ( retval != SUCCESSFUL_RETURN )
		{
			if ( retval == RET_KKT_MATRIX_SINGULAR && options.printLevel == PL_HIGH )
				MyPrintf( "In addBound: KKT matrix singular when resetting Schur complement\n" );
			else if ( options.printLevel == PL_HIGH )
				MyPrintf( "In addBound, resetSchurComplement failed with retval = %d\n", retval);
			return THROWERROR( RET_ADDCONSTRAINT_FAILED );
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	a d d B o u n d _ c h e c k L I
 */
returnValue SQProblemSchur::addBound_checkLI( int_t number )
{
	/* Get space for the multipliers xi in linear independence test */
	int_t nAC = getNAC();
	int_t nFX = getNFX();
	real_t *xiC = new real_t[nAC];
	real_t *xiB = new real_t[nFX];

	/* I) Check if new constraint is linearly independent from the active ones. */
	returnValue returnvalueCheckLI = addBound_checkLISchur( number, xiC, xiB );

	delete [] xiB;
	delete [] xiC;

	return returnvalueCheckLI;
}

/*
 *	a d d B o u n d _ c h e c k L I S c h u r
 */
returnValue SQProblemSchur::addBound_checkLISchur( int_t number, real_t* xiC, real_t* xiB )
{
	returnValue returnvalue = RET_LINEARLY_DEPENDENT;


	int_t ii;
	int_t nV  = getNV( );
	int_t nFR = getNFR( );
	int_t nC  = getNC( );
	int_t nAC = getNAC();
	int_t nFX = getNFX();
	int_t *FR_idx;

	bounds.getFree( )->getNumberArray( &FR_idx );

	/* For the Schur complement version we only use options.enableFullLITests = TRUE */
	{
		/*
		 * expensive LI test. Backsolve with refinement using special right
		 * hand side. This gives an estimate for what should be considered
		 * "zero". We then check linear independence relative to this estimate.
		 */

		real_t *delta_g   = new real_t[nV];
		real_t *delta_xFX = new real_t[nFX];
		real_t *delta_xFR = new real_t[nFR];
		real_t *delta_yAC = xiC;
		real_t *delta_yFX = xiB;

		for (ii = 0; ii < nV; ++ii)
			delta_g[ii] = 0.0;
		delta_g[number] = 1.0;

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
		delete[] delta_xFR;
		delete[] delta_xFX;
		delete[] delta_g;

	}
	return THROWINFO( returnvalue );
}


/*
 *	a d d B o u n d _ e n s u r e L I
 */
returnValue SQProblemSchur::addBound_ensureLI( int_t number, SubjectToStatus B_status )
{
	/* Get space for the multipliers xi in linear independence test */
	int_t nAC = getNAC();
	int_t nFX = getNFX();
	real_t *xiC = new real_t[nAC];
	real_t *xiB = new real_t[nFX];

	/* I) Check if new constraint is linearly independent from the active ones. */
	returnValue returnvalueCheckLI = addBound_checkLISchur( number, xiC, xiB );

	if ( returnvalueCheckLI == RET_INDEXLIST_CORRUPTED )
	{
		delete [] xiB;
		delete [] xiC;
		return THROWERROR( RET_ENSURELI_FAILED );
	}

	if ( returnvalueCheckLI == RET_LINEARLY_INDEPENDENT )
	{
		delete [] xiB;
		delete [] xiC;
		return SUCCESSFUL_RETURN;
	}

 	/* II) NEW BOUND IS LINEARLY DEPENDENT: */
	/* 1) Coefficients of linear combination, have already been computed, but we need to correct the sign.  */
	int_t i, ii;

	if ( B_status != ST_LOWER )
	{
		for( i=0; i<nAC; ++i )
			xiC[i] = -xiC[i];
		for( i=0; i<nFX; ++i )
			xiB[i] = -xiB[i];
	}

	int_t nV  = getNV( );

	int_t* FX_idx;
	bounds.getFixed( )->getNumberArray( &FX_idx );

	int_t* AC_idx;
	constraints.getActive( )->getNumberArray( &AC_idx );

	real_t* num = new real_t[nV];

	real_t y_min = options.maxDualJump;
	int_t y_min_number = -1;
	int_t y_min_number_bound = -1;
	BooleanType y_min_isBound = BT_FALSE;

	returnValue returnvalue = SUCCESSFUL_RETURN;

	/* III) DETERMINE CONSTRAINT/BOUND TO BE REMOVED. */

	/* 1) Constraints. */
	for( i=0; i<nAC; ++i )
	{
		ii = AC_idx[i];
		num[i] = y[nV+ii];
	}

	performRatioTest( nAC,AC_idx,&constraints, num,xiC, options.epsNum, options.epsDen, y_min,y_min_number );

	/* 2) Bounds. */
	for( i=0; i<nFX; ++i )
	{
		ii = FX_idx[i];
		num[i] = y[ii];
	}

	performRatioTest( nFX,FX_idx,&bounds, num,xiB, options.epsNum, options.epsDen, y_min,y_min_number_bound );

	if ( y_min_number_bound >= 0 )
	{
		y_min_number = y_min_number_bound;
		y_min_isBound = BT_TRUE;
	}

	/* IV) REMOVE CONSTRAINT/BOUND FOR RESOLVING LINEAR DEPENDENCE: */
	char messageString[80];

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
			#ifndef __XPCTARGET__
			snprintf( messageString,80,"bound no. %d.",(int)y_min_number );
			getGlobalMessageHandler( )->throwInfo( RET_REMOVE_FROM_ACTIVESET,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
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
			#ifndef __XPCTARGET__
			snprintf( messageString,80,"constraint no. %d.",(int)y_min_number );
			getGlobalMessageHandler( )->throwInfo( RET_REMOVE_FROM_ACTIVESET,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
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
	delete[] xiC;

	getGlobalMessageHandler( )->throwInfo( RET_LI_RESOLVED,0,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );

	return (returnvalue != SUCCESSFUL_RETURN) ? THROWERROR (returnvalue) : returnvalue;
}



/*
 *	r e m o v e C o n s t r a i n t
 */
returnValue SQProblemSchur::removeConstraint(	int_t number,
												BooleanType updateCholesky,
												BooleanType allowFlipping,
												BooleanType ensureNZC
												)
{
	returnValue returnvalue = SUCCESSFUL_RETURN;

	int_t sModType = 0;
	int_t idxDeleted = -1;
	SubjectToStatus oldStatus;
	real_t oldDet, newDet;

	/* consistency check */
	if ( ( getStatus( ) == QPS_NOTINITIALISED )    ||
		 ( getStatus( ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( getStatus( ) == QPS_HOMOTOPYQPSOLVED )  ||
 		 ( getStatus( ) == QPS_SOLVED )            )
	{
		return THROWERROR( RET_UNKNOWN_BUG );
	}

	/* some definitions */
	int_t nAC = getNAC( );
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

	/* N) PERFORM ZERO CURVATURE TEST. */
	if (ensureNZC == BT_TRUE)
	{
		returnvalue = ensureNonzeroCurvature(BT_FALSE, number, exchangeHappened, addBoundNotConstraint, addIdx, addStatus);

		if (returnvalue != SUCCESSFUL_RETURN)
			return returnvalue;
	}

	/* save old constraint status and determinant of old S for flipping strategy */
	oldStatus = constraints.getStatus( number );
	oldDet = detS;

	/* I) UPDATE INDICES */
	tabularOutput.idxRemC = number;
	if ( constraints.moveActiveToInactive( number ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_REMOVECONSTRAINT_FAILED );

	/* Also update the Schur complement. */

	/* First check if this constraint had been added before. In that
	   case delete this constraint from the Schur complement. */
	bool found = false;
	for ( int_t i=0; i<nS; i++ )
	{
		if ( schurUpdate[i] == SUT_ConAdded && number == schurUpdateIndex[i] )
		{
			if ( deleteFromSchurComplement( i, BT_TRUE ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_REMOVECONSTRAINT_FAILED );
			found = true;
			idxDeleted = i;
			sModType = 2;
			break;
		}
	}
	if ( !found )
	{
		if ( nS < 0 || nS==nSmax )
		{
			/* The schur complement has become too large, reset. */
			/* Don't check inertia here, may be corrected later! */
			returnValue retval = resetSchurComplement( BT_FALSE );
			if ( retval != SUCCESSFUL_RETURN )
			{
				if ( retval == RET_KKT_MATRIX_SINGULAR && options.printLevel == PL_HIGH )
					MyPrintf( "In removeConstraint: KKT matrix singular when resetting Schur complement\n" );
				else if ( options.printLevel == PL_HIGH )
					MyPrintf( "In removeConstraint, resetSchurComplement failed with retval = %d\n", retval);
				return THROWERROR( RET_REMOVECONSTRAINT_FAILED );
			}
			found = true;
			sModType = 3;
		}
		else
		{
			/* If the constraint was not yet in Schur complement, add it now. */
			int_t nFRStart = boundsFreeStart.getLength();
			int_t nACStart = constraintsActiveStart.getLength();
			int_t* AC_idxStart;
			constraintsActiveStart.getNumberArray( &AC_idxStart );

			for ( int_t i=0; i<nACStart; i++ )
				if ( AC_idxStart[i] == number )
				{
					real_t one = 1.0;
					sparse_int_t pos = nFRStart+i;
					if ( addToSchurComplement( number, SUT_ConRemoved, 1, &pos, &one, 0, 0, 0, 0.0 ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_REMOVECONSTRAINT_FAILED );
					found = true;
					break;
				}

			sModType = 1;
		}
	}

	if ( !found )
		return THROWERROR( RET_REMOVECONSTRAINT_FAILED );

	/* Now we have a new Schur complement (might be smaller, larger, or empty). Update QR factorization. */

	/* Flipping bounds strategy */
	if ( ( options.enableFlippingBounds == BT_TRUE ) && ( allowFlipping == BT_TRUE ) && ( exchangeHappened == BT_FALSE ) )
	{
		if ( sModType == 1 )
		{/* Case 1: We added a row and column to S. */

			/* Check if a direction of negative curvature showed up, i.e. determinants have THE SAME sign */
			newDet = calcDetSchur( idxDeleted );

			if ( oldDet * newDet > 0 )
			{
				hessianType = HST_SEMIDEF;

				/* Restore old S */
				nS--;

				/* Flip bounds */
				tabularOutput.idxAddC = number;
				tabularOutput.excAddC = 2;
				switch ( oldStatus )
				{
					case ST_LOWER:
						constraints.moveInactiveToActive( number, ST_UPPER );
						ubA[number] = lbA[number];
						Ax_l[number] = -Ax_u[number];
						break;
					case ST_UPPER:
						constraints.moveInactiveToActive( number, ST_LOWER );
						lbA[number] = ubA[number];
						Ax_u[number] = -Ax_l[number];
						break;
					default:
						return THROWERROR( RET_REMOVECONSTRAINT_FAILED );
				}
			}
			else
			{/* Determinants have the correct sign, compute QR of new (larger) S */
				updateSchurQR( idxDeleted );
			}
		}
		else if ( sModType == 2 )
		{/* Case 2: We deleted a row and column of S. */

			/* Check if a direction of negative curvature showed up, i.e. determinants have DIFFERENT signs */
			newDet = calcDetSchur( idxDeleted );

			if ( oldDet * newDet < 0.0 )
			{
				hessianType = HST_SEMIDEF;

				/* Restore old S */
				undoDeleteFromSchurComplement( idxDeleted );

				/* Flip bounds */
				tabularOutput.idxAddC = number;
				tabularOutput.excAddC = 2;
				switch ( oldStatus )
				{
					case ST_LOWER:
						constraints.moveInactiveToActive( number, ST_UPPER );
						ubA[number] = lbA[number];
						Ax_l[number] = -Ax_u[number];
						break;
					case ST_UPPER:
						constraints.moveInactiveToActive( number, ST_LOWER );
						lbA[number] = ubA[number];
						Ax_u[number] = -Ax_l[number];
						break;
					default:
						return THROWERROR( RET_REMOVECONSTRAINT_FAILED );
				}
			}
			else
			{/* Determinants have the correct sign, compute QR of new (smaller) S */
				updateSchurQR( idxDeleted );
			}
		}
		else if ( sModType == 3 )
		{/* Case 3: S was reset. */

			/* Check inertia of new factorization given by the sparse solver: must be ( nFR, nAC, 0 ) */
			int_t neig = sparseSolver->getNegativeEigenvalues( );
			if( neig > getNAC( ) ) // Wrong inertia!
			{
				/* Flip bounds and update Schur complement */
				tabularOutput.idxAddC = number;
				tabularOutput.excAddC = 2;
				switch ( oldStatus )
				{
					case ST_LOWER:
						ubA[number] = lbA[number];
						Ax_l[number] = -Ax_u[number];
						addConstraint( number, ST_UPPER, BT_TRUE, BT_FALSE );
						break;
					case ST_UPPER:
						lbA[number] = ubA[number];
						Ax_u[number] = -Ax_l[number];
						addConstraint( number, ST_LOWER, BT_TRUE, BT_FALSE );
						break;
					default:
						return THROWERROR( RET_REMOVECONSTRAINT_FAILED );
				}
			}

			/* Check if flipping deleted the negative eigenvalue */
			if( correctInertia( ) )
				return THROWERROR( RET_REMOVECONSTRAINT_FAILED );
		}
		else
		{/* None of the three cases happened */
			return THROWERROR( RET_REMOVECONSTRAINT_FAILED );
		}
	}
	else
	{/* No flipping strategy, update QR factorization of S */
		updateSchurQR( idxDeleted );
	}

	/* If reciprocal of condition number becomes to small, refactorize KKT matrix */
	if( rcondS < options.rcondSMin )
	{
		returnValue retval = resetSchurComplement( BT_TRUE );
		if ( retval != SUCCESSFUL_RETURN )
		{
			if ( retval == RET_KKT_MATRIX_SINGULAR && options.printLevel == PL_HIGH )
				MyPrintf( "In removeConstraint: KKT matrix singular when resetting Schur complement\n" );
			else if ( options.printLevel == PL_HIGH )
				MyPrintf( "In removeConstraint, resetSchurComplement failed with retval = %d\n", retval);
			return THROWERROR( RET_ADDCONSTRAINT_FAILED );
		}
	}

	if ( exchangeHappened == BT_TRUE )
	{
		/* add bound or constraint */

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
returnValue SQProblemSchur::removeBound(	int_t number,
											BooleanType updateCholesky,
											BooleanType allowFlipping,
											BooleanType ensureNZC
											)
{
	returnValue returnvalue = SUCCESSFUL_RETURN;
	int_t addIdx;
	BooleanType addBoundNotConstraint;
	SubjectToStatus addStatus;
	BooleanType exchangeHappened = BT_FALSE;

	int_t sModType = 0;
	int_t idxDeleted = -1;
	SubjectToStatus oldStatus;
	real_t oldDet, newDet;

	/* consistency checks */
	if ( bounds.getStatus( number ) == ST_INACTIVE )
		return THROWERROR( RET_BOUND_NOT_ACTIVE );

	if ( ( getStatus( ) == QPS_NOTINITIALISED )    ||
		 ( getStatus( ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( getStatus( ) == QPS_HOMOTOPYQPSOLVED )  ||
 		 ( getStatus( ) == QPS_SOLVED )           )
	{
		return THROWERROR( RET_UNKNOWN_BUG );
	}

	/* N) PERFORM ZERO CURVATURE TEST. */
	if (ensureNZC == BT_TRUE)
	{
		returnvalue = ensureNonzeroCurvature(BT_TRUE, number, exchangeHappened, addBoundNotConstraint, addIdx, addStatus);

		if (returnvalue != SUCCESSFUL_RETURN)
			return returnvalue;
	}

	/* save old bound status and determinant of old S for flipping strategy */
	oldStatus = bounds.getStatus( number );
	oldDet = detS;

	/* I) UPDATE INDICES */
	tabularOutput.idxRemB = number;
	if ( bounds.moveFixedToFree( number ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_REMOVEBOUND_FAILED );

	/* Also update the Schur complement. */

	/* First check if this variable had been fixed before. In that
	   case delete this variable from the Schur complement. */
	bool found = false;
	for ( int_t i=0; i<nS; i++ )
	{
		if ( schurUpdate[i] == SUT_VarFixed && number == schurUpdateIndex[i] )
		{
			if ( deleteFromSchurComplement( i, BT_TRUE ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_REMOVEBOUND_FAILED );
			found = true;
			idxDeleted = i;
			sModType = 2;
			break;
		}
	}
	if ( !found )
	{
		if ( nS < 0 || nS==nSmax )
		{
			/* The schur complement has become too large, reset. */
			/* Don't correct inertia here, may be corrected by flipping bounds! */
			returnValue retval = resetSchurComplement( BT_FALSE );
			if ( retval != SUCCESSFUL_RETURN )
			{
				if ( retval == RET_KKT_MATRIX_SINGULAR && options.printLevel == PL_HIGH  )
					MyPrintf( "In removeBound: KKT matrix singular when resetting Schur complement\n" );
				else if ( options.printLevel == PL_HIGH )
					MyPrintf( "In removeBound, resetSchurComplement failed with retval = %d\n", retval);
				return THROWERROR( RET_REMOVEBOUND_FAILED );
			}
			found = true;
			sModType = 3;
		}
		else
		{
			/* If the variable was not yet in Schur complement, add it now. */
			int_t nFRStart = boundsFreeStart.getLength();
			int_t nACStart = constraintsActiveStart.getLength();
			int_t* FR_idxStart;
			boundsFreeStart.getNumberArray( &FR_idxStart );
			int_t* AC_idxStart;
			constraintsActiveStart.getNumberArray( &AC_idxStart );

			int_t numNonzerosM = 0;
			sparse_int_t* Mpos = new sparse_int_t[nFRStart+nACStart+nS]; // This is an overestimate
			real_t* Mvals = new real_t[nFRStart+nACStart+nS];
			int_t numNonzerosN = 0;
			sparse_int_t* Npos = new sparse_int_t[nFRStart+nACStart+nS]; // This is an overestimate
			real_t* Nvals = new real_t[nFRStart+nACStart+nS];
			real_t N_diag;

			int_t* irn = new int_t[nFRStart+nACStart+nS+1];
			int_t* jcn = new int_t[nFRStart+nACStart+nS+1];
			real_t* vals = new real_t[nFRStart+nACStart+nS+1];
			int_t iLength;
			int_t* iNumber = new int_t[nFRStart+nACStart+nS+1];
			int_t numNonzeros;
			int_t* iSIdx = new int_t[nS];

			/* First the Hessian part. */
			real_t regularisation = options.epsRegularisation;
			switch ( hessianType )
			{
				case HST_ZERO:
					N_diag = regularisation;
					break;

				case HST_IDENTITY:
					N_diag = 1.0 + regularisation;
					break;

				default:
					N_diag = regularisation;
					for ( int_t i=0; i<nFRStart; i++ )
						iNumber[i] = FR_idxStart[i];
					iLength = nFRStart;
					for ( int_t i=0; i<nS; i++ )
						if ( schurUpdate[i] == SUT_VarFreed )
						{
							iNumber[iLength] = schurUpdateIndex[i];
							iSIdx[iLength-nFRStart] = i;
							iLength++;
						}
					iNumber[iLength++] = number;

					H->getSparseSubmatrix( iLength, iNumber, 1, &number, 0, 0, numNonzeros, irn, jcn, vals );

					for ( int_t i=0; i<numNonzeros; i++ )
					{
						if ( irn[i] < nFRStart )
						{
							Mpos[numNonzerosM] = irn[i];
							Mvals[numNonzerosM] = vals[i];
							numNonzerosM++;
						}
						else if ( irn[i] != iLength-1 )
						{
							Npos[numNonzerosN] = iSIdx[irn[i] - nFRStart];
							Nvals[numNonzerosN] = vals[i];
							numNonzerosN++;
						}
						else
							N_diag += vals[i];
					}
					break;
			}

			if ( constraintProduct != 0 )
			{
				MyPrintf( "In SQProblemSchur::removeBound, constraintProduct not yet implemented.\n");
				return THROWERROR(RET_NOT_YET_IMPLEMENTED);
			}

			for ( int_t i=0; i<nACStart; i++ )
				iNumber[i] = AC_idxStart[i];
			iLength = nACStart;
			for ( int_t i=0; i<nS; i++ )
				if ( schurUpdate[i] == SUT_ConAdded )
				{
					iNumber[iLength] = schurUpdateIndex[i];
					iSIdx[iLength-nACStart] = i;
					iLength++;
				}

			A->getSparseSubmatrix( iLength, iNumber, 1, &number, 0, 0, numNonzeros, irn, jcn, vals );

			for ( int_t i=0; i<numNonzeros; i++ )
			{
				if ( irn[i] < nACStart )
				{
					Mpos[numNonzerosM] = irn[i] + nFRStart;
					Mvals[numNonzerosM] = vals[i];
					numNonzerosM++;
				}
				else
				{
					Npos[numNonzerosN] = iSIdx[irn[i] - nACStart];
					Nvals[numNonzerosN] = vals[i];
					numNonzerosN++;
				}
			}

			delete [] iSIdx;
			delete [] iNumber;
			delete [] vals;
			delete [] jcn;
			delete [] irn;

			returnvalue = addToSchurComplement( number, SUT_VarFreed, numNonzerosM, Mpos, Mvals, numNonzerosN, Npos, Nvals, N_diag );

			delete [] Mvals;
			delete [] Mpos;
			delete [] Nvals;
			delete [] Npos;

			if ( returnvalue != SUCCESSFUL_RETURN )
				return THROWERROR( RET_REMOVEBOUND_FAILED );

			found = true;
			sModType = 1;
		}
	}

	if ( !found )
		return THROWERROR( RET_REMOVEBOUND_FAILED );

	/* Now we have a new Schur complement (might be smaller, larger, or empty). Update QR factorization. */

	/* Flipping bounds strategy */
	if ( ( options.enableFlippingBounds == BT_TRUE ) && ( allowFlipping == BT_TRUE ) && ( exchangeHappened == BT_FALSE ) )
	{
		if ( sModType == 1 )
		{/* Case 1: We added a row and column to S. */

			/* Check if a direction of negative curvature showed up, i.e. determinants have THE SAME sign */
			newDet = calcDetSchur( idxDeleted );

			if ( oldDet * newDet > 0.0 )
			{
				hessianType = HST_SEMIDEF;

				/* Restore old S */
				nS--;

				/* Flip bounds */
				tabularOutput.idxAddB = number;
				tabularOutput.excAddB = 2;
				switch ( oldStatus )
				{
					case ST_LOWER:
						bounds.moveFreeToFixed( number, ST_UPPER );
						ub[number] = lb[number];
						break;
					case ST_UPPER:
						bounds.moveFreeToFixed( number, ST_LOWER );
						lb[number] = ub[number];
						break;
					default:
						return THROWERROR( RET_MOVING_BOUND_FAILED );
				}
			}
			else
			{/* Determinants have the correct sign, compute QR of new (larger) S */
				updateSchurQR( idxDeleted );
			}
		}
		else if ( sModType == 2 )
		{/* Case 2: We deleted a row and column of S. */

			/* Check if a direction of negative curvature showed up, i.e. determinants have DIFFERENT signs */
			newDet = calcDetSchur( idxDeleted );

			if ( oldDet * newDet < 0.0 )
			{
				hessianType = HST_SEMIDEF;

				/* Restore old S */
				undoDeleteFromSchurComplement( idxDeleted );

				/* Flip bounds */
				tabularOutput.idxAddB = number;
				tabularOutput.excAddB = 2;
				switch ( oldStatus )
				{
					case ST_LOWER:
						bounds.moveFreeToFixed( number, ST_UPPER );
						ub[number] = lb[number];
						break;
					case ST_UPPER:
						bounds.moveFreeToFixed( number, ST_LOWER );
						lb[number] = ub[number];
						break;
					default:
						return THROWERROR( RET_MOVING_BOUND_FAILED );
				}
			}
			else
			{/* Determinants have the correct sign, compute QR of new (smaller) S */
				updateSchurQR( idxDeleted );
			}
		}
		else if ( sModType == 3 )
		{/* Case 3: S was reset. */

			/* Check inertia of new factorization given by the sparse solver: must be ( nFR, nAC, 0 ) */
			int_t neig = sparseSolver->getNegativeEigenvalues( );
			if( neig > getNAC( ) ) // Wrong inertia, flip bounds!
			{
				/* Flip bounds and update Schur complement */
				tabularOutput.idxAddB = number;
				tabularOutput.excAddB = 2;
				switch ( oldStatus )
				{
					case ST_LOWER:
						ub[number] = lb[number];
						addBound( number, ST_UPPER, BT_TRUE, BT_FALSE );
						break;
					case ST_UPPER:
						lb[number] = ub[number];
						addBound( number, ST_LOWER, BT_TRUE, BT_FALSE );
						break;
					default:
						return THROWERROR( RET_MOVING_BOUND_FAILED );
				}
			}

			/* Check if flipping deleted the negative eigenvalue */
			if( correctInertia( ) )
				return THROWERROR( RET_REMOVEBOUND_FAILED );
		}
		else
		{/* None of the three cases happened */
			return THROWERROR( RET_REMOVEBOUND_FAILED );
		}
	}
	else
	{/* No flipping strategy, update QR factorization of S */
		updateSchurQR( idxDeleted );
	}

	/* If reciprocal of condition number becomes to small, refactorize KKT matrix */
	if( rcondS < options.rcondSMin )
	{
		returnValue retval = resetSchurComplement( BT_TRUE );
		if ( retval != SUCCESSFUL_RETURN )
		{
			if ( retval == RET_KKT_MATRIX_SINGULAR && options.printLevel == PL_HIGH )
				MyPrintf( "In removeBound: KKT matrix singular when resetting Schur complement\n" );
			else if ( options.printLevel == PL_HIGH )
				MyPrintf( "In removeBound, resetSchurComplement failed with retval = %d\n", retval);
			return THROWERROR( RET_ADDCONSTRAINT_FAILED );
		}
	}

	if ( exchangeHappened == BT_TRUE )
	{
		/* add bound or constraint */

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
 *	s e t u p T Q f a c t o r i s a t i o n
 */
returnValue SQProblemSchur::backsolveT( const real_t* const b, BooleanType transposed, real_t* const a ) const
{
	return THROWERROR( RET_UNKNOWN_BUG );
}


/*
 *	b a c k s o l v e R
 */
returnValue SQProblemSchur::backsolveR(	const real_t* const b, BooleanType transposed, real_t* const a 	) const
{
	return THROWERROR( RET_UNKNOWN_BUG );
}


/*
 *	b a c k s o l v e R
 */
returnValue SQProblemSchur::backsolveR(	const real_t* const b, BooleanType transposed, BooleanType removingBound, real_t* const a ) const
{
	return THROWERROR( RET_UNKNOWN_BUG );
}


/*
 *	c a l c D e t S c h u r
 */
real_t SQProblemSchur::calcDetSchur( int_t idxDel )
{
	if ( nS <= 0 )
		return 1.0;

	real_t newDet;
	int_t i, j;
	real_t c, s, nu;

	/* Case 1: S has been bordered by one row and column */
	if( idxDel < 0 )
	{
		/* Do a solve with the old S to check determinant of new (bordered) S */
		real_t *temp1 = new real_t[nS-1];
		real_t *temp2 = new real_t[nS-1];
		for( i=0; i<nS-1; i++ )
			temp1[i] = S[i + (nS-1)*nSmax];
		backsolveSchurQR( nS-1, temp1, 1, temp2 );

		newDet = S[(nS-1) + (nS-1)*nSmax];
		for( i=0; i<nS-1; i++ )
			newDet -= temp1[i]*temp2[i];
		newDet *= detS;
		delete [] temp1;
		delete [] temp2;
	}
	/* Case 2: row and column idxDel have been deleted from S */
	else
	{
		const int_t dim = nS+1;
		real_t *tempR = new real_t[dim*(dim-1)];
		real_t *tempColQ = new real_t[dim];

		/* Copy current R without column idxDel*/
		for( j=0; j<idxDel; j++ )
			for( i=0; i<dim; i++ )
				tempR[i+j*dim] = R_[i+j*nSmax];
		for( j=idxDel; j<dim-1; j++ )
			for( i=0; i<dim; i++ )
				tempR[i+j*dim] = R_[i+(j+1)*nSmax];
		/* Copy row idxDel of Q */
		for( j=0; j<dim; j++ )
			tempColQ[j] = Q_[idxDel+j*nSmax];

		/* Bring tempR to triangular form with nS-idxDel Givens rotations */
		for ( i=idxDel; i<nS; i++ )
		{
			computeGivens( tempR[i+i*dim], tempR[(i+1)+i*dim], tempR[i+i*dim], tempR[(i+1)+i*dim], c, s );
			nu = s/(1.0+c);
			/// \todo I think we do not need to transform all columns of R, i+3 or so should be sufficient
			for ( j=i+1; j<nS; j++ )
				applyGivens( c, s, nu, tempR[i+j*dim], tempR[(i+1)+j*dim], tempR[i+j*dim], tempR[(i+1)+j*dim] );

			/* Simultaneously transform relevant column of Q**T */
			applyGivens( c, s, nu, tempColQ[i], tempColQ[i+1], tempColQ[i], tempColQ[i+1] );
		}

		/* Delete row: nS Givens rotations to transform last column (and row!) of (old) Q**T to getAbs((nS+1)-th unity vector) */
		for ( i=nS; i>0; i-- )
		{
			computeGivens( tempColQ[nS], tempColQ[i-1], tempColQ[nS], tempColQ[i-1], c, s );
			nu = s/(1.0+c);

			/* Simultaneously transform diagonal elements of R (coldim is already one less than Q) */
			applyGivens( c, s, nu, tempR[nS+(i-1)*dim], tempR[(i-1)+(i-1)*dim], tempR[nS+(i-1)*dim], tempR[(i-1)+(i-1)*dim] );
		}

		/* Note that we implicitly did a row permutation of Q.
		 * If we did an  odd permutation AND deleted a positive unity vector or
		 * if we did an even permutation AND deleted a negative unity vector, then det(Q)=-1
		 * ->Change signs of first column of Q and first row of R */
		if ( (( (nS - idxDel) % 2 == 1 ) && ( tempColQ[nS] > 0.0 )) ||
			(( (nS - idxDel) % 2 == 0 ) && ( tempColQ[nS] < 0.0 )) )
		{
			tempR[0] = -tempR[0];
		}

		newDet = 1.0;
		//for( i=0; i<nS; i++ )
			//newDet *= tempR[i+i*dim];
		for( i=0; i<nS; i++ )
			if( tempR[i+i*dim] < 0.0 ) newDet = -newDet;
		delete [] tempR;
		delete [] tempColQ;
	}

	return newDet;
}


/*
 *	u p d a t e S c h u r Q R
 */
returnValue SQProblemSchur::updateSchurQR( int_t idxDel )
{
	int_t i, j;
	real_t c, s, nu;

	if ( nS <= 0 )
	{
		detS = 1.0;
		rcondS = 1.0;
		return SUCCESSFUL_RETURN;
	}

	/* Case 1: S has been bordered by one row and column */
	if ( idxDel < 0 )
	{
		/* I: Augment Q**T by nS-th unity vector (row and column) */
		for ( i=0; i<nS; i++ )
		{
			Q_[i+(nS-1)*nSmax] = 0.0;
			Q_[(nS-1)+i*nSmax] = 0.0;
		}
		Q_[(nS-1)+(nS-1)*nSmax] = 1.0;

		/* IIa: Augment rows of R by last row of S */
		for ( i=0; i<nS; i++ )
			R_[(nS-1)+i*nSmax] = S[(nS-1)+i*nSmax];

		/* IIb: Augment columns of R by Q**T * S[nS,:] */
		for ( i=0; i<nS; i++ )
		{
			R_[i+(nS-1)*nSmax] = 0.0;
			for ( j=0; j<nS; j++ )
				R_[i+(nS-1)*nSmax] += Q_[j+i*nSmax] * S[j+(nS-1)*nSmax];
		}

		/* III: Restore triangular form of R by nS-1 Givens rotations */
		for ( i=0; i<nS-1; i++ )
		{
			computeGivens( R_[i+i*nSmax], R_[(nS-1)+i*nSmax], R_[i+i*nSmax], R_[(nS-1)+i*nSmax], c, s );
			nu = s/(1.0+c);
			for ( j=i+1; j<nS; j++ )
				applyGivens( c, s, nu, R_[i+j*nSmax], R_[(nS-1)+j*nSmax], R_[i+j*nSmax], R_[(nS-1)+j*nSmax] );

			/* Simultaneously transform Q**T */
			for ( j=0; j<nS; j++ )
				applyGivens( c, s, nu, Q_[j+i*nSmax], Q_[j+(nS-1)*nSmax], Q_[j+i*nSmax], Q_[j+(nS-1)*nSmax] );
		}
	}
	/* Case 2: row and column idxDel have been deleted from S */
	else
	{
		/* I: Delete column idxDel of R */
		for ( j=idxDel; j<nS; j++ )
			for ( i=0; i<nS+1; i++ )
				R_[i+j*nSmax] = R_[i+(j+1)*nSmax];

		/* II: Bring R back to triangular form with nS-idxDel Givens rotations */
		for ( i=idxDel; i<nS; i++ )
		{
			computeGivens( R_[i+i*nSmax], R_[(i+1)+i*nSmax], R_[i+i*nSmax], R_[(i+1)+i*nSmax], c, s );
			nu = s/(1.0+c);
			for ( j=i+1; j<nS; j++ )
				applyGivens( c, s, nu, R_[i+j*nSmax], R_[(i+1)+j*nSmax], R_[i+j*nSmax], R_[(i+1)+j*nSmax] );

			/* Simultaneously transform (old) Q**T (coldim is one larger)*/
			for ( j=0; j<nS+1; j++ )
				applyGivens( c, s, nu, Q_[j+i*nSmax], Q_[j+(i+1)*nSmax], Q_[j+i*nSmax], Q_[j+(i+1)*nSmax] );
		}

		/* III: Permute rows of Q: move row idxDel to position nS */
		real_t temp;
		for ( j=0; j<nS+1; j++ )
		{
			temp = Q_[idxDel+j*nSmax];
			for ( i=idxDel; i<nS; i++ )
				Q_[i+j*nSmax] = Q_[(i+1)+j*nSmax];
			Q_[nS+j*nSmax] = temp;
		}

		/* IV: Delete row: nS Givens rotations to transform last column (and row!) of (old) Q**T to getAbs((nS+1)-th unity vector) */
		for ( i=nS; i>0; i-- )
		{
			computeGivens( Q_[nS+nS*nSmax], Q_[nS+(i-1)*nSmax], Q_[nS+nS*nSmax], Q_[nS+(i-1)*nSmax], c, s );
			nu = s/(1.0+c);
			for ( j=0; j<nS; j++ )
				applyGivens( c, s, nu, Q_[j+nS*nSmax], Q_[j+(i-1)*nSmax], Q_[j+nS*nSmax], Q_[j+(i-1)*nSmax] );

			/* Simultaneously transform R (coldim is already one less than Q) */
			for ( j=i-1; j<nS; j++ )
				applyGivens( c, s, nu, R_[nS+j*nSmax], R_[(i-1)+j*nSmax], R_[nS+j*nSmax], R_[(i-1)+j*nSmax] );
		}

		/* If we did an  odd permutation AND deleted a positive unity vector or
		 * if we did an even permutation AND deleted a negative unity vector, then det(Q)=-1
		 * ->Change signs of first column of Q and first row of R s.t. we always maintain det(Q)=1 */
		if ( (( (nS - idxDel) % 2 == 1 ) && ( Q_[nS+nS*nSmax] > 0.0 )) ||
			(( (nS - idxDel) % 2 == 0 ) && ( Q_[nS+nS*nSmax] < 0.0 )) )
		{
			for ( i=0; i<nS+1; i++ )
				Q_[i] = -Q_[i];
			for ( i=0; i<nS; i++ )
				R_[i*nSmax] = -R_[i*nSmax];
		}
	}

	/* Compute determinant */
	detS = 1.0;
	//for ( i=0; i<nS; i++ )
		//detS *= R_[i+i*nSmax];
	for ( i=0; i<nS; i++ )
		if( R_[i+i*nSmax] < 0.0 ) detS = -detS;

	/* Estimate condition number of R (= condition number of S)*/
	real_t *WORK;
	unsigned long N = (unsigned long)nS;
	unsigned long LDA = (unsigned long)nSmax;
	unsigned long *IWORK;
	long INFO = 0;
	IWORK = new unsigned long[N];
	WORK = new real_t[3*N];
	TRCON( "1", "U", "N", &N, R_, &LDA, &rcondS, WORK, IWORK, &INFO );
	if ( INFO != 0 )
	{
		MyPrintf( "TRCON returns INFO = %d\n",(int)INFO );
	}

	if ( options.printLevel == PL_HIGH )
		MyPrintf( "1/cond(S) = %23.16e.\n", rcondS );

	delete[] IWORK;
	delete[] WORK;

	return SUCCESSFUL_RETURN;
}


/*
 *	b a c k s o l v e S c h u r Q R
 */
returnValue SQProblemSchur::backsolveSchurQR( int_t dimS, const real_t* const rhs, int_t dimRhs, real_t* const sol )
{
	if( dimS < 1 || dimRhs < 1 )
		return SUCCESSFUL_RETURN;

	if( dimRhs > 1 )
	{
		MyPrintf("backsolve not implemented for dimRhs = %d\n", dimRhs);
		return RET_QR_FACTORISATION_FAILED;
	}

	int_t i, j;
	long INFO = 0;
	unsigned long NRHS = 1;
	unsigned long M = (unsigned long)dimS;
	unsigned long LDA = (unsigned long)nSmax;
	unsigned long LDC = (unsigned long)dimS;

	for( i=0; i<dimS; i++ )
		sol[i] = 0.0;

	/* Compute sol = Q**T * rhs */
	for( i=0; i<dimS; i++ )
		for( j=0; j<dimS; j++ )
			sol[i] += Q_[j+i*nSmax] * rhs[j];

	/* Solve Rx = sol */
	TRTRS( "U", "N", "N", &M, &NRHS, R_, &LDA, sol, &LDC, &INFO );
	if ( INFO != 0 )
	{
		MyPrintf("TRTRS returns INFO = %d\n", INFO);
		if ( INFO == ((long)0xDEADBEEF) )
			MyPrintf( "If SQProblemSchur is to be used, system LAPACK must be used instead of the qpOASES LAPACK replacement" );
		return RET_QR_FACTORISATION_FAILED;
	}

	return SUCCESSFUL_RETURN;
}


returnValue SQProblemSchur::stepCalcRhs(	int_t nFR, int_t nFX, int_t nAC, int_t* FR_idx, int_t* FX_idx, int_t* AC_idx, real_t& rhs_max, 
											const real_t* const delta_g, const real_t* const delta_lbA, const real_t* const delta_ubA,
											const real_t* const delta_lb, const real_t* const delta_ub,
											BooleanType Delta_bC_isZero, BooleanType Delta_bB_isZero,
											real_t* const delta_xFX, real_t* const delta_xFR,
											real_t* const delta_yAC, real_t* const delta_yFX
											)
{
	int_t i, ii;
	returnValue retval;

	if ( nS < 0 )
	{
		retval = resetSchurComplement( BT_FALSE );
		if (retval != SUCCESSFUL_RETURN)
		{
			MyPrintf( "In SQProblemSchur::stepCalcRhs, resetSchurComplement returns %d\n", retval);
			return THROWERROR( retval );
		}
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
	if ( ( hessianType != HST_IDENTITY ) && ( hessianType != HST_ZERO ) )
	{
		/* tempA becomes RHS for reduced augmented system, gFR+H_FX*delta_xFR */
		H->times(bounds.getFree(), bounds.getFixed(), 1, 1.0, delta_xFX, nFX, 1.0, tempA, nFR);
	}
	/* tempB becomes RHS for reduced augmented system, bA-A_CX*delta_xFR */
	A->times(constraints.getActive(), bounds.getFixed(), 1, -1.0, delta_xFX, nFX, 1.0, tempB, nAC);

	/* If iterative refinement is requested, compute max-norm of RHS for termination test. */
	rhs_max = 0.0;
	if ( options.numRefinementSteps > 0 )
	{
		for ( i=0; i<nFR; i++ )
			rhs_max = getMax(rhs_max, getAbs(tempA[i]));
		for ( i=0; i<nAC; i++ )
			rhs_max = getMax(rhs_max, getAbs(tempB[i]));
	}
	return SUCCESSFUL_RETURN;
}

returnValue SQProblemSchur::stepCalcReorder(int_t nFR, int_t nAC, int_t* FR_idx, int_t* AC_idx, int_t nFRStart, int_t nACStart, int_t* FR_idxStart, int_t* AC_idxStart, int_t* FR_iSort, int_t* FR_iSortStart, int_t* AC_iSort, int_t* AC_iSortStart, real_t* rhs)
{
	int_t i, ii;
	/* Reorder information for the new to the old free variables. */
	i = 0;
	ii = 0;
	while ( ii < nFRStart )
	{
		if ( i == nFR )
			rhs[FR_iSortStart[ii++]] = 0.0;
		else
		{
			int_t idx = FR_idx[FR_iSort[i]];
			int_t idxStart = FR_idxStart[FR_iSortStart[ii]];

			if ( idx == idxStart )
				rhs[FR_iSortStart[ii++]] = -tempA[FR_iSort[i++]];
			else if ( idx < idxStart )
				i++;
			else
				rhs[FR_iSortStart[ii++]] = 0.0;
		}
	}
	/* Reorder information for the new to the old active constraints. */
	i = 0;
	ii = 0;
	while ( ii < nACStart )
	{
		if ( i == nAC )
			rhs[nFRStart+AC_iSortStart[ii++]] = 0.0;
		else
		{
			int_t idx = AC_idx[AC_iSort[i]];
			int_t idxStart = AC_idxStart[AC_iSortStart[ii]];

			if ( idx == idxStart )
				rhs[nFRStart+AC_iSortStart[ii++]] = tempB[AC_iSort[i++]];
			else if ( idx < idxStart )
				i++;
			else
				rhs[nFRStart+AC_iSortStart[ii++]] = 0.0;
		}
	}
	return SUCCESSFUL_RETURN;
}

returnValue SQProblemSchur::stepCalcBacksolveSchur( int_t nFR, int_t nFX, int_t nAC, int_t* FR_idx, int_t* FX_idx, int_t* AC_idx, int_t dim, real_t* rhs, real_t* sol )
{
	returnValue retval;
	int_t i, ii;

	real_t* q = new real_t[nS];

	/* Compute extra compoments of the RHS */
	for ( ii=0; ii<nS; ii++ )
	{
		int_t idx = schurUpdateIndex[ii];
		switch ( schurUpdate[ii] ) // TODO: All the loops below could be done faster by binary search or so
		{
			case SUT_VarFixed:
				q[ii] = 0.0;
				break;

			case SUT_VarFreed:
				/* Find index of freed variable */
				for( i=0; i<nFR; ++i )
					if ( FR_idx[i] == idx )
					{
						q[ii] = -tempA[i];
						break;
					}
				break;

			case SUT_ConAdded:
				/* Find index of added constraint */
				for( i=0; i<nAC; ++i )
					if ( AC_idx[i] == idx )
					{
						q[ii] = tempB[i];
						break;
					}
				break;

			case SUT_ConRemoved:
				q[ii] = 0.0;
				break;

			default:
				return THROWERROR( RET_UNKNOWN_BUG );
		}
	}

	/* compute q = M^T K^{-1} r - q */
	computeMTransTimes(1.0, sol, -1.0, q);

	/* Solve linear system with Schur complement. */
	real_t* p = new real_t[nS];
	backsolveSchurQR( nS, q, 1, p );

	computeMTimes(-1.0, p, 1.0, rhs);

	retval = sparseSolver->solve(dim, rhs, sol);
	if (retval != SUCCESSFUL_RETURN)
	{
		MyPrintf( "sparseSolver->solve (second time) failed.\n");
		return THROWERROR(RET_MATRIX_FACTORISATION_FAILED); // TODO: Different return code
	}

	/* Transfer extra compoments of the Schur complement solution to the correct place. */
	for ( ii=0; ii<nS; ii++ )
	{
		int_t idx = schurUpdateIndex[ii];
		switch ( schurUpdate[ii] ) // TODO: All the loops below could be done faster by binary search or so
		{
			case SUT_VarFixed:
				break;

			case SUT_VarFreed:
				/* Find index of freed variable */
				for( i=0; i<nFR; ++i )
					if ( FR_idx[i] == idx )
					{
						delta_xFR_TMP[i] = p[ii];
						break;
					}
				break;

			case SUT_ConAdded:
				/* Find index of added constraint */
				for( i=0; i<nAC; ++i )
					if ( AC_idx[i] == idx )
					{
						delta_yAC_TMP[i] = -p[ii];
						break;
					}
				break;

			case SUT_ConRemoved:
				break;

			default:
				return THROWERROR( RET_UNKNOWN_BUG );
		}
	}

	delete [] p;
	delete [] q;
	return SUCCESSFUL_RETURN;
}

returnValue SQProblemSchur::stepCalcReorder2(int_t nFR, int_t nAC, int_t* FR_idx, int_t* AC_idx, int_t nFRStart, int_t nACStart, int_t* FR_idxStart, int_t* AC_idxStart, int_t* FR_iSort, int_t* FR_iSortStart, int_t* AC_iSort, int_t* AC_iSortStart, real_t* sol, real_t* const delta_xFR, real_t* const delta_yAC)
{
	int_t i, ii;
			i = 0;
			ii = 0;
			while ( ii < nFRStart && i < nFR )
			{
				int_t idx = FR_idx[FR_iSort[i]];
				int_t idxStart = FR_idxStart[FR_iSortStart[ii]];

				if ( idx == idxStart )
					delta_xFR_TMP[FR_iSort[i++]] = sol[FR_iSortStart[ii++]];
				else if ( idx < idxStart )
					i++;
				else
					ii++;
			}
			/* Transfer Schur complement solution for the active constraint multipliers to the correct places */
			i = 0;
			ii = 0;
			while ( ii < nACStart && i < nAC )
			{
				int_t idx = AC_idx[AC_iSort[i]];
				int_t idxStart = AC_idxStart[AC_iSortStart[ii]];

				if ( idx == idxStart )
					delta_yAC_TMP[AC_iSort[i++]] = -sol[nFRStart+AC_iSortStart[ii++]];
				else if ( idx < idxStart )
					i++;
				else
					ii++;
			}

			/* refine the solution found so far */
			for ( i=0; i<nFR; ++i )
				delta_xFR[i] += delta_xFR_TMP[i];
			for ( i=0; i<nAC; ++i )
				delta_yAC[i] += delta_yAC_TMP[i];
	return SUCCESSFUL_RETURN;
}

returnValue SQProblemSchur::stepCalcResid(int_t nFR, int_t nFX, int_t nAC, int_t* FR_idx, int_t* FX_idx, int_t* AC_idx, BooleanType Delta_bC_isZero, real_t* const delta_xFX, real_t* const delta_xFR, real_t* const delta_yAC, const real_t* const delta_g, const real_t* const delta_lbA, const real_t* const delta_ubA, real_t& rnrm)
{
	int_t i, ii;
				/* compute residuals in tempA and tempB, and max-norm */
				for ( i=0; i<nFR; ++i )
				{
					ii = FR_idx[i];
					tempA[i] = delta_g[ii];
				}

				switch ( hessianType )
				{
					case HST_ZERO:
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

				for ( i=0; i<nFR; ++i )
					tempA[i] += options.epsRegularisation*delta_xFR[i];

				A->transTimes(constraints.getActive(), bounds.getFree(), 1, -1.0, delta_yAC, nAC, 1.0, tempA, nFR);
				rnrm = 0.0;
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

	return SUCCESSFUL_RETURN;
}

returnValue SQProblemSchur::stepCalcDeltayFx(int_t nFR, int_t nFX, int_t nAC, int_t* FX_idx, const real_t* const delta_g, real_t* const delta_xFX, real_t* const delta_xFR, real_t* const delta_yAC, real_t* const delta_yFX)
{
	int_t i;
		for( i=0; i<nFX; ++i )
			delta_yFX[i] = delta_g[FX_idx[i]];

		A->transTimes(constraints.getActive(), bounds.getFixed(), 1, -1.0, delta_yAC, nAC, 1.0, delta_yFX, nFX);

		if ( hessianType == HST_ZERO )
		{
		  // TODO: if ( usingRegularisation( ) == BT_TRUE )
				for( i=0; i<nFX; ++i )
					delta_yFX[i] += options.epsRegularisation*delta_xFX[i];
		}
		else if ( hessianType == HST_IDENTITY )
		{
			for( i=0; i<nFX; ++i )
				delta_yFX[i] += delta_xFX[i];
		}
		else
		{
			H->times(bounds.getFixed(), bounds.getFree(), 1, 1.0, delta_xFR, nFR, 1.0, delta_yFX, nFX);
			H->times(bounds.getFixed(), bounds.getFixed(), 1, 1.0, delta_xFX, nFX, 1.0, delta_yFX, nFX);
		}
	return SUCCESSFUL_RETURN;
}

returnValue SQProblemSchur::determineStepDirection(	const real_t* const delta_g, const real_t* const delta_lbA, const real_t* const delta_ubA,
												const real_t* const delta_lb, const real_t* const delta_ub,
												BooleanType Delta_bC_isZero, BooleanType Delta_bB_isZero,
												real_t* const delta_xFX, real_t* const delta_xFR,
												real_t* const delta_yAC, real_t* const delta_yFX
												)
{
	returnValue retval = determineStepDirection2(	delta_g, delta_lbA, delta_ubA, delta_lb, delta_ub,
													Delta_bC_isZero, Delta_bB_isZero, delta_xFX, delta_xFR,
													delta_yAC, delta_yFX
													);

	if ( retval == RET_QR_FACTORISATION_FAILED )
	{
		retval = resetSchurComplement( BT_FALSE );
		if (retval != SUCCESSFUL_RETURN)
		{
			MyPrintf( "In SQProblem::determineStepDirection, resetSchurComplement returns %d\n", retval);
			return THROWERROR( retval );
		}
		retval = determineStepDirection2(	delta_g, delta_lbA, delta_ubA, delta_lb, delta_ub,
													Delta_bC_isZero, Delta_bB_isZero, delta_xFX, delta_xFR,
													delta_yAC, delta_yFX
													);
	}
	return retval;
}

/*
 *	d e t e r m i n e S t e p D i r e c t i o n
 */
returnValue SQProblemSchur::determineStepDirection2(	const real_t* const delta_g, const real_t* const delta_lbA, const real_t* const delta_ubA,
												const real_t* const delta_lb, const real_t* const delta_ub,
												BooleanType Delta_bC_isZero, BooleanType Delta_bB_isZero,
												real_t* const delta_xFX, real_t* const delta_xFR,
												real_t* const delta_yAC, real_t* const delta_yFX
												)
{
	/* The linear system to be solved here is this:

	   / H_FF  H_FX  A_CF^T  0 \ /  delta_xFR \   / -delta_g_F \
	   | H_XF  H_XX  A_CX^T  I | |  delta_xFX |   | -delta_g_X |
	   | A_CF  A_CX    0     0 | | -delta_yAC | = |  delta_bA  |  <-- active entries of delta_lbA and delta_ubA with corresponding sign
	   \  0     I      0     0 / \ -delta_yFX /   \  delta_bX  /  <-- fixed entries of delta_lb and delta_ub with corresponding sign

	*/


	int_t i, ii, r;

	returnValue retval;

  //int_t nV  = getNV( );
	int_t nFR = getNFR( );
	int_t nFX = getNFX( );
	int_t nAC = getNAC( );

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

	if ( nFR+nAC>0 ) {
		real_t rhs_max = 0.0;
		retval = stepCalcRhs( nFR, nFX, nAC, FR_idx, FX_idx, AC_idx, rhs_max, delta_g, delta_lbA, delta_ubA,
							  delta_lb, delta_ub, Delta_bC_isZero, Delta_bB_isZero, delta_xFX, delta_xFR,
							  delta_yAC, delta_yFX );

		if (retval != SUCCESSFUL_RETURN)
			return retval;
		int_t nFRStart = boundsFreeStart.getLength();
		int_t nACStart = constraintsActiveStart.getLength();

		int_t* FR_iSort;
		int_t* AC_iSort;
		bounds.getFree( )->getISortArray( &FR_iSort );
		constraints.getActive( )->getISortArray( &AC_iSort );

		int_t* FR_idxStart;
		int_t* AC_idxStart;
		boundsFreeStart.getNumberArray( &FR_idxStart );
		constraintsActiveStart.getNumberArray( &AC_idxStart );

		int_t* FR_iSortStart;
		int_t* AC_iSortStart;
		boundsFreeStart.getISortArray( &FR_iSortStart );
		constraintsActiveStart.getISortArray( &AC_iSortStart );

		int_t dim = nFRStart + nACStart;
		real_t* rhs = new real_t[dim];
		real_t* sol = new real_t[dim];

		/* Iterative refinement loop for delta_xFR, delta_yAC */
		for ( r=0; r<=options.numRefinementSteps; ++r )
		{
		  retval = stepCalcReorder(nFR, nAC, FR_idx, AC_idx, nFRStart, nACStart, FR_idxStart, AC_idxStart, FR_iSort, FR_iSortStart, AC_iSort, AC_iSortStart, rhs);
			if (retval != SUCCESSFUL_RETURN)
				return retval;

			retval = sparseSolver->solve(dim, rhs, sol);

			if (retval != SUCCESSFUL_RETURN)
			{
				MyPrintf( "sparseSolver->solve (first time) failed.\n");
				return THROWERROR(RET_MATRIX_FACTORISATION_FAILED); // TODO: Different return code
			}

			if ( nS > 0 )
			{
				retval = stepCalcBacksolveSchur( nFR, nFX, nAC, FR_idx, FX_idx, AC_idx, dim, rhs, sol );
				if (retval != SUCCESSFUL_RETURN)
					return retval;
			}

			/* Transfer Schur complement solution for the free variables to the correct places */
			retval = stepCalcReorder2(nFR, nAC, FR_idx, AC_idx, nFRStart, nACStart, FR_idxStart, AC_idxStart, FR_iSort, FR_iSortStart, AC_iSort, AC_iSortStart, sol, delta_xFR, delta_yAC);
			if (retval != SUCCESSFUL_RETURN)
				return retval;

			if ( r < options.numRefinementSteps ) // TODO: use "<" to avoid computation in last round
			{
				real_t rnrm;
				retval = stepCalcResid(nFR, nFX, nAC, FR_idx, FX_idx, AC_idx, Delta_bC_isZero, delta_xFX, delta_xFR, delta_yAC, delta_g, delta_lbA, delta_ubA, rnrm);
				if (retval != SUCCESSFUL_RETURN)
					return retval;

				/* early termination of residual norm small enough */
				if ( options.printLevel == PL_HIGH )
					MyPrintf( "In iterative refinement (iter %d) rnrm = %e and epsIterRef*rhs_max = %e.\n", r, rnrm, options.epsIterRef*rhs_max);

				if ( rnrm <= options.epsIterRef*rhs_max )
					break;
			}

		}

		delete [] sol;
		delete [] rhs;
	}

	/* IV) DETERMINE delta_yFX */
	if ( nFX > 0 )
	{
		retval = stepCalcDeltayFx(nFR, nFX, nAC, FX_idx, delta_g, delta_xFX, delta_xFR, delta_yAC, delta_yFX);
		if (retval != SUCCESSFUL_RETURN)
			return retval;
	}

	return SUCCESSFUL_RETURN;
}

returnValue SQProblemSchur::resetSchurComplement( BooleanType allowInertiaCorrection )
{
	int_t j;
	int_t nFR = getNFR( );
	int_t nAC = getNAC( );

	if ( options.printLevel == PL_HIGH )
		MyPrintf( "Resetting Schur complement.\n");

	nS = 0;
	detS = 1.0;
	rcondS = 1.0;
	boundsFreeStart = *bounds.getFree();
	constraintsActiveStart = *constraints.getActive();

	if ( nSmax > 0 )
		M_jc[0] = 0;

	int_t dim = nFR+nAC;
	// Count the number of nonzeros
	int_t numNonzeros;
	switch ( hessianType )
	{
		case HST_ZERO:
			numNonzeros = 0;
			break;

		case HST_IDENTITY:
			numNonzeros = nFR;
			break;

		default:
			H->getSparseSubmatrix( bounds.getFree(), bounds.getFree(), 1, 1, numNonzeros, 0, 0, 0, BT_TRUE);
			break;
	}
	// TODO: For now, we regularize every time
	if (options.epsRegularisation > 0.0)
		numNonzeros += nFR;

	int_t numNonzerosA;

	if ( constraintProduct != 0 )
	{
		MyPrintf( "In SQProblemSchur::determineStepDirection, constraintProduct not yet implemented.\n");
		return THROWERROR(RET_NOT_YET_IMPLEMENTED);
	}
	A->getSparseSubmatrix( constraints.getActive(), bounds.getFree(), nFR+1, 1, numNonzerosA, 0, 0, 0, BT_FALSE);
	numNonzeros += numNonzerosA;

	// Get the values
	real_t* avals = new real_t[numNonzeros];
	int_t* irn = new int_t[numNonzeros];
	int_t* jcn = new int_t[numNonzeros];
	numNonzeros = 0;
	switch ( hessianType )
	{
		case HST_ZERO:
			break;

		case HST_IDENTITY:
			numNonzeros += nFR;
			for (j = 0; j<nFR; j++)
			{
				irn[j] = j+1;
				jcn[j] = j+1;
				avals[j] = 1.0;
			}
			break;

		default:
			H->getSparseSubmatrix( bounds.getFree(), bounds.getFree(), 1, 1, numNonzeros, irn, jcn, avals, BT_TRUE);
			break;
	}

	// For now, we regularize every time
	if (options.epsRegularisation > 0.0)
	{
		for (j = 0; j<nFR; j++)
		{
			irn[numNonzeros] = j+1;
			jcn[numNonzeros] = j+1;
			avals[numNonzeros++] = options.epsRegularisation;
		}
	}

	A->getSparseSubmatrix( constraints.getActive(), bounds.getFree(), nFR+1, 1, numNonzerosA, irn+numNonzeros, jcn+numNonzeros, avals+numNonzeros, BT_FALSE);
	numNonzeros += numNonzerosA;

	// Call the linear solver
	sparseSolver->reset();
	returnValue retval = sparseSolver->setMatrixData(dim, numNonzeros, irn, jcn, avals);
	delete [] jcn;
	delete [] irn;
	delete [] avals;
	
	if (retval != SUCCESSFUL_RETURN)
		return THROWERROR(RET_NO_SPARSE_SOLVER);

	// Factorize the matrix for later backsolves
	retval = sparseSolver->factorize();
	numFactorizations++;

	// If matrix is singular, add bounds/remove constraints according to zero pivots
	if (retval == RET_KKT_MATRIX_SINGULAR)
	{
		if( repairSingularWorkingSet( ) == SUCCESSFUL_RETURN )
			return resetSchurComplement( allowInertiaCorrection );
		else
			return RET_KKT_MATRIX_SINGULAR;
	}

	// If matrix has wrong inertia, add bounds until inertia is correct
	if (retval == SUCCESSFUL_RETURN && allowInertiaCorrection)
	{
		int_t neig = sparseSolver->getNegativeEigenvalues( );
		if( neig > getNAC( ) )
		{
			if ( options.printLevel == PL_HIGH )
				MyPrintf( "WARNING: After new factorization, reduced Hessian has %i negative eigenvalues, should be %i.\n", neig, getNAC( ) );

			retval = correctInertia();
		}
	}

	if (retval != SUCCESSFUL_RETURN)
		return THROWERROR(RET_MATRIX_FACTORISATION_FAILED);

	nS = 0;

	return SUCCESSFUL_RETURN;
}

returnValue SQProblemSchur::computeMTimes( real_t alpha, const real_t* const x_, real_t beta, real_t* const y_ )
{
	if ( isEqual( alpha, -1.0 ) == BT_FALSE || isEqual( beta, 1.0 ) == BT_FALSE )
		return THROWERROR(RET_NOT_YET_IMPLEMENTED);

	int_t i, j;

	for ( j=0; j<nS; j++ )
	{
		const real_t xval = x_[j];
		for ( i=M_jc[j]; i<M_jc[j+1]; i++)
			y_[M_ir[i]] -= M_vals[i]*xval;
	}

	return SUCCESSFUL_RETURN;
}

returnValue SQProblemSchur::computeMTransTimes( real_t alpha, const real_t* const x_, real_t beta, real_t* const y_ )
{
	if ( isEqual( alpha, 1.0 ) == BT_FALSE || ( isZero( beta ) == BT_FALSE && isEqual( beta, -1.0 ) == BT_FALSE ) )
		return THROWERROR(RET_NOT_YET_IMPLEMENTED);

	int_t i, j;

	if ( isZero( beta ) == BT_TRUE )
	{
		for ( j=0; j<nS; j++ )
		{
			y_[j] = 0.0;
			for ( i=M_jc[j]; i<M_jc[j+1]; i++)
				y_[j] += M_vals[i]*x_[M_ir[i]];
		}
	}
	else
	{
		for ( j=0; j<nS; j++ )
		{
			y_[j] = -y_[j];
			for ( i=M_jc[j]; i<M_jc[j+1]; i++)
				y_[j] += M_vals[i]*x_[M_ir[i]];
		}
	}

	return SUCCESSFUL_RETURN;
}

returnValue SQProblemSchur::addToSchurComplement( int_t number, SchurUpdateType update, int_t numNonzerosM, const sparse_int_t* Mpos, const real_t* const Mvals, int_t numNonzerosN, const sparse_int_t* Npos, const real_t* const Nvals, real_t N_diag )
{
	int_t i;

	int_t nFRStart = boundsFreeStart.getLength();
	int_t nACStart = constraintsActiveStart.getLength();

	real_t* new_Scol = new real_t[nS];

	int_t dim = nFRStart + nACStart;
	real_t* rhs = new real_t[dim];
	real_t* sol = new real_t[dim];

	for ( i=0; i<dim; i++ )
		rhs[i] = 0.0;

	for ( i=0; i<numNonzerosM; i++ )
		rhs[Mpos[i]] = Mvals[i];

	returnValue retval = sparseSolver->solve(dim, rhs, sol);
	if (retval != SUCCESSFUL_RETURN)
	{
		MyPrintf( "sparseSolver->solve in SQProblemSchur::addToSchurComplement failed.\n");
		return THROWERROR(RET_MATRIX_FACTORISATION_FAILED); // TODO: Different return code
	}

	computeMTransTimes(1.0, sol, 0.0, new_Scol);

	/* Take care of off-diagonal elements in N. */
	for ( i=0; i<numNonzerosN; i++ )
		new_Scol[Npos[i]] -= Nvals[i];

	real_t sdiag = -N_diag;
	for ( i=0; i<numNonzerosM; i++ )
		sdiag += Mvals[i] * sol[Mpos[i]];

	/* Now augment S */
	for ( i=0; i<nS; i++)
		S[nS*nSmax + i] = new_Scol[i];
	for ( i=0; i<nS; i++)
		S[i*nSmax + nS] = new_Scol[i];
	S[nS*nSmax + nS] = sdiag;

	schurUpdateIndex[nS] = number;
	schurUpdate[nS] = update;

	/* Augment M matrix.  */
	if ( M_physicallength < M_jc[nS] + numNonzerosM )
	{
		/* If necessary, allocate more memory for M. */
		int_t M_physicallength_new = getMax(2*M_physicallength, M_physicallength + 2*numNonzerosM);
		real_t* M_vals_new = new real_t[M_physicallength_new];
		sparse_int_t* M_ir_new = new sparse_int_t[M_physicallength_new];
		memcpy( M_vals_new, M_vals, ((unsigned int)(M_jc[nS]))*sizeof(real_t) );
		memcpy( M_ir_new, M_ir, ((unsigned int)(M_jc[nS]))*sizeof(sparse_int_t) );
		M_physicallength = M_physicallength_new;
		delete [] M_vals;
		delete [] M_ir;
		M_vals = M_vals_new;
		M_ir = M_ir_new;
	}

	for ( i=0; i<numNonzerosM; i++ )
	{
		M_vals[M_jc[nS] + i] = Mvals[i];
		M_ir[M_jc[nS] + i] = Mpos[i];
	}
	M_jc[nS+1] = M_jc[nS] + numNonzerosM;

	nS++;

	delete [] sol;
	delete [] rhs;
	delete [] new_Scol;

	if ( options.printLevel == PL_HIGH )
		MyPrintf( "added index %d with update type %d to Schur complement.  nS = %d\n", number, update, nS);

	return SUCCESSFUL_RETURN;
}


returnValue SQProblemSchur::deleteFromSchurComplement( int_t idx, BooleanType allowUndo )
{
	if ( options.printLevel == PL_HIGH )
		MyPrintf( "deleting entry %d with idx = %d and type %d from Schur complement.", idx, schurUpdateIndex[idx], schurUpdate[idx]);

	if ( idx != nS-1 )
	{
		real_t *temp_vals = NULL;
		int_t *temp_ir = NULL;
		int_t schurUpdateIndexTemp = -1;
		SchurUpdateType schurUpdateTemp = SUT_UNDEFINED;

		/* temporarily save the column of S to be deleted */
		if( allowUndo == BT_TRUE )
		{
			temp_vals = new real_t[nS];
			for ( int_t i=0; i<nS; i++ )
				temp_vals[i] = S[idx*nSmax + i];

			schurUpdateIndexTemp = schurUpdateIndex[idx];
			schurUpdateTemp = schurUpdate[idx];
		}

		/* Shift rows and columns >idx of S by one to the upper left */
		for ( int_t i=0; i<idx; i++ )
			for ( int_t j=idx+1; j<nS; j++ )
				S[i*nSmax + j-1] = S[i*nSmax + j];
		for ( int_t i=idx+1; i<nS; i++ )
		{
			for ( int_t j=0; j<idx; j++ )
				S[(i-1)*nSmax + j] = S[i*nSmax + j];
			for ( int_t j=idx+1; j<nS; j++ )
				S[(i-1)*nSmax + j-1] = S[i*nSmax + j];
		}
		for ( int_t i=idx+1; i<nS; i++ )
		{
			schurUpdateIndex[i-1] = schurUpdateIndex[i];
			schurUpdate[i-1] = schurUpdate[i];
		}

		/* Store deleted row/column in the last row/column of S, can retrieve it from there later */
		if( allowUndo == BT_TRUE )
		{
			for ( int_t i=0; i<nS; i++ )
			{
				S[(nS-1)*nSmax + i] = temp_vals[i];
				S[i*nSmax + (nS-1)] = temp_vals[i];
			}
			schurUpdateIndex[nS-1] = schurUpdateIndexTemp;
			schurUpdate[nS-1] = schurUpdateTemp;
			delete[] temp_vals;
		}

		/* temporarily save the (sparse) column of M to be deleted */
		int_t numEntries = M_jc[idx+1] - M_jc[idx];
		if( allowUndo == BT_TRUE )
		{
			temp_ir = new int_t[numEntries];
			temp_vals = new real_t[numEntries];

			for ( int_t i=M_jc[idx]; i<M_jc[idx+1]; i++ )
			{
				temp_ir[i-M_jc[idx]] = M_ir[i];
				temp_vals[i-M_jc[idx]] = M_vals[i];
			}
		}

		/* Shift all columns >idx one to the left */
		for ( int_t i=M_jc[idx+1]; i<M_jc[nS]; i++ )
		{
			M_ir[i-numEntries] = M_ir[i];
			M_vals[i-numEntries] = M_vals[i];
		}
		for ( int_t i=idx; i<nS; i++ )
			M_jc[i] = M_jc[i+1] - numEntries;

		/* Store deleted column of M in the last column, can retrieve it from there later */
		if( allowUndo == BT_TRUE )
		{
			for ( int_t i=M_jc[nS-1]; i<M_jc[nS]; i++ )
			{
				M_ir[i] = temp_ir[i-M_jc[nS-1]];
				M_vals[i] = temp_vals[i-M_jc[nS-1]];
			}

			delete[] temp_ir;
			delete[] temp_vals;
		}
	}

	nS--;

	if ( options.printLevel == PL_HIGH )
		MyPrintf( "  nS = %d\n", nS);

	return SUCCESSFUL_RETURN;
}


returnValue SQProblemSchur::undoDeleteFromSchurComplement( int_t idx )
{
	if ( options.printLevel == PL_HIGH )
		MyPrintf( "undo deletion of entry %d with idx = %d and type %d from Schur complement. nS = %i\n", idx, schurUpdateIndex[nS-1], schurUpdate[nS-1],nS+1);

	if ( idx != nS )
	{
		real_t *temp_vals;
		int_t *temp_ir;
		int_t schurUpdateIndexTemp = -1;
		SchurUpdateType schurUpdateTemp = SUT_UNDEFINED;

		/* temporarily save the last column of S */
		temp_vals = new real_t[nS+1];
		for ( int_t i=0; i<nS+1; i++ )
			temp_vals[i] = S[i+nS*nSmax];

		schurUpdateIndexTemp = schurUpdateIndex[nS];
		schurUpdateTemp = schurUpdate[nS];

		/* Shift rows and columns =>idx of S by one to the lower right */
		for ( int_t i=idx-1; i>-1; i-- )
			for ( int_t j=nS-1; j>idx-1; j-- )
				S[(j+1)+i*nSmax] = S[j+i*nSmax];
		for ( int_t i=nS-1; i>idx-1; i-- )
		{
			for ( int_t j=idx-1; j>-1; j-- )
				S[j+(i+1)*nSmax] = S[j+i*nSmax];
			for ( int_t j=nS-1; j>idx-1; j-- )
				S[(j+1)+(i+1)*nSmax] = S[j+i*nSmax];
		}
		for ( int_t i=nS-1; i>idx-1; i-- )
		{
			schurUpdateIndex[i+1] = schurUpdateIndex[i];
			schurUpdate[i+1] = schurUpdate[i];
		}

		/* Insert stored row/column of S at position idx */
		for ( int_t i=0; i<nS+1; i++ )
		{
			S[idx*nSmax + i] = temp_vals[i];
			S[i*nSmax + idx] = temp_vals[i];
		}
		schurUpdateIndex[idx] = schurUpdateIndexTemp;
		schurUpdate[idx] = schurUpdateTemp;
		delete[] temp_vals;

		/* temporarily save the last (sparse) column of M */
		int_t numEntries = M_jc[nS+1] - M_jc[nS];
		temp_ir = new int_t[numEntries];
		temp_vals = new real_t[numEntries];
		for ( int_t i=M_jc[nS]; i<M_jc[nS+1]; i++ )
		{
			temp_ir[i-M_jc[nS]] = M_ir[i];
			temp_vals[i-M_jc[nS]] = M_vals[i];
		}

		/* Shift all columns =>idx one to the right */
		for ( int_t i=M_jc[nS]-1; i>M_jc[idx]-1; i-- )
		{
			M_ir[i+numEntries] = M_ir[i];
			M_vals[i+numEntries] = M_vals[i];
		}
		for ( int_t i=nS; i>idx-1; i-- )
			M_jc[i+1] = M_jc[i] + numEntries;

		/* Insert stored column of M at position idx */
		for ( int_t i=M_jc[idx]; i<M_jc[idx+1]; i++ )
		{
			M_ir[i] = temp_ir[i-M_jc[idx]];
			M_vals[i] = temp_vals[i-M_jc[idx]];
		}

		delete[] temp_ir;
		delete[] temp_vals;
	}

	nS++;

	if ( options.printLevel == PL_HIGH )
		MyPrintf( "  nS = %d\n", nS);

	return SUCCESSFUL_RETURN;
}

returnValue SQProblemSchur::correctInertia( )
{
	SubjectToStatus B_status;
	real_t oldDetS;
	int_t nFR = getNFR( );
	int_t k, number, neig, nAdded;
	int_t *freeBoundIdx = new int_t[nFR];
	int_t *numberarray;

	/* method may only be called after refactorization or if one bound/constraint
	 * has been added, i.e. when a bound has flipped after refactorization */
	if( nS != 0 && nS != 1 )
		return THROWERROR( RET_INERTIA_CORRECTION_FAILED );
	neig = sparseSolver->getNegativeEigenvalues( );

	/* if a bound flipped, check if it did in fact remove a negative eigenvalue */
	if( nS == 1 && detS < 0 )
		neig--;

	/* if this method is triggered after flipping bounds, inertia is now probably correct */
	if( neig == getNAC( ) )
		return SUCCESSFUL_RETURN;

	/* get bound numbers in the order in which they are in the non-basis */
	bounds.getFree()->getNumberArray( &numberarray );
	for( k=0; k<nFR; k++ )
		freeBoundIdx[k] = numberarray[k];

	k = 0;
	nAdded = getNFR( );
	while ( neig > getNAC( ) && k < nFR )
	{
		oldDetS = detS;

		/* If it's linearly independent, fix the next free variable at the nearest bound */
		number = freeBoundIdx[k];
		if( addBound_checkLI( number ) == RET_LINEARLY_INDEPENDENT )
		{
			/* This is just heuristics: we need the bound which gives correct multiplier sign */
			if ( x[number] - lb[number] < ub[number] - x[number] )
				B_status = ST_LOWER;
			else
				B_status = ST_UPPER;

			/* Update Schur complement */
			if( addBound( number, B_status, BT_TRUE, BT_FALSE ) != SUCCESSFUL_RETURN )
			{
				if ( options.printLevel == PL_HIGH )
					MyPrintf("In correctInertia: Adding bound[%i] = %i failed!\n", k, number );
				return THROWERROR( RET_INERTIA_CORRECTION_FAILED );
			}

			/* Adjust bounds */
			if ( B_status == ST_LOWER )
				lb[number] = x[number];
			else
				ub[number] = x[number];
		}
		else
		{
			if ( options.printLevel == PL_HIGH )
				MyPrintf("bound[%i] = %i is linearly dependent. Do not add.\n", k, number );
			k++;
			continue;
		}

		/* Case 1: Schur complement has been reset, check inertia of new factorization */
		if( nS == 0 )
			neig = sparseSolver->getNegativeEigenvalues( );
		/* Case 2: Schur complement has grown, check if determinant changed sign */
		else if( oldDetS * detS < 0 )
			neig--;
		/* NB: Case 3: (Schur complement has shrunk) cannot happen here:
		 * This method is called after a factorization reset or after ONE bound has been added */

		k++;
	}
	nAdded -= getNFR( );

	delete[] freeBoundIdx;

	/* if there are still too many negative eigenvalues, exit */
	if( neig > getNAC( ) )
	{
		if ( options.printLevel == PL_HIGH )
			MyPrintf( "Added %i bounds but KKT matrix still has %i negative eigenvalues, should be %i.\n", nAdded, neig, getNAC( ) );
		return THROWERROR( RET_INERTIA_CORRECTION_FAILED );
	}
	else
	{
		if ( options.printLevel == PL_HIGH )
			MyPrintf( "After adding %i bounds, reduced Hessian has correct inertia.\n", nAdded, neig );
		return SUCCESSFUL_RETURN;
	}
}


returnValue SQProblemSchur::repairSingularWorkingSet( )
{
	int_t k, number;
	SubjectToStatus B_status;
	int_t rank = sparseSolver->getRank( );
	int_t nFR = getNFR( );
	int_t defect = nFR + getNAC( ) - rank;

	/* Rank detection not supported by linear solver */
	if ( rank < 0 )
		return RET_KKT_MATRIX_SINGULAR;

	/* Consistency check */
	if ( defect <= 0 )
		return RET_UNKNOWN_BUG;

	/* Determine zero pivots */
	int_t *zeroPivots = new int_t[defect];
	sparseSolver->getZeroPivots( zeroPivots );

	/* Determination of zero pivots not supported by linear solver */
	if ( zeroPivots == 0 )
		return RET_KKT_MATRIX_SINGULAR;

	/* We assume implicitly that pivots are sorted in ascending order */
	/// \todo make sure that this is so.
	/* Remove the one with the highest index first so not to mess up index lists */
	int_t bndsAdded = 0;
	for ( k=defect-1; k>-1; k-- )
	{
		/* Zero curvature in the Hessian: add a bound */
		if ( zeroPivots[k] < nFR )
		{
			number = bounds.getFree()->getNumber( zeroPivots[k] );

			if ( options.printLevel == PL_HIGH )
				MyPrintf( "WARNING: KKT matrix singular! Add bound %i before refactorization.\n", number);

			/* This is just heuristics: we need the bound which gives correct multiplier sign */
			if ( x[number] - lb[number] < ub[number] - x[number] )
				B_status = ST_LOWER;
			else
				B_status = ST_UPPER;

			/* Here we do not need to update the Schur complement because KKT matrix is factorized afterwards */
			if ( bounds.moveFreeToFixed( number, B_status ) != SUCCESSFUL_RETURN )
				return RET_ADDBOUND_FAILED;

			/* Adjust bounds */
			if ( B_status == ST_LOWER )
				lb[number] = x[number];
			else
				ub[number] = x[number];

			bndsAdded++;
		}
		/* Linearly dependent row in the Jacobian: remove a constraint */
		else
		{
			number = constraints.getActive()->getNumber( zeroPivots[k]-nFR );
			if ( options.printLevel == PL_HIGH )
				MyPrintf( "WARNING: KKT matrix singular! Removing constraint %i before refactorization.\n", number);

			if ( constraints.moveActiveToInactive( number ) != SUCCESSFUL_RETURN )
				return RET_REMOVECONSTRAINT_FAILED;

			// AW: If this is an equality constraint, it is now inactive and
			// will not be considered in the step computation which leads to
			// violation of that constraint in the future. Here, I try to
			// fix this by simply making this constraint no longer an
			// equality.
			// TODO: This is probably also necessary for bound constraints
			if ( constraints.getType(number) == ST_EQUALITY )
			{
				if ( options.printLevel == PL_HIGH )
					MyPrintf( "WARNING: Making this constraint no longer an equality.\n");
				constraints.setType( number, ST_BOUNDED );
			}

			/* Adjust dual variable */
			y[number] = 0.0;
		}
	}

	if ( options.printLevel == PL_HIGH )
		MyPrintf( "WARNING: KKT matrix singular! Removed %i constraints and added %i bounds before refactorization.\n",
					defect-bndsAdded, bndsAdded );

	delete[] zeroPivots;

	return SUCCESSFUL_RETURN;
}

END_NAMESPACE_QPOASES


/*
 *	end of file
 */
