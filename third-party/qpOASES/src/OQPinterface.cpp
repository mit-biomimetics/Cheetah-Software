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
 *	\file src/OQPinterface.cpp
 *	\author Hans Joachim Ferreau
 *	\version 3.2
 *	\date 2008-2017
 *
 *	Implementation of an interface comprising several utility functions
 *	for solving test problems from the Online QP Benchmark Collection
 *	(This collection is no longer maintained, see 
 *	http://www.qpOASES.org/onlineQP for a backup).
 *
 */


#include <qpOASES/extras/OQPinterface.hpp>
#include <qpOASES/QProblem.hpp>


BEGIN_NAMESPACE_QPOASES


/*
 *	r e a d O q p D i m e n s i o n s
 */
returnValue readOqpDimensions(	const char* path,
								int_t& nQP, int_t& nV, int_t& nC, int_t& nEC
								)
{
	/* 1) Setup file name where dimensions are stored. */
	char filename[MAX_STRING_LENGTH];
	snprintf( filename,MAX_STRING_LENGTH,"%sdims.oqp",path );

	/* 2) Load dimensions from file. */
	int_t dims[4];
	if ( readFromFile( dims,4,filename ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_UNABLE_TO_READ_FILE );

	nQP = dims[0];
	nV  = dims[1];
	nC  = dims[2];
	nEC = dims[3];


	/* consistency check */
	if ( ( nQP <= 0 ) || ( nV <= 0 ) || ( nC < 0 ) || ( nEC < 0 ) )
		return THROWERROR( RET_FILEDATA_INCONSISTENT );

	return SUCCESSFUL_RETURN;
}


/*
 *	r e a d O q p D a t a
 */
returnValue readOqpData(	const char* path,
							int_t& nQP, int_t& nV, int_t& nC, int_t& nEC,
							real_t** H, real_t** g, real_t** A, real_t** lb, real_t** ub, real_t** lbA, real_t** ubA,
							real_t** xOpt, real_t** yOpt, real_t** objOpt
							)
{
	char filename[MAX_STRING_LENGTH];

	/* consistency check */
	if ( ( H == 0 ) || ( g == 0 ) || ( lb == 0 ) || ( ub == 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );


	/* 1) Obtain OQP dimensions. */
	if ( readOqpDimensions( path, nQP,nV,nC,nEC ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_UNABLE_TO_READ_FILE );


	/* another consistency check */
	if ( ( nC > 0 ) && ( ( A == 0 ) || ( lbA == 0 ) || ( ubA == 0 ) ) )
		return THROWERROR( RET_FILEDATA_INCONSISTENT );


	/* 2) Allocate memory and load OQP data: */
	/* Hessian matrix */
	*H  = new real_t[nV*nV];
	snprintf( filename,MAX_STRING_LENGTH,"%sH.oqp",path );
	if ( readFromFile( *H,nV,nV,filename ) != SUCCESSFUL_RETURN )
	{
		delete[] *H;
		return THROWERROR( RET_UNABLE_TO_READ_FILE );
	}

	/* gradient vector sequence */
	*g  = new real_t[nQP*nV];
	snprintf( filename,MAX_STRING_LENGTH,"%sg.oqp",path );
	if ( readFromFile( *g,nQP,nV,filename ) != SUCCESSFUL_RETURN )
	{
		delete[] *g; delete[] *H;
		return THROWERROR( RET_UNABLE_TO_READ_FILE );
	}

	/* lower bound vector sequence */
	*lb  = new real_t[nQP*nV];
	snprintf( filename,MAX_STRING_LENGTH,"%slb.oqp",path );
	if ( readFromFile( *lb,nQP,nV,filename ) != SUCCESSFUL_RETURN )
	{
		delete[] *lb; delete[] *g; delete[] *H;
		return THROWERROR( RET_UNABLE_TO_READ_FILE );
	}

	/* upper bound vector sequence */
	*ub  = new real_t[nQP*nV];
	snprintf( filename,MAX_STRING_LENGTH,"%sub.oqp",path );
	if ( readFromFile( *ub,nQP,nV,filename ) != SUCCESSFUL_RETURN )
	{
		delete[] *ub; delete[] *lb; delete[] *g; delete[] *H;
		return THROWERROR( RET_UNABLE_TO_READ_FILE );
	}

	if ( nC > 0 )
	{
		/* Constraint matrix */
		*A   = new real_t[nC*nV];
		snprintf( filename,MAX_STRING_LENGTH,"%sA.oqp",path );
		if ( readFromFile( *A,nC,nV,filename ) != SUCCESSFUL_RETURN )
		{
			delete[] *A;
			delete[] *ub; delete[] *lb; delete[] *g; delete[] *H;
			return THROWERROR( RET_UNABLE_TO_READ_FILE );
		}

		/* lower constraints' bound vector sequence */
		*lbA = new real_t[nQP*nC];
		snprintf( filename,MAX_STRING_LENGTH,"%slbA.oqp",path );
		if ( readFromFile( *lbA,nQP,nC,filename ) != SUCCESSFUL_RETURN )
		{
			delete[] *lbA; delete[] *A;
			delete[] *ub; delete[] *lb; delete[] *g; delete[] *H;
			return THROWERROR( RET_UNABLE_TO_READ_FILE );
		}

		/* upper constraints' bound vector sequence */
		*ubA = new real_t[nQP*nC];
		snprintf( filename,MAX_STRING_LENGTH,"%subA.oqp",path );
		if ( readFromFile( *ubA,nQP,nC,filename ) != SUCCESSFUL_RETURN )
		{
			delete[] *ubA; delete[] *lbA; delete[] *A;
			delete[] *ub; delete[] *lb; delete[] *g; delete[] *H;
			return THROWERROR( RET_UNABLE_TO_READ_FILE );
		}
	}
	else
	{
		*A = 0;
		*lbA = 0;
		*ubA = 0;
	}

	if ( xOpt != 0 )
	{
		/* primal solution vector sequence */
		*xOpt = new real_t[nQP*nV];
		snprintf( filename,MAX_STRING_LENGTH,"%sx_opt.oqp",path );
		if ( readFromFile( *xOpt,nQP,nV,filename ) != SUCCESSFUL_RETURN )
		{
			delete[] xOpt;
			if ( nC > 0 ) { delete[] *ubA; delete[] *lbA; delete[] *A; };
			delete[] *ub; delete[] *lb; delete[] *g; delete[] *H;
			return THROWERROR( RET_UNABLE_TO_READ_FILE );
		}
	}

	if ( yOpt != 0 )
	{
		/* dual solution vector sequence */
		*yOpt = new real_t[nQP*(nV+nC)];
		snprintf( filename,MAX_STRING_LENGTH,"%sy_opt.oqp",path );
		if ( readFromFile( *yOpt,nQP,nV+nC,filename ) != SUCCESSFUL_RETURN )
		{
			delete[] yOpt;
			if ( xOpt != 0 ) { delete[] xOpt; };
			if ( nC > 0 ) { delete[] *ubA; delete[] *lbA; delete[] *A; };
			delete[] *ub; delete[] *lb; delete[] *g; delete[] *H;
			return THROWERROR( RET_UNABLE_TO_READ_FILE );
		}
	}

	if ( objOpt != 0 )
	{
		/* dual solution vector sequence */
		*objOpt = new real_t[nQP];
		snprintf( filename,MAX_STRING_LENGTH,"%sobj_opt.oqp",path );
		if ( readFromFile( *objOpt,nQP,1,filename ) != SUCCESSFUL_RETURN )
		{
			delete[] objOpt;
			if ( yOpt != 0 ) { delete[] yOpt; };
			if ( xOpt != 0 ) { delete[] xOpt; };
			if ( nC > 0 ) { delete[] *ubA; delete[] *lbA; delete[] *A; };
			delete[] *ub; delete[] *lb; delete[] *g; delete[] *H;
			return THROWERROR( RET_UNABLE_TO_READ_FILE );
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s o l v e O q p B e n c h m a r k
 */
returnValue solveOqpBenchmark(	int_t nQP, int_t nV, int_t nC, int_t nEC,
								const real_t* const _H, const real_t* const g, const real_t* const _A,
								const real_t* const lb, const real_t* const ub,
								const real_t* const lbA, const real_t* const ubA,
								BooleanType isSparse, 
								const Options& options, int_t& nWSR, real_t& maxCPUtime,
								real_t& maxStationarity, real_t& maxFeasibility, real_t& maxComplementarity
								)
{
	real_t maxNWSR = 0.0;
	real_t avgNWSR = 0.0;
	real_t avgCPUtime = 0.0;

	returnValue returnvalue = solveOqpBenchmark(	nQP,nV,nC,nEC,
													_H,g,_A,lb,ub,lbA,ubA,
													isSparse,BT_TRUE,
													options,nWSR,
													maxNWSR,avgNWSR,maxCPUtime,avgCPUtime, 
													maxStationarity,maxFeasibility,maxComplementarity
													);
	nWSR = (int_t)maxNWSR;

	return returnvalue;
}



/*
 *	s o l v e O q p B e n c h m a r k
 */
returnValue solveOqpBenchmark(	int_t nQP, int_t nV, int_t nC, int_t nEC,
								const real_t* const _H, const real_t* const g, const real_t* const _A,
								const real_t* const lb, const real_t* const ub,
								const real_t* const lbA, const real_t* const ubA,
								BooleanType isSparse, BooleanType useHotstarts, 
								const Options& options, int_t maxAllowedNWSR,
								real_t& maxNWSR, real_t& avgNWSR, real_t& maxCPUtime, real_t& avgCPUtime,
								real_t& maxStationarity, real_t& maxFeasibility, real_t& maxComplementarity
								)
{
	int_t k;

	/* I) SETUP AUXILIARY VARIABLES: */
	/* 1) Keep nWSR and store current and maximum number of
	 *    working set recalculations in temporary variables */
	int_t nWSRcur;

	real_t CPUtimeLimit = maxCPUtime;
	real_t CPUtimeCur = CPUtimeLimit;
	maxNWSR = 0.0;
	avgNWSR = 0.0;
	maxCPUtime = 0.0;
	avgCPUtime = 0.0;
	maxStationarity    = 0.0;
	maxFeasibility     = 0.0;
	maxComplementarity = 0.0;
	real_t stat, feas, cmpl;

	/* 2) Pointers to data of current QP ... */
	const real_t* gCur;
	const real_t* lbCur;
	const real_t* ubCur;
	const real_t* lbACur;
	const real_t* ubACur;

	/* 3) Vectors for solution obtained by qpOASES. */
	real_t* x = new real_t[nV];
	real_t* y = new real_t[nV+nC];
	//real_t obj;

	/* 4) Prepare matrix objects */
	SymmetricMatrix *H; 
	Matrix *A;

	real_t* H_cpy = new real_t[nV*nV];
	memcpy( H_cpy,_H, ((uint_t)(nV*nV))*sizeof(real_t) );
	real_t* A_cpy = new real_t[nC*nV];
	memcpy( A_cpy,_A, ((uint_t)(nC*nV))*sizeof(real_t) );

	if ( isSparse == BT_TRUE )
	{
		SymSparseMat *Hs;
		H = Hs = new SymSparseMat(nV, nV, nV, H_cpy);
		A = new SparseMatrixRow(nC, nV, nV, A_cpy);
		Hs->createDiagInfo();
		delete[] A_cpy; delete[] H_cpy;
	}
	else
	{
		H = new SymDenseMat(nV, nV, nV, const_cast<real_t *>(H_cpy));
		A = new DenseMatrix(nC, nV, nV, const_cast<real_t *>(A_cpy));
	}

	H->doFreeMemory( );
	A->doFreeMemory( );

	/* II) SETUP QPROBLEM OBJECT */
	QProblem qp( nV,nC );
	qp.setOptions( options );
	//qp.setPrintLevel( PL_LOW );

	//qp.printOptions();

	/* III) RUN BENCHMARK SEQUENCE: */
	returnValue returnvalue;

	for( k=0; k<nQP; ++k )
	{
		//if ( k%50 == 0 )
		//	printf( "%d\n",k );

		/* 1) Update pointers to current QP data. */
		gCur   = &( g[k*nV] );
		lbCur  = &( lb[k*nV] );
		ubCur  = &( ub[k*nV] );
		lbACur = &( lbA[k*nC] );
		ubACur = &( ubA[k*nC] );

		/* 2) Set nWSR and maximum CPU time. */
		nWSRcur = maxAllowedNWSR;
		CPUtimeCur = CPUtimeLimit;

		/* 3) Solve current QP. */
		if ( ( k == 0 ) || ( useHotstarts == BT_FALSE ) )
		{
			/* initialise */
			returnvalue = qp.init( H,gCur,A,lbCur,ubCur,lbACur,ubACur, nWSRcur,&CPUtimeCur );
			if ( ( returnvalue != SUCCESSFUL_RETURN ) && ( returnvalue != RET_MAX_NWSR_REACHED ) )
			{
				delete A; delete H; delete[] y; delete[] x;
				return THROWERROR( returnvalue );
			}
		}
		else
		{
			/* hotstart */
			returnvalue = qp.hotstart( gCur,lbCur,ubCur,lbACur,ubACur, nWSRcur,&CPUtimeCur );
			if ( ( returnvalue != SUCCESSFUL_RETURN ) && ( returnvalue != RET_MAX_NWSR_REACHED ) )
			{
				delete A; delete H; delete[] y; delete[] x;
				return THROWERROR( returnvalue );
			}
		}

		/* 4) Obtain solution vectors and objective function value */
		qp.getPrimalSolution( x );
		qp.getDualSolution( y );
		//obj = qp.getObjVal( );

		/* 5) Compute KKT residuals */
		getKktViolation( nV,nC, _H,gCur,_A,lbCur,ubCur,lbACur,ubACur, x,y, stat,feas,cmpl );
		
		/* 6) Update maximum and average values. */
		if ( ((double)nWSRcur) > maxNWSR )
			maxNWSR = ((double)nWSRcur);
		if (stat > maxStationarity) maxStationarity = stat;
		if (feas > maxFeasibility) maxFeasibility = feas;
		if (cmpl > maxComplementarity) maxComplementarity = cmpl;

		if ( CPUtimeCur > maxCPUtime )
			maxCPUtime = CPUtimeCur;

		avgNWSR += ((double)nWSRcur);
		avgCPUtime += CPUtimeCur;
	}
	avgNWSR /= ((double)nQP);
	avgCPUtime /= ((double)nQP);

	delete A; delete H; delete[] y; delete[] x;

	return SUCCESSFUL_RETURN;
}


/*
 *	s o l v e O q p B e n c h m a r k
 */
returnValue solveOqpBenchmark(	int_t nQP, int_t nV,
								const real_t* const _H, const real_t* const g,
								const real_t* const lb, const real_t* const ub,
								BooleanType isSparse, 
								const Options& options, int_t& nWSR, real_t& maxCPUtime,
								real_t& maxStationarity, real_t& maxFeasibility, real_t& maxComplementarity
								)
{
	real_t maxNWSR = 0.0;
	real_t avgNWSR = 0.0;
	real_t avgCPUtime = 0.0;

	returnValue returnvalue = solveOqpBenchmark(	nQP,nV,
													_H,g,lb,ub,
													isSparse,BT_TRUE,
													options,nWSR,
													maxNWSR,avgNWSR,maxCPUtime,avgCPUtime, 
													maxStationarity,maxFeasibility,maxComplementarity
													);
	nWSR = (int_t)maxNWSR;

	return returnvalue;
}


/*
 *	s o l v e O q p B e n c h m a r k
 */
returnValue solveOqpBenchmark(	int_t nQP, int_t nV,
								const real_t* const _H, const real_t* const g,
								const real_t* const lb, const real_t* const ub,
								BooleanType isSparse, BooleanType useHotstarts, 
								const Options& options, int_t maxAllowedNWSR,
								real_t& maxNWSR, real_t& avgNWSR, real_t& maxCPUtime, real_t& avgCPUtime,
								real_t& maxStationarity, real_t& maxFeasibility, real_t& maxComplementarity
								)
{
	int_t k;

	/* I) SETUP AUXILIARY VARIABLES: */
	/* 1) Keep nWSR and store current and maximum number of
	 *    working set recalculations in temporary variables */
	int_t nWSRcur;

	real_t CPUtimeLimit = maxCPUtime;
	real_t CPUtimeCur = CPUtimeLimit;
	real_t stat, feas, cmpl;
	maxNWSR = 0;
	avgNWSR = 0;
	maxCPUtime = 0.0;
	avgCPUtime = 0.0;
	maxStationarity = 0.0;
	maxFeasibility = 0.0;
	maxComplementarity = 0.0;

	/* 2) Pointers to data of current QP ... */
	const real_t* gCur;
	const real_t* lbCur;
	const real_t* ubCur;

	/* 3) Vectors for solution obtained by qpOASES. */
	real_t* x = new real_t[nV];
	real_t* y = new real_t[nV];
	//real_t  obj;

	/* 4) Prepare matrix objects */
	SymmetricMatrix *H; 
	real_t* H_cpy = new real_t[nV*nV];
	memcpy( H_cpy,_H, ((uint_t)(nV*nV))*sizeof(real_t) );

	if ( isSparse == BT_TRUE )
	{
		SymSparseMat *Hs;
		H = Hs = new SymSparseMat(nV, nV, nV, H_cpy);
		Hs->createDiagInfo();
		delete[] H_cpy;
	}
	else
	{
		H = new SymDenseMat(nV, nV, nV, const_cast<real_t *>(H_cpy));
	}
	
	H->doFreeMemory( );

	/* II) SETUP QPROBLEM OBJECT */
	QProblemB qp( nV );
	qp.setOptions( options );
	//qp.setPrintLevel( PL_LOW );


	/* III) RUN BENCHMARK SEQUENCE: */
	returnValue returnvalue;

	for( k=0; k<nQP; ++k )
	{
		//if ( k%50 == 0 )
		//	printf( "%d\n",k );

		/* 1) Update pointers to current QP data. */
		gCur   = &( g[k*nV] );
		lbCur  = &( lb[k*nV] );
		ubCur  = &( ub[k*nV] );

		/* 2) Set nWSR and maximum CPU time. */
		nWSRcur = maxAllowedNWSR;
		CPUtimeCur = CPUtimeLimit;

		/* 3) Solve current QP. */
		if ( ( k == 0 ) || ( useHotstarts == BT_FALSE ) )
		{
			/* initialise */
			returnvalue = qp.init( H,gCur,lbCur,ubCur, nWSRcur,&CPUtimeCur );
			if ( ( returnvalue != SUCCESSFUL_RETURN ) && ( returnvalue != RET_MAX_NWSR_REACHED ) )
			{
				delete H; delete[] y; delete[] x;
				return THROWERROR( returnvalue );
			}
		}
		else
		{
			/* hotstart */
			returnvalue = qp.hotstart( gCur,lbCur,ubCur, nWSRcur,&CPUtimeCur );
			if ( ( returnvalue != SUCCESSFUL_RETURN ) && ( returnvalue != RET_MAX_NWSR_REACHED ) )
			{
				delete H; delete[] y; delete[] x;
				return THROWERROR( returnvalue );
			}
		}

		/* 4) Obtain solution vectors and objective function value ... */
		qp.getPrimalSolution( x );
		qp.getDualSolution( y );
		//obj = qp.getObjVal( );

		/* 5) Compute KKT residuals */
		getKktViolation( nV, _H,gCur,lbCur,ubCur, x,y, stat,feas,cmpl );

		/* 6) update maximum values. */
		if ( ((real_t)nWSRcur) > maxNWSR )
			maxNWSR = ((real_t)nWSRcur);
		if (stat > maxStationarity) maxStationarity = stat;
		if (feas > maxFeasibility) maxFeasibility = feas;
		if (cmpl > maxComplementarity) maxComplementarity = cmpl;

		if ( CPUtimeCur > maxCPUtime )
			maxCPUtime = CPUtimeCur;

		avgNWSR += (real_t)nWSRcur;
		avgCPUtime += CPUtimeCur;
	}
	avgNWSR /= ((real_t)nQP);
	avgCPUtime /= ((real_t)nQP);

	delete H; delete[] y; delete[] x;

	return SUCCESSFUL_RETURN;
}


/*
 *	r u n O q p B e n c h m a r k
 */
returnValue runOqpBenchmark(	const char* path, BooleanType isSparse, const Options& options,
								int_t& nWSR, real_t& maxCPUtime,
								real_t& maxStationarity, real_t& maxFeasibility, real_t& maxComplementarity
								)
{
	real_t maxNWSR = 0.0;
	real_t avgNWSR = 0.0;
	real_t avgCPUtime = 0.0;

	returnValue returnvalue = runOqpBenchmark(	path,isSparse,BT_TRUE,
												options,nWSR,
												maxNWSR,avgNWSR,maxCPUtime,avgCPUtime, 
												maxStationarity,maxFeasibility,maxComplementarity
												);
	nWSR = (int_t)maxNWSR;

	return returnvalue;
}


/*
 *	r u n O q p B e n c h m a r k
 */
returnValue runOqpBenchmark(	const char* path, BooleanType isSparse, BooleanType useHotstarts, 
								const Options& options, int_t maxAllowedNWSR,
								real_t& maxNWSR, real_t& avgNWSR, real_t& maxCPUtime, real_t& avgCPUtime,
								real_t& maxStationarity, real_t& maxFeasibility, real_t& maxComplementarity
								)
{
	int_t nQP=0, nV=0, nC=0, nEC=0;

	real_t *H=0, *g=0, *A=0, *lb=0, *ub=0, *lbA=0, *ubA=0;


	returnValue returnvalue;

	/* I) SETUP BENCHMARK: */
	/* 1) Obtain QP sequence dimensions. */
	if ( readOqpDimensions( path, nQP,nV,nC,nEC ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_UNABLE_TO_READ_BENCHMARK );

	/* 2) Read OQP benchmark data. */
	if ( readOqpData(	path,
						nQP,nV,nC,nEC,
						&H,&g,&A,&lb,&ub,&lbA,&ubA,
						0,0,0
						) != SUCCESSFUL_RETURN )
	{
		return THROWERROR( RET_UNABLE_TO_READ_BENCHMARK );
	}

	// normaliseConstraints( nV,nC,A,lbA,ubA ); //only works when nP==1

	/* II) SOLVE BENCHMARK */
	if ( nC > 0 )
	{
		returnvalue = solveOqpBenchmark(	nQP,nV,nC,nEC,
											H,g,A,lb,ub,lbA,ubA,
											isSparse,useHotstarts,
											options,maxAllowedNWSR,
											maxNWSR,avgNWSR,maxCPUtime,avgCPUtime,
											maxStationarity,maxFeasibility,maxComplementarity
											);

		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			if ( H != 0 )   delete[] H;
			if ( g != 0 )   delete[] g;
			if ( A != 0 )   delete[] A;
			if ( lb != 0 )  delete[] lb;
			if ( ub != 0 )  delete[] ub;
			if ( lbA != 0 ) delete[] lbA;
			if ( ubA != 0 ) delete[] ubA;
			return THROWERROR( returnvalue );
		}
	}
	else
	{
		returnvalue = solveOqpBenchmark(	nQP,nV,
											H,g,lb,ub,
											isSparse,useHotstarts,
											options,maxAllowedNWSR,
											maxNWSR,avgNWSR,maxCPUtime,avgCPUtime,
											maxStationarity,maxFeasibility,maxComplementarity
											);

		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			if ( H != 0 )   delete[] H;
			if ( g != 0 )   delete[] g;
			if ( A != 0 )   delete[] A;
			if ( lb != 0 )  delete[] lb;
			if ( ub != 0 )  delete[] ub;
			return THROWERROR( returnvalue );
		}
	}

	if ( H != 0 )   delete[] H;
	if ( g != 0 )   delete[] g;
	if ( A != 0 )   delete[] A;
	if ( lb != 0 )  delete[] lb;
	if ( ub != 0 )  delete[] ub;
	if ( lbA != 0 ) delete[] lbA;
	if ( ubA != 0 ) delete[] ubA;

	return SUCCESSFUL_RETURN;
}


END_NAMESPACE_QPOASES


/*
 *	end of file
 */
