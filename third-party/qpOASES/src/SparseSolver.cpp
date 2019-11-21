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
 *	\file src/SparseSolver.cpp
 *	\author Andreas Waechter, Dennis Janka
 *	\version 3.2
 *	\date 2012-2017
 *
 *	Interfaces to sparse linear solvers that are used in a Schur-complement
 *	implementation in qpOASES.
 */


#include <qpOASES/SparseSolver.hpp>

#ifndef __MATLAB__
# include <cstdarg>
void MyPrintf(const char* pformat, ... );
#else
# include <mex.h>
# define MyPrintf mexPrintf
#endif

BEGIN_NAMESPACE_QPOASES

/*****************************************************************************
 *  P U B L I C                                                              *
 ****************************************************************************/

/*
 *	S p a r s e S o l v e r
 */
SparseSolver::SparseSolver( )
{
}


/*
 *	S p a r s e S o l v e r
 */
SparseSolver::SparseSolver( const SparseSolver& rhs )
{
	copy( rhs );
}


/*
 *	~ S p a r s e S o l v e r
 */
SparseSolver::~SparseSolver( )
{
	clear( );
}


/*
 *	o p e r a t o r =
 */
SparseSolver& SparseSolver::operator=( const SparseSolver& rhs )
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
returnValue SparseSolver::reset( )
{
	return SUCCESSFUL_RETURN;
}

/*
 *	g e t N e g a t i v e E i g e n v a l u e s
 */
int_t SparseSolver::getNegativeEigenvalues( )
{
	return -1;
}

/*
 *	g e t R a n k
 */
int_t SparseSolver::getRank( )
{
	return -1;
}

/*
 *	g e t Z e r o P i v o t s
 */
returnValue SparseSolver::getZeroPivots( int_t *&zeroPivots )
{
	if ( zeroPivots ) delete[] zeroPivots;
	zeroPivots = 0;
	return SUCCESSFUL_RETURN;
}


/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/

/*
 *	c l e a r
 */
returnValue SparseSolver::clear( )
{
	return SUCCESSFUL_RETURN;
}


/*
 *	c o p y
 */
returnValue SparseSolver::copy( 	const SparseSolver& rhs
									)
{
	return SUCCESSFUL_RETURN;
}

#ifdef SOLVER_MA27

/*****************************************************************************
 *****************************************************************************
 *****************************************************************************
 *  M A 2 7 S P A R E S E S O L V E R                                        *
 *****************************************************************************
 *****************************************************************************
 *****************************************************************************/

#define MA27ID ma27id_
#define MA27AD ma27ad_
#define MA27BD ma27bd_
#define MA27CD ma27cd_

extern "C" {
  void MA27ID( fint_t* ICNTL, double* CNTL );
  void MA27AD( fint_t *N, fint_t *NZ, const fint_t *IRN, const fint_t* ICN,
               fint_t *IW, fint_t* LIW, fint_t* IKEEP, fint_t *IW1,
               fint_t* NSTEPS, fint_t* IFLAG, fint_t* ICNTL,
               double* CNTL, fint_t *INFO, double* OPS);
  void MA27BD( fint_t *N, fint_t *NZ, const fint_t *IRN, const fint_t* ICN,
               double* A, fint_t* LA, fint_t* IW, fint_t* LIW,
               fint_t* IKEEP, fint_t* NSTEPS, fint_t* MAXFRT,
               fint_t* IW1, fint_t* ICNTL, double* CNTL,
               fint_t* INFO);
  void MA27CD( fint_t *N, double* A, fint_t* LA, fint_t* IW,
               fint_t* LIW, double* W, fint_t* MAXFRT,
               double* RHS, fint_t* IW1, fint_t* NSTEPS,
               fint_t* ICNTL, double* CNTL);
}

/*****************************************************************************
 *  P U B L I C                                                              *
 ****************************************************************************/

/*
 *	M a 2 7 S p a r s e S o l v e r
 */
Ma27SparseSolver::Ma27SparseSolver( ) : SparseSolver()
{
	a_ma27 = 0;
	irn_ma27 = 0;
	jcn_ma27 = 0;
	iw_ma27 = 0;
	ikeep_ma27 = 0;
	clear( );

	/* Set default options for MA27 */
	MA27ID(icntl_ma27, cntl_ma27);
	icntl_ma27[0] = 0;       /* Suppress error messages */
	icntl_ma27[1] = 0;       /* Suppress diagnostic messages */
	cntl_ma27[0] = 1e-8;     /* Set pivot tolerance */
}


/*
 *	M a 2 7 S p a r s e S o l v e r
 */
Ma27SparseSolver::Ma27SparseSolver( const Ma27SparseSolver& rhs )
{
	copy( rhs );
}


/*
 *	~ M a 2 7 S p a r s e S o l v e r
 */
Ma27SparseSolver::~Ma27SparseSolver( )
{
	clear( );
}


/*
 *	o p e r a t o r =
 */
Ma27SparseSolver& Ma27SparseSolver::operator=( const SparseSolver& rhs )
{
	const Ma27SparseSolver* ma27_rhs = dynamic_cast<const Ma27SparseSolver*>(&rhs);
	if (!ma27_rhs)
	{
		fprintf(getGlobalMessageHandler()->getOutputFile(),"Error in Ma27SparseSolver& Ma27SparseSolver::operator=( const SparseSolver& rhs )\n");
		throw; /* TODO: More elegant exit? */
	}
	if ( this != ma27_rhs )
	{
		clear( );
		SparseSolver::operator=( rhs );
		copy( *ma27_rhs );
	}

	return *this;
}

/*
 *	s e t M a t r i x D a t a
 */
returnValue Ma27SparseSolver::setMatrixData( int_t dim_,
						   int_t numNonzeros_,
						   const int_t* const irn,
						   const int_t* const jcn,
						   const real_t* const avals
						   )
{
	reset( );
	dim = dim_;
	numNonzeros = numNonzeros_;

	if ( numNonzeros_ > 0 )
	{
		a_ma27 = new double[numNonzeros_];
		irn_ma27 = new fint_t[numNonzeros_];
		jcn_ma27 = new fint_t[numNonzeros_];

		numNonzeros=0;
		for (int_t i=0; i<numNonzeros_; ++i)
			if ( avals[i] != 0 )
			{
				a_ma27[numNonzeros] = avals[i];
				irn_ma27[numNonzeros] = irn[i];
				jcn_ma27[numNonzeros] = jcn[i];
				numNonzeros++;
			}
	}
	else
	{
		numNonzeros = 0;
		a_ma27 = 0;
		irn_ma27 = 0;
		jcn_ma27 = 0;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	f a c t o r i z e
 */
returnValue Ma27SparseSolver::factorize( )
{
	if ( dim == 0 )
	{
		have_factorization = true;
		neig = 0;
		rank = 0;
		return SUCCESSFUL_RETURN;
	}

	/******************************************
	 * Call MA27AD for symbolic factorization *
	******************************************/

	/* Overstimation factor for LIW (20% recommended in MA27 documentation) */
	const double LiwFact = 2.0;   /* This is 200% overestimation */
	liw_ma27 = (fint_t)(LiwFact*(double(2*numNonzeros+3*dim+1)));
	iw_ma27 = new fint_t[liw_ma27];

	ikeep_ma27 = new fint_t[3*dim];

	fint_t iflag_ma27 = 0;
	double ops_ma27;
	fint_t info_ma27[20];
	fint_t* iw1_ma27 = new fint_t[2*dim];
	MA27AD(&dim, &numNonzeros, irn_ma27, jcn_ma27, iw_ma27, &liw_ma27, ikeep_ma27,
                            iw1_ma27, &nsteps_ma27, &iflag_ma27, icntl_ma27, cntl_ma27,
                            info_ma27, &ops_ma27);

	/* Receive some information from MA27AD */
	fint_t iflag  = info_ma27[0];  /* Information flag */
	fint_t ierror = info_ma27[1];  /* Error flag */
	fint_t nrlnec = info_ma27[4];  /* recommended value for la */
	fint_t nirnec = info_ma27[5];  /* recommended value for liw */
	if (iflag != 0)
	{
		MyPrintf("MA27AD returns iflag = %d with ierror = %d\n", iflag, ierror);
		delete [] iw1_ma27;
		clear( );
		return THROWERROR(RET_MATRIX_FACTORISATION_FAILED);
	}

	/* Allocate memory for actual factorization */
	delete [] iw_ma27;
	double liw_init_factor = 5.0; /* This could be an option. */
	liw_ma27 = (fint_t)(liw_init_factor * (double)(nirnec));
	iw_ma27 = new fint_t[liw_ma27];

	double la_init_factor = 20.0; /* This could be an option. */
	la_ma27 = getMax(numNonzeros,(fint_t)(la_init_factor * (double)(nrlnec)));
	double* a_new = new double[la_ma27];
	for (int_t i=0; i<numNonzeros; ++i)
		a_new[i] = a_ma27[i];
	delete [] a_ma27;
	a_ma27 = a_new;

    /*******************************************
	 * Call MA27BD for numerical factorization *
     *******************************************/

	MA27BD(&dim, &numNonzeros, irn_ma27, jcn_ma27, a_ma27,
		   &la_ma27, iw_ma27, &liw_ma27, ikeep_ma27, &nsteps_ma27,
		   &maxfrt_ma27, iw1_ma27, icntl_ma27, cntl_ma27, info_ma27);

	delete [] iw1_ma27;
	/* Receive some information from MA27BD */
	iflag = info_ma27[0];   /* Information flag */
	ierror = info_ma27[1];  /* Error flag */
	neig = info_ma27[14];   /* Number of negative eigenvalues */
	if (iflag == 3)
	{
		rank = info_ma27[1];
		return RET_KKT_MATRIX_SINGULAR;
	}
	else if (iflag == -5)
	{ /* DJ: I think this is more severe. Can this actually happen? */
		rank = -1;
		return RET_KKT_MATRIX_SINGULAR;
	}
	else if (iflag != 0)
	{
		MyPrintf("MA27BD returns iflag = %d with ierror = %d\n", iflag, ierror);
		clear( );
		return THROWERROR(RET_MATRIX_FACTORISATION_FAILED);
	}
	else
		rank = dim;

	have_factorization = true;
	return SUCCESSFUL_RETURN;
}


/*
 *	s o l v e
 */
returnValue Ma27SparseSolver::solve( int_t dim_,
					   const real_t* const rhs,
					   real_t* const sol
					   )
{
	/* consistency check */
	if ( dim_ != dim )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	if ( !have_factorization )
	{
	  MyPrintf("Factorization not called before solve in Ma27SparseSolver::solve.\n");
	  return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	if ( dim == 0 )
		return SUCCESSFUL_RETURN;

	/* Call MA27CD to solve the system */
	double* w_ma27 = new double[maxfrt_ma27];
	fint_t* iw1_ma27 = new fint_t[nsteps_ma27];

	/* MA27CD overwrites rhs */
	for (int_t i=0; i<dim; ++i) sol[i] = rhs[i];
	MA27CD(&dim, a_ma27, &la_ma27, iw_ma27, &liw_ma27, w_ma27, &maxfrt_ma27,
		   sol, iw1_ma27, &nsteps_ma27, icntl_ma27, cntl_ma27);

	delete [] iw1_ma27;
	delete [] w_ma27;

	return SUCCESSFUL_RETURN;
}

/*
 *	r e s e t
 */
returnValue Ma27SparseSolver::reset( )
{
	/* AW: We probably want to avoid resetting factorization in QProblem */
	if ( SparseSolver::reset( ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_RESET_FAILED );

	clear( );
	return SUCCESSFUL_RETURN;
}

/*
 *	g e t N e g a t i v e E i g e n v a l u e s
 */
int_t Ma27SparseSolver::getNegativeEigenvalues( )
{
	if( !have_factorization )
		return -1;
	else
		return neig;
}

/*
 *	g e t R a n k
 */
int_t Ma27SparseSolver::getRank( )
{
	return rank;
}


/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/

/*
 *	c l e a r
 */
returnValue Ma27SparseSolver::clear( )
{
	delete [] a_ma27;
	delete [] irn_ma27;
	delete [] jcn_ma27;
	delete [] iw_ma27;
	delete [] ikeep_ma27;

	dim = -1;
	numNonzeros = -1;
	neig = -1;
	rank = -1;
	la_ma27 = -1;
	a_ma27 = 0;
	irn_ma27 = 0;
	jcn_ma27 = 0;

	liw_ma27 = -1;
	iw_ma27 = 0;
	ikeep_ma27 = 0;
	nsteps_ma27 = -1;
	maxfrt_ma27 = -1;

	have_factorization = false;
	return SUCCESSFUL_RETURN;
}


/*
 *	c o p y
 */
returnValue Ma27SparseSolver::copy( 	const Ma27SparseSolver& rhs
										)
{
	dim = rhs.dim;
	numNonzeros = rhs.numNonzeros;
	la_ma27 = rhs.la_ma27;
	if ( rhs.a_ma27 != 0 )
	{
	  if (rhs.have_factorization)
		{
		  a_ma27 = new double[la_ma27];
		  memcpy( a_ma27,rhs.a_ma27,la_ma27*sizeof(double) );
		}
	  else
		{
		  a_ma27 = new double[numNonzeros];
		  memcpy( a_ma27,rhs.a_ma27,numNonzeros*sizeof(double) );
		}
	}
	else
		a_ma27 = 0;

	if ( rhs.irn_ma27 != 0 )
	{
		irn_ma27 = new fint_t[numNonzeros];
		memcpy( irn_ma27,rhs.irn_ma27,numNonzeros*sizeof(fint_t) );
	}
	else
		irn_ma27 = 0;

	if ( rhs.jcn_ma27 != 0 )
	{
		jcn_ma27 = new fint_t[numNonzeros];
		memcpy( jcn_ma27,rhs.jcn_ma27,numNonzeros*sizeof(fint_t) );
	}
	else
		jcn_ma27 = 0;

	for ( int_t i=0; i<30; ++i)
		icntl_ma27[i] = rhs.icntl_ma27[i];

	for ( int_t i=0; i<5; ++i)
		cntl_ma27[i] = rhs.cntl_ma27[i];

	liw_ma27 = rhs.liw_ma27;

	if ( rhs.iw_ma27 != 0 )
	{
		iw_ma27 = new fint_t[liw_ma27];
		memcpy( iw_ma27,rhs.iw_ma27,liw_ma27*sizeof(fint_t) );
	}
	else
		iw_ma27 = 0;

	if ( rhs.ikeep_ma27 != 0 )
	{
		ikeep_ma27 = new fint_t[3*dim];
		memcpy( ikeep_ma27,rhs.ikeep_ma27,3*dim*sizeof(fint_t) );
	}
	else
		ikeep_ma27 = 0;

	nsteps_ma27 = rhs.nsteps_ma27;
	maxfrt_ma27 = rhs.maxfrt_ma27;

	have_factorization = rhs.have_factorization;
	neig = rhs.neig;
	rank = rhs.rank;

	return SUCCESSFUL_RETURN;
}

#endif /* SOLVER_MA27 */

#ifdef SOLVER_MA57

/*****************************************************************************
 *****************************************************************************
 *****************************************************************************
 *  M A 5 7 S P A R E S E S O L V E R                                        *
 *****************************************************************************
 *****************************************************************************
 *****************************************************************************/

#define MA57ID ma57id_
#define MA57AD ma57ad_
#define MA57BD ma57bd_
#define MA57CD ma57cd_

extern "C"
{
	/*
	*  MA57ID -- Initialize solver.
	*/
	extern void  MA57ID (
		double	*cntl,
		fint_t		*icntl);

	/*
	*  MA57AD -- Symbolic Factorization.
	*/
	extern void  MA57AD (
		fint_t			*n,		/* Order of matrix. */
		fint_t			*ne,	/* Number of entries. */
		const fint_t	*irn,	/* Matrix nonzero row structure */
		const fint_t	*jcn,	/* Matrix nonzero column structure */
		fint_t			*lkeep,	/* Workspace for the pivot order of lenght 3*n */
		fint_t			*keep,	/* Workspace for the pivot order of lenght 3*n */
								/* Automatically iflag = 0; ikeep pivot order iflag = 1 */
		fint_t			*iwork,	/* Integer work space. */
		fint_t			*icntl,	/* Integer Control parameter of length 30*/
		fint_t			*info,	/* Statistical Information; Integer array of length 20 */
		double			*rinfo);/* Double Control parameter of length 5 */

	/*
	* MA57BD -- Numerical Factorization.
	*/
	extern void  MA57BD (
		fint_t	*n,			/* Order of matrix. */
		fint_t	*ne,		/* Number of entries. */
		double	*a,			/* Numerical values. */
		double	*fact,		/* Entries of factors. */
		fint_t	*lfact,		/* Length of array `fact'. */
		fint_t	*ifact,		/* Indexing info for factors. */
		fint_t	*lifact,	/* Length of array `ifact'. */
		fint_t	*lkeep,		/* Length of array `keep'. */
		fint_t	*keep,		/* Integer array. */
		fint_t	*iwork,		/* Workspace of length `n'. */
		fint_t	*icntl,		/* Integer Control parameter of length 20. */
		double	*cntl,		/* Double Control parameter of length 5. */
		fint_t	*info,		/* Statistical Information; Integer array of length 40. */
		double	*rinfo);	/* Statistical Information; Real array of length 20. */

	/*
	* MA57CD -- Solution.
	*/
	extern void  MA57CD (
		fint_t	*job,		/* Solution job.  Solve for... */
							/* JOB <= 1:  A */
							/* JOB == 2:  PLP^t */
							/* JOB == 3:  PDP^t */
							/* JOB >= 4:  PL^t P^t */
		fint_t	*n,			/* Order of matrix. */
		double	*fact,		/* Entries of factors. */
		fint_t	*lfact,		/* Length of array `fact'. */
		fint_t	*ifact,		/* Indexing info for factors. */
		fint_t	*lifact,	/* Length of array `ifact'. */
		fint_t	*nrhs,		/* Number of right hand sides. */
		double	*rhs,		/* Numerical Values. */
		fint_t	*lrhs,		/* Leading dimensions of `rhs'. */
		double	*work,		/* Real workspace. */
		fint_t	*lwork,		/* Length of `work', >= N*NRHS. */
		fint_t	*iwork,		/* Integer array of length `n'. */
		fint_t	*icntl,		/* Integer Control parameter array of length 20. */
		fint_t	*info);		/* Statistical Information; Integer array of length 40. */
}

/*****************************************************************************
 *  P U B L I C                                                              *
 ****************************************************************************/


/*
 *	M a 5 7 S p a r s e S o l v e r
 */
Ma57SparseSolver::Ma57SparseSolver( ) : SparseSolver()
{
	a_ma57 = 0;
	irn_ma57 = 0;
	jcn_ma57 = 0;
	fact_ma57 = 0;
	ifact_ma57 = 0;
	pivots = 0;
	clear( );

	/* Set default options for MA57 */
	MA57ID( cntl_ma57, icntl_ma57 );

	icntl_ma57[0] = -1;			/* Suppress error messages */
	icntl_ma57[1] = -1;			/* Suppress warning messages */
	icntl_ma57[2] = -1;			/* Suppress monitoring messages */
	/*icntl_ma57[4] = 4;		// Print everything (for debugging) */
	icntl_ma57[15] = 1;			/* Place small pivots at the end of the factorization (default: 0) */

	/* \todo good default values?
	cntl_ma57[1] = 5.0e-16;		// Pivots smaller than this are treated as zero and are placed at the end of the factorization (default: 1e-20)
	cntl_ma57[0] = 0.5;			// Set pivot tolerance: Higher values = more stable but slower/less sparse (default: 0.01, max 0.5) */
}


/*
 *	M a 5 7 S p a r s e S o l v e r
 */
Ma57SparseSolver::Ma57SparseSolver( const Ma57SparseSolver& rhs )
{
	copy( rhs );
}


/*
 *	~ M a 5 7 S p a r s e S o l v e r
 */
Ma57SparseSolver::~Ma57SparseSolver( )
{
	clear( );
}


/*
 *	o p e r a t o r =
 */
Ma57SparseSolver& Ma57SparseSolver::operator=( const SparseSolver& rhs )
{
	const Ma57SparseSolver* ma57_rhs = dynamic_cast<const Ma57SparseSolver*>(&rhs);
	if (!ma57_rhs)
	{
		fprintf(getGlobalMessageHandler()->getOutputFile(),"Error in Ma57SparseSolver& Ma57SparseSolver::operator=( const SparseSolver& rhs )\n");
		throw; /* TODO: More elegant exit? */
	}
	if ( this != ma57_rhs )
	{
		clear( );
		SparseSolver::operator=( rhs );
		copy( *ma57_rhs );
	}

	return *this;
}

/*
 *	s e t M a t r i x D a t a
 */
returnValue Ma57SparseSolver::setMatrixData(	int_t dim_,
												int_t numNonzeros_,
												const int_t* const irn,
												const int_t* const jcn,
												const real_t* const avals
												)
{
	reset( );
	dim = dim_;
	numNonzeros = numNonzeros_;

	if ( numNonzeros_ > 0 )
	{
		a_ma57 = new double[numNonzeros_];
		irn_ma57 = new fint_t[numNonzeros_];
		jcn_ma57 = new fint_t[numNonzeros_];

		numNonzeros=0;
		for (int_t i=0; i<numNonzeros_; ++i)
			if ( isZero(avals[i]) == BT_FALSE )
			{
				a_ma57[numNonzeros] = avals[i];
				irn_ma57[numNonzeros] = irn[i];
				jcn_ma57[numNonzeros] = jcn[i];
				numNonzeros++;
			}
	}
	else
	{
		numNonzeros = 0;
		a_ma57 = 0;
		irn_ma57 = 0;
		jcn_ma57 = 0;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	f a c t o r i z e
 */
returnValue Ma57SparseSolver::factorize( )
{
	if ( dim == 0 )
	{
		have_factorization = true;
		neig = 0;
		rank = 0;
		return SUCCESSFUL_RETURN;
	}

	/******************************************
	 * Call MA57AD for symbolic factorization *
	******************************************/

	fint_t lkeep_ma57 = 5*dim + numNonzeros + getMax(numNonzeros,dim) + 42;
	fint_t *keep_ma57 = new fint_t[lkeep_ma57];
	fint_t *iwork_ma57 = new fint_t[5*dim];

	/* Set initial pivot sequence. */
	for (fint_t i = 0; i < lkeep_ma57; i++) keep_ma57[i] = i+1;

	fint_t info_ma57[40];
	double rinfo_ma57[20];

	MA57AD(&dim, &numNonzeros, irn_ma57, jcn_ma57, &lkeep_ma57, keep_ma57,
			iwork_ma57, icntl_ma57, info_ma57, rinfo_ma57);

	/* Receive some information from MA57AD */
	fint_t iflag = info_ma57[0];   /* Information flag */
	fint_t ierror = info_ma57[1];  /* Error flag */
	if (iflag != 0)
	{
		MyPrintf("MA57AD returns iflag = %d with ierror = %d\n", iflag, ierror);
		delete [] keep_ma57;
		delete [] iwork_ma57;
		clear( );
		return THROWERROR(RET_MATRIX_FACTORISATION_FAILED);
	}

	/* Allocate memory for actual factorization */
	double lfact_factor = 10.0; /* This could be an option */

	lfact_ma57 = (fint_t)(lfact_factor * (double)(info_ma57[8]));
	fact_ma57 = new double[lfact_ma57];

	lifact_ma57 = (fint_t)(lfact_factor * (double)(info_ma57[9]));
	ifact_ma57 = new int_t[lifact_ma57];

    /*******************************************
	 * Call MA57BD for numerical factorization *
     *******************************************/

	MA57BD(	&dim, &numNonzeros, a_ma57, fact_ma57, &lfact_ma57,
			ifact_ma57, &lifact_ma57, &lkeep_ma57, keep_ma57,
			iwork_ma57, icntl_ma57, cntl_ma57, info_ma57, rinfo_ma57 );

	delete [] iwork_ma57;
	delete [] keep_ma57;

	/* Receive some information from MA57BD */
	iflag = info_ma57[0];   /* Information flag */
	ierror = info_ma57[1];  /* Error flag */
	neig = info_ma57[23];   /* Number of negative eigenvalues */
	rank = info_ma57[24];   /* Rank of matrix */

	/* Read pivot sequence (see MA57UD source code) */
	pivots = new fint_t[dim];
	fint_t nrows, ncols;
	fint_t nblk = ifact_ma57[2];
	int_t iwpos = 3; /* = 4-1 */
    int_t k, kk, count = 0;
    for ( k=0; k<nblk; k++ )
    {
        ncols = ifact_ma57[iwpos];
        nrows = ifact_ma57[iwpos+1];

        for ( kk=0; kk<nrows; kk++ )
            pivots[count++] = ifact_ma57[iwpos+2+kk]-1; /* convert Fortran to C indices! */

        iwpos = iwpos+ncols+2;
    }

	if (iflag == 4)
	{
		/*MyPrintf("dim = %i, rank = %i. Pivots: ", dim, rank);
		for( k=rank; k<dim; k++ )
			MyPrintf("%i ", pivots[k]);
		MyPrintf("\n");*/

		return RET_KKT_MATRIX_SINGULAR;
	}
	else if (iflag != 0)
	{
		MyPrintf("MA57BD returns iflag = %d with ierror = %d\n", iflag, ierror);
		clear( );
		return THROWERROR(RET_MATRIX_FACTORISATION_FAILED);
	}

	have_factorization = true;
	return SUCCESSFUL_RETURN;
}


/*
 *	s o l v e
 */
returnValue Ma57SparseSolver::solve(	int_t dim_,
										const real_t* const rhs,
										real_t* const sol
										)
{
	/* consistency check */
	if ( dim_ != dim )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	if ( !have_factorization )
	{
	  MyPrintf("Factorization not called before solve in Ma57SparseSolver::solve.\n");
	  return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	if ( dim == 0 )
		return SUCCESSFUL_RETURN;

	/* Call MA57CD to solve the system */
	fint_t job_ma57 = 1;
	fint_t nrhs_ma57 = 1;
	fint_t lrhs_ma57 = dim;
	fint_t info_ma57[40];

	fint_t lwork_ma57 = dim*nrhs_ma57;
	double* work_ma57 = new double[lwork_ma57];
	fint_t* iwork_ma57 = new fint_t[dim];

	/* MA57CD overwrites rhs */
	for (int_t i=0; i<dim; ++i) sol[i] = rhs[i];
	MA57CD(&job_ma57, &dim, fact_ma57, &lfact_ma57, ifact_ma57, &lifact_ma57,
			&nrhs_ma57, sol, &lrhs_ma57, work_ma57, &lwork_ma57, iwork_ma57,
			icntl_ma57, info_ma57);

	delete [] work_ma57;
	delete [] iwork_ma57;

	fint_t iflag = info_ma57[0];   /* Information flag */
	fint_t ierror = info_ma57[1];  /* Error flag */
	if (iflag != 0)
	{
		MyPrintf("MA57CD returns iflag = %d with ierror = %d\n", iflag, ierror);
		clear( );
		return THROWERROR(RET_MATRIX_FACTORISATION_FAILED);
	}

	return SUCCESSFUL_RETURN;
}

/*
 *	r e s e t
 */
returnValue Ma57SparseSolver::reset( )
{
	/* AW: We probably want to avoid resetting factorization in QProblem */
	if ( SparseSolver::reset( ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_RESET_FAILED );

	clear( );
	return SUCCESSFUL_RETURN;
}

/*
 *	g e t N e g a t i v e E i g e n v a l u e s
 */
int_t Ma57SparseSolver::getNegativeEigenvalues( )
{
	if( !have_factorization )
		return -1;
	else
		return neig;
}

/*
 *	g e t R a n k
 */
int_t Ma57SparseSolver::getRank( )
{
	return rank;
}

/*
 *	g e t Z e r o P i v o t s
 */
returnValue Ma57SparseSolver::getZeroPivots( int_t *&zeroPivots )
{
	for ( int_t k=0; k<dim-rank; k++ )
		zeroPivots[k] = pivots[rank+k];

	return SUCCESSFUL_RETURN;
}


/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/

/*
 *	c l e a r
 */
returnValue Ma57SparseSolver::clear( )
{
	delete [] a_ma57;
	delete [] irn_ma57;
	delete [] jcn_ma57;
	delete [] fact_ma57;
	delete [] ifact_ma57;
	delete [] pivots;

	dim = -1;
	numNonzeros = -1;
	neig = -1;
	rank = -1;
	pivots = 0;

	a_ma57 = 0;
	irn_ma57 = 0;
	jcn_ma57 = 0;

	fact_ma57 = 0;
	lfact_ma57 = -1;
	ifact_ma57 = 0;
	lifact_ma57 = -1;

	have_factorization = false;
	return SUCCESSFUL_RETURN;
}


/*
 *	c o p y
 */
returnValue Ma57SparseSolver::copy( 	const Ma57SparseSolver& rhs
										)
{
	dim = rhs.dim;
	numNonzeros = rhs.numNonzeros;
	neig = rhs.neig;
	rank = rhs.rank;
	have_factorization = rhs.have_factorization;

	if ( rhs.a_ma57 != 0 )
	{
		a_ma57 = new double[numNonzeros];
		memcpy( a_ma57,rhs.a_ma57,numNonzeros*sizeof(double) );
	}
	else
		a_ma57 = 0;

	if ( rhs.irn_ma57 != 0 )
	{
		irn_ma57 = new fint_t[numNonzeros];
		memcpy( irn_ma57,rhs.irn_ma57,numNonzeros*sizeof(fint_t) );
	}
	else
		irn_ma57 = 0;

	if ( rhs.jcn_ma57 != 0 )
	{
		jcn_ma57 = new fint_t[numNonzeros];
		memcpy( jcn_ma57,rhs.jcn_ma57,numNonzeros*sizeof(fint_t) );
	}
	else
		jcn_ma57 = 0;

	for ( int_t i=0; i<30; ++i)
		icntl_ma57[i] = rhs.icntl_ma57[i];

	for ( int_t i=0; i<5; ++i)
		cntl_ma57[i] = rhs.cntl_ma57[i];

	lfact_ma57 = rhs.lfact_ma57;
	if ( rhs.fact_ma57 != 0 )
	{
		fact_ma57 = new double[lfact_ma57];
		memcpy( fact_ma57,rhs.fact_ma57,lfact_ma57*sizeof(double) );
	}
	else
		fact_ma57 = 0;

	lifact_ma57 = rhs.lifact_ma57;
	if ( rhs.ifact_ma57 != 0 )
	{
		ifact_ma57 = new fint_t[lifact_ma57];
		memcpy( ifact_ma57,rhs.ifact_ma57,lifact_ma57*sizeof(fint_t) );
	}
	else
		ifact_ma57 = 0;

	if ( have_factorization )
	{
		pivots = new fint_t[dim];
		memcpy( pivots, rhs.pivots, dim*sizeof(fint_t) );
	}
	else
		pivots = 0;

	return SUCCESSFUL_RETURN;
}

#endif /* SOLVER_MA57 */


#ifdef SOLVER_NONE

returnValue DummySparseSolver::setMatrixData( 	int_t dim, /**< Dimension of the linear system. */
												int_t numNonzeros, /**< Number of nonzeros in the matrix. */
												const int_t* const airn, /**< Row indices for each matrix entry. */
												const int_t* const acjn, /**< Column indices for each matrix entry. */
												const real_t* const avals /**< Values for each matrix entry. */
												)
{
	return THROWERROR(RET_NO_SPARSE_SOLVER);
}

returnValue DummySparseSolver::factorize( )
{
	return THROWERROR(RET_NO_SPARSE_SOLVER);
}

returnValue DummySparseSolver::solve(	int_t dim, /**< Dimension of the linear system. */
										const real_t* const rhs, /**< Values for the right hand side. */
										real_t* const sol /**< Solution of the linear system. */
										)
{
	return THROWERROR(RET_NO_SPARSE_SOLVER);
}

#endif /* SOLVER_NONE */

END_NAMESPACE_QPOASES


/*
 *	end of file
 */
