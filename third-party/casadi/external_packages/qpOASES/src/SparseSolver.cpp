/*
 *  This file is part of qpOASES.
 *
 *  qpOASES -- An Implementation of the Online Active Set Strategy.
 *
 *  Copyright (C) 2012 by Andreas Waechter. All rights reserved.
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
 *  \file src/SparseSolver.cpp
 *  \author Andreas Waechter, Dennis Janka
 *  \version 3.2
 *  \date 2012-2015
 *
 *  Interfaces to sparse linear solvers that are used in a Schur-complement
 *  implementation in qpOASES.
 */


#include <qpOASES/SparseSolver.hpp>

#ifndef __MATLAB__
# include <cstdarg>
void MyPrintf(const char* pformat, ... );
#else
# include <mex.h>
# define MyPrintf mexPrintf
#endif

//#define __DEBUG_ITER__

BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 ****************************************************************************/


/*
 *  S p a r s e S o l v e r B a s e
 */
SparseSolver::SparseSolver( )
{
}


/*
 *  S p a r s e S o l v e r B a s e
 */
SparseSolver::SparseSolver( const SparseSolver& rhs )
{
    copy( rhs );
}


/*
 *  ~ S p a r s e S o l v e r B a s e
 */
SparseSolver::~SparseSolver( )
{
    clear( );
}


/*
 *  o p e r a t o r =
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
 *  r e s e t
 */
returnValue SparseSolver::reset( )
{
    return SUCCESSFUL_RETURN;
}

/*
 *  g e t N e g a t i v e E i g e n v a l u e s
 */
int_t SparseSolver::getNegativeEigenvalues( )
{
    return -1;
}

/*
 *  g e t R a n k
 */
int_t SparseSolver::getRank( )
{
    return -1;
}

/*
 *  g e t Z e r o P i v o t s
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
 *  c l e a r
 */
returnValue SparseSolver::clear( )
{
    return SUCCESSFUL_RETURN;
}


/*
 *  c o p y
 */
returnValue SparseSolver::copy(     const SparseSolver& rhs
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
  void MA27ID(fint* ICNTL, double* CNTL);
  void MA27AD(fint *N, fint *NZ, const fint *IRN, const fint* ICN,
                               fint *IW, fint* LIW, fint* IKEEP, fint *IW1,
                               fint* NSTEPS, fint* IFLAG, fint* ICNTL,
                               double* CNTL, fint *INFO, double* OPS);
  void MA27BD(fint *N, fint *NZ, const fint *IRN, const fint* ICN,
                               double* A, fint* LA, fint* IW, fint* LIW,
                               fint* IKEEP, fint* NSTEPS, fint* MAXFRT,
                               fint* IW1, fint* ICNTL, double* CNTL,
                               fint* INFO);
  void MA27CD(fint *N, double* A, fint* LA, fint* IW,
                               fint* LIW, double* W, fint* MAXFRT,
                               double* RHS, fint* IW1, fint* NSTEPS,
                               fint* ICNTL, double* CNTL);
}

/*****************************************************************************
 *  P U B L I C                                                              *
 ****************************************************************************/


/*
 *  S p a r s e S o l v e r B a s e
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
    icntl_ma27[0] = 0;       // Suppress error messages
    icntl_ma27[1] = 0;       // Suppress diagnostic messages
    cntl_ma27[0] = 1e-8;     // Set pivot tolerance
}


/*
 *  S p a r s e S o l v e r B a s e
 */
Ma27SparseSolver::Ma27SparseSolver( const Ma27SparseSolver& rhs )
{
    copy( rhs );
}


/*
 *  ~ S p a r s e S o l v e r B a s e
 */
Ma27SparseSolver::~Ma27SparseSolver( )
{
    clear( );
}


/*
 *  o p e r a t o r =
 */
Ma27SparseSolver& Ma27SparseSolver::operator=( const SparseSolver& rhs )
{
    const Ma27SparseSolver* ma27_rhs = dynamic_cast<const Ma27SparseSolver*>(&rhs);
    if (!ma27_rhs)
    {
        fprintf(getGlobalMessageHandler()->getOutputFile(),"Error in Ma27SparseSolver& Ma27SparseSolver::operator=( const SparseSolver& rhs )\n");
        throw; // TODO: More elegant exit?
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
 *  s e t M a t r i x D a t a
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
        irn_ma27 = new fint[numNonzeros_];
        jcn_ma27 = new fint[numNonzeros_];

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
 *  f a c t o r i z e
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

    // Overstimation factor for LIW (20% recommended in MA27 documentation)
    const double LiwFact = 2.0;   // This is 200% overestimation
    liw_ma27 = (fint)(LiwFact*(double(2*numNonzeros+3*dim+1)));
    iw_ma27 = new fint[liw_ma27];

    ikeep_ma27 = new fint[3*dim];

    fint iflag_ma27 = 0;
    double ops_ma27;
    fint info_ma27[20];
    fint* iw1_ma27 = new fint[2*dim];
    MA27AD(&dim, &numNonzeros, irn_ma27, jcn_ma27, iw_ma27, &liw_ma27, ikeep_ma27,
                            iw1_ma27, &nsteps_ma27, &iflag_ma27, icntl_ma27, cntl_ma27,
                            info_ma27, &ops_ma27);

    /* Receive some information from MA27AD */
    fint iflag = info_ma27[0];   // Information flag
    fint ierror = info_ma27[1];  // Error flag
    fint nrlnec = info_ma27[4];  // recommended value for la
    fint nirnec = info_ma27[5];  // recommended value for liw
    if (iflag != 0)
    {
        MyPrintf("MA27AD returns iflag = %d with ierror = %d\n", iflag, ierror);
        delete [] iw1_ma27;
        clear( );
        return THROWERROR(RET_MATRIX_FACTORISATION_FAILED);
    }

    /* Allocate memory for actual factorization */
    delete [] iw_ma27;
    double liw_init_factor = 5.0; // This could be an option.
    liw_ma27 = (fint)(liw_init_factor * (double)(nirnec));
    iw_ma27 = new fint[liw_ma27];

    double la_init_factor = 20.0; // This could be an option.
    la_ma27 = getMax(numNonzeros,(fint)(la_init_factor * (double)(nrlnec)));
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
    iflag = info_ma27[0];   // Information flag
    ierror = info_ma27[1];  // Error flag
    neig = info_ma27[14];   // Number of negative eigenvalues
    if (iflag == 3)
    {
        rank = info_ma27[1];
        return RET_KKT_MATRIX_SINGULAR;
    }
    else if (iflag == -5)
    { //DJ: I think this is more severe. Can this actually happen?
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
 *  s o l v e
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
    fint* iw1_ma27 = new fint[nsteps_ma27];

    /* MA27CD overwrites rhs */
    for (int_t i=0; i<dim; ++i) sol[i] = rhs[i];
    MA27CD(&dim, a_ma27, &la_ma27, iw_ma27, &liw_ma27, w_ma27, &maxfrt_ma27,
           sol, iw1_ma27, &nsteps_ma27, icntl_ma27, cntl_ma27);

    delete [] iw1_ma27;
    delete [] w_ma27;

    return SUCCESSFUL_RETURN;
}

/*
 *  r e s e t
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
 *  g e t N e g a t i v e E i g e n v a l u e s
 */
int_t Ma27SparseSolver::getNegativeEigenvalues( )
{
    if( !have_factorization )
        return -1;
    else
        return neig;
}

/*
 *  g e t R a n k
 */
int_t Ma27SparseSolver::getRank( )
{
    return rank;
}

/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/

/*
 *  c l e a r
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
 *  c o p y
 */
returnValue Ma27SparseSolver::copy(     const Ma27SparseSolver& rhs
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
        irn_ma27 = new fint[numNonzeros];
        memcpy( irn_ma27,rhs.irn_ma27,numNonzeros*sizeof(fint) );
    }
    else
        irn_ma27 = 0;

    if ( rhs.jcn_ma27 != 0 )
    {
        jcn_ma27 = new fint[numNonzeros];
        memcpy( jcn_ma27,rhs.jcn_ma27,numNonzeros*sizeof(fint) );
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
        iw_ma27 = new fint[liw_ma27];
        memcpy( iw_ma27,rhs.iw_ma27,liw_ma27*sizeof(fint) );
    }
    else
        iw_ma27 = 0;

    if ( rhs.ikeep_ma27 != 0 )
    {
        ikeep_ma27 = new fint[3*dim];
        memcpy( ikeep_ma27,rhs.ikeep_ma27,3*dim*sizeof(fint) );
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

#endif // SOLVER_MA27

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
        double  *cntl,
        fint        *icntl);

    /*
    *  MA57AD -- Symbolic Factorization.
    */
    extern void  MA57AD (
        fint        *n,     /* Order of matrix. */
        fint        *ne,    /* Number of entries. */
        const fint  *irn,   /* Matrix nonzero row structure */
        const fint  *jcn,   /* Matrix nonzero column structure */
        fint        *lkeep, /* Workspace for the pivot order of lenght 3*n */
        fint        *keep,  /* Workspace for the pivot order of lenght 3*n */
                            /* Automatically iflag = 0; ikeep pivot order iflag = 1 */
        fint        *iwork, /* Integer work space. */
        fint        *icntl, /* Integer Control parameter of length 30*/
        fint        *info,  /* Statistical Information; Integer array of length 20 */
        double      *rinfo);/* Double Control parameter of length 5 */

    /*
    * MA57BD -- Numerical Factorization.
    */
    extern void  MA57BD (
        fint    *n,         /* Order of matrix. */
        fint    *ne,        /* Number of entries. */
        double  *a,         /* Numerical values. */
        double  *fact,      /* Entries of factors. */
        fint    *lfact,     /* Length of array `fact'. */
        fint    *ifact,     /* Indexing info for factors. */
        fint    *lifact,    /* Length of array `ifact'. */
        fint    *lkeep,     /* Length of array `keep'. */
        fint    *keep,      /* Integer array. */
        fint    *iwork,     /* Workspace of length `n'. */
        fint    *icntl,     /* Integer Control parameter of length 20. */
        double  *cntl,      /* Double Control parameter of length 5. */
        fint    *info,      /* Statistical Information; Integer array of length 40. */
        double  *rinfo);    /* Statistical Information; Real array of length 20. */

    /*
    * MA57CD -- Solution.
    */
    extern void  MA57CD (
        fint    *job,       /* Solution job.  Solve for... */
                            /* JOB <= 1:  A */
                            /* JOB == 2:  PLP^t */
                            /* JOB == 3:  PDP^t */
                            /* JOB >= 4:  PL^t P^t */
        fint    *n,         /* Order of matrix. */
        double  *fact,      /* Entries of factors. */
        fint    *lfact,     /* Length of array `fact'. */
        fint    *ifact,     /* Indexing info for factors. */
        fint    *lifact,    /* Length of array `ifact'. */
        fint    *nrhs,      /* Number of right hand sides. */
        double  *rhs,       /* Numerical Values. */
        fint    *lrhs,      /* Leading dimensions of `rhs'. */
        double  *work,      /* Real workspace. */
        fint    *lwork,     /* Length of `work', >= N*NRHS. */
        fint    *iwork,     /* Integer array of length `n'. */
        fint    *icntl,     /* Integer Control parameter array of length 20. */
        fint    *info);     /* Statistical Information; Integer array of length 40. */
}

/*****************************************************************************
 *  P U B L I C                                                              *
 ****************************************************************************/


/*
 *  S p a r s e S o l v e r B a s e
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

    icntl_ma57[0] = -1;         // Suppress error messages
    icntl_ma57[1] = -1;         // Suppress warning messages
    icntl_ma57[2] = -1;         // Suppress monitoring messages
    //icntl_ma57[4] = 4;        // Print everything (for debugging)
    icntl_ma57[15] = 1;         // Place small pivots at the end of the factorization (default: 0)

    /// \todo good default values?
    //cntl_ma57[1] = 5.0e-16;   // Pivots smaller than this are treated as zero and are placed at the end of the factorization (default: 1e-20)
    //cntl_ma57[0] = 0.5;       // Set pivot tolerance: Higher values = more stable but slower/less sparse (default: 0.01, max 0.5)
}


/*
 *  S p a r s e S o l v e r B a s e
 */
Ma57SparseSolver::Ma57SparseSolver( const Ma57SparseSolver& rhs )
{
    copy( rhs );
}


/*
 *  ~ S p a r s e S o l v e r B a s e
 */
Ma57SparseSolver::~Ma57SparseSolver( )
{
    clear( );
}


/*
 *  o p e r a t o r =
 */
Ma57SparseSolver& Ma57SparseSolver::operator=( const SparseSolver& rhs )
{
    const Ma57SparseSolver* ma57_rhs = dynamic_cast<const Ma57SparseSolver*>(&rhs);
    if (!ma57_rhs)
    {
        fprintf(getGlobalMessageHandler()->getOutputFile(),"Error in Ma57SparseSolver& Ma57SparseSolver::operator=( const SparseSolver& rhs )\n");
        throw; // TODO: More elegant exit?
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
 *  s e t M a t r i x D a t a
 */
returnValue Ma57SparseSolver::setMatrixData(    int_t dim_,
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
        irn_ma57 = new fint[numNonzeros_];
        jcn_ma57 = new fint[numNonzeros_];

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
 *  f a c t o r i z e
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

    fint lkeep_ma57 = 5*dim + numNonzeros + getMax(numNonzeros,dim) + 42;
    fint *keep_ma57 = new fint[lkeep_ma57];
    fint *iwork_ma57 = new fint[5*dim];

    fint info_ma57[40];
    double rinfo_ma57[20];

    MA57AD(&dim, &numNonzeros, irn_ma57, jcn_ma57, &lkeep_ma57, keep_ma57,
            iwork_ma57, icntl_ma57, info_ma57, rinfo_ma57);

    /* Receive some information from MA57AD */
    fint iflag = info_ma57[0];   // Information flag
    fint ierror = info_ma57[1];  // Error flag
    if (iflag != 0)
    {
        MyPrintf("MA57AD returns iflag = %d with ierror = %d\n", iflag, ierror);
        delete [] keep_ma57;
        delete [] iwork_ma57;
        clear( );
        return THROWERROR(RET_MATRIX_FACTORISATION_FAILED);
    }

    /* Allocate memory for actual factorization */
    double lfact_factor = 10.0; // This could be an option

    lfact_ma57 = (fint)(lfact_factor * (double)(info_ma57[8]));
    fact_ma57 = new double[lfact_ma57];

    lifact_ma57 = (fint)(lfact_factor * (double)(info_ma57[9]));
    ifact_ma57 = new int_t[lifact_ma57];

    /*******************************************
     * Call MA57BD for numerical factorization *
     *******************************************/

    MA57BD( &dim, &numNonzeros, a_ma57, fact_ma57, &lfact_ma57,
            ifact_ma57, &lifact_ma57, &lkeep_ma57, keep_ma57,
            iwork_ma57, icntl_ma57, cntl_ma57, info_ma57, rinfo_ma57 );

    delete [] iwork_ma57;
    delete [] keep_ma57;

    /* Receive some information from MA57BD */
    iflag = info_ma57[0];   // Information flag
    ierror = info_ma57[1];  // Error flag
    neig = info_ma57[23];   // Number of negative eigenvalues
    rank = info_ma57[24];   // Rank of matrix

    /* Read pivot sequence (see MA57UD source code) */
    pivots = new fint[dim];
    fint nrows, ncols;
    fint nblk = ifact_ma57[2];
    int_t iwpos = 3; // = 4-1
    int_t k, kk, count = 0;
    for ( k=0; k<nblk; k++ )
    {
        ncols = ifact_ma57[iwpos];
        nrows = ifact_ma57[iwpos+1];

        for ( kk=0; kk<nrows; kk++ )
            pivots[count++] = ifact_ma57[iwpos+2+kk]-1; //convert Fortran to C indices!

        iwpos = iwpos+ncols+2;
    }

    if (iflag == 4)
    {
        //MyPrintf("dim = %i, rank = %i. Pivots: ", dim, rank);
        //for( k=rank; k<dim; k++ )
            //MyPrintf("%i ", pivots[k]);
        //MyPrintf("\n");

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
 *  s o l v e
 */
returnValue Ma57SparseSolver::solve(    int_t dim_,
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
    fint job_ma57 = 1;
    fint nrhs_ma57 = 1;
    fint lrhs_ma57 = dim;
    fint info_ma57[40];

    fint lwork_ma57 = dim*nrhs_ma57;
    double* work_ma57 = new double[lwork_ma57];
    fint* iwork_ma57 = new fint[dim];

    /* MA57CD overwrites rhs */
    for (int_t i=0; i<dim; ++i) sol[i] = rhs[i];
    MA57CD(&job_ma57, &dim, fact_ma57, &lfact_ma57, ifact_ma57, &lifact_ma57,
            &nrhs_ma57, sol, &lrhs_ma57, work_ma57, &lwork_ma57, iwork_ma57,
            icntl_ma57, info_ma57);

    delete [] work_ma57;
    delete [] iwork_ma57;

    fint iflag = info_ma57[0];   // Information flag
    fint ierror = info_ma57[1];  // Error flag
    if (iflag != 0)
    {
        MyPrintf("MA57CD returns iflag = %d with ierror = %d\n", iflag, ierror);
        clear( );
        return THROWERROR(RET_MATRIX_FACTORISATION_FAILED);
    }

    return SUCCESSFUL_RETURN;
}

/*
 *  r e s e t
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
 *  g e t N e g a t i v e E i g e n v a l u e s
 */
int_t Ma57SparseSolver::getNegativeEigenvalues( )
{
    if( !have_factorization )
        return -1;
    else
        return neig;
}

/*
 *  g e t R a n k
 */
int_t Ma57SparseSolver::getRank( )
{
    return rank;
}

/*
 *  g e t Z e r o P i v o t s
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
 *  c l e a r
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
 *  c o p y
 */
returnValue Ma57SparseSolver::copy(     const Ma57SparseSolver& rhs
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
        irn_ma57 = new fint[numNonzeros];
        memcpy( irn_ma57,rhs.irn_ma57,numNonzeros*sizeof(fint) );
    }
    else
        irn_ma57 = 0;

    if ( rhs.jcn_ma57 != 0 )
    {
        jcn_ma57 = new fint[numNonzeros];
        memcpy( jcn_ma57,rhs.jcn_ma57,numNonzeros*sizeof(fint) );
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
        ifact_ma57 = new fint[lifact_ma57];
        memcpy( ifact_ma57,rhs.ifact_ma57,lifact_ma57*sizeof(fint) );
    }
    else
        ifact_ma57 = 0;

    if ( have_factorization )
    {
        pivots = new fint[dim];
        memcpy( pivots, rhs.pivots, dim*sizeof(fint) );
    }
    else
        pivots = 0;

    return SUCCESSFUL_RETURN;
}

#endif // SOLVER_MA57


#ifdef SOLVER_NONE

UserSparseSolver::UserSparseSolver(linsol_memory_t _linsol_data,
                                     linsol_init_t _linsol_init,
                                     linsol_sfact_t _linsol_sfact,
                                     linsol_nfact_t _linsol_nfact,
                                     linsol_solve_t _linsol_solve) :
  linsol_data(_linsol_data),
  linsol_init(_linsol_init),
  linsol_sfact(_linsol_sfact),
  linsol_nfact(_linsol_nfact),
  linsol_solve(_linsol_solve)
{
  dim = 0;
  nnz = 0;
  allocated_nnz = 0;
  val = 0;
  row = 0;
  col = 0;
  neig = -1;
  rank = 0;
}

UserSparseSolver::~UserSparseSolver()
{
  if (row) delete[] row;
  if (col) delete[] col;
  if (val) delete[] val;
}

returnValue UserSparseSolver::setMatrixData(   int_t dim, /**< Dimension of the linear system. */
                                                int_t nnz, /**< Number of nonzeros in the matrix. */
                                                const int_t* const airn, /**< Row indices for each matrix entry. */
                                                const int_t* const acjn, /**< Column indices for each matrix entry. */
                                                const real_t* const avals /**< Values for each matrix entry. */
                                                )
{
  // Reset linear solver
  reset();

  // Trivial return
  this->dim = dim;
  if (dim==0) return SUCCESSFUL_RETURN;

  // No user-defined linear solver
  if (linsol_init==0) return THROWERROR(RET_NO_SPARSE_SOLVER);

  // Count actual nonzeros
  this->nnz = 0;
  for (int_t i=0; i<nnz; ++i) if (avals[i]!=0) this->nnz++;

  // Allocate more memory if needed
  if (this->nnz > allocated_nnz) {
    // Free existing memory
    if (row) delete[] row;
    if (col) delete[] col;
    if (val) delete[] val;

    // Allocate new memory (2x factor to avoid frequent reallocation)
    allocated_nnz = 2*this->nnz;
    row = new int_t[allocated_nnz];
    col = new int_t[allocated_nnz];
    val = new double[allocated_nnz];
  }

  // Save nonzeros
  int_t k = 0;
  for (int_t i=0; i<nnz; ++i) {
    if (avals[i]!=0) {
      row[k] = airn[i];
      col[k] = acjn[i];
      val[k] = avals[i];
      k++;
    }
  }

  // Call initialization function
  if (linsol_init(linsol_data, dim, this->nnz, row, col)) {
    return THROWERROR(RET_MATRIX_FACTORISATION_FAILED);
  }

  // Number of eigenvalues and rank not available
  neig = -1;
  rank = 0;

  return SUCCESSFUL_RETURN;
}

returnValue UserSparseSolver::factorize( )
{
  // Trivial return
  if (dim==0) {
    neig = 0;
    rank = 0;
    return SUCCESSFUL_RETURN;
  }

  // Symbolic factorization (if any)
  if (linsol_sfact) {
    if (linsol_sfact(linsol_data, val)) {
      return THROWERROR(RET_MATRIX_FACTORISATION_FAILED);
    }
  }

  // No user-defined linear solver
  if (linsol_nfact==0) return THROWERROR(RET_NO_SPARSE_SOLVER);

  // Numerical factorization
  if (linsol_nfact(linsol_data, val, &neig, &rank)) {
    rank = 0;
    neig = -1;
    return THROWERROR(RET_MATRIX_FACTORISATION_FAILED);
  }

  // Is the matrix singular?
  if (rank<dim) return RET_KKT_MATRIX_SINGULAR;

  return SUCCESSFUL_RETURN;
}

int_t UserSparseSolver::getNegativeEigenvalues( ) {
  return neig;
}

int UserSparseSolver::getRank( ) {
  return rank;
}

returnValue UserSparseSolver::solve(   int_t dim, /**< Dimension of the linear system. */
                                        const real_t* const rhs, /**< Values for the right hand side. */
                                        real_t* const sol /**< Solution of the linear system. */
                                        )
{
  // Consistency check
  if (dim!=this->dim) return THROWERROR( RET_INVALID_ARGUMENTS );

  // Trivial return
  if (dim==0) return SUCCESSFUL_RETURN;

  // No user-defined linear solver
  if (linsol_solve==0) return THROWERROR(RET_NO_SPARSE_SOLVER);

  // Solution overwrites the right-hand-side
  for (int_t i=0; i<dim; ++i) sol[i] = rhs[i];

  // Linear system solve
  if (linsol_solve(linsol_data, 1, sol)) {
    return THROWERROR(RET_MATRIX_FACTORISATION_FAILED);
  }

  return SUCCESSFUL_RETURN;
}

#endif // SOLVER_NONE

END_NAMESPACE_QPOASES


/*
 *  end of file
 */
