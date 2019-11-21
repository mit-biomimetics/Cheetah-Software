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
 *	\file src/LAPACKReplacement.cpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.2
 *	\date 2007-2017
 *
 *  LAPACK replacement routines.
 */


#include <qpOASES/Utils.hpp>


extern "C" void dpotrf_(	const char* uplo, const la_uint_t* _n, double* a,
							const la_uint_t* _lda, la_int_t* info
							)
{
	double sum;
	la_int_t i, j, k;
	la_int_t n = (la_int_t)(*_n);
	la_int_t lda = (la_int_t)(*_lda);

	for( i=0; i<n; ++i )
	{
		/* j == i */
		sum = a[i + lda*i];

		for( k=(i-1); k>=0; --k )
			sum -= a[k+lda*i] * a[k+lda*i];

		if ( sum > 0.0 )
			a[i+lda*i] = REFER_NAMESPACE_QPOASES getSqrt( sum );
		else
		{
			a[0] = sum; /* tunnel negative diagonal element to caller */
			if (info != 0)
				*info = (la_int_t)i+1;
			return;
		}

		for( j=(i+1); j<n; ++j )
		{
			sum = a[j*lda + i];

			for( k=(i-1); k>=0; --k )
				sum -= a[k+lda*i] * a[k+lda*j];

			a[i+lda*j] = sum / a[i+lda*i];
		}
	}
	if (info != 0)
		*info = 0;
}


extern "C" void spotrf_(	const char* uplo, const la_uint_t* _n, float* a,
							const la_uint_t* _lda, la_int_t* info
							)
{
	float sum;
	la_int_t i, j, k;
	la_int_t n = (la_int_t)(*_n);
	la_int_t lda = (la_int_t)(*_lda);

	for( i=0; i<n; ++i )
	{
		/* j == i */
		sum = a[i + lda*i];

		for( k=(i-1); k>=0; --k )
			sum -= a[k+lda*i] * a[k+lda*i];

		if ( sum > 0.0 )
			a[i+lda*i] = (float)(REFER_NAMESPACE_QPOASES getSqrt( sum ));
		else
		{
			a[0] = sum; /* tunnel negative diagonal element to caller */
			if (info != 0)
				*info = (la_int_t)i+1;
			return;
		}

		for( j=(i+1); j<n; ++j )
		{
			sum = a[j*lda + i];

			for( k=(i-1); k>=0; --k )
				sum -= a[k+lda*i] * a[k+lda*j];

			a[i+lda*j] = sum / a[i+lda*i];
		}
	}
	if (info != 0)
		*info = 0;
}

extern "C" void dtrtrs_(	const char* UPLO, const char* TRANS, const char* DIAG,
							const la_uint_t* N, const la_uint_t* NRHS,
							double* A, const la_uint_t* LDA, double* B, const la_uint_t* LDB, la_int_t* INFO
							)
{
	INFO[0] = ((la_int_t)0xDEADBEEF); /* Dummy. If SQProblemSchur is to be used, system LAPACK must be used */
}

extern "C" void strtrs_(	const char* UPLO, const char* TRANS, const char* DIAG,
							const la_uint_t* N, const la_uint_t* NRHS,
							float* A, const la_uint_t* LDA, float* B, const la_uint_t* LDB, la_int_t* INFO
							)
{
	INFO[0] = ((la_int_t)0xDEADBEEF); /* Dummy. If SQProblemSchur is to be used, system LAPACK must be used */
}

extern "C" void dtrcon_(	const char* NORM, const char* UPLO, const char* DIAG,
							const la_uint_t* N, double* A, const la_uint_t*LDA,
							double* RCOND, double* WORK, const la_uint_t* IWORK, la_int_t* INFO
							)
{
	INFO[0] = ((la_int_t)0xDEADBEEF); /* Dummy. If SQProblemSchur is to be used, system LAPACK must be used */
}

extern "C" void strcon_(	const char* NORM, const char* UPLO, const char* DIAG,
							const la_uint_t* N, float* A, const la_uint_t* LDA,
							float* RCOND, float* WORK, const la_uint_t* IWORK, la_int_t* INFO
							)
{
	INFO[0] = ((la_int_t)0xDEADBEEF); /* Dummy. If SQProblemSchur is to be used, system LAPACK must be used */
}

