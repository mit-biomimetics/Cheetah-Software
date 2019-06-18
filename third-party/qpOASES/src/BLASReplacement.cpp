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
 *	\file src/BLASReplacement.cpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.2
 *	\date 2007-2017
 *
 *	BLAS Level 3 replacement routines.
 */


#include <qpOASES/Utils.hpp>


extern "C" void dgemm_(	const char* TRANSA, const char* TRANSB,
						const la_uint_t* M, const la_uint_t* N, const la_uint_t* K,
						const double* ALPHA, const double* A, const la_uint_t* LDA, const double* B, const la_uint_t* LDB,
						const double* BETA, double* C, const la_uint_t* LDC
						)
{
	la_uint_t i, j, k;

	if ( REFER_NAMESPACE_QPOASES isZero(*BETA) == REFER_NAMESPACE_QPOASES BT_TRUE )
		for (k = 0; k < *N; k++)
			for (j = 0; j < *M; j++)
				C[j+(*LDC)*k] = 0.0;
	else if ( REFER_NAMESPACE_QPOASES isEqual(*BETA,-1.0) == REFER_NAMESPACE_QPOASES BT_TRUE )
		for (k = 0; k < *N; k++)
			for (j = 0; j < *M; j++)
				C[j+(*LDC)*k] = -C[j+(*LDC)*k];
	else if ( REFER_NAMESPACE_QPOASES isEqual(*BETA,1.0) == REFER_NAMESPACE_QPOASES BT_FALSE )
		for (k = 0; k < *N; k++)
			for (j = 0; j < *M; j++)
				C[j+(*LDC)*k] *= *BETA;

	if (TRANSA[0] == 'N')
		if ( REFER_NAMESPACE_QPOASES isEqual(*ALPHA,1.0) == REFER_NAMESPACE_QPOASES BT_TRUE )
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] += A[j+(*LDA)*i] * B[i+(*LDB)*k];
		else if ( REFER_NAMESPACE_QPOASES isEqual(*ALPHA,-1.0) == REFER_NAMESPACE_QPOASES BT_TRUE )
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] -= A[j+(*LDA)*i] * B[i+(*LDB)*k];
		else
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] += *ALPHA * A[j+(*LDA)*i] * B[i+(*LDB)*k];
	else
		if ( REFER_NAMESPACE_QPOASES isEqual(*ALPHA,1.0) == REFER_NAMESPACE_QPOASES BT_TRUE )
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] += A[i+(*LDA)*j] * B[i+(*LDB)*k];
		else if ( REFER_NAMESPACE_QPOASES isEqual(*ALPHA,-1.0) == REFER_NAMESPACE_QPOASES BT_TRUE )
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] -= A[i+(*LDA)*j] * B[i+(*LDB)*k];
		else
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] += *ALPHA * A[i+(*LDA)*j] * B[i+(*LDB)*k];
}

extern "C" void sgemm_(	const char* TRANSA, const char* TRANSB,
						const la_uint_t* M, const la_uint_t* N, const la_uint_t* K,
						const float* ALPHA, const float* A, const la_uint_t* LDA, const float* B, const la_uint_t* LDB,
						const float* BETA, float* C, const la_uint_t* LDC
						)
{
	la_uint_t i, j, k;

	if ( REFER_NAMESPACE_QPOASES isZero(*BETA) == REFER_NAMESPACE_QPOASES BT_TRUE )
		for (k = 0; k < *N; k++)
			for (j = 0; j < *M; j++)
				C[j+(*LDC)*k] = 0.0;
	else if ( REFER_NAMESPACE_QPOASES isEqual(*BETA,-1.0) == REFER_NAMESPACE_QPOASES BT_TRUE )
		for (k = 0; k < *N; k++)
			for (j = 0; j < *M; j++)
				C[j+(*LDC)*k] = -C[j+(*LDC)*k];
	else if ( REFER_NAMESPACE_QPOASES isEqual(*BETA,1.0) == REFER_NAMESPACE_QPOASES BT_FALSE )
		for (k = 0; k < *N; k++)
			for (j = 0; j < *M; j++)
				C[j+(*LDC)*k] *= *BETA;

	if (TRANSA[0] == 'N')
		if ( REFER_NAMESPACE_QPOASES isEqual(*ALPHA,1.0) == REFER_NAMESPACE_QPOASES BT_TRUE )
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] += A[j+(*LDA)*i] * B[i+(*LDB)*k];
		else if ( REFER_NAMESPACE_QPOASES isEqual(*ALPHA,-1.0) == REFER_NAMESPACE_QPOASES BT_TRUE )
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] -= A[j+(*LDA)*i] * B[i+(*LDB)*k];
		else
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] += *ALPHA * A[j+(*LDA)*i] * B[i+(*LDB)*k];
	else
		if ( REFER_NAMESPACE_QPOASES isEqual(*ALPHA,1.0) == REFER_NAMESPACE_QPOASES BT_TRUE )
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] += A[i+(*LDA)*j] * B[i+(*LDB)*k];
		else if ( REFER_NAMESPACE_QPOASES isEqual(*ALPHA,-1.0) == REFER_NAMESPACE_QPOASES BT_TRUE )
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] -= A[i+(*LDA)*j] * B[i+(*LDB)*k];
		else
			for (k = 0; k < *N; k++)
				for (j = 0; j < *M; j++)
					for (i = 0; i < *K; i++)
						C[j+(*LDC)*k] += *ALPHA * A[i+(*LDA)*j] * B[i+(*LDB)*k];
}
