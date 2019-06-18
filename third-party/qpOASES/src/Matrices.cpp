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
 *	\file src/Matrices.cpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.2
 *	\date 2007-2017
 *
 *	Implementation of the matrix classes.
 */


#include <qpOASES/Matrices.hpp>
#include <qpOASES/LapackBlasReplacement.hpp>


BEGIN_NAMESPACE_QPOASES



/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/

returnValue Matrix::getSparseSubmatrix(
				const Indexlist* const irows,
				const Indexlist* const icols,
				int_t rowoffset,
				int_t coloffset,
				int_t& numNonzeros,
				int_t* irn,
				int_t* jcn,
				real_t* avals,
				BooleanType only_lower_triangular) const
{
	int_t* rowsNumbers;
	irows->getNumberArray( &rowsNumbers );
	int_t* colsNumbers;
	icols->getNumberArray( &colsNumbers );

	return getSparseSubmatrix(irows->getLength(), rowsNumbers, icols->getLength(), colsNumbers, rowoffset, coloffset, numNonzeros, irn, jcn, avals, only_lower_triangular);
}

returnValue Matrix::getSparseSubmatrix(
				const Indexlist* const irows,
				int_t idx_icol,
				int_t rowoffset,
				int_t coloffset,
				int_t& numNonzeros,
				int_t* irn,
				int_t* jcn,
				real_t* avals,
				BooleanType only_lower_triangular) const
{
	int_t* rowsNumbers;
	irows->getNumberArray( &rowsNumbers );

	return getSparseSubmatrix(irows->getLength(), rowsNumbers, 1, &idx_icol, rowoffset, coloffset, numNonzeros, irn, jcn, avals, only_lower_triangular);
}

returnValue Matrix::getSparseSubmatrix(
				int_t idx_row,
				const Indexlist* const icols,
				int_t rowoffset,
				int_t coloffset,
				int_t& numNonzeros,
				int_t* irn,
				int_t* jcn,
				real_t* avals,
				BooleanType only_lower_triangular) const
{
	int_t* colsNumbers;
	icols->getNumberArray( &colsNumbers );

	return getSparseSubmatrix(1, &idx_row, icols->getLength(), colsNumbers, rowoffset, coloffset, numNonzeros, irn, jcn, avals, only_lower_triangular);
}


DenseMatrix::~DenseMatrix()
{
	if ( needToFreeMemory( ) == BT_TRUE )
		free( );
}

void DenseMatrix::free( )
{
	if (val != 0)
		delete[] val;
	val = 0;
}

Matrix *DenseMatrix::duplicate( ) const
{
	DenseMatrix *dupl = 0;

	if ( needToFreeMemory( ) == BT_TRUE )
	{
		real_t* val_new = new real_t[nRows*nCols];
		memcpy( val_new,val, ((uint_t)(nRows*nCols))*sizeof(real_t) );
		dupl = new DenseMatrix(nRows, nCols, nCols, val_new);
		dupl->doFreeMemory( );
	}
	else
	{
		dupl = new DenseMatrix(nRows, nCols, nCols, val);
	}

	return dupl;
}

real_t DenseMatrix::diag(	int_t i
							) const
{
	return val[i*(leaDim+1)];
}

BooleanType DenseMatrix::isDiag( ) const
{
	int_t i, j;

	if (nRows != nCols)
		return BT_FALSE;

	for ( i=0; i<nRows; ++i )
		for ( j=0; j<i; ++j )
			if ( ( getAbs( val[i*leaDim+j] ) > EPS ) || ( getAbs( val[j*leaDim+i] ) > EPS ) )
				return BT_FALSE;

	return BT_TRUE;
}


real_t DenseMatrix::getNorm(	int_t type
								) const
{
	return REFER_NAMESPACE_QPOASES getNorm( val,nCols*nRows,type );
}


real_t DenseMatrix::getRowNorm( int_t rNum, int_t type ) const
{
	return REFER_NAMESPACE_QPOASES getNorm( &(val[rNum*leaDim]),nCols,type );
}

returnValue DenseMatrix::getRowNorm( real_t* norm, int_t type ) const
{
	int_t i;
	for (i = 0; i < nRows; ++i)
	{
		norm[i] = REFER_NAMESPACE_QPOASES getNorm( &(val[i*leaDim]),nCols,type );
	}
	return SUCCESSFUL_RETURN;
}

returnValue DenseMatrix::getRow(int_t rNum, const Indexlist* const icols, real_t alpha, real_t* row) const
{
	int_t i;
	if (icols != 0)
    {
	    if ( isEqual(alpha,1.0) == BT_TRUE )
		    for (i = 0; i < icols->length; i++)
			    row[i] = val[rNum*leaDim+icols->number[i]];
	    else if ( isEqual(alpha,-1.0) == BT_TRUE )
			for (i = 0; i < icols->length; i++)
				row[i] = -val[rNum*leaDim+icols->number[i]];
		else
			for (i = 0; i < icols->length; i++)
				row[i] = alpha*val[rNum*leaDim+icols->number[i]];
	}
	else
	{
		if ( isEqual(alpha,1.0) == BT_TRUE )
			for (i = 0; i < nCols; i++)
				row[i] = val[rNum*leaDim+i];
		else if ( isEqual(alpha,-1.0) == BT_TRUE )
			for (i = 0; i < nCols; i++)
				row[i] = -val[rNum*leaDim+i];
		else
			for (i = 0; i < nCols; i++)
				row[i] = alpha*val[rNum*leaDim+i];
	}
	return SUCCESSFUL_RETURN;
}


returnValue DenseMatrix::getCol(int_t cNum, const Indexlist* const irows, real_t alpha, real_t* col) const
{
	int_t i;

	if ( isEqual(alpha,1.0) == BT_TRUE )
		for (i = 0; i < irows->length; i++)
			col[i] = val[irows->number[i]*leaDim+cNum];
	else if ( isEqual(alpha,-1.0) == BT_TRUE )
		for (i = 0; i < irows->length; i++)
			col[i] = -val[irows->number[i]*leaDim+cNum];
	else
		for (i = 0; i < irows->length; i++)
			col[i] = alpha*val[irows->number[i]*leaDim+cNum];

	return SUCCESSFUL_RETURN;
}

returnValue DenseMatrix::getSparseSubmatrix (int_t irowsLength, const int_t* const irowsNumber,
											 int_t icolsLength, const int_t* const icolsNumber,
											 int_t rowoffset, int_t coloffset, int_t& numNonzeros,	int_t* irn,
											 int_t* jcn, real_t* avals,
											 BooleanType only_lower_triangular /*= BT_FALSE */) const
{
	int_t i, j, irA;
	real_t v;
	numNonzeros = 0;
	if ( only_lower_triangular == BT_FALSE )
	{
		if (irn == 0)
		{
			if (jcn != 0 || avals != 0)
				return THROWERROR( RET_INVALID_ARGUMENTS );
			for (j = 0; j<irowsLength; j++)
			{
				irA = irowsNumber[j] * leaDim;
				for (i = 0; i<icolsLength; i++)
					if (isZero( val[irA+icolsNumber[i]] ) == BT_FALSE)
						numNonzeros++;
			}
		}
		else
		{
			for (j = 0; j<irowsLength; j++)
			{
				irA = irowsNumber[j] * leaDim;
				for (i = 0; i<icolsLength; i++)
				{
					v = val[irA+icolsNumber[i]];
					if (isZero( v ) == BT_FALSE)
					{
						irn[numNonzeros] = j+rowoffset;
						jcn[numNonzeros] = i+coloffset;
						avals[numNonzeros] = v;
						numNonzeros++;
					}
				}
			}
		}
	}
	else
	{
		if (irn == 0)
		{
			if (jcn != 0 || avals != 0)
				return THROWERROR( RET_INVALID_ARGUMENTS );
			for (j = 0; j<irowsLength; j++)
			{
				irA = irowsNumber[j] * leaDim;
				for (i = 0; i<=j; i++)
					if (isZero( val[irA+irowsNumber[i]] ) == BT_FALSE)
						numNonzeros++;
			}
		}
		else
		{
			for (j = 0; j<irowsLength; j++)
			{
				irA = irowsNumber[j] * leaDim;
				for (i = 0; i<=j; i++)
				{
					v = val[irA+irowsNumber[i]];
					if (isZero( v ) == BT_FALSE)
					{
						irn[numNonzeros] = j+rowoffset;
						jcn[numNonzeros] = i+coloffset;
						avals[numNonzeros] = v;
						numNonzeros++;
					}
				}
			}
		}
	}

	return SUCCESSFUL_RETURN;
}

returnValue DenseMatrix::times(	int_t xN, real_t alpha, const real_t* x, int_t xLD, real_t beta, real_t* y, int_t yLD ) const
{
	la_uint_t _xN     = (la_uint_t)xN;
	la_uint_t _nRows  = (la_uint_t)nRows;
	la_uint_t _nCols  = (la_uint_t)nCols;
	la_uint_t _leaDim = (la_uint_t)getMax(1,nCols);
	la_uint_t _xLD    = (la_uint_t)getMax(1,xLD);
	la_uint_t _yLD    = (la_uint_t)getMax(1,yLD);

	/* Call BLAS. Mind row major format! */
	GEMM( "TRANS", "NOTRANS", &_nRows, &_xN, &_nCols, &alpha, val, &_leaDim, x, &_xLD, &beta, y, &_yLD );
	return SUCCESSFUL_RETURN;
}

returnValue DenseMatrix::transTimes( int_t xN, real_t alpha, const real_t* x, int_t xLD, real_t beta, real_t* y, int_t yLD ) const
{
	la_uint_t _xN     = (la_uint_t)xN;
	la_uint_t _nRows  = (la_uint_t)nRows;
	la_uint_t _nCols  = (la_uint_t)nCols;
	la_uint_t _leaDim = (la_uint_t)getMax(1,nCols);
	la_uint_t _xLD    = (la_uint_t)getMax(1,xLD);
	la_uint_t _yLD    = (la_uint_t)getMax(1,yLD);

	/* Call BLAS. Mind row major format! */
	GEMM( "NOTRANS", "NOTRANS", &_nCols, &_xN, &_nRows, &alpha, val, &_leaDim, x, &_xLD, &beta, y, &_yLD );
	return SUCCESSFUL_RETURN;
}

returnValue DenseMatrix::times(	const Indexlist* const irows, const Indexlist* const icols,
								int_t xN, real_t alpha, const real_t* x, int_t xLD, real_t beta, real_t* y, int_t yLD,
								BooleanType yCompr ) const
{
	int_t i, j, k, row, col, iy, irA;

	if (yCompr == BT_TRUE)
	{
		if ( isZero(beta) == BT_TRUE )
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[j+k*yLD] = 0.0;
		else if ( isEqual(beta,-1.0) == BT_TRUE )
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[j+k*yLD] = -y[j+k*yLD];
		else if ( isEqual(beta,1.0) == BT_FALSE )
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[j+k*yLD] *= beta;

		if (icols == 0)
			if ( isEqual(alpha,1.0) == BT_TRUE )
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->iSort[j];
						iy = row + k * yLD;
						irA = irows->number[row] * leaDim;
						for (i = 0; i < nCols; i++)
							y[iy] += val[irA+i] * x[k*xLD+i];
					}
			else if ( isEqual(alpha,-1.0) == BT_TRUE )
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->iSort[j];
						iy = row + k * yLD;
						irA = irows->number[row] * leaDim;
						for (i = 0; i < nCols; i++)
							y[iy] -= val[irA+i] * x[k*xLD+i];
					}
			else
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->iSort[j];
						iy = row + k * yLD;
						irA = irows->number[row] * leaDim;
						for (i = 0; i < nCols; i++)
							y[iy] += alpha * val[irA+i] * x[k*xLD+i];
					}
		else /* icols != 0 */
			if ( isEqual(alpha,1.0) == BT_TRUE )
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->iSort[j];
						iy = row + k * yLD;
						irA = irows->number[row] * leaDim;
						for (i = 0; i < icols->length; i++)
						{
							col = icols->iSort[i];
							y[iy] += val[irA+icols->number[col]] * x[k*xLD+col];
						}
					}
			else if ( isEqual(alpha,-1.0) == BT_TRUE )
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->iSort[j];
						iy = row + k * yLD;
						irA = irows->number[row] * leaDim;
						for (i = 0; i < icols->length; i++)
						{
							col = icols->iSort[i];
							y[iy] -= val[irA+icols->number[col]] * x[k*xLD+col];
						}
					}
			else
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->iSort[j];
						iy = row + k * yLD;
						irA = irows->number[row] * leaDim;
						for (i = 0; i < icols->length; i++)
						{
							col = icols->iSort[i];
							y[iy] += alpha * val[irA+icols->number[col]] * x[k*xLD+col];
						}
					}
	}
	else /* y not compressed */
	{
		if ( isZero(beta) == BT_TRUE )
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[irows->number[j]+k*yLD] = 0.0;
		else if ( isEqual(beta,-1.0) == BT_TRUE )
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[irows->number[j]+k*yLD] = -y[j+k*yLD];
		else if ( isEqual(beta,1.0) == BT_FALSE )
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[irows->number[j]+k*yLD] *= beta;

		if (icols == 0)
			if ( isEqual(alpha,1.0) == BT_TRUE )
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->number[irows->iSort[j]];
						iy = row + k * yLD;
						irA = row * leaDim;
						for (i = 0; i < nCols; i++)
							y[iy] += val[irA+i] * x[k*xLD+i];
					}
			else if ( isEqual(alpha,-1.0) == BT_TRUE )
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->number[irows->iSort[j]];
						iy = row + k * yLD;
						irA = row * leaDim;
						for (i = 0; i < nCols; i++)
							y[iy] -= val[irA+i] * x[k*xLD+i];
					}
			else
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->number[irows->iSort[j]];
						iy = row + k * yLD;
						irA = row * leaDim;
						for (i = 0; i < nCols; i++)
							y[iy] += alpha * val[irA+i] * x[k*xLD+i];
					}
		else /* icols != 0 */
			if ( isEqual(alpha,1.0) == BT_TRUE )
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->number[irows->iSort[j]];
						iy = row + k * yLD;
						irA = row * leaDim;
						for (i = 0; i < icols->length; i++)
						{
							col = icols->iSort[i];
							y[iy] += val[irA+icols->number[col]] * x[k*xLD+col];
						}
					}
			else if ( isEqual(alpha,-1.0) == BT_TRUE )
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->number[irows->iSort[j]];
						iy = row + k * yLD;
						irA = row * leaDim;
						for (i = 0; i < icols->length; i++)
						{
							col = icols->iSort[i];
							y[iy] -= val[irA+icols->number[col]] * x[k*xLD+col];
						}
					}
			else
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
					{
						row = irows->number[irows->iSort[j]];
						iy = row + k * yLD;
						irA = row * leaDim;
						for (i = 0; i < icols->length; i++)
						{
							col = icols->iSort[i];
							y[iy] += alpha * val[irA+icols->number[col]] * x[k*xLD+col];
						}
					}
	}

	return SUCCESSFUL_RETURN;
}

returnValue DenseMatrix::transTimes(	const Indexlist* const irows, const Indexlist* const icols,
										int_t xN, real_t alpha, const real_t* x, int_t xLD, real_t beta, real_t* y, int_t yLD ) const
{
	int_t i, j, k, row, col;

	if ( isZero(beta) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < icols->length; j++)
				y[j+k*yLD] = 0.0;
	else if ( isEqual(beta,-1.0) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < icols->length; j++)
				y[j+k*yLD] = -y[j+k*yLD];
	else if ( isEqual(beta,1.0) == BT_FALSE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < icols->length; j++)
				y[j+k*yLD] *= beta;

	if ( isEqual(alpha,1.0) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < irows->length; j++)
			{
				row = irows->iSort[j];
				for (i = 0; i < icols->length; i++)
				{
					col = icols->iSort[i];
					y[col+k*yLD] += val[irows->number[row]*leaDim+icols->number[col]] * x[row+k*xLD];
				}
			}
	else if ( isEqual(alpha,-1.0) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < irows->length; j++)
			{
				row = irows->iSort[j];
				for (i = 0; i < icols->length; i++)
				{
					col = icols->iSort[i];
					y[col+k*yLD] -= val[irows->number[row]*leaDim+icols->number[col]] * x[row+k*xLD];
				}
			}
	else
		for (k = 0; k < xN; k++)
			for (j = 0; j < irows->length; j++)
			{
				row = irows->iSort[j];
				for (i = 0; i < icols->length; i++)
				{
					col = icols->iSort[i];
					y[col+k*yLD] += alpha * val[irows->number[row]*leaDim+icols->number[col]] * x[row+k*xLD];
				}
			}

	return SUCCESSFUL_RETURN;
}


returnValue DenseMatrix::addToDiag( real_t alpha )
{
	int_t i;
	for (i = 0; i < nRows && i < nCols; i++)
		val[i*(leaDim+1)] += alpha;

	return SUCCESSFUL_RETURN;
}


returnValue DenseMatrix::writeToFile( FILE* output_file, const char* prefix ) const
{
	return THROWERROR( RET_NOT_YET_IMPLEMENTED );
}


real_t* DenseMatrix::full() const
{
	real_t* v = new real_t[nRows*nCols];
	memcpy( v,val, ((uint_t)(nRows*nCols))*sizeof(real_t) );
	return v;
}


returnValue DenseMatrix::print( const char* name ) const
{
	return REFER_NAMESPACE_QPOASES print( val,nRows,nCols,name );
}



Matrix *SymDenseMat::duplicate( ) const
{
	return duplicateSym();
}


SymmetricMatrix *SymDenseMat::duplicateSym( ) const
{
	/* "same" as duplicate() in DenseMatrix */
	SymDenseMat *dupl = 0;

	if ( needToFreeMemory( ) == BT_TRUE )
	{
		real_t* val_new = new real_t[nRows*nCols];
		memcpy( val_new,val, ((uint_t)(nRows*nCols))*sizeof(real_t) );
		dupl = new SymDenseMat(nRows, nCols, nCols, val_new);
		dupl->doFreeMemory( );
	}
	else
	{
		dupl = new SymDenseMat(nRows, nCols, nCols, val);
	}

	return dupl;
}


returnValue SymDenseMat::bilinear(	const Indexlist* const icols,
									int_t xN, const real_t* x, int_t xLD, real_t* y, int_t yLD ) const
{
	int_t ii, jj, kk, col;
	int_t i,j,k,irA;

	for (ii = 0; ii < xN; ii++)
		for (jj = 0; jj < xN; jj++)
			y[ii*yLD+jj] = 0.0;

	real_t* Ax = new real_t[icols->length * xN];

	for (i=0;i<icols->length * xN;++i)
		Ax[i]=0.0;

	/* exploit symmetry of A ! */
	for (j = 0; j < icols->length; j++) {
		irA = icols->number[j] * leaDim;
		for (i = 0; i < icols->length; i++)
		{
			real_t h = val[irA+icols->number[i]];
			for (k = 0; k < xN; k++)
				Ax[j + k * icols->length] += h * x[k*xLD+icols->number[i]];
		}
	}

	for (ii = 0; ii < icols->length; ++ii) {
		col = icols->number[ii];
		for (jj = 0; jj < xN; ++jj) {
			for (kk = 0; kk < xN; ++kk) {
				y[kk + jj*yLD] += x[col + jj*xLD] * Ax[ii + kk*icols->length];
			}
		}
	}
	delete[] Ax;

	return SUCCESSFUL_RETURN;
}



SparseMatrix::SparseMatrix() : nRows(0), nCols(0), ir(0), jc(0), jd(0), val(0) {}

SparseMatrix::SparseMatrix(	int_t nr, int_t nc, sparse_int_t* r, sparse_int_t* c, real_t* v )
								: nRows(nr), nCols(nc), ir(r), jc(c), jd(0), val(v) { doNotFreeMemory(); }

SparseMatrix::SparseMatrix(	int_t nr, int_t nc, int_t ld, const real_t*  const v )
								: nRows(nr), nCols(nc), jd(0)
{
	int_t i, j, nnz;

	jc = new sparse_int_t[nc+1];
	ir = new sparse_int_t[nr*nc];
	val = new real_t[nr*nc];

	nnz = 0;
	for (j = 0; j < nCols; j++)
	{
		jc[j] = nnz;
		for (i = 0; i < nRows; i++)
			if ( ( isZero( v[i*ld+j],0.0 ) == BT_FALSE ) || ( i == j ) ) /* also include zero diagonal elemets! */
			{
				ir[nnz] = i;
				val[nnz++] = v[i*ld+j];
			}
	}
	jc[nCols] = nnz;

	doFreeMemory( );
}


SparseMatrix::~SparseMatrix()
{
	if (jd != 0)
	{
		delete[] jd;
		jd = 0;
	}

	if ( needToFreeMemory() == BT_TRUE )
		free( );
}


void SparseMatrix::free( )
{
	if (ir != 0) delete[] ir;
	ir = 0;
	if (jc != 0) delete[] jc;
	jc = 0;
	if (val != 0) delete[] val;
	val = 0;

	doNotFreeMemory( );
}


Matrix *SparseMatrix::duplicate( ) const
{
	long i, length = jc[nCols];
	SparseMatrix *dupl = new SparseMatrix;

	dupl->nRows = nRows;
	dupl->nCols = nCols;
	dupl->ir = new sparse_int_t[length];
	dupl->jc = new sparse_int_t[nCols+1];
	dupl->val = new real_t[length];

	for (i = 0; i < length; i++) dupl->ir[i] = ir[i];
	for (i = 0; i <= nCols; i++) dupl->jc[i] = jc[i];
	for (i = 0; i < length; i++) dupl->val[i] = val[i];

	if ( jd != 0 )
	{
		dupl->jd = new sparse_int_t[nCols];
		for (i = 0; i < nCols; i++) dupl->jd[i] = jd[i];
	}
	else
		dupl->jd = 0;

	dupl->doFreeMemory( );

	return dupl;
}


void SparseMatrix::setVal( const real_t* newVal )
{
	long length = jc[nCols]; /* [nCols+1] ?? */
	for (long index = 0; index < length; ++index)
	{
		val[index] = newVal[index];
	}
}


real_t SparseMatrix::diag( int_t i ) const
{
	if ( jd == 0 )
	{
		THROWERROR( RET_DIAGONAL_NOT_INITIALISED );
		return INFTY;
	}

	int_t entry = jd[i];
	return (entry < jc[i+1] && ir[entry] == i) ? val[entry] : 0.0;
}


BooleanType SparseMatrix::isDiag( ) const
{
	int_t j;

	if ( nCols != nRows )
		return BT_FALSE;

	for (j = 0; j < nCols; ++j)
	{
		if ( jc[j+1] > jc[j]+1 )
			return BT_FALSE;

		if ( ( jc[j+1] == jc[j]+1 ) && ( ir[jc[j]] != j ) )
			return BT_FALSE;
	}

	return BT_TRUE;
}



real_t SparseMatrix::getNorm(	int_t type
								) const
{
	int_t length = jc[nCols];
	return REFER_NAMESPACE_QPOASES getNorm( val,length,type );
}


real_t SparseMatrix::getRowNorm( int_t rNum, int_t type ) const
{
	int_t i,j;
	real_t norm = 0.0;

	switch( type )
	{
		case 2:
			for ( j=0; j < nCols; ++j ) {
				for (i = jc[j]; i < jc[j+1] && ir[i] < rNum; i++) {};
				if (i < jc[j+1] && ir[i] == rNum) norm += val[i]*val[i];
			}
			return getSqrt(norm);

		case 1:
			for ( j=0; j < nCols; ++j ) {
				for (i = jc[j]; i < jc[j+1] && ir[i] < rNum; i++) {};
				if (i < jc[j+1] && ir[i] == rNum) norm += getAbs( val[i] );
			}
			return norm;

		default:
			THROWERROR( RET_INVALID_ARGUMENTS );
			return -INFTY;
	}
}

returnValue SparseMatrix::getRowNorm( real_t* norm, int_t type ) const
{
	int_t i,j;

	for ( j=0; j < nRows; ++j ) norm[j] = 0.0;

	switch( type )
	{
		case 2:
			for ( j=0; j < nCols; ++j ) {
				for (i = jc[j]; i < jc[j+1]; i++)
				  norm[ir[i]] += val[i]*val[i];
			}
			for ( j=0; j < nRows; ++j ) norm[j] = getSqrt(norm[j]);
			break;
		case 1:
			for ( j=0; j < nCols; ++j ) {
				for (i = jc[j]; i < jc[j+1]; i++);
				  norm[ir[i]] += getAbs( val[i] );
			}
			break;
		default:
			return RET_INVALID_ARGUMENTS;
	}
	return SUCCESSFUL_RETURN;
}


returnValue SparseMatrix::getRow( int_t rNum, const Indexlist* const icols, real_t alpha, real_t* row ) const
{
	long i, j, k;

	if (icols != 0)
	{
		if ( isEqual(alpha,1.0) == BT_TRUE )
			for (k = 0; k < icols->length; k++)
			{
				j = icols->number[icols->iSort[k]];
				for (i = jc[j]; i < jc[j+1] && ir[i] < rNum; i++);
				row[icols->iSort[k]] = (i < jc[j+1] && ir[i] == rNum) ? val[i] : 0.0;
			}
		else if ( isEqual(alpha,-1.0) == BT_TRUE )
			for (k = 0; k < icols->length; k++)
			{
				j = icols->number[icols->iSort[k]];
				for (i = jc[j]; i < jc[j+1] && ir[i] < rNum; i++);
				row[icols->iSort[k]] = (i < jc[j+1] && ir[i] == rNum) ? -val[i] : 0.0;
			}
		else
			for (k = 0; k < icols->length; k++)
			{
				j = icols->number[icols->iSort[k]];
				for (i = jc[j]; i < jc[j+1] && ir[i] < rNum; i++);
				row[icols->iSort[k]] = (i < jc[j+1] && ir[i] == rNum) ? alpha*val[i] : 0.0;
			}
	}
	else
	{
		if ( isEqual(alpha,1.0) == BT_TRUE )
			for (j = 0; j < nCols; j++)
			{
				for (i = jc[j]; i < jc[j+1] && ir[i] < rNum; i++);
				row[j] = (i < jc[j+1] && ir[i] == rNum) ? val[i] : 0.0;
			}
		else if ( isEqual(alpha,-1.0) == BT_TRUE )
			for (j = 0; j < icols->length; j++)
			{
				for (i = jc[j]; i < jc[j+1] && ir[i] < rNum; i++);
				row[j] = (i < jc[j+1] && ir[i] == rNum) ? -val[i] : 0.0;
			}
		else
			for (j = 0; j < icols->length; j++)
			{
				for (i = jc[j]; i < jc[j+1] && ir[i] < rNum; i++);
				row[j] = (i < jc[j+1] && ir[i] == rNum) ? alpha*val[i] : 0.0;
			}
	}
	return SUCCESSFUL_RETURN;
}


returnValue SparseMatrix::getCol( int_t cNum, const Indexlist* const irows, real_t alpha, real_t* col ) const
{
	long i, j;

	i = jc[cNum];
	j = 0;
	if ( isEqual(alpha,1.0) == BT_TRUE )
		while (i < jc[cNum+1] && j < irows->length)
			if (ir[i] == irows->number[irows->iSort[j]])
				col[irows->iSort[j++]] = val[i++];
			else if (ir[i] > irows->number[irows->iSort[j]])
				col[irows->iSort[j++]] = 0.0;
			else
				i++;
	else if ( isEqual(alpha,-1.0) == BT_TRUE )
		while (i < jc[cNum+1] && j < irows->length)
			if (ir[i] == irows->number[irows->iSort[j]])
				col[irows->iSort[j++]] = -val[i++];
			else if (ir[i] > irows->number[irows->iSort[j]])
				col[irows->iSort[j++]] = 0.0;
			else
				i++;
	else
		while (i < jc[cNum+1] && j < irows->length)
			if (ir[i] == irows->number[irows->iSort[j]])
				col[irows->iSort[j++]] = alpha * val[i++];
			else if (ir[i] > irows->number[irows->iSort[j]])
				col[irows->iSort[j++]] = 0.0;
			else
				i++;

	/* fill in remaining zeros */
	while (j < irows->length)
		col[irows->iSort[j++]] = 0.0;

	return SUCCESSFUL_RETURN;
}

returnValue SparseMatrix::getSparseSubmatrix(	int_t irowsLength, const int_t* const irowsNumber,
												int_t icolsLength, const int_t* const icolsNumber,
												int_t rowoffset, int_t coloffset, int_t& numNonzeros,	int_t* irn,
												int_t* jcn, real_t* avals,
												BooleanType only_lower_triangular /*= BT_FALSE */ ) const
{
	int_t i, j, k, l;

	// Compute the "inverse" of the irows->number array
	// TODO: Ideally this should be a part of Indexlist
	int_t* rowNumberInv = new int_t[nRows];
	for (i=0; i<nRows; i++)
		rowNumberInv[i] = -1;
	for (i=0; i<irowsLength; i++)
		rowNumberInv[irowsNumber[i]] = i;

	numNonzeros = 0;
	if ( only_lower_triangular == BT_FALSE )
	{
		if (irn == 0)
		{
			if (jcn != 0 || avals != 0)
				return THROWERROR( RET_INVALID_ARGUMENTS );
			for (k = 0; k < icolsLength; k++)
			{
				j = icolsNumber[k];
				for (i = jc[j]; i < jc[j+1]; i++)
				{
					l = rowNumberInv[ir[i]];
					if (l >= 0)
						numNonzeros++;
				}
			}
		}
		else
		{
			for (k = 0; k < icolsLength; k++)
			{
				j = icolsNumber[k];
				for (i = jc[j]; i < jc[j+1]; i++)
				{
					l = rowNumberInv[ir[i]];
					if (l >= 0)
					{
						irn[numNonzeros] = l+rowoffset;
						jcn[numNonzeros] = k+coloffset;
						avals[numNonzeros] = val[i];
						numNonzeros++;
					}
				}
			}
		}
	}
	else
	{
		if (irn == 0)
		{
			if (jcn != 0 || avals != 0)
				return THROWERROR( RET_INVALID_ARGUMENTS );
			for (k = 0; k < icolsLength; k++)
			{
				j = icolsNumber[k];
				for (i = jc[j]; i < jc[j+1]; i++)
				{
					l = rowNumberInv[ir[i]];
					if (l >= k)
						numNonzeros++;
				}
			}
		}
		else
		{
			for (k = 0; k < icolsLength; k++)
			{
				j = icolsNumber[k];
				for (i = jc[j]; i < jc[j+1]; i++)
				{
					l = rowNumberInv[ir[i]];
					if (l >= k)
					{
						irn[numNonzeros] = l+rowoffset;
						jcn[numNonzeros] = k+coloffset;
						avals[numNonzeros] = val[i];
						numNonzeros++;
					}
				}
			}
		}
	}
	delete [] rowNumberInv;

	return SUCCESSFUL_RETURN;
}

returnValue SparseMatrix::times( int_t xN, real_t alpha, const real_t* x, int_t xLD,
		real_t beta, real_t* y, int_t yLD ) const
{
	long i, j, k;

	if ( isZero(beta) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < nRows; j++)
				y[j+k*yLD] = 0.0;
	else if ( isEqual(beta,-1.0) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < nRows; j++)
				y[j+k*yLD] = -y[j+k*yLD];
	else if ( isEqual(beta,1.0) == BT_FALSE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < nRows; j++)
				y[j+k*yLD] *= beta;

	if ( isEqual(alpha,1.0) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < nCols; j++)
				for (i = jc[j]; i < jc[j+1]; i++)
					y[ir[i]+k*yLD] += val[i] * x[j+k*xLD];
	else if ( isEqual(alpha,-1.0) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < nCols; j++)
				for (i = jc[j]; i < jc[j+1]; i++)
					y[ir[i]+k*yLD] -= val[i] * x[j+k*xLD];
	else
		for (k = 0; k < xN; k++)
			for (j = 0; j < nCols; j++)
				for (i = jc[j]; i < jc[j+1]; i++)
					y[ir[i]+k*yLD] += alpha * val[i] * x[j+k*xLD];

	return SUCCESSFUL_RETURN;
}


returnValue SparseMatrix::transTimes( int_t xN, real_t alpha, const real_t* x, int_t xLD,
		real_t beta, real_t* y, int_t yLD ) const
{
	long i, j, k;

	if ( isZero(beta) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < nCols; j++)
				y[j+k*yLD] = 0.0;
	else if ( isEqual(beta,-1.0) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < nCols; j++)
				y[j+k*yLD] = -y[j+k*yLD];
	else if ( isEqual(beta,1.0) == BT_FALSE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < nCols; j++)
				y[j+k*yLD] *= beta;

	if ( isEqual(alpha,1.0) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < nCols; j++)
				for (i = jc[j]; i < jc[j+1]; i++)
					y[j+k*yLD] += val[i] * x[ir[i]+k*xLD];
	else if ( isEqual(alpha,-1.0) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < nCols; j++)
				for (i = jc[j]; i < jc[j+1]; i++)
					y[j+k*yLD] -= val[i] * x[ir[i]+k*xLD];
	else
		for (k = 0; k < xN; k++)
			for (j = 0; j < nCols; j++)
				for (i = jc[j]; i < jc[j+1]; i++)
					y[j+k*yLD] += alpha * val[i] * x[ir[i]+k*xLD];

	return SUCCESSFUL_RETURN;
}


returnValue SparseMatrix::times( const Indexlist* const irows, const Indexlist* const icols,
		int_t xN, real_t alpha, const real_t* x, int_t xLD, real_t beta, real_t* y, int_t yLD,
		BooleanType yCompr ) const
{
	long i, j, k, l, col;
	real_t xcol;

	if ( isEqual(alpha,0.0) == BT_TRUE )
	{
		if (yCompr == BT_TRUE)
		{
			if ( isZero(beta) == BT_TRUE )
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
						y[j+k*yLD] = 0.0;
			else if ( isEqual(beta,-1.0) == BT_TRUE )
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
						y[j+k*yLD] = -y[j+k*yLD];
			else if ( isEqual(beta,1.0) == BT_FALSE )
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
						y[j+k*yLD] *= beta;
		}
		else
		{
			if (isZero( beta ) == BT_TRUE)
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
						y[irows->number[j]+k*yLD] = 0.0;
			else if (isEqual( beta, -1.0 ) == BT_TRUE)
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
						y[irows->number[j]+k*yLD] = -y[irows->number[j]+k*yLD];
			else if (isEqual( beta, 1.0 ) == BT_FALSE)
				for (k = 0; k < xN; k++)
					for (j = 0; j < irows->length; j++)
						y[irows->number[j]+k*yLD] *= beta;
		}
		return SUCCESSFUL_RETURN;
	}

	// First, work with full, unordered copy of y and store matrix times x in there
	const int_t yfullLength = nRows;
	real_t* ytmp = new real_t[xN*yfullLength];
	for (k = 0; k < xN*yfullLength; k++)
		ytmp[k] = 0.0;

	if (icols!=0)
	{
		if (xN==1)
		{
			for (l = 0; l < icols->length; l++)
			{
				col = icols->iSort[l];
				xcol = x[col];
				if (isZero( xcol ) == BT_FALSE)
				{
					j = icols->number[col];
					for (i = jc[j]; i < jc[j+1]; i++)
						ytmp[ir[i]] += val[i] * xcol;
				}
			}
		}
		else
		{
			// AW: I didn't test the case xN>1, but I hope it is working
			real_t* xcols = new real_t[xN];
			for (l = 0; l < icols->length; l++)
			{
				col = icols->iSort[l];
				real_t xmax = 0.0;
				for (k=0; k<xN; k++)
				{
					xcols[k] = x[k*xLD+col];
					xmax = getMax(xmax,getAbs(xcols[k]));
				}
				if (isZero( xmax ) == BT_FALSE)
				{
					j = icols->number[col];
					for (i = jc[j]; i < jc[j+1]; i++)
						for (k=0; k<xN; k++)
						  // AW: Maybe it makes more sense to order ytmp by vectors, not vector entries, for better cache peformance?
							ytmp[k*yfullLength+ir[i]] += val[i] * xcols[k];
				}
			}
			delete [] xcols;
		}
	}
	else /* icols == 0 */
	{
		if (xN==1)
		{
			for (col = 0; col < nCols; col++)
			{
				xcol = x[col];
				if (isZero( xcol ) == BT_FALSE)
					for (i = jc[col]; i < jc[col+1]; i++)
						ytmp[ir[i]] += val[i] * xcol;
			}
		}
		else
		{
			// AW: I didn't test the case xN>1, but I hope it is working
			real_t* xcols = new real_t[xN];
			for (col = 0; col < nCols; col++)
			{
				real_t xmax = 0.0;
				for (k=0; k<xN; k++)
				{
					xcols[k] = x[k*xLD+col];
					xmax = getMax(xmax,getAbs(xcols[k]));
				}
				if (isZero( xmax ) == BT_FALSE)
					for (i = jc[col]; i < jc[col+1]; i++)
						for (k=0; k<xN; k++)
							ytmp[k*yfullLength+ir[i]] += val[i] * xcols[k];
				delete [] xcols;
			}
		}
	}

	if (yCompr == BT_TRUE)
	{
		if ( isZero(beta) == BT_TRUE )
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[j+k*yLD] = alpha*ytmp[irows->number[j]+k*yfullLength];
		else if (isEqual( beta, 1.0 ) == BT_TRUE)
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[j+k*yLD] += alpha*ytmp[irows->number[j]+k*yfullLength];
		else if (isEqual( beta, -1.0 ) == BT_TRUE)
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[j+k*yLD] = alpha*ytmp[irows->number[j]+k*yfullLength]-y[j+k*yLD];
		else if (isEqual( beta, 1.0 ) == BT_FALSE)
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[j+k*yLD] = alpha*ytmp[irows->number[j]+k*yfullLength]+beta*y[j+k*yLD];
	}
	else
	{
		if (isZero( beta ) == BT_TRUE)
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[irows->number[j]+k*yLD] = alpha*ytmp[irows->number[j]+k*yfullLength];
		else if (isEqual( beta, 1.0 ) == BT_TRUE)
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[irows->number[j]+k*yLD] = alpha*ytmp[irows->number[j]+k*yfullLength]+y[j+k*yLD];
		else if (isEqual( beta, -1.0 ) == BT_TRUE)
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[irows->number[j]+k*yLD] = alpha*ytmp[irows->number[j]+k*yfullLength]-y[j+k*yLD];
		else if (isEqual( beta, 1.0 ) == BT_FALSE)
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[irows->number[j]+k*yLD] = alpha*ytmp[irows->number[j]+k*yfullLength]+beta*y[j+k*yLD];
	}

	delete [] ytmp;
	return SUCCESSFUL_RETURN;
}


returnValue SparseMatrix::transTimes( const Indexlist* const irows, const Indexlist* const icols,
		int_t xN, real_t alpha, const real_t* x, int_t xLD, real_t beta, real_t* y, int_t yLD ) const
{
	long i, j, k, l, col;
	real_t yadd;

	if ( isZero(beta) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < icols->length; j++)
				y[j+k*yLD] = 0.0;
	else if ( isEqual(beta,-1.0) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < icols->length; j++)
				y[j+k*yLD] = -y[j+k*yLD];
	else if ( isEqual(beta,1.0) == BT_FALSE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < icols->length; j++)
				y[j+k*yLD] *= beta;
	if ( isEqual(alpha,0.0) == BT_TRUE )
		return SUCCESSFUL_RETURN;

	// work with full, unordered copy of x
	const int_t xfullLength = nRows;
	real_t* xtmp = new real_t[xfullLength];
	for (k = 0; k < xN; k++)
	{
		for (i = 0; i < xfullLength; i++)
			xtmp[i] = 0.0;
		for (i = 0; i < irows->length; i++)
			xtmp[irows->number[i]] = x[k*xLD+i];
		for (l = 0; l < icols->length; l++)
		{
			col = icols->iSort[l];
			yadd = 0.0;
			j = icols->number[col];
			for (i = jc[j]; i < jc[j+1]; i++)
				yadd += val[i] * xtmp[ir[i]];
			y[col] += alpha*yadd;
		}
		y += yLD; // move on to next RHS
	}

	delete [] xtmp;

	return SUCCESSFUL_RETURN;
}



returnValue SparseMatrix::addToDiag( real_t alpha )
{
	long i;

	if ( jd == 0 )
		return THROWERROR( RET_DIAGONAL_NOT_INITIALISED );

	if ( isZero( alpha ) == BT_FALSE )
	{
		for (i = 0; i < nRows && i < nCols; i++)
		{
			if (ir[jd[i]] == i)
				val[jd[i]] += alpha;
			else
				return RET_NO_DIAGONAL_AVAILABLE;
		}
	}

	return SUCCESSFUL_RETURN;
}


sparse_int_t* SparseMatrix::createDiagInfo( )
{
	sparse_int_t i, j;

	if (jd == 0) {
		jd = new sparse_int_t[nCols];

		for (j = 0; j < nCols; j++)
		{
			for (i = jc[j]; i < jc[j+1] && ir[i] < j; i++);
			jd[j] = i;
		}
	}

	return jd;
}



real_t* SparseMatrix::full( ) const
{
	sparse_int_t i, j;
	real_t* v = new real_t[nRows*nCols];

	for (i = 0; i < nCols*nRows; i++)
		v[i] = 0.0;

	for (j = 0; j < nCols; j++)
		for (i = jc[j]; i < jc[j+1]; i++)
			v[ir[i] * nCols + j] = val[i];

	return v;
}


returnValue SparseMatrix::print( const char* name ) const
{
	real_t* tmp = this->full();
	returnValue retVal = REFER_NAMESPACE_QPOASES print( tmp,nRows,nCols,name );
	delete[] tmp;

	return retVal;
}

returnValue SparseMatrix::writeToFile( FILE* output_file, const char* prefix ) const
{
	for (int_t i=0; i<=nCols; i++) {
		fprintf( output_file,"%sjc[%d] = %d\n",prefix,(int)i,(int)(jc[i]) );
	}
	for (int_t i=0; i<jc[nCols]; i++) {
		fprintf( output_file,"%sir[%d] = %d\n",prefix,(int)i,(int)(ir[i]) );
	}
	for (int_t i=0; i<jc[nCols]; i++) {
		fprintf( output_file,"%sval[%d] = %23.16e\n",prefix,(int)i,val[i] );
	}

	return SUCCESSFUL_RETURN;
}


SparseMatrixRow::SparseMatrixRow( ) : nRows(0), nCols(0), jr(0), ic(0), jd(0), val(0) {}

SparseMatrixRow::SparseMatrixRow( int_t nr, int_t nc, sparse_int_t* r, sparse_int_t* c, real_t* v )
	: nRows(nr), nCols(nc), jr(r), ic(c), jd(0), val(v) { doNotFreeMemory(); }

SparseMatrixRow::SparseMatrixRow( int_t nr, int_t nc, int_t ld, const real_t*  const v ) : nRows(nr), nCols(nc), jd(0)
{
	int_t i, j, nnz;

	jr = new sparse_int_t[nr+1];
	ic = new sparse_int_t[nr*nc];
	val = new real_t[nr*nc];

	nnz = 0;
	for (j = 0; j < nRows; j++)
	{
		jr[j] = nnz;
		for (i = 0; i < nCols; i++)
			if ( ( isZero( v[j*ld+i],0.0 ) == BT_FALSE ) || ( j == i ) )
			{
				ic[nnz] = i;
				val[nnz++] = v[j*ld+i];
			}
	}
	jr[nRows] = nnz;

	doFreeMemory( );
}


SparseMatrixRow::~SparseMatrixRow( )
{
	if (jd != 0)
	{
		delete[] jd;
		jd = 0;
	}

	if ( needToFreeMemory() == BT_TRUE )
		free( );
}


void SparseMatrixRow::free( )
{
	if (jr != 0) delete[] jr;
	jr = 0;
	if (ic != 0) delete[] ic;
	ic = 0;
	if (val != 0) delete[] val;
	val = 0;

	doNotFreeMemory( );
}


Matrix *SparseMatrixRow::duplicate( ) const
{
	long i, length = jr[nRows];
	SparseMatrixRow *dupl = new SparseMatrixRow;

	dupl->nRows = nRows;
	dupl->nCols = nCols;
	dupl->jr = new sparse_int_t[nRows+1];
	dupl->ic = new sparse_int_t[length];
	dupl->val = new real_t[length];

	for (i = 0; i < length; i++) dupl->jr[i] = jr[i];
	for (i = 0; i <= nCols; i++) dupl->ic[i] = ic[i];
	for (i = 0; i < length; i++) dupl->val[i] = val[i];

	if ( jd != 0 )
	{
		dupl->jd = new sparse_int_t[nRows];
		for (i = 0; i < nCols; i++) dupl->jd[i] = jd[i];
	}
	else
		dupl->jd = 0;

	dupl->doFreeMemory( );

	return dupl;
}



real_t SparseMatrixRow::diag( int_t i ) const
{
	if ( jd == 0 )
	{
		THROWERROR( RET_DIAGONAL_NOT_INITIALISED );
		return INFTY;
	}

	int_t entry = jd[i];
	return (entry < jr[i+1] && ic[entry] == i) ? val[entry] : 0.0;
}


BooleanType SparseMatrixRow::isDiag( ) const
{
	int_t i;

	if ( nCols != nRows )
		return BT_FALSE;

	for (i = 0; i < nRows; ++i)
	{
		if ( jr[i+1] > jr[i]+1 )
			return BT_FALSE;

		if ( ( jr[i+1] == jr[i]+1 ) && ( ic[jr[i]] != i ) )
			return BT_FALSE;
	}

	return BT_TRUE;
}



real_t SparseMatrixRow::getNorm(	int_t type
									) const
{
	int_t length = jr[nRows];
	return REFER_NAMESPACE_QPOASES getNorm( val,length,type );

}


real_t SparseMatrixRow::getRowNorm( int_t rNum, int_t type ) const
{
	int_t length = jr[rNum+1] - jr[rNum];
	return REFER_NAMESPACE_QPOASES getNorm( &(val[jr[rNum]]),length,type );
}


returnValue SparseMatrixRow::getRowNorm( real_t* norm, int_t type ) const
{
	int_t i;
	for (i = 0; i < nRows; ++i)
	{
	  int_t length = jr[i+1] - jr[i];
	  norm[i] = REFER_NAMESPACE_QPOASES getNorm( &(val[jr[i]]),length,type );
	}
	return SUCCESSFUL_RETURN;
}


returnValue SparseMatrixRow::getRow( int_t rNum, const Indexlist* const icols, real_t alpha, real_t* row ) const
{
	long i, j;

	if (icols != 0)
	{
		j = jr[rNum];
		i = 0;
		if ( isEqual(alpha,1.0) == BT_TRUE )
			while (j < jr[rNum+1] && i < icols->length)
				if (ic[j] == icols->number[icols->iSort[i]])
					row[icols->iSort[i++]] = val[j++];
				else if (ic[j] > icols->number[icols->iSort[i]])
					row[icols->iSort[i++]] = 0.0;
				else
					j++;
		else if ( isEqual(alpha,-1.0) == BT_TRUE )
			while (j < jr[rNum+1] && i < icols->length)
				if (ic[j] == icols->number[icols->iSort[i]])
					row[icols->iSort[i++]] = -val[j++];
				else if (ic[j] > icols->number[icols->iSort[i]])
					row[icols->iSort[i++]] = 0.0;
				else
					j++;
		else
			while (j < jr[rNum+1] && i < icols->length)
				if (ic[j] == icols->number[icols->iSort[i]])
					row[icols->iSort[i++]] = alpha * val[j++];
				else if (ic[j] > icols->number[icols->iSort[i]])
					row[icols->iSort[i++]] = 0.0;
				else
					j++;

		/* fill in remaining zeros */
		while (i < icols->length)
			row[icols->iSort[i++]] = 0.0;
	}
	else
	{
		for (i = 0; i < nCols; i++)
			row[i] = 0;

		if ( isEqual(alpha,1.0) == BT_TRUE )
			for (j = jr[rNum]; j < jr[rNum+1]; j++)
				row[ic[j]] = val[j];
		else if ( isEqual(alpha,-1.0) == BT_TRUE )
			for (j = jr[rNum]; j < jr[rNum+1]; j++)
				row[ic[j]] = -val[j];
		else
			for (j = jr[rNum]; j < jr[rNum+1]; j++)
				row[ic[j]] = alpha * val[j];
	}

	return SUCCESSFUL_RETURN;
}


returnValue SparseMatrixRow::getCol( int_t cNum, const Indexlist* const irows, real_t alpha, real_t* col ) const
{
	long i, j, k, srt;

	if (irows != 0)
	{
		if ( isEqual(alpha,1.0) == BT_TRUE )
			for (k = 0; k < irows->length; k++)
			{
				srt = irows->iSort[k];
				j = irows->number[srt];
				for (i = jr[j]; i < jr[j+1] && ic[i] < cNum; i++);
				col[srt] = (i < jr[j+1] && ic[i] == cNum) ? val[i] : 0.0;
			}
		else if ( isEqual(alpha,-1.0) == BT_TRUE )
			for (k = 0; k < irows->length; k++)
			{
				srt = irows->iSort[k];
				j = irows->number[srt];
				for (i = jr[j]; i < jr[j+1] && ic[i] < cNum; i++);
				col[srt] = (i < jr[j+1] && ic[i] == cNum) ? -val[i] : 0.0;
			}
		else
			for (k = 0; k < irows->length; k++)
			{
				srt = irows->iSort[k];
				j = irows->number[srt];
				for (i = jr[j]; i < jr[j+1] && ic[i] < cNum; i++);
				col[srt] = (i < jr[j+1] && ic[i] == cNum) ? alpha*val[i] : 0.0;
			}
	}
	else
	{
		if ( isEqual(alpha,1.0) == BT_TRUE )
			for (j = 0; j < nCols; j++)
			{
				for (i = jr[j]; i < jr[j+1] && ic[i] < cNum; i++);
				col[j] = (i < jr[j+1] && ic[i] == cNum) ? val[i] : 0.0;
			}
		else if ( isEqual(alpha,-1.0) == BT_TRUE )
			for (j = 0; j < irows->length; j++)
			{
				for (i = jr[j]; i < jr[j+1] && ic[i] < cNum; i++);
				col[j] = (i < jr[j+1] && ic[i] == cNum) ? -val[i] : 0.0;
			}
		else
			for (j = 0; j < irows->length; j++)
			{
				for (i = jr[j]; i < jr[j+1] && ic[i] < cNum; i++);
				col[j] = (i < jr[j+1] && ic[i] == cNum) ? alpha*val[i] : 0.0;
			}
	}
	return SUCCESSFUL_RETURN;
}

returnValue SparseMatrixRow::getSparseSubmatrix (
				int_t irowsLength, const int_t* const irowsNumber,
				int_t icolsLength, const int_t* const icolsNumber,
				int_t rowoffset, int_t coloffset, int_t& numNonzeros,	int_t* irn,
				int_t* jcn, real_t* avals, BooleanType only_lower_triangular /*= BT_FALSE */) const
{
	fprintf(stderr, "SparseMatrixRow::getSparseSubmatrix not implemented!\n");

	return THROWERROR(RET_NOT_YET_IMPLEMENTED);
}

returnValue SparseMatrixRow::times( int_t xN, real_t alpha, const real_t* x, int_t xLD,
		real_t beta, real_t* y, int_t yLD ) const
{
	long i, j, k;

	if ( isZero(beta) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < nRows; j++)
				y[j+k*yLD] = 0.0;
	else if ( isEqual(beta,-1.0) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < nRows; j++)
				y[j+k*yLD] = -y[j+k*yLD];
	else if ( isEqual(beta,1.0) == BT_FALSE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < nRows; j++)
				y[j+k*yLD] *= beta;

	if ( isEqual(alpha,1.0) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < nRows; j++)
				for (i = jr[j]; i < jr[j+1]; i++)
					y[j+k*yLD] += val[i] * x[ic[i]+k*xLD];
	else if ( isEqual(alpha,-1.0) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < nRows; j++)
				for (i = jr[j]; i < jr[j+1]; i++)
					y[j+k*yLD] -= val[i] * x[ic[i]+k*xLD];
	else
		for (k = 0; k < xN; k++)
			for (j = 0; j < nRows; j++)
				for (i = jr[j]; i < jr[j+1]; i++)
					y[j+k*yLD] += alpha * val[i] * x[ic[i]+k*xLD];

	return SUCCESSFUL_RETURN;
}


returnValue SparseMatrixRow::transTimes( int_t xN, real_t alpha, const real_t* x, int_t xLD,
		real_t beta, real_t* y, int_t yLD ) const
{
	long i, j, k;

	if ( isZero(beta) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < nCols; j++)
				y[j+k*yLD] = 0.0;
	else if ( isEqual(beta,-1.0) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < nCols; j++)
				y[j+k*yLD] = -y[j+k*yLD];
	else if ( isEqual(beta,1.0) == BT_FALSE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < nCols; j++)
				y[j+k*yLD] *= beta;

	if ( isEqual(alpha,1.0) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (i = 0; i < nRows; i++)
				for (j = jr[i]; j < jr[i+1]; j++)
					y[ic[j]+k*yLD] += val[j] * x[i+k*xLD];
	else if ( isEqual(alpha,-1.0) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (i = 0; i < nRows; i++)
				for (j = jr[i]; j < jr[i+1]; j++)
					y[ic[j]+k*yLD] -= val[j] * x[i+k*xLD];
	else
		for (k = 0; k < xN; k++)
			for (i = 0; i < nRows; i++)
				for (j = jr[i]; j < jr[i+1]; j++)
					y[ic[j]+k*yLD] += alpha * val[j] * x[i+k*xLD];

	return SUCCESSFUL_RETURN;
}


returnValue SparseMatrixRow::times( const Indexlist* const irows, const Indexlist* const icols,
		int_t xN, real_t alpha, const real_t* x, int_t xLD, real_t beta, real_t* y, int_t yLD,
		BooleanType yCompr ) const
{
	long i, j, k, l, srt, row;

	if (yCompr == BT_TRUE)
	{
		if ( isZero(beta) == BT_TRUE )
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[j+k*yLD] = 0.0;
		else if ( isEqual(beta,-1.0) == BT_TRUE )
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[j+k*yLD] = -y[j+k*yLD];
		else if ( isEqual(beta,1.0) == BT_FALSE )
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[j+k*yLD] *= beta;

		if (icols == 0)
			if ( isEqual(alpha,1.0) == BT_TRUE )
				for (l = 0; l < irows->length; l++)
				{
					srt = irows->iSort[l];
					row = irows->number[srt];
					for (j = jr[row]; j < jr[row+1]; j++)
						for (k = 0; k < xN; k++)
							y[k*yLD+srt] += val[j] * x[k*xLD+ic[j]];
				}
			else if ( isEqual(alpha,-1.0) == BT_TRUE )
				for (l = 0; l < irows->length; l++)
				{
					srt = irows->iSort[l];
					row = irows->number[srt];
					for (j = jr[row]; j < jr[row+1]; j++)
						for (k = 0; k < xN; k++)
							y[k*yLD+srt] -= val[j] * x[k*xLD+ic[j]];
				}
			else
				for (l = 0; l < irows->length; l++)
				{
					srt = irows->iSort[l];
					row = irows->number[srt];
					for (j = jr[row]; j < jr[row+1]; j++)
						for (k = 0; k < xN; k++)
							y[k*yLD+srt] += alpha * val[j] * x[k*xLD+ic[j]];
				}
		else /* icols != 0 */
			if ( isEqual(alpha,1.0) == BT_TRUE )
				for (l = 0; l < irows->length; l++)
				{
					srt = irows->iSort[l];
					row = irows->number[srt];
					j = jr[row];
					i = 0;
					while (j < jr[row+1] && i < icols->length)
						if (ic[j] == icols->number[icols->iSort[i]])
						{
							for (k = 0; k < xN; k++)
								y[k*yLD+srt] += val[j] * x[k*xLD+icols->iSort[i]];
							j++, i++;
						}
						else if (ic[j] > icols->number[icols->iSort[i]]) i++;
						else j++;
				}
			else if ( isEqual(alpha,-1.0) == BT_TRUE )
				for (l = 0; l < irows->length; l++)
				{
					srt = irows->iSort[l];
					row = irows->number[srt];
					j = jr[row];
					i = 0;
					while (j < jr[row+1] && i < icols->length)
						if (ic[j] == icols->number[icols->iSort[i]])
						{
							for (k = 0; k < xN; k++)
								y[k*yLD+srt] -= val[j] * x[k*xLD+icols->iSort[i]];
							j++, i++;
						}
						else if (ic[j] > icols->number[icols->iSort[i]]) i++;
						else j++;
				}
			else
				for (l = 0; l < irows->length; l++)
				{
					srt = irows->iSort[l];
					row = irows->number[srt];
					j = jr[row];
					i = 0;
					while (j < jr[row+1] && i < icols->length)
						if (ic[j] == icols->number[icols->iSort[i]])
						{
							for (k = 0; k < xN; k++)
								y[k*yLD+srt] += alpha * val[j] * x[k*xLD+icols->iSort[i]];
							j++, i++;
						}
						else if (ic[j] > icols->number[icols->iSort[i]]) i++;
						else j++;
				}
	}
	else /* y not compressed */
	{
		if ( isZero(beta) == BT_TRUE )
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[irows->number[j]+k*yLD] = 0.0;
		else if ( isEqual(beta,-1.0) == BT_TRUE )
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[irows->number[j]+k*yLD] = -y[j+k*yLD];
		else if ( isEqual(beta,1.0) == BT_FALSE )
			for (k = 0; k < xN; k++)
				for (j = 0; j < irows->length; j++)
					y[irows->number[j]+k*yLD] *= beta;

		if (icols == 0)
			if ( isEqual(alpha,1.0) == BT_TRUE )
				for (l = 0; l < irows->length; l++)
				{
					row = irows->number[irows->iSort[l]];
					for (j = jr[row]; j < jr[row+1]; j++)
						for (k = 0; k < xN; k++)
							y[k*yLD+row] += val[j] * x[k*xLD+ic[j]];
				}
			else if ( isEqual(alpha,-1.0) == BT_TRUE )
				for (l = 0; l < irows->length; l++)
				{
					row = irows->number[irows->iSort[l]];
					for (j = jr[row]; j < jr[row+1]; j++)
						for (k = 0; k < xN; k++)
							y[k*yLD+row] -= val[j] * x[k*xLD+ic[j]];
				}
			else
				for (l = 0; l < irows->length; l++)
				{
					row = irows->number[irows->iSort[l]];
					for (j = jr[row]; j < jr[row+1]; j++)
						for (k = 0; k < xN; k++)
							y[k*yLD+row] += alpha * val[j] * x[k*xLD+ic[j]];
				}
		else /* icols != 0 */
			if ( isEqual(alpha,1.0) == BT_TRUE )
				for (l = 0; l < irows->length; l++)
				{
					row = irows->iSort[l];
					j = jr[irows->number[row]];
					i = 0;
					while (j < jr[irows->number[row]+1] && i < icols->length)
						if (ic[j] == icols->number[icols->iSort[i]])
						{
							for (k = 0; k < xN; k++)
								y[k*yLD+row] += val[j] * x[k*xLD+icols->iSort[i]];
							j++, i++;
						}
						else if (ic[j] > icols->number[icols->iSort[i]]) i++;
						else j++;
				}
			else if ( isEqual(alpha,-1.0) == BT_TRUE )
				for (l = 0; l < irows->length; l++)
				{
					row = irows->iSort[l];
					j = jr[irows->number[row]];
					i = 0;
					while (j < jr[irows->number[row]+1] && i < icols->length)
						if (ic[j] == icols->number[icols->iSort[i]])
						{
							for (k = 0; k < xN; k++)
								y[k*yLD+row] -= val[j] * x[k*xLD+icols->iSort[i]];
							j++, i++;
						}
						else if (ic[j] > icols->number[icols->iSort[i]]) i++;
						else j++;
				}
			else
				for (l = 0; l < irows->length; l++)
				{
					row = irows->iSort[l];
					j = jr[irows->number[row]];
					i = 0;
					while (j < jr[irows->number[row]+1] && i < icols->length)
						if (ic[j] == icols->number[icols->iSort[i]])
						{
							for (k = 0; k < xN; k++)
								y[k*yLD+row] += alpha * val[j] * x[k*xLD+icols->iSort[i]];
							j++, i++;
						}
						else if (ic[j] > icols->number[icols->iSort[i]]) i++;
						else j++;
				}
	}
	return SUCCESSFUL_RETURN;
}


returnValue SparseMatrixRow::transTimes( const Indexlist* const irows, const Indexlist* const icols,
		int_t xN, real_t alpha, const real_t* x, int_t xLD, real_t beta, real_t* y, int_t yLD ) const
{
	long i, j, k, l, row, srt;

	if ( isZero(beta) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < icols->length; j++)
				y[j+k*yLD] = 0.0;
	else if ( isEqual(beta,-1.0) == BT_TRUE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < icols->length; j++)
				y[j+k*yLD] = -y[j+k*yLD];
	else if ( isEqual(beta,1.0) == BT_FALSE )
		for (k = 0; k < xN; k++)
			for (j = 0; j < icols->length; j++)
				y[j+k*yLD] *= beta;

	if ( isEqual(alpha,1.0) == BT_TRUE )
		for (l = 0; l < irows->length; l++)
		{
			srt = irows->iSort[l];
			row = irows->number[srt];
			j = jr[row];
			i = 0;
			while (j < jr[row+1] && i < icols->length)
				if (ic[j] == icols->number[icols->iSort[i]])
				{
					for (k = 0; k < xN; k++)
						y[k*yLD+icols->iSort[i]] += val[j] * x[k*xLD+srt];
					j++, i++;
				}
				else if (ic[j] > icols->number[icols->iSort[i]]) i++;
				else j++;
		}
	else if ( isEqual(alpha,-1.0) == BT_TRUE )
		for (l = 0; l < irows->length; l++)
		{
			srt = irows->iSort[l];
			row = irows->number[srt];
			j = jr[row];
			i = 0;
			while (j < jr[row+1] && i < icols->length)
				if (ic[j] == icols->number[icols->iSort[i]])
				{
					for (k = 0; k < xN; k++)
						y[k*yLD+icols->iSort[i]] -= val[j] * x[k*xLD+srt];
					j++, i++;
				}
				else if (ic[j] > icols->number[icols->iSort[i]]) i++;
				else j++;
		}
	else
		for (l = 0; l < irows->length; l++)
		{
			srt = irows->iSort[l];
			row = irows->number[srt];
			j = jr[row];
			i = 0;
			while (j < jr[row+1] && i < icols->length)
				if (ic[j] == icols->number[icols->iSort[i]])
				{
					for (k = 0; k < xN; k++)
						y[k*yLD+icols->iSort[i]] += alpha * val[j] * x[k*xLD+srt];
					j++, i++;
				}
				else if (ic[j] > icols->number[icols->iSort[i]]) i++;
				else j++;
		}

	return SUCCESSFUL_RETURN;
}



returnValue SparseMatrixRow::addToDiag( real_t alpha )
{
	long i;

	if ( jd == 0 )
		return THROWERROR( RET_DIAGONAL_NOT_INITIALISED );

	if ( isZero(alpha) == BT_FALSE )
	{
		for (i = 0; i < nRows && i < nCols; i++)
		{
			if (ic[jd[i]] == i)
				val[jd[i]] += alpha;
			else
				return RET_NO_DIAGONAL_AVAILABLE;
		}
	}

	return SUCCESSFUL_RETURN;
}


sparse_int_t* SparseMatrixRow::createDiagInfo( )
{
	sparse_int_t i, j;

	if (jd == 0) {
		jd = new sparse_int_t[nRows];

		for (i = 0; i < nRows; i++)
		{
			for (j = jr[i]; j < jr[i+1] && ic[j] < i; j++);
			jd[i] = j;
		}
	}

	return jd;
}


real_t* SparseMatrixRow::full( ) const
{
	sparse_int_t i, j;
	real_t* v = new real_t[nRows*nCols];

	for (i = 0; i < nCols*nRows; i++)
		v[i] = 0.0;

	for (i = 0; i < nRows; i++)
		for (j = jr[i]; j < jr[i+1]; j++)
			v[ic[j] + i * nCols] = val[j];

	return v;
}


returnValue SparseMatrixRow::print( const char* name ) const
{
	real_t* tmp = this->full();
	returnValue retVal = REFER_NAMESPACE_QPOASES print( tmp,nRows,nCols,name );
	delete[] tmp;

	return retVal;
}

returnValue SparseMatrixRow::writeToFile( FILE* output_file, const char* prefix ) const
{
	return THROWERROR( RET_NOT_YET_IMPLEMENTED );
}

Matrix *SymSparseMat::duplicate() const
{
	return duplicateSym();
}


SymmetricMatrix *SymSparseMat::duplicateSym( ) const
{
	/* "same" as duplicate() in SparseMatrix */
	long i, length = jc[nCols];
	SymSparseMat *dupl = new SymSparseMat;

	dupl->nRows = nRows;
	dupl->nCols = nCols;
	dupl->ir = new sparse_int_t[length];
	dupl->jc = new sparse_int_t[nCols+1];
	dupl->val = new real_t[length];

	for (i = 0; i < length; i++) dupl->ir[i] = ir[i];
	for (i = 0; i <= nCols; i++) dupl->jc[i] = jc[i];
	for (i = 0; i < length; i++) dupl->val[i] = val[i];

	if ( jd != 0 )
	{
		dupl->jd = new sparse_int_t[nCols];
		for (i = 0; i < nCols; i++) dupl->jd[i] = jd[i];
	}
	else
		dupl->jd = 0;

	dupl->doFreeMemory( );

	return dupl;
}


returnValue SymSparseMat::bilinear( const Indexlist* const icols,
		int_t xN, const real_t* x, int_t xLD, real_t* y, int_t yLD ) const
{
	int_t i, j, k, l, idx, row, col;

	if ( jd == 0 )
		return THROWERROR( RET_DIAGONAL_NOT_INITIALISED );

	/* clear output */
	for (i = 0; i < xN*xN; i++)
		y[i] = 0.0;

	/* compute lower triangle */
	for (l = 0; l < icols->length; l++)
	{
		col = icols->number[icols->iSort[l]];
		idx = jd[col];
		k = 0;
		while (idx < jc[col+1] && k < icols->length)
		{
			row = icols->number[icols->iSort[k]];
			if (ir[idx] == row)
			{
				/* TODO: It is possible to formulate this as DSYR and DSYR2
				 * operations. */
				if (row == col) /* diagonal element */
					for (i = 0; i < xN; i++)
						for (j = i; j < xN; j++)
							y[i*yLD+j] += val[idx] * x[i*xLD+col] * x[j*xLD+col];
				else /* subdiagonal elements */
					for (i = 0; i < xN; i++)
						for (j = i; j < xN; j++)
							y[i*yLD+j] += val[idx] * (x[i*xLD+col] * x[j*xLD+row] + x[i*xLD+row] * x[j*xLD+col]);
				idx++, k++;
			}
			else if (ir[idx] > row) k++;
			else idx++;
		}
	}

	/* fill upper triangle */
	for (i = 0; i < xN; i++)
		for (j = i; j < xN; j++)
			y[j*yLD+i] = y[i*yLD+j];

	return SUCCESSFUL_RETURN;
}


END_NAMESPACE_QPOASES


/*
 *	end of file
 */
