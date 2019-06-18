/*
 * -----------------------------------------------------------------
 * $Revision: 4378 $
 * $Date: 2015-02-19 10:55:14 -0800 (Thu, 19 Feb 2015) $
 * -----------------------------------------------------------------
 * Programmer: Carol Woodward @ LLNL,
 *             Daniel R. Reynolds @ SMU.
 * -----------------------------------------------------------------
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This header file contains definitions and declarations for use by
 * sparse linear solvers for Ax = b. 
 * -----------------------------------------------------------------
 */

#ifndef _SUNDIALS_SPARSE_H
#define _SUNDIALS_SPARSE_H

#include <sundials/sundials_types.h>
#include <sundials/sundials_direct.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * ==================================================================
 * Type definitions
 * ==================================================================
 */

/*
 * -----------------------------------------------------------------
 * Type : SlsMat
 * -----------------------------------------------------------------
 * The type SlsMat is defined to be a pointer to a structure
 * with various sizes, a data field, and arrays for the row and 
 * column information for the sparse matrix entries.
 * The M and N fields indicates the number 
 * of rows and columns, respectively. The data field is a one 
 * dimensional array used for component storage. The NNZ field indicates
 * the number of nonzero entries in the matrix. The integer array, asub, 
 * holds the row index for each of the matrix entries.  The integer
 * array, xa, holds the index entry for the starting value of each column.
 * -----------------------------------------------------------------
 * The relevant fields in DlsMat are:
 *    M     - number of rows
 *    N     - number of columns
 *    NNZ   - the number of nonzero entries in the matrix
 *    data  - pointer to a contiguous block of realtype variables
 *    rowvals - row indices of each nonzero entry
 *    colptrs - starting index of the first entry in data in each column
 *
 * The nonzero entries of the matrix are stored in
 * compressed column format.  Row indices of entries in 
 * column j are stored in rowvals[colptrs[j]] through rowvals[colptrs[j+i]-1]
 * and corresponding numerical values of the matrix are stored 
 * in the same entries of data.
 * -----------------------------------------------------------------
 */

typedef struct _SlsMat {
  int M;
  int N;
  int NNZ;
  realtype *data;
  int *rowvals;
  int *colptrs;
} *SlsMat;

/*
 * ==================================================================
 * Exported function prototypes (functions working on SlsMat)
 * ==================================================================
 */

/*
 * -----------------------------------------------------------------
 * Function: NewSparseMat
 * -----------------------------------------------------------------
 * NewSparseMat allocates memory for a compressed column sparse
 * matrix with M rows, N columns, and NNZ nonzeros. NewSparseMat
 * returns NULL if the request for matrix storage cannot be
 * satisfied. See the above documentation for the type SlsMat
 * for matrix storage details.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SlsMat NewSparseMat(int M, int N, int NNZ);

/*
 * -----------------------------------------------------------------
 * Function: SlsConvertDls
 * -----------------------------------------------------------------
 * SlsConvertDense creates a new sparse matrix from an existing
 * dense/band matrix by copying all nonzero values into the sparse 
 * matrix structure.  SlsConvertDense returns NULL if the request 
 * for matrix storage cannot be satisfied. 
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SlsMat SlsConvertDls(DlsMat A);

/*
 * -----------------------------------------------------------------
 * Functions: DestroySparseMat
 * -----------------------------------------------------------------
 * DestroySparseMat frees the memory allocated by NewSparseMat
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void DestroySparseMat(SlsMat A);

/*
 * -----------------------------------------------------------------
 * Function : SlsSetToZero
 * -----------------------------------------------------------------
 * SetToZero sets all the elements of the sparse matrix A to 0.0.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void SlsSetToZero(SlsMat A);

/*
 * -----------------------------------------------------------------
 * Functions: CopySparseMat
 * -----------------------------------------------------------------
 * This function copies sparse matrix A into sparse matrix B.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void CopySparseMat(SlsMat A, SlsMat B);

/*
 * -----------------------------------------------------------------
 * Functions: ScaleSparseMat
 * -----------------------------------------------------------------
 * This function scales all data entries of a sparse matrix A 
 * by the realtype number in b.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void ScaleSparseMat(realtype b, SlsMat A);

/*
 * -----------------------------------------------------------------
 * Functions: AddIdentitySparseMat
 * -----------------------------------------------------------------
 * This function adds 1 to every diagonal entry of A.
 * Note that the resulting matrix may have more nonzero entries than 
 * the original.  This is accounted for, so that the return matrix 
 * may be larger than the one sent in.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void AddIdentitySparseMat(SlsMat A);

/*
 * -----------------------------------------------------------------
 * Functions: SlsAddMat
 * -----------------------------------------------------------------
 * This function adds two sparse matrices: A = A+B.
 * Note that the resulting matrix may have more nonzero entries than
 * either of the original matrices.  This is accounted for, so that 
 * the return matrix may be larger than the ones sent in.  Upon 
 * successful completion, the return value is zero; otherwise 1 is 
 * returned.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int SlsAddMat(SlsMat A, SlsMat B);

/*
 * -----------------------------------------------------------------
 * Functions: ReallocSparseMat
 * -----------------------------------------------------------------
 * This function reallocs internal arrays so that the resulting matrix 
 * holds colptrs[N] nonzeros.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void ReallocSparseMat(SlsMat A);

/*
 * -----------------------------------------------------------------
 * Functions: SlsMatvec
 * -----------------------------------------------------------------
 * This function computes the matrix-vector product, y=A*x, where A
 * is a sparse matrix of dimension MxN, x is a realtype array of 
 * length N, and y is a realtype array of length M. Upon successful
 * completion, the return value is zero; otherwise 1 is returned.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int SlsMatvec(SlsMat A, realtype *x, realtype *y);

/*
 * -----------------------------------------------------------------
 * Functions: PrintSparseMat
 * -----------------------------------------------------------------
 * This function prints the compressed column matrix information for 
 * matrix A to standard output.
 * It is intended as a debugging tool with small values of NNZ.
 * The elements are printed using the %g/%lg/%Lg option. 
 * A blank line is printed before and after the matrix.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void PrintSparseMat(SlsMat A);


#ifdef __cplusplus
}
#endif

#endif
