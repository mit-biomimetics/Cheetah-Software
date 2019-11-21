#ifndef QDLDL_H
#define QDLDL_H

// Include qdldl type options
#include "qdldl_types.h"

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

/**
  * Compute the elimination tree for a quasidefinite matrix
  * in compressed sparse column form, where the input matrix is
  * assumed to contain data for the upper triangular part of A only,
  * and there are no duplicate indices.
  *
  * Returns an elimination tree for the factorization A = LDL^T and a
  * count of the nonzeros in each column of L that are strictly below the
  * diagonal.
  *
  * Does not use MALLOC.  It is assumed that the arrays work, Lnz, and
  * etree will be allocated with a number of elements equal to n.
  *
  * The data in (n,Ap,Ai) are from a square matrix A in CSC format, and
  * should include the upper triangular part of A only.
  *
  * This function is only intended for factorisation of QD matrices specified
  * by their upper triangular part.  An error is returned if any column has
  * data below the diagonal or s completely empty.
  *
  * For matrices with a non-empty column but a zero on the corresponding diagonal,
  * this function will *not* return an error, as it may still be possible to factor
  * such a matrix in LDL form.   No promises are made in this case though...
  *
  * @param   n     number of columns in CSC matrix A (assumed square)
  * @param  Ap     column pointers (size n+1) for columns of A
  * @param  Ai     row indices of A.  Has Ap[n] elements
  * @param  work   work vector (size n) (no meaning on return)
  * @param  Lnz    count of nonzeros in each column of L (size n) below diagonal
  * @param  etree  elimination tree (size n)
  * @return total  sum of Lnz (i.e. total nonzeros in L below diagonal). Returns
  *                -1 if the input does not have triu structure or has an empty
  *                column.
  *
  *
*/

 QDLDL_int QDLDL_etree(const QDLDL_int   n,
                       const QDLDL_int* Ap,
                       const QDLDL_int* Ai,
                       QDLDL_int* work,
                       QDLDL_int* Lnz,
                       QDLDL_int* etree);

/**
  * Compute an LDL decomposition for a quasidefinite matrix
  * in compressed sparse column form, where the input matrix is
  * assumed to contain data for the upper triangular part of A only,
  * and there are no duplicate indices.
  *
  * Returns factors L, D and Dinv = 1./D.
  *
  * Does not use MALLOC.  It is assumed that L will be a compressed
  * sparse column matrix with data (Ln,Lp,Li)  with sufficient space
  * allocated, with a number of nonzeros equal to the count given
  * as a return value by osqp_ldl_etree
  *
  * @param   n     number of columns in L and A (both square)
  * @param  Ap     column pointers (size n+1) for columns of A
  * @param  Ai     row indices of A.  Has Ap[n] elements
  * @param  Ln     number of columns in CSC matrix L
  * @param  Lp     column pointers (size Ln+1) for columns of L
  * @param  Li     row indices of L.  Has Lp[Ln] elements
  * @param  D      vectorized factor D.  Length is n
  * @param  Dinv   reciprocal of D.  Length is n
  * @param  Lnz    count of nonzeros in each column of L below diagonal,
  *                as given by osqp_ldl_etree (not modified)
  * @param  etree  elimination tree as as given by osqp_ldl_etree (not modified)
  * @param  bwork  working array of bools. Length is n
  * @param  iwork  working array of integers. Length is 3*n
  * @param  fwork  working array of floats. Length is n
  * @return        Returns a count of the number of positive elements
  *                in D.  Returns -1 and exits immediately if any element
  *                of D evaluates exactly to zero (matrix is not quasidefinite
  *                or otherwise LDL factorisable)
  *
*/


QDLDL_int QDLDL_factor(const QDLDL_int    n,
                  const QDLDL_int*   Ap,
                  const QDLDL_int*   Ai,
                  const QDLDL_float* Ax,
                  QDLDL_int*   Lp,
                  QDLDL_int*   Li,
                  QDLDL_float* Lx,
                  QDLDL_float* D,
                  QDLDL_float* Dinv,
                  const QDLDL_int* Lnz,
                  const QDLDL_int* etree,
                  QDLDL_bool* bwork,
                  QDLDL_int* iwork,
                  QDLDL_float* fwork);


/**
  * Solves LDL'x = b
  *
  * It is assumed that L will be a compressed
  * sparse column matrix with data (Ln,Lp,Li).
  *
  * @param   n     number of columns in L (both square)
  * @param  Ln     number of columns in CSC matrix L
  * @param  Lp     column pointers (size Ln+1) for columns of L
  * @param  Li     row indices of L.  Has Lp[Ln] elements
  * @param  Dinv   reciprocal of D.  Length is n
  * @param  x      initialized to b.  Equal to x on return
  *
  *
*/
void QDLDL_solve(const QDLDL_int    n,
                 const QDLDL_int*   Lp,
                 const QDLDL_int*   Li,
                 const QDLDL_float* Lx,
                 const QDLDL_float* Dinv,
                 QDLDL_float* x);


/**
 * Solves (L+I)x = b
 *
 * It is assumed that L will be a compressed
 * sparse column matrix with data (Ln,Lp,Li).
 *
 * @param   n     number of columns in L (both square)
 * @param  Ln     number of columns in CSC matrix L
 * @param  Lp     column pointers (size Ln+1) for columns of L
 * @param  Li     row indices of L.  Has Lp[Ln] elements
 * @param  Dinv   reciprocal of D.  Length is n
 * @param  x      initialized to b.  Equal to x on return
 *
 *
*/

void QDLDL_Lsolve(const QDLDL_int    n,
                  const QDLDL_int*   Lp,
                  const QDLDL_int*   Li,
                  const QDLDL_float* Lx,
                  QDLDL_float* x);

/**
 * Solves (L+I)'x = b
 *
 * It is assumed that L will be a compressed
 * sparse column matrix with data (Ln,Lp,Li).
 *
 * @param   n     number of columns in L (both square)
 * @param  Ln     number of columns in CSC matrix L
 * @param  Lp     column pointers (size Ln+1) for columns of L
 * @param  Li     row indices of L.  Has Lp[Ln] elements
 * @param  Dinv   reciprocal of D.  Length is n
 * @param  x      initialized to b.  Equal to x on return
 *
 *
*/

void QDLDL_Ltsolve(const QDLDL_int    n,
                   const QDLDL_int*   Lp,
                   const QDLDL_int*   Li,
                   const QDLDL_float* Lx,
                   QDLDL_float* x);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef QDLDL_H
