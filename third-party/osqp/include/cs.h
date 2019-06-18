#ifndef CS_H
# define CS_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

# include "types.h"   // CSC matrix type
# include "lin_alg.h" // Vector copy operations

/*****************************************************************************
* Create and free CSC Matrices                                              *
*****************************************************************************/

/**
 * Create Compressed-Column-Sparse matrix from existing arrays
    (no MALLOC to create inner arrays x, i, p)
 * @param  m     First dimension
 * @param  n     Second dimension
 * @param  nzmax Maximum number of nonzero elements
 * @param  x     Vector of data
 * @param  i     Vector of row indices
 * @param  p     Vector of column pointers
 * @return       New matrix pointer
 */
csc* csc_matrix(c_int    m,
                c_int    n,
                c_int    nzmax,
                c_float *x,
                c_int   *i,
                c_int   *p);


/**
 * Create uninitialized CSC matrix atricture
    (uses MALLOC to create inner arrays x, i, p)
 * @param  m       First dimension
 * @param  n       Second dimension
 * @param  nzmax   Maximum number of nonzero elements
 * @param  values  Allocate values (0/1)
 * @param  triplet Allocate CSC or triplet format matrix (1/0)
 * @return         Matrix pointer
 */
csc* csc_spalloc(c_int m,
                 c_int n,
                 c_int nzmax,
                 c_int values,
                 c_int triplet);


/**
 * Free sparse matrix
    (uses FREE to free inner arrays x, i, p)
 * @param  A Matrix in CSC format
 */
void csc_spfree(csc *A);


/**
 * free workspace and return a sparse matrix result
 * @param  C  CSC matrix
 * @param  w  Workspace vector
 * @param  x  Workspace vector
 * @param  ok flag
 * @return    Return result C if OK, otherwise free it
 */
csc* csc_done(csc  *C,
              void *w,
              void *x,
              c_int ok);

/*****************************************************************************
* Copy Matrices                                                             *
*****************************************************************************/

/**
 *  Copy sparse CSC matrix A to output.
 *  output is allocated by this function (uses MALLOC)
 */
csc* copy_csc_mat(const csc *A);


/**
 *  Copy sparse CSC matrix A to B (B is preallocated, NO MALOC)
 */
void prea_copy_csc_mat(const csc *A,
                       csc       *B);


/*****************************************************************************
* Matrices Conversion                                                       *
*****************************************************************************/


/**
 * C = compressed-column CSC from matrix T in triplet form
 *
 * TtoC stores the vector of indices from T to C
 *  -> C[TtoC[i]] = T[i]
 *
 * @param  T    matrix in triplet format
 * @param  TtoC vector of indices from triplet to CSC format
 * @return      matrix in CSC format
 */
csc* triplet_to_csc(const csc *T,
                    c_int     *TtoC);


/**
 * C = compressed-row CSR from matrix T in triplet form
 *
 * TtoC stores the vector of indices from T to C
 *  -> C[TtoC[i]] = T[i]
 *
 * @param  T    matrix in triplet format
 * @param  TtoC vector of indices from triplet to CSR format
 * @return      matrix in CSR format
 */
csc* triplet_to_csr(const csc *T,
                    c_int     *TtoC);


/**
 * Convert sparse to dense
 */
//c_float* csc_to_dns(csc *M);


/**
 * Convert square CSC matrix into upper triangular one
 *
 * @param  M         Matrix to be converted
 * @return           Upper triangular matrix in CSC format
 */
csc* csc_to_triu(csc *M);


/*****************************************************************************
* Extra operations                                                          *
*****************************************************************************/

/**
 * p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c
 *
 * @param  p Create cumulative sum into p
 * @param  c Vector of which we compute cumulative sum
 * @param  n Number of elements
 * @return   Exitflag
 */
c_int csc_cumsum(c_int *p,
                 c_int *c,
                 c_int  n);

/**
 * Compute inverse of permuation matrix stored in the vector p.
 * The computed inverse is also stored in a vector.
 */
c_int* csc_pinv(c_int const *p,
                c_int        n);

/**
 * C = A(p,p)= PAP' where A and C are symmetric the upper part stored;
 *  NB: pinv not p!
 * @param  A      Original matrix (upper-triangular)
 * @param  pinv   Inverse of permutation vector
 * @param  AtoC   Mapping from indices of A-x to C->x
 * @param  values Are values of A allocated?
 * @return        New matrix (allocated)
 */
csc* csc_symperm(const csc   *A,
                 const c_int *pinv,
                 c_int       *AtoC,
                 c_int        values);


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef CS_H
