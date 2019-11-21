#ifndef LIN_ALG_H
# define LIN_ALG_H


# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

# include "types.h"


/* VECTOR FUNCTIONS ----------------------------------------------------------*/

# ifndef EMBEDDED

/* copy vector a into output (Uses MALLOC)*/
c_float* vec_copy(c_float *a,
                  c_int    n);
# endif // ifndef EMBEDDED

/* copy vector a into preallocated vector b */
void prea_vec_copy(const c_float *a,
                   c_float       *b,
                   c_int          n);

/* copy integer vector a into preallocated vector b */
void prea_int_vec_copy(const c_int *a,
                       c_int       *b,
                       c_int        n);

/* set float vector to scalar */
void vec_set_scalar(c_float *a,
                    c_float  sc,
                    c_int    n);

/* set integer vector to scalar */
void int_vec_set_scalar(c_int *a,
                        c_int  sc,
                        c_int  n);

/* add scalar to vector*/
void vec_add_scalar(c_float *a,
                    c_float  sc,
                    c_int    n);

/* multiply scalar to vector */
void vec_mult_scalar(c_float *a,
                     c_float  sc,
                     c_int    n);

/* c = a + sc*b */
void vec_add_scaled(c_float       *c,
                    const c_float *a,
                    const c_float *b,
                    c_int          n,
                    c_float        sc);

/* ||v||_inf */
c_float vec_norm_inf(const c_float *v,
                     c_int          l);

/* ||Sv||_inf */
c_float vec_scaled_norm_inf(const c_float *S,
                            const c_float *v,
                            c_int          l);

/* ||a - b||_inf */
c_float vec_norm_inf_diff(const c_float *a,
                          const c_float *b,
                          c_int          l);

/* mean of vector elements */
c_float vec_mean(const c_float *a,
                 c_int          n);

# if EMBEDDED != 1

/* Vector elementwise reciprocal b = 1./a (needed for scaling)*/
void vec_ew_recipr(const c_float *a,
                   c_float       *b,
                   c_int          n);
# endif // if EMBEDDED != 1

/* Inner product a'b */
c_float vec_prod(const c_float *a,
                 const c_float *b,
                 c_int          n);

/* elementwse product a.*b stored in c*/
void vec_ew_prod(const c_float *a,
                 const c_float *b,
                 c_float       *c,
                 c_int          n);

# if EMBEDDED != 1

/* elementwise sqrt of the vector elements */
void vec_ew_sqrt(c_float *a,
                 c_int    n);

/* elementwise max between each vector component and max_val */
void vec_ew_max(c_float *a,
                c_int    n,
                c_float  max_val);

/* elementwise min between each vector component and max_val */
void vec_ew_min(c_float *a,
                c_int    n,
                c_float  min_val);

/* Elementwise maximum between vectors c = max(a, b) */
void vec_ew_max_vec(const c_float *a,
                    const c_float *b,
                    c_float       *c,
                    c_int          n);

/* Elementwise minimum between vectors c = min(a, b) */
void vec_ew_min_vec(const c_float *a,
                    const c_float *b,
                    c_float       *c,
                    c_int          n);

# endif // if EMBEDDED != 1


/* MATRIX FUNCTIONS ----------------------------------------------------------*/

/* multiply scalar to matrix */
void mat_mult_scalar(csc    *A,
                     c_float sc);

/* Premultiply matrix A by diagonal matrix with diagonal d,
   i.e. scale the rows of A by d
 */
void mat_premult_diag(csc           *A,
                      const c_float *d);

/* Premultiply matrix A by diagonal matrix with diagonal d,
   i.e. scale the columns of A by d
 */
void mat_postmult_diag(csc           *A,
                       const c_float *d);


/* Matrix-vector multiplication
 *    y  =  A*x  (if plus_eq == 0)
 *    y +=  A*x  (if plus_eq == 1)
 *    y -=  A*x  (if plus_eq == -1)
 */
void mat_vec(const csc     *A,
             const c_float *x,
             c_float       *y,
             c_int          plus_eq);


/* Matrix-transpose-vector multiplication
 *    y  =  A'*x  (if plus_eq == 0)
 *    y +=  A'*x  (if plus_eq == 1)
 *    y -=  A'*x  (if plus_eq == -1)
 * If skip_diag == 1, then diagonal elements of A are assumed to be zero.
 */
void mat_tpose_vec(const csc     *A,
                   const c_float *x,
                   c_float       *y,
                   c_int          plus_eq,
                   c_int          skip_diag);


# if EMBEDDED != 1

/**
 * Infinity norm of each matrix column
 * @param M	Input matrix
 * @param E     Vector of infinity norms
 *
 */
void mat_inf_norm_cols(const csc *M,
                       c_float   *E);

/**
 * Infinity norm of each matrix row
 * @param M	Input matrix
 * @param E     Vector of infinity norms
 *
 */
void mat_inf_norm_rows(const csc *M,
                       c_float   *E);

/**
 * Infinity norm of each matrix column
 * Matrix M is symmetric upper-triangular
 *
 * @param M	Input matrix (symmetric, upper-triangular)
 * @param E     Vector of infinity norms
 *
 */
void mat_inf_norm_cols_sym_triu(const csc *M,
                                c_float   *E);

# endif // EMBEDDED != 1

/**
 * Compute quadratic form f(x) = 1/2 x' P x
 * @param  P quadratic matrix in CSC form (only upper triangular)
 * @param  x argument float vector
 * @return   quadratic form value
 */
c_float quad_form(const csc     *P,
                  const c_float *x);


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef LIN_ALG_H
