#ifndef QDLDL_INTERFACE_H
#define QDLDL_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "types.h"
#include "qdldl_types.h"

/**
 * QDLDL solver structure
 */
typedef struct qdldl qdldl_solver;

struct qdldl {
    enum linsys_solver_type type;

    /**
     * @name Functions
     * @{
     */
    c_int (*solve)(struct qdldl * self, c_float * b, const OSQPSettings * settings);

    #ifndef EMBEDDED
    void (*free)(struct qdldl * self); ///< Free workspace (only if desktop)
    #endif

    // This used only in non embedded or embedded 2 version
    #if EMBEDDED != 1
    c_int (*update_matrices)(struct qdldl * self, const csc *P, const csc *A, const OSQPSettings *settings); ///< Update solver matrices
    c_int (*update_rho_vec)(struct qdldl * self, const c_float * rho_vec, const c_int m); ///< Update solver matrices
    #endif

    #ifndef EMBEDDED
    c_int nthreads;
    #endif
    /** @} */

    /**
     * @name Attributes
     * @{
     */
    csc *L;         ///< lower triangular matrix in LDL factorization
    c_float *Dinv;  ///< inverse of diag matrix in LDL (as a vector)
    c_int   *P;     ///< permutation of KKT matrix for factorization
    c_float *bp;    ///< workspace memory for solves


    #if EMBEDDED != 1
    // These are required for matrix updates
    c_int * Pdiag_idx, Pdiag_n;  ///< index and number of diagonal elements in P
    csc   * KKT;                 ///< Permuted KKT matrix in sparse form (used to update P and A matrices)
    c_int * PtoKKT, * AtoKKT;    ///< Index of elements from P and A to KKT matrix
    c_int * rhotoKKT;            ///< Index of rho places in KKT matrix
    // QDLDL Numeric workspace
    QDLDL_float *D;
    QDLDL_int   *etree;
    QDLDL_int   *Lnz;
    QDLDL_int   *iwork;
    QDLDL_bool  *bwork;
    QDLDL_float *fwork;
    #endif

    /** @} */
};



/**
 * Initialize QDLDL Solver
 *
 * @param  P      Cost function matrix (upper triangular form)
 * @param  A      Constraints matrix
 * @param	sigma   Algorithm parameter. If polish, then sigma = delta.
 * @param	rho_vec Algorithm parameter. If polish, then rho_vec = OSQP_NULL.
 * @param  polish Flag whether we are initializing for polish or not
 * @return        Initialized private structure
 */
qdldl_solver *init_linsys_solver_qdldl(const csc * P, const csc * A, c_float sigma, c_float * rho_vec, c_int polish);

/**
 * Solve linear system and store result in b
 * @param  s        Linear system solver structure
 * @param  b        Right-hand side
 * @param  settings OSQP solver settings
 * @return          Exitflag
 */
c_int solve_linsys_qdldl(qdldl_solver * s, c_float * b, const OSQPSettings * settings);


#if EMBEDDED != 1
/**
 * Update linear system solver matrices
 * @param  s        Linear system solver structure
 * @param  P        Matrix P
 * @param  A        Matrix A
 * @param  settings Solver settings
 * @return          Exitflag
 */
c_int update_linsys_solver_matrices_qdldl(qdldl_solver * s,
		const csc *P, const csc *A, const OSQPSettings *settings);




/**
 * Update rho parameter in linear system solver structure
 * @param  s   Linear system solver structure
 * @param  rho new rho value
 * @param  m   number of constraints
 * @return     exitflag
 */
c_int update_linsys_solver_rho_vec_qdldl(qdldl_solver * s, const c_float * rho_vec, const c_int m);

#endif

#ifndef EMBEDDED
/**
 * Free linear system solver
 * @param s linear system solver object
 */
void free_linsys_solver_qdldl(qdldl_solver * s);
#endif

#ifdef __cplusplus
}
#endif

#endif
