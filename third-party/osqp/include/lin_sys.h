#ifndef LIN_SYS_H
# define LIN_SYS_H

/* KKT linear system definition and solution */

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

# include "types.h"

# ifndef EMBEDDED


/**
 * Load linear system solver shared library
 * @param	linsys_solver  Linear system solver
 * @return Zero on success, nonzero on failure.
 */
c_int load_linsys_solver(enum linsys_solver_type linsys_solver);


/**
 * Unload linear system solver shared library
 * @param	linsys_solver  Linear system solver
 * @return Zero on success, nonzero on failure.
 */
c_int unload_linsys_solver(enum linsys_solver_type linsys_solver);


// NB: Only the upper triangular part of P is stuffed!

/**
 * Initialize linear system solver structure
 * @param P             Cost function matrix
 * @param	A             Constraint matrix
 * @param	sigma         Algorithm parameter
 * @param	rho_vec       Algorithm parameter
 * @param	polish        0/1 depending whether we are allocating for
 *polishing or not
 * @param	linsys_solver Linear system solver
 * @return Pointer to linear system solver structure on success, OSQP_NULL on
 *failure.
 */
LinSysSolver* init_linsys_solver(const csc              *P,
                                 const csc              *A,
                                 c_float                 sigma,
                                 c_float                *rho_vec,
                                 enum linsys_solver_type linsys_solver,
                                 c_int                   polish);


# endif // EMBEDDED


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef LIN_SYS_H
