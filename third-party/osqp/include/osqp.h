#ifndef OSQP_H
# define OSQP_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

/* Includes */
# include "types.h"
# include "util.h" // Needed for osqp_set_default_settings functions


// Library to deal with sparse matrices enabled only if embedded not defined
# ifndef EMBEDDED
#  include "cs.h"
# endif // ifndef EMBEDDED

/********************
* Main Solver API  *
********************/

/**
 * @name Main solver API
 * @{
 */

/**
 * Set default settings from constants.h file
 * assumes settings already allocated in memory
 * @param settings settings structure
 */
void osqp_set_default_settings(OSQPSettings *settings);


# ifndef EMBEDDED

/**
 * Initialize OSQP solver allocating memory.
 *
 * All the inputs must be already allocated in memory before calling.
 *
 * It performs:
 * - data and settings validation
 * - problem data scaling
 * - automatic parameters tuning (if enabled)
 * - setup linear system solver:
 *      - direct solver: KKT matrix factorization is performed here
 *      - indirect solver: KKT matrix preconditioning is performed here
 *
 * NB: This is the only function that allocates dynamic memory and is not used
 *during code generation
 *
 * @param  data         Problem data
 * @param  settings     Solver settings
 * @return              Solver environment
 */
OSQPWorkspace* osqp_setup(const OSQPData *data,
                          OSQPSettings   *settings);

# endif // #ifndef EMBEDDED

/**
 * Solve quadratic program
 *
 * The final solver information is stored in the \a work->info  structure
 *
 * The solution is stored in the  \a work->solution  structure
 *
 * If the problem is primal infeasible, the certificate is stored
 * in \a work->delta_y
 *
 * If the problem is dual infeasible, the certificate is stored in \a
 * work->delta_x
 *
 * @param  work Workspace allocated
 * @return      Exitflag for errors
 */
c_int osqp_solve(OSQPWorkspace *work);


# ifndef EMBEDDED

/**
 * Cleanup workspace by deallocating memory
 *
 * This function is not used in code generation
 * @param  work Workspace
 * @return      Exitflag for errors
 */
c_int osqp_cleanup(OSQPWorkspace *work);

# endif // ifndef EMBEDDED

/** @} */


/********************************************
* Sublevel API                             *
*                                          *
* Edit data without performing setup again *
********************************************/

/**
 * @name Sublevel API
 * @{
 */

/**
 * Update linear cost in the problem
 * @param  work  Workspace
 * @param  q_new New linear cost
 * @return       Exitflag for errors and warnings
 */
c_int osqp_update_lin_cost(OSQPWorkspace *work,
                           const c_float *q_new);


/**
 * Update lower and upper bounds in the problem constraints
 * @param  work   Workspace
 * @param  l_new New lower bound
 * @param  u_new New upper bound
 * @return        Exitflag: 1 if new lower bound is not <= than new upper bound
 */
c_int osqp_update_bounds(OSQPWorkspace *work,
                         const c_float *l_new,
                         const c_float *u_new);


/**
 * Update lower bound in the problem constraints
 * @param  work   Workspace
 * @param  l_new New lower bound
 * @return        Exitflag: 1 if new lower bound is not <= than upper bound
 */
c_int osqp_update_lower_bound(OSQPWorkspace *work,
                              const c_float *l_new);


/**
 * Update upper bound in the problem constraints
 * @param  work   Workspace
 * @param  u_new New upper bound
 * @return        Exitflag: 1 if new upper bound is not >= than lower bound
 */
c_int osqp_update_upper_bound(OSQPWorkspace *work,
                              const c_float *u_new);


/**
 * Warm start primal and dual variables
 * @param  work Workspace structure
 * @param  x    Primal variable
 * @param  y    Dual variable
 * @return      Exitflag
 */
c_int osqp_warm_start(OSQPWorkspace *work,
                      const c_float *x,
                      const c_float *y);


/**
 * Warm start primal variable
 * @param  work Workspace structure
 * @param  x    Primal variable
 * @return      Exitflag
 */
c_int osqp_warm_start_x(OSQPWorkspace *work,
                        const c_float *x);


/**
 * Warm start dual variable
 * @param  work Workspace structure
 * @param  y    Dual variable
 * @return      Exitflag
 */
c_int osqp_warm_start_y(OSQPWorkspace *work,
                        const c_float *y);


# if EMBEDDED != 1

/**
 * Update elements of matrix P (upper triangular)
 * without changing sparsity structure.
 *
 *
 *  If Px_new_idx is OSQP_NULL, Px_new is assumed to be as long as P->x
 *  and the whole P->x is replaced.
 *
 * @param  work       Workspace structure
 * @param  Px_new     Vector of new elements in P->x (upper triangular)
 * @param  Px_new_idx Index mapping new elements to positions in P->x
 * @param  P_new_n    Number of new elements to be changed
 * @return            output flag:  0: OK
 *                                  1: P_new_n > nnzP
 *                                 <0: error in the update
 */
c_int osqp_update_P(OSQPWorkspace *work,
                    const c_float *Px_new,
                    const c_int   *Px_new_idx,
                    c_int          P_new_n);


/**
 * Update elements of matrix A without changing sparsity structure.
 *
 *
 *  If Ax_new_idx is OSQP_NULL, Ax_new is assumed to be as long as A->x
 *  and the whole P->x is replaced.
 *
 * @param  work       Workspace structure
 * @param  Ax_new     Vector of new elements in A->x
 * @param  Ax_new_idx Index mapping new elements to positions in A->x
 * @param  A_new_n    Number of new elements to be changed
 * @return            output flag:  0: OK
 *                                  1: A_new_n > nnzA
 *                                 <0: error in the update
 */
c_int osqp_update_A(OSQPWorkspace *work,
                    const c_float *Ax_new,
                    const c_int   *Ax_new_idx,
                    c_int          A_new_n);


/**
 * Update elements of matrix P (upper triangular) and elements of matrix A
 * without changing sparsity structure.
 *
 *
 *  If Px_new_idx is OSQP_NULL, Px_new is assumed to be as long as P->x
 *  and the whole P->x is replaced.
 *
 *  If Ax_new_idx is OSQP_NULL, Ax_new is assumed to be as long as A->x
 *  and the whole P->x is replaced.
 *
 * @param  work       Workspace structure
 * @param  Px_new     Vector of new elements in P->x (upper triangular)
 * @param  Px_new_idx Index mapping new elements to positions in P->x
 * @param  P_new_n    Number of new elements to be changed
 * @param  Ax_new     Vector of new elements in A->x
 * @param  Ax_new_idx Index mapping new elements to positions in A->x
 * @param  A_new_n    Number of new elements to be changed
 * @return            output flag:  0: OK
 *                                  1: P_new_n > nnzP
 *                                  2: A_new_n > nnzA
 *                                 <0: error in the update
 */
c_int osqp_update_P_A(OSQPWorkspace *work,
                      const c_float *Px_new,
                      const c_int   *Px_new_idx,
                      c_int          P_new_n,
                      const c_float *Ax_new,
                      const c_int   *Ax_new_idx,
                      c_int          A_new_n);

/**
 * Update rho. Limit it between RHO_MIN and RHO_MAX.
 * @param  work         Workspace
 * @param  rho_new      New rho setting
 * @return              Exitflag
 */
c_int osqp_update_rho(OSQPWorkspace *work,
                      c_float        rho_new);

# endif // if EMBEDDED != 1

/** @} */


/**
 * @name Update settings
 * @{
 */


/**
 * Update max_iter setting
 * @param  work         Workspace
 * @param  max_iter_new New max_iter setting
 * @return              Exitflag
 */
c_int osqp_update_max_iter(OSQPWorkspace *work,
                           c_int          max_iter_new);


/**
 * Update absolute tolernace value
 * @param  work        Workspace
 * @param  eps_abs_new New absolute tolerance value
 * @return             Exitflag
 */
c_int osqp_update_eps_abs(OSQPWorkspace *work,
                          c_float        eps_abs_new);


/**
 * Update relative tolernace value
 * @param  work        Workspace
 * @param  eps_rel_new New relative tolerance value
 * @return             Exitflag
 */
c_int osqp_update_eps_rel(OSQPWorkspace *work,
                          c_float        eps_rel_new);


/**
 * Update primal infeasibility tolerance
 * @param  work          Workspace
 * @param  eps_prim_inf_new  New primal infeasibility tolerance
 * @return               Exitflag
 */
c_int osqp_update_eps_prim_inf(OSQPWorkspace *work,
                               c_float        eps_prim_inf_new);


/**
 * Update dual infeasibility tolerance
 * @param  work          Workspace
 * @param  eps_dual_inf_new  New dual infeasibility tolerance
 * @return               Exitflag
 */
c_int osqp_update_eps_dual_inf(OSQPWorkspace *work,
                               c_float        eps_dual_inf_new);


/**
 * Update relaxation parameter alpha
 * @param  work  Workspace
 * @param  alpha_new New relaxation parameter value
 * @return       Exitflag
 */
c_int osqp_update_alpha(OSQPWorkspace *work,
                        c_float        alpha_new);


/**
 * Update warm_start setting
 * @param  work           Workspace
 * @param  warm_start_new New warm_start setting
 * @return                Exitflag
 */
c_int osqp_update_warm_start(OSQPWorkspace *work,
                             c_int          warm_start_new);


/**
 * Update scaled_termination setting
 * @param  work                 Workspace
 * @param  scaled_termination_new  New scaled_termination setting
 * @return                      Exitflag
 */
c_int osqp_update_scaled_termination(OSQPWorkspace *work,
                                     c_int          scaled_termination_new);

/**
 * Update check_termination setting
 * @param  work                   Workspace
 * @param  check_termination_new  New check_termination setting
 * @return                        Exitflag
 */
c_int osqp_update_check_termination(OSQPWorkspace *work,
                                    c_int          check_termination_new);


# ifndef EMBEDDED

/**
 * Update regularization parameter in polish
 * @param  work      Workspace
 * @param  delta_new New regularization parameter
 * @return           Exitflag
 */
c_int osqp_update_delta(OSQPWorkspace *work,
                        c_float        delta_new);


/**
 * Update polish setting
 * @param  work          Workspace
 * @param  polish_new New polish setting
 * @return               Exitflag
 */
c_int osqp_update_polish(OSQPWorkspace *work,
                         c_int          polish_new);


/**
 * Update number of iterative refinement steps in polish
 * @param  work                Workspace
 * @param  polish_refine_iter_new New iterative reginement steps
 * @return                     Exitflag
 */
c_int osqp_update_polish_refine_iter(OSQPWorkspace *work,
                                     c_int          polish_refine_iter_new);


/**
 * Update verbose setting
 * @param  work        Workspace
 * @param  verbose_new New verbose setting
 * @return             Exitflag
 */
c_int osqp_update_verbose(OSQPWorkspace *work,
                          c_int          verbose_new);


# endif // #ifndef EMBEDDED

# ifdef PROFILING

/**
 * Update time_limit setting
 * @param  work            Workspace
 * @param  time_limit_new  New time_limit setting
 * @return                 Exitflag
 */
c_int osqp_update_time_limit(OSQPWorkspace *work,
                             c_float        time_limit_new);
# endif // ifdef PROFILING

/** @} */


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef OSQP_H
