#ifndef CONSTANTS_H
# define CONSTANTS_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


/*******************
* OSQP Versioning *
*******************/
# define OSQP_VERSION ("0.5.0") /* string literals automatically null-terminated
                                   */


/******************
* Solver Status  *
******************/
# define OSQP_DUAL_INFEASIBLE_INACCURATE (4)
# define OSQP_PRIMAL_INFEASIBLE_INACCURATE (3)
# define OSQP_SOLVED_INACCURATE (2)
# define OSQP_SOLVED (1)
# define OSQP_MAX_ITER_REACHED (-2)
# define OSQP_PRIMAL_INFEASIBLE (-3) /* primal infeasible  */
# define OSQP_DUAL_INFEASIBLE (-4)   /* dual infeasible */
# define OSQP_SIGINT (-5)            /* interrupted by user */
# ifdef PROFILING
#  define OSQP_TIME_LIMIT_REACHED (-6)
# endif // ifdef PROFILING
# define OSQP_NON_CVX (-7)           /* problem non convex */
# define OSQP_UNSOLVED (-10) /* Unsolved. Only setup function has been called */


/*************************
* Linear System Solvers *
*************************/
enum linsys_solver_type { QDLDL_SOLVER, MKL_PARDISO_SOLVER };
__attribute__((unused))
static const char *LINSYS_SOLVER_NAME[] = {
  "qdldl", "mkl pardiso"
};

/**********************************
* Solver Parameters and Settings *
**********************************/

# define RHO (0.1)
# define SIGMA (1E-06)
# define MAX_ITER (4000)
# define EPS_ABS (1E-3)
# define EPS_REL (1E-3)
# define EPS_PRIM_INF (1E-4)
# define EPS_DUAL_INF (1E-4)
# define ALPHA (1.6)
# define LINSYS_SOLVER (QDLDL_SOLVER)

# define RHO_MIN (1e-06)
# define RHO_MAX (1e06)
# define RHO_EQ_OVER_RHO_INEQ (1e03)
# define RHO_TOL (1e-04) ///< Tolerance for detecting if an inequality is set to
                         // equality


# ifndef EMBEDDED
#  define DELTA (1E-6)
#  define POLISH (0)
#  define POLISH_REFINE_ITER (3)
#  define VERBOSE (1)
# endif // ifndef EMBEDDED

# define SCALED_TERMINATION (0)
# define CHECK_TERMINATION (25)
# define WARM_START (1)
# define SCALING (10)

# define MIN_SCALING (1e-04) ///< Minimum scaling value
# define MAX_SCALING (1e+04) ///< Maximum scaling value


# ifndef OSQP_NULL
#  define OSQP_NULL 0
# endif /* ifndef OSQP_NULL */

# ifndef OSQP_NAN
#  define OSQP_NAN ((c_float)0x7fc00000UL) // Not a Number
# endif /* ifndef OSQP_NAN */

# ifndef OSQP_INFTY
#  define OSQP_INFTY ((c_float)1e20) // Infinity
# endif /* ifndef OSQP_INFTY */


# if EMBEDDED != 1
#  define ADAPTIVE_RHO (1)
#  define ADAPTIVE_RHO_INTERVAL (0)
#  define ADAPTIVE_RHO_FRACTION (0.4)           ///< Fraction of setup time
                                                // after which we update rho
#  define ADAPTIVE_RHO_MULTIPLE_TERMINATION (4) ///< Multiple of termination
                                                // check time we update rho if
                                                // profiling disabled
#  define ADAPTIVE_RHO_FIXED (100)              ///< Number of iterations after
                                                // which we update rho if
                                                // termination check if disabled
                                                // and profiling disabled
#  define ADAPTIVE_RHO_TOLERANCE (5)            ///< Tolerance for adopting new
                                                // rho. Number of times the new
                                                // rho is larger or smaller than
                                                // the current one
# endif // if EMBEDDED != 1

# ifdef PROFILING
#  define TIME_LIMIT (0) ///< Disable time limit as default
# endif // ifdef PROFILING

/* Printing */
# define PRINT_INTERVAL 200


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef CONSTANTS_H
