#ifndef OSQP_TYPES_H
# define OSQP_TYPES_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

# include "glob_opts.h"
# include "constants.h"


/******************
* Internal types *
******************/

/**
 *  Matrix in compressed-column or triplet form
 */
typedef struct {
  c_int    nzmax; ///< maximum number of entries.
  c_int    m;     ///< number of rows
  c_int    n;     ///< number of columns
  c_int   *p;     ///< column pointers (size n+1) (col indices (size nzmax)
                  // start from 0 when using triplet format (direct KKT matrix
                  // formation))
  c_int   *i;     ///< row indices, size nzmax starting from 0
  c_float *x;     ///< numerical values, size nzmax
  c_int    nz;    ///< # of entries in triplet matrix, -1 for csc
} csc;

/**
 * Linear system solver structure (sublevel objects initialize it differently)
 */

typedef struct linsys_solver LinSysSolver;

/**
 * OSQP Timer for statistics
 */
typedef struct OSQP_TIMER OSQPTimer;

/**
 * Problem scaling matrices stored as vectors
 */
typedef struct {
  c_float  c;    ///< cost function scaling
  c_float *D;    ///< primal variable scaling
  c_float *E;    ///< dual variable scaling
  c_float  cinv; ///< cost function rescaling
  c_float *Dinv; ///< primal variable rescaling
  c_float *Einv; ///< dual variable rescaling
} OSQPScaling;

/**
 * Solution structure
 */
typedef struct {
  c_float *x; ///< Primal solution
  c_float *y; ///< Lagrange multiplier associated to \f$l <= Ax <= u\f$
} OSQPSolution;


/**
 * Solver return information
 */
typedef struct {
  c_int iter;          ///< number of iterations taken
  char  status[32];    ///< status string, e.g. 'solved'
  c_int status_val;    ///< status as c_int, defined in constants.h

# ifndef EMBEDDED
  c_int status_polish; ///< polish status: successful (1), unperformed (0), (-1)
                       // unsuccessful
# endif // ifndef EMBEDDED

  c_float obj_val;     ///< primal objective
  c_float pri_res;     ///< norm of primal residual
  c_float dua_res;     ///< norm of dual residual

# ifdef PROFILING
  c_float setup_time;  ///< time taken for setup phase (seconds)
  c_float solve_time;  ///< time taken for solve phase (seconds)
  c_float update_time; ///< time taken for update phase (seconds)
  c_float polish_time; ///< time taken for polish phase (seconds)
  c_float run_time;    ///< total time  (seconds)
# endif // ifdef PROFILING

# if EMBEDDED != 1
  c_int   rho_updates;  ///< number of rho updates
  c_float rho_estimate; ///< best rho estimate so far from residuals
# endif // if EMBEDDED != 1
} OSQPInfo;


# ifndef EMBEDDED

/**
 * Polish structure
 */
typedef struct {
  csc *Ared;          ///< Active rows of A.
  ///<    Ared = vstack[Alow, Aupp]
  c_int    n_low;     ///< number of lower-active rows
  c_int    n_upp;     ///< number of upper-active rows
  c_int   *A_to_Alow; ///< Maps indices in A to indices in Alow
  c_int   *A_to_Aupp; ///< Maps indices in A to indices in Aupp
  c_int   *Alow_to_A; ///< Maps indices in Alow to indices in A
  c_int   *Aupp_to_A; ///< Maps indices in Aupp to indices in A
  c_float *x;         ///< optimal x-solution obtained by polish
  c_float *z;         ///< optimal z-solution obtained by polish
  c_float *y;         ///< optimal y-solution obtained by polish
  c_float  obj_val;   ///< objective value at polished solution
  c_float  pri_res;   ///< primal residual at polished solution
  c_float  dua_res;   ///< dual residual at polished solution
} OSQPPolish;
# endif // ifndef EMBEDDED


/**********************************
* Main structures and Data Types *
**********************************/

/**
 * Data structure
 */
typedef struct {
  c_int    n; ///< number of variables n
  c_int    m; ///< number of constraints m
  csc     *P; ///< quadratic part of the cost P in csc format (size n x n). It
              ///  can be either the full P or only the upper triangular part. The
              ///  workspace stores only the upper triangular part
  csc     *A; ///< linear constraints matrix A in csc format (size m x n)
  c_float *q; ///< dense array for linear part of cost function (size n)
  c_float *l; ///< dense array for lower bound (size m)
  c_float *u; ///< dense array for upper bound (size m)
} OSQPData;


/**
 * Settings struct
 */
typedef struct {
  c_float rho;                    ///< ADMM step rho
  c_float sigma;                  ///< ADMM step sigma
  c_int   scaling;                ///< heuristic data scaling iterations. If 0,
                                  // scaling disabled

# if EMBEDDED != 1
  c_int   adaptive_rho;           ///< boolean, is rho step size adaptive?
  c_int   adaptive_rho_interval;  ///< Number of iterations between rho
                                  // adaptations rho. If 0, it is automatic
  c_float adaptive_rho_tolerance; ///< Tolerance X for adapting rho. The new rho
                                  // has to be X times larger or 1/X times
                                  // smaller than the current one to trigger a
                                  // new factorization.
#  ifdef PROFILING
  c_float adaptive_rho_fraction;  ///< Interval for adapting rho (fraction of
                                  // the setup time)
#  endif // Profiling
# endif // EMBEDDED != 1

  c_int                   max_iter;      ///< maximum iterations
  c_float                 eps_abs;       ///< absolute convergence tolerance
  c_float                 eps_rel;       ///< relative convergence tolerance
  c_float                 eps_prim_inf;  ///< primal infeasibility tolerance
  c_float                 eps_dual_inf;  ///< dual infeasibility tolerance
  c_float                 alpha;         ///< relaxation parameter
  enum linsys_solver_type linsys_solver; ///< linear system solver to use

# ifndef EMBEDDED
  c_float delta;                         ///< regularization parameter for
                                         // polish
  c_int   polish;                        ///< boolean, polish ADMM solution
  c_int   polish_refine_iter;            ///< iterative refinement steps in
                                         // polish

  c_int verbose;                         ///< boolean, write out progres
# endif // ifndef EMBEDDED

  c_int scaled_termination;              ///< boolean, use scaled termination
                                         // criteria
  c_int check_termination;               ///< integer, check termination
                                         // interval. If 0, termination checking
                                         // is disabled
  c_int warm_start;                      ///< boolean, warm start

# ifdef PROFILING
  c_float time_limit;                    ///< maximum seconds allowed to solve
                                         // the problem
# endif // ifdef PROFILING
} OSQPSettings;


/**
 * OSQP Workspace
 */
typedef struct {
  /// Problem data to work on (possibly scaled)
  OSQPData *data;

  /// Linear System solver structure
  LinSysSolver *linsys_solver;

        # ifndef EMBEDDED

  /// Polish structure
  OSQPPolish *pol;
        # endif // ifndef EMBEDDED

  /**
   * @name Vector used to store a vectorized rho parameter
   * @{
   */
  c_float *rho_vec;     ///< vector of rho values
  c_float *rho_inv_vec; ///< vector of inv rho values

  /** @} */

        # if EMBEDDED != 1
  c_int *constr_type; ///< Type of constraints: loose (-1), equality (1),
                      // inequality (0)
        # endif // if EMBEDDED != 1

  /**
   * @name Iterates
   * @{
   */
  c_float *x;        ///< Iterate x
  c_float *y;        ///< Iterate y
  c_float *z;        ///< Iterate z
  c_float *xz_tilde; ///< Iterate xz_tilde

  c_float *x_prev;   ///< Previous x

  /**< NB: Used also as workspace vector for dual residual */
  c_float *z_prev;   ///< Previous z

  /**< NB: Used also as workspace vector for primal residual */

  /**
   * @name Primal and dual residuals workspace variables
   *
   * Needed for residuals computation, tolerances computation,
   * approximate tolerances computation and adapting rho
   * @{
   */
  c_float *Ax;  ///< Scaled A * x
  c_float *Px;  ///< Scaled P * x
  c_float *Aty; ///< Scaled A * x

  /** @} */

  /**
   * @name Primal infeasibility variables
   * @{
   */
  c_float *delta_y;   ///< Difference of consecutive dual iterates
  c_float *Atdelta_y; ///< A' * delta_y

  /** @} */

  /**
   * @name Dual infeasibility variables
   * @{
   */
  c_float *delta_x;  ///< Difference of consecutive primal iterates
  c_float *Pdelta_x; ///< P * delta_x
  c_float *Adelta_x; ///< A * delta_x

  /** @} */

  /**
   * @name Temporary vectors used in scaling
   * @{
   */

  c_float *D_temp;   ///< temporary primal variable scaling vectors
  c_float *D_temp_A; ///< temporary primal variable scaling vectors storing
                     // norms of A columns
  c_float *E_temp;   ///< temporary constraints scaling vectors storing norms of
                     // A' columns


  /** @} */

  OSQPSettings *settings; ///< Problem settings
  OSQPScaling  *scaling;  ///< Scaling vectors
  OSQPSolution *solution; ///< Problem solution
  OSQPInfo     *info;     ///< Solver information

# ifdef PROFILING
  OSQPTimer *timer;       ///< Timer object

  /// flag indicating whether the solve function has been run before
  c_int first_run;

  /// flag indicating whether the update_time should be cleared
  c_int clear_update_time;
# endif // ifdef PROFILING

# ifdef PRINTING
  c_int summary_printed; ///< Has last summary been printed? (true/false)
# endif // ifdef PRINTING

} OSQPWorkspace;


/**
 * Define linsys_solver prototype structure
 *
 * NB: The details are defined when the linear solver is initialized depending
 *      on the choice
 */
struct linsys_solver {
  enum linsys_solver_type type; ///< Linear system solver type (see type.h)
  // Functions
  c_int (*solve)(LinSysSolver       *self,
                 c_float            *b,
                 const OSQPSettings *settings); ///< Solve linear system

    # ifndef EMBEDDED
  void (*free)(LinSysSolver *self);             ///< Free linear system solver
                                                // (only in desktop version)
    # endif // ifndef EMBEDDED

    # if EMBEDDED != 1
  c_int (*update_matrices)(LinSysSolver *self, const csc *P, const csc *A,
                           const OSQPSettings *settings); ///< Update matrices P
                                                          // and A in the solver
  c_int (*update_rho_vec)(LinSysSolver  *s,
                          const c_float *rho_vec,
                          const c_int    m);              ///< Update rho
    # endif // if EMBEDDED != 1

# ifndef EMBEDDED
  c_int nthreads; ///< Number of threads active
# endif // ifndef EMBEDDED
};


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef OSQP_TYPES_H
