#include "osqp.h" // For OSQP rho update
#include "auxil.h"
#include "proj.h"
#include "lin_alg.h"
#include "constants.h"
#include "scaling.h"
#include "util.h"

/***********************************************************
* Auxiliary functions needed to compute ADMM iterations * *
***********************************************************/
#if EMBEDDED != 1
c_float compute_rho_estimate(OSQPWorkspace *work) {
  c_int   n, m;                       // Dimensions
  c_float pri_res, dua_res;           // Primal and dual residuals
  c_float pri_res_norm, dua_res_norm; // Normalization for the residuals
  c_float temp_res_norm;              // Temporary residual norm
  c_float rho_estimate;               // Rho estimate value

  // Get problem dimensions
  n = work->data->n;
  m = work->data->m;

  // Get primal and dual residuals
  pri_res = vec_norm_inf(work->z_prev, m);
  dua_res = vec_norm_inf(work->x_prev, n);

  // Normalize primal residual
  pri_res_norm  = vec_norm_inf(work->z, m);           // ||z||
  temp_res_norm = vec_norm_inf(work->Ax, m);          // ||Ax||
  pri_res_norm  = c_max(pri_res_norm, temp_res_norm); // max (||z||,||Ax||)
  pri_res      /= (pri_res_norm + 1e-10);             // Normalize primal
                                                      // residual (prevent 0
                                                      // division)

  // Normalize dual residual
  dua_res_norm  = vec_norm_inf(work->data->q, n);     // ||q||
  temp_res_norm = vec_norm_inf(work->Aty, n);         // ||A' y||
  dua_res_norm  = c_max(dua_res_norm, temp_res_norm);
  temp_res_norm = vec_norm_inf(work->Px, n);          //  ||P x||
  dua_res_norm  = c_max(dua_res_norm, temp_res_norm); // max(||q||,||A' y||,||P
                                                      // x||)
  dua_res      /= (dua_res_norm + 1e-10);             // Normalize dual residual
                                                      // (prevent 0 division)


  // Return rho estimate
  rho_estimate = work->settings->rho * c_sqrt(pri_res / (dua_res + 1e-10)); // (prevent
                                                                            // 0
                                                                            // division)
  rho_estimate = c_min(c_max(rho_estimate, RHO_MIN), RHO_MAX);              // Constrain
                                                                            // rho
                                                                            // values
  return rho_estimate;
}

c_int adapt_rho(OSQPWorkspace *work) {
  c_int   exitflag; // Exitflag
  c_float rho_new;  // New rho value

  exitflag = 0;     // Initialize exitflag to 0

  // Compute new rho
  rho_new = compute_rho_estimate(work);

  // Set rho estimate in info
  work->info->rho_estimate = rho_new;

  // Check if the new rho is large or small enough and update it in case
  if ((rho_new > work->settings->rho * ADAPTIVE_RHO_TOLERANCE) ||
      (rho_new < work->settings->rho /  ADAPTIVE_RHO_TOLERANCE)) {
    exitflag                 = osqp_update_rho(work, rho_new);
    work->info->rho_updates += 1;
  }

  return exitflag;
}

void set_rho_vec(OSQPWorkspace *work) {
  c_int i;

  work->settings->rho = c_min(c_max(work->settings->rho, RHO_MIN), RHO_MAX);

  for (i = 0; i < work->data->m; i++) {
    if ((work->data->l[i] < -OSQP_INFTY * MIN_SCALING) &&
        (work->data->u[i] > OSQP_INFTY * MIN_SCALING)) {
      // Loose bounds
      work->constr_type[i] = -1;
      work->rho_vec[i]     = RHO_MIN;
    } else if (work->data->u[i] - work->data->l[i] < RHO_TOL) {
      // Equality constraints
      work->constr_type[i] = 1;
      work->rho_vec[i]     = RHO_EQ_OVER_RHO_INEQ * work->settings->rho;
    } else {
      // Inequality constraints
      work->constr_type[i] = 0;
      work->rho_vec[i]     = work->settings->rho;
    }
    work->rho_inv_vec[i] = 1. / work->rho_vec[i];
  }
}

c_int update_rho_vec(OSQPWorkspace *work) {
  c_int i, exitflag, constr_type_changed;

  exitflag            = 0;
  constr_type_changed = 0;

  for (i = 0; i < work->data->m; i++) {
    if ((work->data->l[i] < -OSQP_INFTY * MIN_SCALING) &&
        (work->data->u[i] > OSQP_INFTY * MIN_SCALING)) {
      // Loose bounds
      if (work->constr_type[i] != -1) {
        work->constr_type[i] = -1;
        work->rho_vec[i]     = RHO_MIN;
        work->rho_inv_vec[i] = 1. / RHO_MIN;
        constr_type_changed  = 1;
      }
    } else if (work->data->u[i] - work->data->l[i] < RHO_TOL) {
      // Equality constraints
      if (work->constr_type[i] != 1) {
        work->constr_type[i] = 1;
        work->rho_vec[i]     = RHO_EQ_OVER_RHO_INEQ * work->settings->rho;
        work->rho_inv_vec[i] = 1. / work->rho_vec[i];
        constr_type_changed  = 1;
      }
    } else {
      // Inequality constraints
      if (work->constr_type[i] != 0) {
        work->constr_type[i] = 0;
        work->rho_vec[i]     = work->settings->rho;
        work->rho_inv_vec[i] = 1. / work->settings->rho;
        constr_type_changed  = 1;
      }
    }
  }

  // Update rho_vec in KKT matrix if constraints type has changed
  if (constr_type_changed == 1) {
    exitflag = work->linsys_solver->update_rho_vec(work->linsys_solver,
                                                   work->rho_vec,
                                                   work->data->m);
  }

  return exitflag;
}

#endif // EMBEDDED != 1


void swap_vectors(c_float **a, c_float **b) {
  c_float *temp;

  temp = *b;
  *b   = *a;
  *a   = temp;
}

void cold_start(OSQPWorkspace *work) {
  vec_set_scalar(work->x, 0., work->data->n);
  vec_set_scalar(work->z, 0., work->data->m);
  vec_set_scalar(work->y, 0., work->data->m);
}

static void compute_rhs(OSQPWorkspace *work) {
  c_int i; // Index

  for (i = 0; i < work->data->n; i++) {
    // Cycle over part related to x variables
    work->xz_tilde[i] = work->settings->sigma * work->x_prev[i] -
                        work->data->q[i];
  }

  for (i = 0; i < work->data->m; i++) {
    // Cycle over dual variable in the first step (nu)
    work->xz_tilde[i + work->data->n] = work->z_prev[i] - work->rho_inv_vec[i] *
                                        work->y[i];
  }
}

static void update_z_tilde(OSQPWorkspace *work) {
  c_int i; // Index

  for (i = 0; i < work->data->m; i++) {
    work->xz_tilde[i + work->data->n] = work->z_prev[i] + work->rho_inv_vec[i] *
                                        (work->xz_tilde[i + work->data->n] -
                                         work->y[i]);
  }
}

void update_xz_tilde(OSQPWorkspace *work) {
  // Compute right-hand side
  compute_rhs(work);

  // Solve linear system
  work->linsys_solver->solve(work->linsys_solver, work->xz_tilde, work->settings);

  // Update z_tilde variable after solving linear system
  update_z_tilde(work);
}

void update_x(OSQPWorkspace *work) {
  c_int i;

  // update x
  for (i = 0; i < work->data->n; i++) {
    work->x[i] = work->settings->alpha * work->xz_tilde[i] +
                 ((c_float)1.0 - work->settings->alpha) * work->x_prev[i];
  }

  // update delta_x
  for (i = 0; i < work->data->n; i++) {
    work->delta_x[i] = work->x[i] - work->x_prev[i];
  }
}

void update_z(OSQPWorkspace *work) {
  c_int i;

  // update z
  for (i = 0; i < work->data->m; i++) {
    work->z[i] = work->settings->alpha * work->xz_tilde[i + work->data->n] +
                 ((c_float)1.0 - work->settings->alpha) * work->z_prev[i] +
                 work->rho_inv_vec[i] * work->y[i];
  }

  // project z
  project(work, work->z);
}

void update_y(OSQPWorkspace *work) {
  c_int i; // Index

  for (i = 0; i < work->data->m; i++) {
    work->delta_y[i] = work->rho_vec[i] *
                       (work->settings->alpha *
                        work->xz_tilde[i + work->data->n] +
                        ((c_float)1.0 - work->settings->alpha) * work->z_prev[i] -
                        work->z[i]);
    work->y[i] += work->delta_y[i];
  }
}

c_float compute_obj_val(OSQPWorkspace *work, c_float *x) {
  c_float obj_val;

  obj_val = quad_form(work->data->P, x) +
            vec_prod(work->data->q, x, work->data->n);

  if (work->settings->scaling) {
    obj_val *= work->scaling->cinv;
  }

  return obj_val;
}

c_float compute_pri_res(OSQPWorkspace *work, c_float *x, c_float *z) {
  // NB: Use z_prev as working vector
  // pr = Ax - z

  mat_vec(work->data->A, x, work->Ax, 0); // Ax
  vec_add_scaled(work->z_prev, work->Ax, z, work->data->m, -1);

  // If scaling active -> rescale residual
  if (work->settings->scaling && !work->settings->scaled_termination) {
    return vec_scaled_norm_inf(work->scaling->Einv, work->z_prev, work->data->m);
  }

  // Return norm of the residual
  return vec_norm_inf(work->z_prev, work->data->m);
}

c_float compute_pri_tol(OSQPWorkspace *work, c_float eps_abs, c_float eps_rel) {
  c_float max_rel_eps, temp_rel_eps;

  // max_rel_eps = max(||z||, ||A x||)
  if (work->settings->scaling && !work->settings->scaled_termination) {
    // ||Einv * z||
    max_rel_eps =
      vec_scaled_norm_inf(work->scaling->Einv, work->z, work->data->m);

    // ||Einv * A * x||
    temp_rel_eps = vec_scaled_norm_inf(work->scaling->Einv,
                                       work->Ax,
                                       work->data->m);

    // Choose maximum
    max_rel_eps = c_max(max_rel_eps, temp_rel_eps);
  } else { // No unscaling required
    // ||z||
    max_rel_eps = vec_norm_inf(work->z, work->data->m);

    // ||A * x||
    temp_rel_eps = vec_norm_inf(work->Ax, work->data->m);

    // Choose maximum
    max_rel_eps = c_max(max_rel_eps, temp_rel_eps);
  }

  // eps_prim
  return eps_abs + eps_rel * max_rel_eps;
}

c_float compute_dua_res(OSQPWorkspace *work, c_float *x, c_float *y) {
  // NB: Use x_prev as temporary vector
  // NB: Only upper triangular part of P is stored.
  // dr = q + A'*y + P*x

  // dr = q
  prea_vec_copy(work->data->q, work->x_prev, work->data->n);

  // P * x (upper triangular part)
  mat_vec(work->data->P, x, work->Px, 0);

  // P' * x (lower triangular part with no diagonal)
  mat_tpose_vec(work->data->P, x, work->Px, 1, 1);

  // dr += P * x (full P matrix)
  vec_add_scaled(work->x_prev, work->x_prev, work->Px, work->data->n, 1);

  // dr += A' * y
  if (work->data->m > 0) {
    mat_tpose_vec(work->data->A, y, work->Aty, 0, 0);
    vec_add_scaled(work->x_prev, work->x_prev, work->Aty, work->data->n, 1);
  }

  // If scaling active -> rescale residual
  if (work->settings->scaling && !work->settings->scaled_termination) {
    return work->scaling->cinv * vec_scaled_norm_inf(work->scaling->Dinv,
                                                     work->x_prev,
                                                     work->data->n);
  }

  return vec_norm_inf(work->x_prev, work->data->n);
}

c_float compute_dua_tol(OSQPWorkspace *work, c_float eps_abs, c_float eps_rel) {
  c_float max_rel_eps, temp_rel_eps;

  // max_rel_eps = max(||q||, ||A' y|, ||P x||)
  if (work->settings->scaling && !work->settings->scaled_termination) {
    // || Dinv q||
    max_rel_eps = vec_scaled_norm_inf(work->scaling->Dinv,
                                      work->data->q,
                                      work->data->n);

    // || Dinv A' y ||
    temp_rel_eps = vec_scaled_norm_inf(work->scaling->Dinv,
                                       work->Aty,
                                       work->data->n);
    max_rel_eps = c_max(max_rel_eps, temp_rel_eps);

    // || Dinv P x||
    temp_rel_eps = vec_scaled_norm_inf(work->scaling->Dinv,
                                       work->Px,
                                       work->data->n);
    max_rel_eps = c_max(max_rel_eps, temp_rel_eps);

    // Multiply by cinv
    max_rel_eps *= work->scaling->cinv;
  } else { // No scaling required
    // ||q||
    max_rel_eps = vec_norm_inf(work->data->q, work->data->n);

    // ||A'*y||
    temp_rel_eps = vec_norm_inf(work->Aty, work->data->n);
    max_rel_eps  = c_max(max_rel_eps, temp_rel_eps);

    // ||P*x||
    temp_rel_eps = vec_norm_inf(work->Px, work->data->n);
    max_rel_eps  = c_max(max_rel_eps, temp_rel_eps);
  }

  // eps_dual
  return eps_abs + eps_rel * max_rel_eps;
}

c_int is_primal_infeasible(OSQPWorkspace *work, c_float eps_prim_inf) {
  // This function checks for the primal infeasibility termination criteria.
  //
  // 1) A' * delta_y < eps * ||delta_y||
  //
  // 2) u'*max(delta_y, 0) + l'*min(delta_y, 0) < -eps * ||delta_y||
  //

  c_int   i; // Index for loops
  c_float norm_delta_y;
  c_float ineq_lhs;

  // Compute infinity norm of delta_y
  if (work->settings->scaling && !work->settings->scaled_termination) { // Unscale
                                                                        // if
                                                                        // necessary
    // Use work->Adelta_x as temporary vector
    vec_ew_prod(work->scaling->E, work->delta_y, work->Adelta_x, work->data->m);
    norm_delta_y = vec_norm_inf(work->Adelta_x, work->data->m);
  } else {
    norm_delta_y = vec_norm_inf(work->delta_y, work->data->m);
  }

  if (norm_delta_y > eps_prim_inf) { // ||delta_y|| > 0
    ineq_lhs = 0;

    for (i = 0; i < work->data->m; i++) {
      ineq_lhs += work->data->u[i] * c_max(work->delta_y[i], 0) + \
                  work->data->l[i] * c_min(work->delta_y[i], 0);
    }

    // Check if the condition is satisfied: ineq_lhs < -eps
    if (ineq_lhs < -eps_prim_inf * norm_delta_y) {
      // Compute and return ||A'delta_y|| < eps_prim_inf
      mat_tpose_vec(work->data->A, work->delta_y, work->Atdelta_y, 0, 0);

      if (work->settings->scaling && !work->settings->scaled_termination) { // Unscale
                                                                            // if
                                                                            // necessary
        vec_ew_prod(work->scaling->Dinv,
                    work->Atdelta_y,
                    work->Atdelta_y,
                    work->data->n);
      }
      return vec_norm_inf(work->Atdelta_y,
                          work->data->n) < eps_prim_inf * norm_delta_y;
    }
  }

  // Conditions not satisfied -> not primal infeasible
  return 0;
}

c_int is_dual_infeasible(OSQPWorkspace *work, c_float eps_dual_inf) {
  // This function checks for the scaled dual infeasibility termination
  // criteria.
  //
  // 1) q * delta_x < - eps * || delta_x ||
  //
  // 2) ||P * delta_x || < eps * || delta_x ||
  //
  // 3) -> (A * delta_x)_i > -eps * || delta_x ||,    l_i != -inf
  //    -> (A * delta_x)_i <  eps * || delta_x ||,    u_i != inf
  //


  c_int   i; // Index for loops
  c_float norm_delta_x;
  c_float cost_scaling;

  // Compute norm of delta_x
  if (work->settings->scaling && !work->settings->scaled_termination) { // Unscale
                                                                        // if
                                                                        // necessary
    norm_delta_x = vec_scaled_norm_inf(work->scaling->D,
                                       work->delta_x,
                                       work->data->n);
    cost_scaling = work->scaling->c;
  } else {
    norm_delta_x = vec_norm_inf(work->delta_x, work->data->n);
    cost_scaling = 1.0;
  }

  // Prevent 0 division || delta_x || > 0
  if (norm_delta_x > eps_dual_inf) {
    // Normalize delta_x by its norm

    /* vec_mult_scalar(work->delta_x, 1./norm_delta_x, work->data->n); */

    // Check first if q'*delta_x < 0
    if (vec_prod(work->data->q, work->delta_x, work->data->n) <
        -cost_scaling * eps_dual_inf * norm_delta_x) {
      // Compute product P * delta_x (NB: P is store in upper triangular form)
      mat_vec(work->data->P, work->delta_x, work->Pdelta_x, 0);
      mat_tpose_vec(work->data->P, work->delta_x, work->Pdelta_x, 1, 1);

      // Scale if necessary
      if (work->settings->scaling && !work->settings->scaled_termination) {
        vec_ew_prod(work->scaling->Dinv,
                    work->Pdelta_x,
                    work->Pdelta_x,
                    work->data->n);
      }

      // Check if || P * delta_x || = 0
      if (vec_norm_inf(work->Pdelta_x, work->data->n) <
          cost_scaling * eps_dual_inf * norm_delta_x) {
        // Compute A * delta_x
        mat_vec(work->data->A, work->delta_x, work->Adelta_x, 0);

        // Scale if necessary
        if (work->settings->scaling && !work->settings->scaled_termination) {
          vec_ew_prod(work->scaling->Einv,
                      work->Adelta_x,
                      work->Adelta_x,
                      work->data->m);
        }

        // De Morgan Law Applied to dual infeasibility conditions for A * x
        // NB: Note that 1e-06 is used to adjust the infinity value
        //      in case the problem is scaled.
        for (i = 0; i < work->data->m; i++) {
          if (((work->data->u[i] < OSQP_INFTY * MIN_SCALING) &&
               (work->Adelta_x[i] >  eps_dual_inf * norm_delta_x)) ||
              ((work->data->l[i] > -OSQP_INFTY * MIN_SCALING) &&
               (work->Adelta_x[i] < -eps_dual_inf * norm_delta_x))) {
            // At least one condition not satisfied -> not dual infeasible
            return 0;
          }
        }

        // All conditions passed -> dual infeasible
        return 1;
      }
    }
  }

  // Conditions not satisfied -> not dual infeasible
  return 0;
}

static c_int has_solution(OSQPInfo * info){

  return ((info->status_val != OSQP_PRIMAL_INFEASIBLE) &&
      (info->status_val != OSQP_PRIMAL_INFEASIBLE_INACCURATE) &&
      (info->status_val != OSQP_DUAL_INFEASIBLE) &&
      (info->status_val != OSQP_DUAL_INFEASIBLE_INACCURATE) &&
      (info->status_val != OSQP_NON_CVX));

}

void store_solution(OSQPWorkspace *work) {
#ifndef EMBEDDED
  c_float norm_vec;
#endif /* ifndef EMBEDDED */

  if (has_solution(work->info)) {
    prea_vec_copy(work->x, work->solution->x, work->data->n); // primal
    prea_vec_copy(work->y, work->solution->y, work->data->m); // dual

    // Unscale solution if scaling has been performed
    if (work->settings->scaling)
      unscale_solution(work);
  } else {
    // No solution present. Solution is NaN
    vec_set_scalar(work->solution->x, OSQP_NAN, work->data->n);
    vec_set_scalar(work->solution->y, OSQP_NAN, work->data->m);

#ifndef EMBEDDED

    // Normalize infeasibility certificates if embedded is off
    // NB: It requires a division
    if ((work->info->status_val == OSQP_PRIMAL_INFEASIBLE) ||
        ((work->info->status_val == OSQP_PRIMAL_INFEASIBLE_INACCURATE))) {
      norm_vec = vec_norm_inf(work->delta_y, work->data->m);
      vec_mult_scalar(work->delta_y, 1. / norm_vec, work->data->m);
    }

    if ((work->info->status_val == OSQP_DUAL_INFEASIBLE) ||
        ((work->info->status_val == OSQP_DUAL_INFEASIBLE_INACCURATE))) {
      norm_vec = vec_norm_inf(work->delta_x, work->data->n);
      vec_mult_scalar(work->delta_x, 1. / norm_vec, work->data->n);
    }

#endif /* ifndef EMBEDDED */

    // Cold start iterates to 0 for next runs (they cannot start from NaN)
    cold_start(work);
  }
}

void update_info(OSQPWorkspace *work,
                 c_int          iter,
                 c_int          compute_objective,
                 c_int          polish) {
  c_float *x, *z, *y;                   // Allocate pointers to variables
  c_float *obj_val, *pri_res, *dua_res; // objective value, residuals

#ifdef PROFILING
  c_float *run_time;                    // Execution time
#endif /* ifdef PROFILING */

#ifndef EMBEDDED

  if (polish) {
    x       = work->pol->x;
    y       = work->pol->y;
    z       = work->pol->z;
    obj_val = &work->pol->obj_val;
    pri_res = &work->pol->pri_res;
    dua_res = &work->pol->dua_res;
# ifdef PROFILING
    run_time = &work->info->polish_time;
# endif /* ifdef PROFILING */
  } else {
#endif // EMBEDDED
  x                = work->x;
  y                = work->y;
  z                = work->z;
  obj_val          = &work->info->obj_val;
  pri_res          = &work->info->pri_res;
  dua_res          = &work->info->dua_res;
  work->info->iter = iter; // Update iteration number
#ifdef PROFILING
  run_time = &work->info->solve_time;
#endif /* ifdef PROFILING */
#ifndef EMBEDDED
}

#endif /* ifndef EMBEDDED */


  // Compute the objective if needed
  if (compute_objective) {
    *obj_val = compute_obj_val(work, x);
  }

  // Compute primal residual
  if (work->data->m == 0) {
    // No constraints -> Always primal feasible
    *pri_res = 0.;
  } else {
    *pri_res = compute_pri_res(work, x, z);
  }

  // Compute dual residual
  *dua_res = compute_dua_res(work, x, y);

  // Update timing
#ifdef PROFILING
  *run_time = osqp_toc(work->timer);
#endif /* ifdef PROFILING */

#ifdef PRINTING
  work->summary_printed = 0; // The just updated info have not been printed
#endif /* ifdef PRINTING */
}


void reset_info(OSQPInfo *info) {
#ifdef PROFILING

  // Initialize info values.
  info->solve_time = 0.0;  // Solve time to zero
# ifndef EMBEDDED
  info->polish_time = 0.0; // Polish time to zero
# endif /* ifndef EMBEDDED */

  // NB: We do not reset the setup_time because it is performed only once
#endif /* ifdef PROFILING */

  update_status(info, OSQP_UNSOLVED); // Problem is unsolved

#if EMBEDDED != 1
  info->rho_updates = 0;              // Rho updates are now 0
#endif /* if EMBEDDED != 1 */
}

void update_status(OSQPInfo *info, c_int status_val) {
  // Update status value
  info->status_val = status_val;

  // Update status string depending on status val
  if (status_val == OSQP_SOLVED) c_strcpy(info->status, "solved");

  if (status_val == OSQP_SOLVED_INACCURATE) c_strcpy(info->status,
                                                     "solved inaccurate");
  else if (status_val == OSQP_PRIMAL_INFEASIBLE) c_strcpy(info->status,
                                                          "primal infeasible");
  else if (status_val == OSQP_PRIMAL_INFEASIBLE_INACCURATE) c_strcpy(info->status,
                                                                     "primal infeasible inaccurate");
  else if (status_val == OSQP_UNSOLVED) c_strcpy(info->status, "unsolved");
  else if (status_val == OSQP_DUAL_INFEASIBLE) c_strcpy(info->status,
                                                        "dual infeasible");
  else if (status_val == OSQP_DUAL_INFEASIBLE_INACCURATE) c_strcpy(info->status,
                                                                   "dual infeasible inaccurate");
  else if (status_val == OSQP_MAX_ITER_REACHED) c_strcpy(info->status,
                                                         "maximum iterations reached");
#ifdef PROFILING
  else if (status_val == OSQP_TIME_LIMIT_REACHED) c_strcpy(info->status,
                                                           "run time limit reached");
#endif /* ifdef PROFILING */
  else if (status_val == OSQP_SIGINT) c_strcpy(info->status, "interrupted");

  else if (status_val == OSQP_NON_CVX) c_strcpy(info->status, "problem non convex");

}

c_int check_termination(OSQPWorkspace *work, c_int approximate) {
  c_float eps_prim, eps_dual, eps_prim_inf, eps_dual_inf;
  c_int   exitflag;
  c_int   prim_res_check, dual_res_check, prim_inf_check, dual_inf_check;
  c_float eps_abs, eps_rel;

  // Initialize variables to 0
  exitflag       = 0;
  prim_res_check = 0; dual_res_check = 0;
  prim_inf_check = 0; dual_inf_check = 0;

  // Initialize tolerances
  eps_abs      = work->settings->eps_abs;
  eps_rel      = work->settings->eps_rel;
  eps_prim_inf = work->settings->eps_prim_inf;
  eps_dual_inf = work->settings->eps_dual_inf;

  // If residuals are too large, the problem is probably non convex
  if ((work->info->pri_res > 2 * OSQP_INFTY) ||
      (work->info->dua_res > 2 * OSQP_INFTY)){
    // Looks like residuals are diverging. Probably the problem is non convex!
    // Terminate and report it
    update_status(work->info, OSQP_NON_CVX);
    work->info->obj_val = OSQP_NAN;
    return 1;
  }

  // If approximate solution required, increase tolerances by 10
  if (approximate) {
    eps_abs      *= 10;
    eps_rel      *= 10;
    eps_prim_inf *= 10;
    eps_dual_inf *= 10;
  }

  // Check residuals
  if (work->data->m == 0) {
    prim_res_check = 1; // No contraints -> Primal feasibility always satisfied
  }
  else {
    // Compute primal tolerance
    eps_prim = compute_pri_tol(work, eps_abs, eps_rel);

    // Primal feasibility check
    if (work->info->pri_res < eps_prim) {
      prim_res_check = 1;
    } else {
      // Primal infeasibility check
      prim_inf_check = is_primal_infeasible(work, eps_prim_inf);
    }
  } // End check if m == 0

  // Compute dual tolerance
  eps_dual = compute_dua_tol(work, eps_abs, eps_rel);

  // Dual feasibility check
  if (work->info->dua_res < eps_dual) {
    dual_res_check = 1;
  } else {
    // Check dual infeasibility
    dual_inf_check = is_dual_infeasible(work, eps_dual_inf);
  }

  // Compare checks to determine solver status
  if (prim_res_check && dual_res_check) {
    // Update final information
    if (approximate) {
      update_status(work->info, OSQP_SOLVED_INACCURATE);
    } else {
      update_status(work->info, OSQP_SOLVED);
    }
    exitflag = 1;
  }
  else if (prim_inf_check) {
    // Update final information
    if (approximate) {
      update_status(work->info, OSQP_PRIMAL_INFEASIBLE_INACCURATE);
    } else {
      update_status(work->info, OSQP_PRIMAL_INFEASIBLE);
    }

    if (work->settings->scaling && !work->settings->scaled_termination) {
      // Update infeasibility certificate
      vec_ew_prod(work->scaling->E, work->delta_y, work->delta_y, work->data->m);
    }
    work->info->obj_val = OSQP_INFTY;
    exitflag            = 1;
  }
  else if (dual_inf_check) {
    // Update final information
    if (approximate) {
      update_status(work->info, OSQP_DUAL_INFEASIBLE_INACCURATE);
    } else {
      update_status(work->info, OSQP_DUAL_INFEASIBLE);
    }

    if (work->settings->scaling && !work->settings->scaled_termination) {
      // Update infeasibility certificate
      vec_ew_prod(work->scaling->D, work->delta_x, work->delta_x, work->data->n);
    }
    work->info->obj_val = -OSQP_INFTY;
    exitflag            = 1;
  }

  return exitflag;
}

#ifndef EMBEDDED


c_int validate_data(const OSQPData *data) {
  c_int j;

  if (!data) {
# ifdef PRINTING
    c_eprint("Missing data");
# endif /* ifdef PRINTING */
    return 1;
  }

  // General dimensions Tests
  if ((data->n <= 0) || (data->m < 0)) {
# ifdef PRINTING
    c_eprint("n must be positive and m nonnegative; n = %i, m = %i",
             (int)data->n, (int)data->m);
# endif /* ifdef PRINTING */
    return 1;
  }

  // Matrix P
  if (data->P->m != data->n) {
# ifdef PRINTING
    c_eprint("P does not have dimension n x n with n = %i", (int)data->n);
# endif /* ifdef PRINTING */
    return 1;
  }

  if (data->P->m != data->P->n) {
# ifdef PRINTING
    c_eprint("P is not square");
# endif /* ifdef PRINTING */
    return 1;
  }

  // Matrix A
  if ((data->A->m != data->m) || (data->A->n != data->n)) {
# ifdef PRINTING
    c_eprint("A does not have dimension m x n with m = %i and n = %i",
             (int)data->m, (int)data->n);
# endif /* ifdef PRINTING */
    return 1;
  }

  // Lower and upper bounds
  for (j = 0; j < data->m; j++) {
    if (data->l[j] > data->u[j]) {
# ifdef PRINTING
      c_eprint("Lower bound at index %d is greater than upper bound: %.4e > %.4e",
               (int)j, data->l[j], data->u[j]);
# endif /* ifdef PRINTING */
      return 1;
    }
  }

  // TODO: Complete with other checks

  return 0;
}

c_int validate_linsys_solver(c_int linsys_solver) {
  if ((linsys_solver != QDLDL_SOLVER) &&
      (linsys_solver != MKL_PARDISO_SOLVER)) {
    return 1;
  }

  // TODO: Add more solvers in case

  // Valid solver
  return 0;
}

c_int validate_settings(const OSQPSettings *settings) {
  if (!settings) {
# ifdef PRINTING
    c_eprint("Missing settings!");
# endif /* ifdef PRINTING */
    return 1;
  }

  if (settings->scaling < 0) {
# ifdef PRINTING
    c_eprint("scaling must be nonnegative");
# endif /* ifdef PRINTING */
    return 1;
  }

  if ((settings->adaptive_rho != 0) && (settings->adaptive_rho != 1)) {
# ifdef PRINTING
    c_eprint("adaptive_rho must be either 0 or 1");
# endif /* ifdef PRINTING */
    return 1;
  }

  if (settings->adaptive_rho_interval < 0) {
# ifdef PRINTING
    c_eprint("adaptive_rho_interval must be nonnegative");
# endif /* ifdef PRINTING */
    return 1;
  }
# ifdef PROFILING

  if (settings->adaptive_rho_fraction <= 0) {
#  ifdef PRINTING
    c_eprint("adaptive_rho_fraction must be positive");
#  endif /* ifdef PRINTING */
    return 1;
  }
# endif /* ifdef PROFILING */

  if (settings->adaptive_rho_tolerance < 1) {
# ifdef PRINTING
    c_eprint("adaptive_rho_tolerance must be >= 1");
# endif /* ifdef PRINTING */
    return 1;
  }

  if (settings->polish_refine_iter < 0) {
# ifdef PRINTING
    c_eprint("polish_refine_iter must be nonnegative");
# endif /* ifdef PRINTING */
    return 1;
  }

  if (settings->rho <= 0) {
# ifdef PRINTING
    c_eprint("rho must be positive");
# endif /* ifdef PRINTING */
    return 1;
  }

  if (settings->delta <= 0) {
# ifdef PRINTING
    c_eprint("delta must be positive");
# endif /* ifdef PRINTING */
    return 1;
  }

  if (settings->max_iter <= 0) {
# ifdef PRINTING
    c_eprint("max_iter must be positive");
# endif /* ifdef PRINTING */
    return 1;
  }

  if (settings->eps_abs < 0) {
# ifdef PRINTING
    c_eprint("eps_abs must be nonnegative");
# endif /* ifdef PRINTING */
    return 1;
  }

  if (settings->eps_rel < 0) {
# ifdef PRINTING
    c_eprint("eps_rel must be nonnegative");
# endif /* ifdef PRINTING */
    return 1;
  }

  if ((settings->eps_rel == 0) && (settings->eps_abs == 0)) {
# ifdef PRINTING
    c_eprint("at least one of eps_abs and eps_rel must be positive");
# endif /* ifdef PRINTING */
    return 1;
  }

  if (settings->eps_prim_inf < 0) {
# ifdef PRINTING
    c_eprint("eps_prim_inf must be nonnegative");
# endif /* ifdef PRINTING */
    return 1;
  }

  if (settings->eps_dual_inf < 0) {
# ifdef PRINTING
    c_eprint("eps_dual_inf must be nonnegative");
# endif /* ifdef PRINTING */
    return 1;
  }

  if ((settings->alpha <= 0) || (settings->alpha >= 2)) {
# ifdef PRINTING
    c_eprint("alpha must be between 0 and 2");
# endif /* ifdef PRINTING */
    return 1;
  }

  if (validate_linsys_solver(settings->linsys_solver)) {
# ifdef PRINTING
    c_eprint("linsys_solver not recognized");
# endif /* ifdef PRINTING */
    return 1;
  }

  if ((settings->verbose != 0) && (settings->verbose != 1)) {
# ifdef PRINTING
    c_eprint("verbose must be either 0 or 1");
# endif /* ifdef PRINTING */
    return 1;
  }

  if ((settings->scaled_termination != 0) &&
      (settings->scaled_termination != 1)) {
# ifdef PRINTING
    c_eprint("scaled_termination must be either 0 or 1");
# endif /* ifdef PRINTING */
    return 1;
  }

  if (settings->check_termination < 0) {
# ifdef PRINTING
    c_eprint("check_termination must be nonnegative");
# endif /* ifdef PRINTING */
    return 1;
  }

  if ((settings->warm_start != 0) && (settings->warm_start != 1)) {
# ifdef PRINTING
    c_eprint("warm_start must be either 0 or 1");
# endif /* ifdef PRINTING */
    return 1;
  }
# ifdef PROFILING

  if (settings->time_limit < 0) {
#  ifdef PRINTING
    c_eprint("time_limit must be nonnegative\n");
#  endif /* ifdef PRINTING */
    return 1;
  }
# endif /* ifdef PROFILING */

  return 0;
}

#endif // #ifndef EMBEDDED
