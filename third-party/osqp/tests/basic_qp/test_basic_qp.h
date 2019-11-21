#include "osqp.h"    // OSQP API
#include "cs.h"      // CSC data structure
#include "util.h"    // Utilities for testing
#include "minunit.h" // Basic testing script header

#include "basic_qp/data.h"


static char* test_basic_qp_solve()
{
  // Problem settings
  OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

  // Structures
  OSQPWorkspace *work; // Workspace
  OSQPData *data;      // Data
  basic_qp_sols_data *sols_data;

  // Populate data
  data      = generate_problem_basic_qp();
  sols_data = generate_problem_basic_qp_sols_data();


  // Define Solver settings as default
  osqp_set_default_settings(settings);
  settings->max_iter   = 2000;
  settings->alpha      = 1.6;
  settings->polish     = 1;
  settings->scaling    = 0;
  settings->verbose    = 1;
  settings->warm_start = 0;

  // Setup workspace
  work = osqp_setup(data, settings);

  // Setup correct
  mu_assert("Basic QP test solve: Setup error!", work != OSQP_NULL);

  // Solve Problem
  osqp_solve(work);

  // Compare solver statuses
  mu_assert("Basic QP test solve: Error in solver status!",
            work->info->status_val == sols_data->status_test);

  // Compare primal solutions
  mu_assert("Basic QP test solve: Error in primal solution!",
            vec_norm_inf_diff(work->solution->x, sols_data->x_test,
                              data->n) < TESTS_TOL);

  // Compare dual solutions
  mu_assert("Basic QP test solve: Error in dual solution!",
            vec_norm_inf_diff(work->solution->y, sols_data->y_test,
                              data->m) < TESTS_TOL);


  // Compare objective values
  mu_assert("Basic QP test solve: Error in objective value!",
            c_absval(work->info->obj_val - sols_data->obj_value_test) <
            TESTS_TOL);

  // Try to set wrong settings
  mu_assert("Basic QP test solve: Wrong value of rho not caught!",
            osqp_update_rho(work, -0.1) == 1);

  mu_assert("Basic QP test solve: Wrong value of max_iter not caught!",
            osqp_update_max_iter(work, -1) == 1);

  mu_assert("Basic QP test solve: Wrong value of eps_abs not caught!",
            osqp_update_eps_abs(work, -1.) == 1);

  mu_assert("Basic QP test solve: Wrong value of eps_rel not caught!",
            osqp_update_eps_rel(work, -1.) == 1);

  mu_assert("Basic QP test solve: Wrong value of eps_prim_inf not caught!",
            osqp_update_eps_prim_inf(work, -0.1) == 1);

  mu_assert("Basic QP test solve: Wrong value of eps_dual_inf not caught!",
            osqp_update_eps_dual_inf(work, -0.1) == 1);

  mu_assert("Basic QP test solve: Wrong value of alpha not caught!",
            osqp_update_alpha(work, 2.0) == 1);

  mu_assert("Basic QP test solve: Wrong value of warm_start not caught!",
            osqp_update_warm_start(work, -1) == 1);

  mu_assert("Basic QP test solve: Wrong value of scaled_termination not caught!",
            osqp_update_scaled_termination(work, 2) == 1);

  mu_assert("Basic QP test solve: Wrong value of check_termination not caught!",
            osqp_update_check_termination(work, -1) == 1);

  mu_assert("Basic QP test solve: Wrong value of delta not caught!",
            osqp_update_delta(work, 0.) == 1);

  mu_assert("Basic QP test solve: Wrong value of polish not caught!",
            osqp_update_polish(work, 2) == 1);

  mu_assert("Basic QP test solve: Wrong value of polish_refine_iter not caught!",
            osqp_update_polish_refine_iter(work, -1) == 1);

  mu_assert("Basic QP test solve: Wrong value of verbose not caught!",
            osqp_update_verbose(work, 2) == 1);


  // Clean workspace
  osqp_cleanup(work);


  // Cleanup data
  clean_problem_basic_qp(data);
  clean_problem_basic_qp_sols_data(sols_data);

  // Cleanup
  c_free(settings);

  return 0;
}

#ifdef ENABLE_MKL_PARDISO
static char* test_basic_qp_solve_pardiso()
{
  // Problem settings
  OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

  // Structures
  OSQPWorkspace *work; // Workspace
  OSQPData *data;      // Data
  basic_qp_sols_data *sols_data;

  // Populate data
  data      = generate_problem_basic_qp();
  sols_data = generate_problem_basic_qp_sols_data();


  // Define Solver settings as default
  osqp_set_default_settings(settings);
  settings->max_iter      = 2000;
  settings->alpha         = 1.6;
  settings->polish        = 1;
  settings->scaling       = 0;
  settings->verbose       = 1;
  settings->warm_start    = 0;
  settings->linsys_solver = MKL_PARDISO_SOLVER;

  // Setup workspace
  work = osqp_setup(data, settings);

  // Setup correct
  mu_assert("Basic QP test solve: Setup error!", work != OSQP_NULL);

  // Solve Problem
  osqp_solve(work);

  // Compare solver statuses
  mu_assert("Basic QP test solve: Error in solver status!",
            work->info->status_val == sols_data->status_test);

  // Compare primal solutions
  mu_assert("Basic QP test solve: Error in primal solution!",
            vec_norm_inf_diff(work->solution->x, sols_data->x_test,
                              data->n) < TESTS_TOL);

  // Compare dual solutions
  mu_assert("Basic QP test solve: Error in dual solution!",
            vec_norm_inf_diff(work->solution->y, sols_data->y_test,
                              data->m) < TESTS_TOL);


  // Compare objective values
  mu_assert("Basic QP test solve: Error in objective value!",
            c_absval(work->info->obj_val - sols_data->obj_value_test) <
            TESTS_TOL);

  // Clean workspace
  osqp_cleanup(work);


  // Cleanup data
  clean_problem_basic_qp(data);
  clean_problem_basic_qp_sols_data(sols_data);

  // Cleanup
  c_free(settings);

  return 0;
}
#endif

static char* test_basic_qp_update()
{
  // Problem settings
  OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

  // Structures
  OSQPWorkspace *work; // Workspace
  OSQPData *data;      // Data
  basic_qp_sols_data *sols_data;

  // Populate data
  data      = generate_problem_basic_qp();
  sols_data = generate_problem_basic_qp_sols_data();


  // Define Solver settings as default
  osqp_set_default_settings(settings);
  settings->max_iter   = 200;
  settings->alpha      = 1.6;
  settings->polish     = 1;
  settings->scaling    = 0;
  settings->verbose    = 1;
  settings->warm_start = 0;

  // Setup workspace
  work = osqp_setup(data, settings);

  // Setup correct
  mu_assert("Basic QP test update: Setup error!", work != OSQP_NULL);


  // ====================================================================
  //  Update data
  // ====================================================================

  // Update linear cost
  osqp_update_lin_cost(work, sols_data->q_new);
  mu_assert("Basic QP test update: Error in updating linear cost!",
            vec_norm_inf_diff(work->data->q, sols_data->q_new,
                              data->n) < TESTS_TOL);

  // UPDATE BOUND
  // Try to update with non-consistent values
  mu_assert("Basic QP test update: Error in bounds update ordering not caught!",
            osqp_update_bounds(work, sols_data->u_new, sols_data->l_new) == 1);

  // Now update with correct values
  mu_assert("Basic QP test update: Error in bounds update ordering!",
            osqp_update_bounds(work, sols_data->l_new, sols_data->u_new) == 0);

  mu_assert("Basic QP test update: Error in bounds update, lower bound!",
            vec_norm_inf_diff(work->data->l, sols_data->l_new,
                              data->m) < TESTS_TOL);

  mu_assert("Basic QP test update: Error in bounds update, upper bound!",
            vec_norm_inf_diff(work->data->u, sols_data->u_new,
                              data->m) < TESTS_TOL);

  // Return original values
  osqp_update_bounds(work, data->l, data->u);


  // UPDATE LOWER BOUND
  // Try to update with non-consistent values
  mu_assert(
    "Basic QP test update: Error in lower bound update ordering not caught!",
    osqp_update_lower_bound(work, sols_data->u_new) == 1);

  // Now update with correct values
  mu_assert("Basic QP test update: Error in lower bound update ordering!",
            osqp_update_lower_bound(work, sols_data->l_new) == 0);

  mu_assert("Basic QP test update: Error in updating lower bound!",
            vec_norm_inf_diff(work->data->l, sols_data->l_new,
                              data->m) < TESTS_TOL);

  // Return original values
  osqp_update_lower_bound(work, data->l);


  // UPDATE UPPER BOUND
  // Try to update with non-consistent values
  mu_assert(
    "Basic QP test update: Error in upper bound update: ordering not caught!",
    osqp_update_upper_bound(work, sols_data->l_new) == 1);

  // Now update with correct values
  mu_assert("Basic QP test update: Error in upper bound update: ordering!",
            osqp_update_upper_bound(work, sols_data->u_new) == 0);

  mu_assert("Basic QP test update: Error in updating upper bound!",
            vec_norm_inf_diff(work->data->u, sols_data->u_new,
                              data->m) < TESTS_TOL);


  // Clean workspace
  osqp_cleanup(work);


  // Cleanup data
  clean_problem_basic_qp(data);
  clean_problem_basic_qp_sols_data(sols_data);

  // Cleanup
  c_free(settings);

  return 0;
}

static char* test_basic_qp_check_termination()
{
  // Problem settings
  OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

  // Structures
  OSQPWorkspace *work; // Workspace
  OSQPData *data;      // Data
  basic_qp_sols_data *sols_data;

  // Populate data
  data      = generate_problem_basic_qp();
  sols_data = generate_problem_basic_qp_sols_data();


  // Define Solver settings as default
  osqp_set_default_settings(settings);
  settings->max_iter          = 200;
  settings->alpha             = 1.6;
  settings->polish            = 0;
  settings->scaling           = 0;
  settings->verbose           = 1;
  settings->check_termination = 0;
  settings->warm_start        = 0;

  // Setup workspace
  work = osqp_setup(data, settings);

  // Setup correct
  mu_assert("Basic QP test solve: Setup error!", work != OSQP_NULL);

  // Solve Problem
  osqp_solve(work);

  // Check if iter == max_iter
  mu_assert(
    "Basic QP test check termination: Error in number of iterations taken!",
    work->info->iter == work->settings->max_iter);

  // Compare solver statuses
  mu_assert("Basic QP test check termination: Error in solver status!",
            work->info->status_val == sols_data->status_test);

  // Compare primal solutions
  mu_assert("Basic QP test check termination: Error in primal solution!",
            vec_norm_inf_diff(work->solution->x, sols_data->x_test,
                              data->n) < TESTS_TOL);

  // Compare dual solutions
  // print_vec(work->solution->y, data->m, "y_sol");
  // print_vec(sols_data->y_test, data->m, "y_test");
  mu_assert("Basic QP test check termination: Error in dual solution!",
            vec_norm_inf_diff(work->solution->y, sols_data->y_test,
                              data->m) < TESTS_TOL);

  // Compare objective values
  mu_assert("Basic QP test check termination: Error in objective value!",
            c_absval(work->info->obj_val - sols_data->obj_value_test) <
            TESTS_TOL);

  // Clean workspace
  osqp_cleanup(work);

  // Cleanup data
  clean_problem_basic_qp(data);
  clean_problem_basic_qp_sols_data(sols_data);

  // Cleanup
  c_free(settings);

  return 0;
}

static char* test_basic_qp_update_rho()
{
  // Problem settings
  OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

  // Structures
  OSQPWorkspace *work; // Workspace
  OSQPData *data;      // Data
  basic_qp_sols_data *sols_data;

  // Exitflag
  c_int exitflag;

  // rho to use
  c_float rho;

  // Define number of iterations to compare
  c_int n_iter_new_solver, n_iter_update_rho;

  // Populate data
  data      = generate_problem_basic_qp();
  sols_data = generate_problem_basic_qp_sols_data();


  // Define Solver settings as default
  rho = 0.7;
  osqp_set_default_settings(settings);
  settings->rho               = rho;
  settings->adaptive_rho      = 0; // Disable adaptive rho for this test
  settings->eps_abs           = 1e-04;
  settings->eps_rel           = 1e-04;
  settings->check_termination = 1;

  // Setup workspace
  work = osqp_setup(data, settings);

  // Setup correct
  mu_assert("Update rho test solve: Setup error!", work != OSQP_NULL);

  // Solve Problem
  osqp_solve(work);

  // Store number of iterations
  n_iter_new_solver = work->info->iter;

  // Compare solver statuses
  mu_assert("Update rho test solve: Error in solver status!",
            work->info->status_val == sols_data->status_test);

  // Compare primal solutions
  mu_assert("Update rho test solve: Error in primal solution!",
            vec_norm_inf_diff(work->solution->x, sols_data->x_test,
                              data->n)/vec_norm_inf(sols_data->x_test, data->n) < TESTS_TOL);

  // Compare dual solutions
  mu_assert("Update rho test solve: Error in dual solution!",
            vec_norm_inf_diff(work->solution->y, sols_data->y_test,
                              data->m)/vec_norm_inf(sols_data->y_test, data->m) < TESTS_TOL);

  // Compare objective values
  mu_assert("Update rho test solve: Error in objective value!",
            c_absval(work->info->obj_val - sols_data->obj_value_test) <
            TESTS_TOL);

  // Clean workspace
  osqp_cleanup(work);


  // Create new problem with different rho and update it
  osqp_set_default_settings(settings);
  settings->rho               = 0.1;
  settings->adaptive_rho      = 0;
  settings->check_termination = 1;
  settings->eps_abs           = 1e-04;
  settings->eps_rel           = 1e-04;

  // Setup workspace
  work = osqp_setup(data, settings);

  // Setup correct
  mu_assert("Update rho test update: Setup error!", work != OSQP_NULL);

  // Update rho
  exitflag = osqp_update_rho(work, rho);
  mu_assert("Update rho test update: Error update rho!", exitflag == 0);

  // Solve Problem
  osqp_solve(work);

  // Compare solver statuses
  mu_assert("Update rho test update: Error in solver status!",
            work->info->status_val == sols_data->status_test);

  // Compare primal solutions
  mu_assert("Update rho test update: Error in primal solution!",
            vec_norm_inf_diff(work->solution->x, sols_data->x_test,
                              data->n)/vec_norm_inf(sols_data->x_test, data->n) < TESTS_TOL);

  // Compare dual solutions
  mu_assert("Update rho test update: Error in dual solution!",
            vec_norm_inf_diff(work->solution->y, sols_data->y_test,
                              data->m)/vec_norm_inf(sols_data->y_test, data->m)< TESTS_TOL);

  // Compare objective values
  mu_assert("Update rho test update: Error in objective value!",
            c_absval(work->info->obj_val - sols_data->obj_value_test) <
            TESTS_TOL);

  // Get number of iterations
  n_iter_update_rho = work->info->iter;

  // Assert same number of iterations
  mu_assert("Update rho test update: Error in number of iterations!",
            n_iter_new_solver == n_iter_update_rho);

  // Cleanup solver
  osqp_cleanup(work);

  // Cleanup data
  clean_problem_basic_qp(data);
  clean_problem_basic_qp_sols_data(sols_data);

  // Cleanup
  c_free(settings);

  return 0;
}

static char* test_basic_qp_time_limit()
{
  // Problem settings
  OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

  // Structures
  OSQPWorkspace *work; // Workspace
  OSQPData *data;      // Data
  basic_qp_sols_data *sols_data;

  // Exitflag
  c_int exitflag;

  // Populate data
  data      = generate_problem_basic_qp();
  sols_data = generate_problem_basic_qp_sols_data();

  // Define Solver settings as default
  osqp_set_default_settings(settings);

  // Check dfault time limit
  mu_assert("Time limit test: Default not correct", settings->time_limit == 0);

  // Setup workspace
  work = osqp_setup(data, settings);

  // Setup correct
  mu_assert("Time limit test: Setup error!", work != OSQP_NULL);

  // Solve Problem
  osqp_solve(work);

  // Compare solver statuses
  mu_assert("Time limit test: Error in no time limit solver status!",
            work->info->status_val == sols_data->status_test);

  // Update time limit
  osqp_update_time_limit(work, 1e-5);
  osqp_update_max_iter(work, 2000000000);
  osqp_update_check_termination(work, 0);

  // Solve Problem
  osqp_solve(work);

  // Compare solver statuses
  mu_assert("Time limit test: Error in timed out solver status!",
            work->info->status_val == OSQP_TIME_LIMIT_REACHED);

  // Cleanup solver
  osqp_cleanup(work);

  // Cleanup data
  clean_problem_basic_qp(data);
  clean_problem_basic_qp_sols_data(sols_data);

  // Cleanup
  c_free(settings);

  return 0;
}

static char* test_basic_qp()
{
  mu_run_test(test_basic_qp_solve);
#ifdef ENABLE_MKL_PARDISO
  mu_run_test(test_basic_qp_solve_pardiso);
#endif
  mu_run_test(test_basic_qp_update);
  mu_run_test(test_basic_qp_check_termination);
  mu_run_test(test_basic_qp_update_rho);
  mu_run_test(test_basic_qp_time_limit);

  return 0;
}
