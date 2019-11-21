#include "osqp.h"
#include "cs.h"
#include "util.h"
#include "minunit.h"

#include "primal_dual_infeasibility/data.h"


static char* test_optimal()
{
  // Structures
  OSQPWorkspace *work;    // Workspace
  OSQPData *problem;      // Problem data
  OSQPSettings *settings; // Settings
  primal_dual_infeasibility_sols_data *data;

  // Load problem data
  data = generate_problem_primal_dual_infeasibility_sols_data();

  // Populate problem data
  problem    = c_malloc(sizeof(OSQPData));
  problem->P = data->P;
  problem->q = data->q;
  problem->A = data->A12;
  problem->l = data->l;
  problem->u = data->u1;
  problem->n = data->P->n;
  problem->m = data->A12->m;

  // Define Solver settings as default
  settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
  osqp_set_default_settings(settings);
  settings->max_iter = 2000;
  settings->alpha    = 1.6;
  settings->polish   = 1;
  settings->scaling  = 0;
  settings->verbose  = 1;

  // Setup workspace
  work = osqp_setup(problem, settings);

  // Setup correct
  mu_assert("Primal dual infeasibility test 1: Setup error!", work != OSQP_NULL);

  // Solve Problem
  osqp_solve(work);

  // Compare solver statuses
  mu_assert("Primal dual infeasibility test 1: Error in solver status!",
            work->info->status_val == OSQP_SOLVED);

  // Compare primal solutions
  mu_assert("Primal dual infeasibility test 1: Error in primal solution!",
            vec_norm_inf_diff(work->solution->x, data->x1,
                              problem->n) < TESTS_TOL);

  // Compare dual solutions
  mu_assert("Primal dual infeasibility test 1: Error in dual solution!",
            vec_norm_inf_diff(work->solution->y, data->y1,
                              problem->m) < TESTS_TOL);


  // Compare objective values
  mu_assert("Primal dual infeasibility test 1: Error in objective value!",
            c_absval(work->info->obj_val - data->obj_value1) < TESTS_TOL);


  // Cleanup
  osqp_cleanup(work);
  clean_problem_primal_dual_infeasibility_sols_data(data);
  c_free(problem);
  c_free(settings);

  return 0;
}

static char* test_prim_infeas()
{
  // Structures
  OSQPWorkspace *work;    // Workspace
  OSQPData *problem;      // Problem data
  OSQPSettings *settings; // Settings
  primal_dual_infeasibility_sols_data *data;

  // Load problem data
  data = generate_problem_primal_dual_infeasibility_sols_data();

  // Populate problem data
  problem    = c_malloc(sizeof(OSQPData));
  problem->P = data->P;
  problem->q = data->q;
  problem->A = data->A12;
  problem->l = data->l;
  problem->u = data->u2;
  problem->n = data->P->n;
  problem->m = data->A12->m;

  // Define Solver settings as default
  settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
  osqp_set_default_settings(settings);
  settings->max_iter = 2000;
  settings->alpha    = 1.6;
  settings->polish   = 0;
  settings->scaling  = 0;
  settings->verbose  = 1;

  // Setup workspace
  work = osqp_setup(problem, settings);

  // Setup correct
  mu_assert("Primal dual infeasibility test 2: Setup error!", work != OSQP_NULL);

  // Solve Problem
  osqp_solve(work);

  // Compare solver statuses
  mu_assert("Primal dual infeasibility test 2: Error in solver status!",
            work->info->status_val == OSQP_PRIMAL_INFEASIBLE);

  // Cleanup
  osqp_cleanup(work);
  clean_problem_primal_dual_infeasibility_sols_data(data);
  c_free(problem);
  c_free(settings);

  return 0;
}

static char* test_dual_infeas()
{
  // Structures
  OSQPWorkspace *work;    // Workspace
  OSQPData *problem;      // Problem data
  OSQPSettings *settings; // Settings
  primal_dual_infeasibility_sols_data *data;

  // Load problem data
  data = generate_problem_primal_dual_infeasibility_sols_data();

  // Populate problem data
  problem    = c_malloc(sizeof(OSQPData));
  problem->P = data->P;
  problem->q = data->q;
  problem->A = data->A34;
  problem->l = data->l;
  problem->u = data->u3;
  problem->n = data->P->n;
  problem->m = data->A34->m;

  // Define Solver settings as default
  settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
  osqp_set_default_settings(settings);
  settings->max_iter = 2000;
  settings->alpha    = 1.6;
  settings->polish   = 0;
  settings->scaling  = 0;
  settings->verbose  = 1;

  // Setup workspace
  work = osqp_setup(problem, settings);

  // Setup correct
  mu_assert("Primal dual infeasibility test 3: Setup error!", work != OSQP_NULL);

  // Solve Problem
  osqp_solve(work);

  // Compare solver statuses
  mu_assert("Primal dual infeasibility test 3: Error in solver status!",
            work->info->status_val == OSQP_DUAL_INFEASIBLE);

  // Cleanup
  osqp_cleanup(work);
  clean_problem_primal_dual_infeasibility_sols_data(data);
  c_free(problem);
  c_free(settings);

  return 0;
}

static char* test_primal_dual_infeas()
{
  // Structures
  OSQPWorkspace *work;    // Workspace
  OSQPData *problem;      // Problem data
  OSQPSettings *settings; // Settings
  primal_dual_infeasibility_sols_data *data;

  // Load problem data
  data = generate_problem_primal_dual_infeasibility_sols_data();

  // Populate problem data
  problem    = c_malloc(sizeof(OSQPData));
  problem->P = data->P;
  problem->q = data->q;
  problem->A = data->A34;
  problem->l = data->l;
  problem->u = data->u4;
  problem->n = data->P->n;
  problem->m = data->A34->m;

  // Define Solver settings as default
  settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
  osqp_set_default_settings(settings);
  settings->max_iter = 2000;
  settings->alpha    = 1.6;
  settings->polish   = 0;
  settings->scaling  = 0;
  settings->verbose  = 1;

  // Setup workspace
  work = osqp_setup(problem, settings);

  // Setup correct
  mu_assert("Primal dual infeasibility test 4: Setup error!", work != OSQP_NULL);

  // Solve Problem
  osqp_solve(work);

  // Compare solver statuses
  mu_assert("Primal dual infeasibility test 4: Error in solver status!",
            ((work->info->status_val == OSQP_PRIMAL_INFEASIBLE) ||
             (work->info->status_val == OSQP_DUAL_INFEASIBLE)));

  // Cleanup
  osqp_cleanup(work);
  clean_problem_primal_dual_infeasibility_sols_data(data);
  c_free(problem);
  c_free(settings);

  return 0;
}

static char* test_primal_dual_infeasibility()
{
  mu_run_test(test_optimal);
  mu_run_test(test_prim_infeas);
  mu_run_test(test_dual_infeas);
  mu_run_test(test_primal_dual_infeas);

  return 0;
}
