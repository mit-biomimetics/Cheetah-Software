#include "osqp.h"    // OSQP API
#include "minunit.h" // Basic testing script header


#include "non_cvx/data.h"


static char* test_non_cvx_solve()
{
  /* local variables */
  c_int exitflag = 0; // No errors

  // Problem settings
  OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

  // Structures
  OSQPWorkspace *work; // Workspace
  OSQPData *data;      // Data
  non_cvx_sols_data *sols_data;


  // Populate data
  data      = generate_problem_non_cvx();
  sols_data = generate_problem_non_cvx_sols_data();


  // Define Solver settings as default
  osqp_set_default_settings(settings);
  settings->verbose = 1;
  settings->sigma = 1e-6;

  // Setup workspace
  work = osqp_setup(data, settings);

  // Setup should fail due to (P + sigma I) having a negative eigenvalue
  mu_assert("Non Convex test solve: Setup should have failed!", work == OSQP_NULL);

  // Update Solver settings
  settings->sigma = sols_data->sigma_new;

  // Setup workspace again
  work = osqp_setup(data, settings);

  // Setup should work this time because (P + sigma I) is positive definite
  mu_assert("Non Convex test solve: Setup error!", work != OSQP_NULL);

  // Solve Problem first time
  osqp_solve(work);

  // Compare solver statuses
  mu_assert("Non Convex test solve: Error in solver status!",
            work->info->status_val == OSQP_NON_CVX);

  // Compare objective values
  mu_assert("Non Convex test solve: Error in objective value!",
            work->info->obj_val == OSQP_NAN);

  // Clean workspace
  osqp_cleanup(work);

  // Cleanup settings and data
  c_free(settings);
  clean_problem_non_cvx(data);
  clean_problem_non_cvx_sols_data(sols_data);

  return 0;
}

static char* test_non_cvx()
{
  mu_run_test(test_non_cvx_solve);

  return 0;
}
