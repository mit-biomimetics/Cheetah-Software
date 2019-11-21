#include "osqp.h"    // OSQP API
#include "minunit.h" // Basic testing script header


#include "unconstrained/data.h"


static char* test_unconstrained_solve()
{
  /* local variables */
  c_int exitflag = 0; // No errors

  // Problem settings
  OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

  // Structures
  OSQPWorkspace *work; // Workspace
  OSQPData *data;      // Data
  unconstrained_sols_data *sols_data;


  // Populate data
  data      = generate_problem_unconstrained();
  sols_data = generate_problem_unconstrained_sols_data();


  // Define Solver settings as default
  osqp_set_default_settings(settings);
  settings->verbose = 1;

  // Setup workspace
  work = osqp_setup(data, settings);

  // Setup correct
  mu_assert("Unconstrained test solve: Setup error!", work != OSQP_NULL);

  // Solve Problem first time
  osqp_solve(work);

  // Compare solver statuses
  mu_assert("Unconstrained test solve: Error in solver status!",
            work->info->status_val == sols_data->status_test);

  // Compare primal solutions
  mu_assert("Unconstrained test solve: Error in primal solution!",
            vec_norm_inf_diff(work->solution->x, sols_data->x_test,
                              data->n) < TESTS_TOL);

  // Compare objective values
  mu_assert("Unconstrained test solve: Error in objective value!",
            c_absval(work->info->obj_val - sols_data->obj_value_test) <
            TESTS_TOL);

  // Clean workspace
  osqp_cleanup(work);

  // Cleanup settings and data
  c_free(settings);
  clean_problem_unconstrained(data);
  clean_problem_unconstrained_sols_data(sols_data);

  return 0;
}

static char* test_unconstrained()
{
  mu_run_test(test_unconstrained_solve);

  return 0;
}
