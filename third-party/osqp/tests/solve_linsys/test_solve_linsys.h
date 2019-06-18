#include <stdio.h>
#include "osqp.h"
#include "cs.h"
#include "util.h"
#include "minunit.h"
#include "lin_sys.h"


#include "solve_linsys/data.h"


static char* test_solveKKT() {
  c_int m, exitflag = 0;
  c_float *rho_vec;
  LinSysSolver *p;  // Private structure to form KKT factorization
  OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings)); // Settings
  solve_linsys_sols_data *data = generate_problem_solve_linsys_sols_data();

  // Settings
  settings->rho   = data->test_solve_KKT_rho;
  settings->sigma = data->test_solve_KKT_sigma;

  // Set rho_vec
  m       = data->test_solve_KKT_A->m;
  rho_vec = c_calloc(m, sizeof(c_float));
  vec_add_scalar(rho_vec, settings->rho, m);

  // Form and factorize KKT matrix
  p = init_linsys_solver(data->test_solve_KKT_Pu, data->test_solve_KKT_A,
                         settings->sigma, rho_vec, LINSYS_SOLVER, 0);

  // // Debug print KKT and LDL
  // print_csc_mat(data->test_solve_KKT_KKT, "KKTpy");
  // // print_csc_mat(p->KKT, "KKT");
  // c_float * KKTdnspy = csc_to_dns(data->test_solve_KKT_KKT);
  // // c_float * KKTdns = csc_to_dns(p->KKT);
  // print_dns_matrix(KKTdnspy, data->test_solve_KKT_KKT->m, data->test_solve_KKT_KKT->n, "KKTdnspy");
  // // print_dns_matrix(KKTdns, p->KKT->m, p->KKT->n, "KKTdns");


  // Solve  KKT x = b via LDL given factorization
  p->solve(p, data->test_solve_KKT_rhs, settings);

  mu_assert(
    "Linear systems solve tests: error in forming and solving KKT system!",
    vec_norm_inf_diff(data->test_solve_KKT_rhs, data->test_solve_KKT_x,
                      data->test_solve_KKT_m + data->test_solve_KKT_n) <
    TESTS_TOL);


  // Cleanup
  p->free(p);
  c_free(settings);
  c_free(rho_vec);
  clean_problem_solve_linsys_sols_data(data);

  return 0;
}

#ifdef ENABLE_MKL_PARDISO
static char* test_solveKKT_pardiso() {
  c_int m, exitflag = 0;
  c_float *rho_vec;
  LinSysSolver *p;                                                         // Private
                                                                           // structure
                                                                           // to
                                                                           // form
                                                                           // KKT
                                                                           // factorization
  OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings)); // Settings

  solve_linsys_sols_data *data = generate_problem_solve_linsys_sols_data();

  // Settings
  settings->rho   = data->test_solve_KKT_rho;
  settings->sigma = data->test_solve_KKT_sigma;

  // Set rho_vec
  m       = data->test_solve_KKT_A->m;
  rho_vec = c_calloc(m, sizeof(c_float));
  vec_add_scalar(rho_vec, settings->rho, m);

  // Load Pardiso shared library
  exitflag = load_linsys_solver(MKL_PARDISO_SOLVER);
  mu_assert("Linear system solve test: error in loading Pardiso shared library",
            exitflag == 0);

  // Form and factorize KKT matrix
  p = init_linsys_solver(data->test_solve_KKT_Pu, data->test_solve_KKT_A,
                         settings->sigma, rho_vec, MKL_PARDISO_SOLVER, 0);

  // Solve  KKT x = b via LDL given factorization
  p->solve(p, data->test_solve_KKT_rhs, settings);

  mu_assert(
    "Linear systems solve tests: error in forming and solving KKT system with PARDISO!",
    vec_norm_inf_diff(data->test_solve_KKT_rhs, data->test_solve_KKT_x,
                      data->test_solve_KKT_m + data->test_solve_KKT_n) <
    TESTS_TOL);


  // Cleanup
  p->free(p);
  c_free(settings);
  c_free(rho_vec);
  clean_problem_solve_linsys_sols_data(data);

  // Unload Pardiso shared library
  exitflag = unload_linsys_solver(MKL_PARDISO_SOLVER);

  return 0;
}
#endif

static char* test_solve_linsys()
{
  mu_run_test(test_solveKKT);
#ifdef ENABLE_MKL_PARDISO
  mu_run_test(test_solveKKT_pardiso);
#endif

  return 0;
}
