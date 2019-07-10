#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Utilities/utilities.h"
#include "osqp.h"

TEST(OSQP, osqp_c_interface) {
  // problem:
  // >> H = [4 1; 1 2]; f = [1;1]; A = [1 1; 1 0; 0 1; -1 -1; -1 0; 0 -1]; b =
  // [1 .7 .7 -1 0 0]';
  // >> x = quadprog(H,f,A,b)
  // x =
  //
  //    0.3000
  //    0.7000

  // hessian
  c_float P_x[4] = {
      4.00,
      1.00,
      1.00,
      2.00,
  };
  c_int P_nnz = 4;
  c_int P_i[4] = {
      0,
      1,
      0,
      1,
  };
  c_int P_p[3] = {
      0,
      2,
      4,
  };

  // gradient
  c_float q[2] = {
      1.00,
      1.00,
  };

  // constraint grad
  c_float A_x[4] = {
      1.00,
      1.00,
      1.00,
      1.00,
  };
  c_int A_nnz = 4;
  c_int A_i[4] = {
      0,
      1,
      0,
      2,
  };
  c_int A_p[3] = {
      0,
      2,
      4,
  };

  // bounds
  c_float l[3] = {
      1.00,
      0.00,
      0.00,
  };
  c_float u[3] = {
      1.00,
      0.70,
      0.70,
  };
  c_int n = 2;  // number of vars
  c_int m = 3;  // number of constraints

  // osqp things
  OSQPSettings* settings = (OSQPSettings*)malloc(sizeof(OSQPSettings));
  OSQPWorkspace* workspace;
  OSQPData* data = (OSQPData*)malloc(sizeof(OSQPData));
  data->n = n;
  data->m = m;
  data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
  data->q = q;
  data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
  data->l = l;
  data->u = u;

  // load defaults
  osqp_set_default_settings(settings);
  // printf("default alpha: %f\n", settings->alpha);
  settings->alpha = 1.0;

  workspace = osqp_setup(data, settings);

  osqp_solve(workspace);

  printf("solution: %6.3f, %6.3f\n", workspace->solution->x[0],
         workspace->solution->x[1]);

  EXPECT_TRUE(fpEqual(workspace->solution->x[0], .3, .0001));
  EXPECT_TRUE(fpEqual(workspace->solution->x[1], .7, .0001));

  osqp_cleanup(workspace);
  c_free(data->A);
  c_free(data->P);
  c_free(data);
  c_free(settings);
}

TEST(OSQP, eigenDataOrder) {
  DMat<double> m(2, 3);
  m << 1, 2, 3, 4, 5, 6;

  double* ptr = m.data();

  double colMajor[6] = {1, 4, 2, 5, 3, 6};

  for (size_t i = 0; i < 6; i++) {
    EXPECT_TRUE(colMajor[i] == ptr[i]);
  }
}
