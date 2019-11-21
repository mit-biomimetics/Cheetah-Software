#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "../third-party/Goldfarb_Optimizer/QuadProg++.hh"
#include "Utilities/utilities.h"

TEST(Goldfarb_Optimizer, Goldfarb_opt_test) {
  // problem:
  // >> H = [4 1; 1 2]; f = [1;1]; A = [1 1; 1 0; 0 1; -1 -1; -1 0; 0 -1]; b =
  // [1 .7 .7 -1 0 0]';
  // >> x = quadprog(H,f,A,b)
  // x =
  //
  //    0.3000
  //    0.7000

  // hessian
  GMatr<double> G(2, 2);
  G[0][0] = 4.0;
  G[0][1] = 1.;
  G[1][0] = 1.;
  G[1][1] = 2.;

  GVect<double> g0(2);
  g0[0] = 1.;
  g0[1] = 1.;

  GMatr<double> CE(0., 2, 1);
  CE[0][0] = 1.;
  CE[1][0] = 1.;
  GVect<double> ce0(1);
  ce0[0] = -1.0;

  GMatr<double> CI(0., 2, 4);
  CI[0][0] = 1.;
  CI[1][1] = 1.;
  CI[0][2] = -1.;
  CI[1][3] = -1.;

  GVect<double> ci0(4);
  ci0[0] = 0.;
  ci0[1] = 0.;
  ci0[2] = 0.7;
  ci0[3] = 0.7;

  GVect<double> x;
  double f = solve_quadprog(G, g0, CE, ce0, CI, ci0, x);

  printf("cost: %f\n", f);
  printf("solution: %6.3f, %6.3f\n", x[0], x[1]);

  EXPECT_TRUE(fpEqual(x[0], .3, .0001));
  EXPECT_TRUE(fpEqual(x[1], .7, .0001));
}
