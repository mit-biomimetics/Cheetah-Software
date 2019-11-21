#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "../third-party/JCQP/ProblemGenerator.h"
#include "../third-party/JCQP/QpProblem.h"
#include "../third-party/JCQP/Timer.h"
#include "Utilities/utilities.h"

TEST(JCQP, test_example) {
  QpProblem<double> problem(2, 3);
  problem.A << 1, 2, 4, 4, 5, 6;
  problem.P << 12, 5, 5, 8;
  problem.u << 1, 2, 3;
  problem.l << -1, -2, -3;
  problem.q << 1, 2;
  problem.runFromDense();

  Vector<double> soln = problem.getSolution();
  Vector<double> qf = 0.5 * soln.transpose() * problem.P * soln;

  double objective = qf(0, 0) + problem.q.transpose() * soln;
  printf("objective: %f\n", objective);

  Vector<double> Ax = problem.A * soln;
  for (int i = 0; i < problem.m; i++) {
    if (Ax[i] < problem.l[i]) {
      printf("lb error (%.3f %.3f, %.3f)!\n", problem.l[i], problem.u[i],
             soln[i]);
      EXPECT_TRUE(false);
    }
    if (Ax[i] > problem.u[i]) {
      printf("ub error (%.3f %.3f, %.3f)!\n", problem.l[i], problem.u[i],
             soln[i]);
      EXPECT_TRUE(false);
    }
  }
}

TEST(JCQP, test_result_dense) {
  QpProblem<double> problem(2, 6);
  problem.A << 1, 1, 1, 0, 0, 1, -1, -1, -1, 0, 0, -1;
  problem.P << 4, 1, 1, 2;
  problem.u << 1, 0.7, 0.7, -1, 0, 0;
  problem.l << -1000, -1000, -1000, -1000, -1000, -1000;
  problem.q << 0, 0;
  problem.runFromDense();

  EXPECT_TRUE(fpEqual(problem.getSolution()[0], 0.3, .0001));
  EXPECT_TRUE(fpEqual(problem.getSolution()[1], 0.7, .0001));
}

TEST(JCQP, test_result_sparse) {
  QpProblem<double> problem(2, 3);
  problem.A << 1, 1, 1., 0., 0., 1.;
  problem.P << 4., 1., 1., 2.;
  problem.u << 1., 0.7, 0.7;
  problem.l << 1., 0., 0.;
  problem.q << 0, 0;
  problem.runFromDense(-1, true);

  EXPECT_TRUE(fpEqual(problem.getSolution()[0], 0.3, .0001));
  EXPECT_TRUE(fpEqual(problem.getSolution()[1], 0.7, .0001));
}

TEST(JCQP, test_sort_triples) {
  std::vector<SparseTriple<double>> tris;
  tris.push_back({1,1,1});
  tris.push_back({0,1,0});
  tris.push_back({0,0,1});
  tris.push_back({0,0,0});

  sortTriples(tris, true);
  EXPECT_TRUE(tris[0].r == 0 && tris[0].c == 0);
  EXPECT_TRUE(tris[1].r == 1 && tris[1].c == 0);
  EXPECT_TRUE(tris[2].r == 0 && tris[2].c == 1);
  EXPECT_TRUE(tris[3].r == 1 && tris[3].c == 1);

  tris.push_back({10,1,1});
  EXPECT_THROW(sortTriples(tris, true), std::runtime_error);

  sortAndSumTriples(tris);
  EXPECT_TRUE(tris[0].value == 10 + 1);

  tris.clear();
  tris.push_back({1,1,1});
  tris.push_back({0,1,0});
  tris.push_back({0,0,1});
  tris.push_back({0,0,0});
  tris.push_back({10,1,1});

  sortAndSumTriples(tris);
  EXPECT_TRUE(tris[0].value == 10 + 1);
}

TEST(JCQP, test_result_sparse_triples) {
  QpProblem<double> problem(2,3);
  problem.A_triples = {{1.,0,0}, {1.,0,1}, {1., 1, 0},
                       {1.,2,1}};
  problem.P_triples = {{4.,0,0}, {1.,0,1},
                       {1.,1,0}, {2.,1,1}};
  problem.u << 1., 0.7, 0.7;
  problem.l << 1., 0., 0.;
  problem.q << 0, 0;
  problem.runFromTriples(-1, true);
  EXPECT_TRUE(fpEqual(problem.getSolution()[0], 0.3, .0001));
  EXPECT_TRUE(fpEqual(problem.getSolution()[1], 0.7, .0001));
}

TEST(JCQP, test_solver_dense_double) {
  DenseMatrix<double> A(100, 100);
  A.setRandom();
  A = A * A.transpose();
  Vector<double> x(100);
  x.setRandom();

  Vector<double> b = x;

  CholeskyDenseSolver<double> solver(false);
  for (s64 i = 0; i < 3; i++) {
    solver.setup(A);
  }

  solver.solve(x);

  Vector<double> diff = A * x - b;
  EXPECT_TRUE(fpEqual(diff.minCoeff(), 0.0, 0.001));
  EXPECT_TRUE(fpEqual(diff.maxCoeff(), 0.0, 0.001));
}

TEST(JCQP, test_solver_dense_float) {
  DenseMatrix<float> A(10, 10);
  A.setRandom();
  A = A * A.transpose();
  Vector<float> x(10);
  x.setRandom();

  Vector<float> b = x;

  CholeskyDenseSolver<float> solver(false);
  for (s64 i = 0; i < 3; i++) {
    solver.setup(A);
  }

  solver.solve(x);

  Vector<float> diff = A * x - b;
  EXPECT_TRUE(fpEqual(diff.minCoeff(), 0.0f, 0.001f));
  EXPECT_TRUE(fpEqual(diff.maxCoeff(), 0.0f, 0.001f));
}

TEST(JCQP, test_solver_sparse) {
  int n = 300;
  Vector<double> x(n);
  x.setRandom();
  Vector<double> b = x;
  DenseMatrix<double> A(n, n);
  A.setRandom();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if ((i + 3 * j) % 7 == 1) {
        A(i, j) = 0;
      }
    }
  }

  A = A * A.transpose();
  SparseMatrix<double> As = A.sparseView();
  CholeskySparseSolver<double> solver;
  solver.preSetup(As);
  solver.setup();
  solver.solve(x);
  Vector<double> diff = As * x - b;
  EXPECT_TRUE(fpEqual(diff.minCoeff(), 0.0, 0.001));
  EXPECT_TRUE(fpEqual(diff.maxCoeff(), 0.0, 0.001));

  printf("diff %g %g\n", diff.minCoeff(), diff.maxCoeff());
}
