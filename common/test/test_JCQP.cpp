#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "../third-party/JCQP/QpProblem.h"
#include "../third-party/JCQP/ProblemGenerator.h"
#include "../third-party/JCQP/Timer.h"
#include "Utilities/utilities.h"


TEST(JCQP, test_example){
  QpProblem<double> problem(2,3);
  problem.A << 1,2,4,4,5,6;
  problem.P << 12,5,5,8;
  problem.u << 1,2,3;
  problem.l << -1,-2,-3;
  problem.q << 1,2;
  problem.run();

  Vector<double> soln = problem.getSolution();
  Vector<double> qf = 0.5 * soln.transpose() * problem.P * soln;

  double objective = qf(0,0) + problem.q.transpose() * soln;
  printf("objective: %f\n", objective);


  Vector<double> Ax = problem.A * soln;
  for(int i = 0; i < problem.m; i++) {
    if(Ax[i] < problem.l[i]) {
      printf("lb error (%.3f %.3f, %.3f)!\n", problem.l[i], problem.u[i], soln[i]);
      EXPECT_TRUE(false);
    }
    if(Ax[i] > problem.u[i]) {
      printf("ub error (%.3f %.3f, %.3f)!\n", problem.l[i], problem.u[i], soln[i]);
      EXPECT_TRUE(false);
    }
  }
}

TEST(JCQP, test_result_dense) {
  QpProblem<double> problem(2,6);
  problem.A << 1, 1,
               1, 0,
               0, 1,
               -1, -1,
               -1, 0,
               0, -1;
  problem.P << 4,1,1,2;
  problem.u << 1, 0.7, 0.7, -1, 0, 0;
  problem.l << -1000, -1000, -1000, -1000, -1000, -1000;
  problem.q << 0, 0;
  problem.run();

  EXPECT_TRUE(fpEqual(problem.getSolution()[0],0.3, .0001));
  EXPECT_TRUE(fpEqual(problem.getSolution()[1],0.7, .0001));
}


TEST(JCQP, test_result_sparse) {
  QpProblem<double> problem(2,3);
  problem.A << 1, 1,
      1., 0.,
      0., 1.;
  problem.P << 4.,1.,1.,2.;
  problem.u << 1., 0.7, 0.7;
  problem.l << 1., 0., 0.;
  problem.q << 0, 0;
  problem.run(-1, true);

  EXPECT_TRUE(fpEqual(problem.getSolution()[0],0.3, .0001));
  EXPECT_TRUE(fpEqual(problem.getSolution()[1],0.7, .0001));
  std::cout << "-----> sol: " << problem.getSolution().transpose() << "\n";
}


TEST(JCQP, test_solver_dense_double){
  DenseMatrix<double> A (100,100);
  A.setRandom();
  A = A * A.transpose();
  Vector<double> x(100);
  x.setRandom();

  Vector<double> b = x;

  CholeskyDenseSolver<double> solver(false);
  for(s64 i = 0; i < 3; i++)
  {
    solver.setup(A);
  }

  solver.solve(x);

  Vector<double> diff = A * x - b;
  EXPECT_TRUE(fpEqual(diff.minCoeff(), 0.0, 0.001));
  EXPECT_TRUE(fpEqual(diff.maxCoeff(), 0.0, 0.001));
}

TEST(JCQP, test_solver_dense_float) {
  DenseMatrix<float> A (10,10);
  A.setRandom();
  A = A * A.transpose();
  Vector<float> x(10);
  x.setRandom();

  Vector<float> b = x;

  CholeskyDenseSolver<float> solver(false);
  for(s64 i = 0; i < 3; i++)
  {
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
  DenseMatrix<double> A (n,n);
  A.setRandom();
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      if((i + 3 * j) % 7 == 1) {
        A(i,j) = 0;
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
