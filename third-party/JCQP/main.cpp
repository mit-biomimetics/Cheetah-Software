#include <iostream>

#include "QpProblem.h"
#include "ProblemGenerator.h"
#include "Timer.h"
#include "OsqpTest.h"
#include "qpOASES.h"

void test_example() {
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
    }
    if(Ax[i] > problem.u[i]) {
      printf("ub error (%.3f %.3f, %.3f)!\n", problem.l[i], problem.u[i], soln[i]);
    }
  }
}

void test_MPC_generator() {
  ProblemGenerator<double> pg;
  auto problem = pg.generateSparseMPC(12,12,300, 0);

  problem.settings.rho = 7;


  //Vector<double> qpOasesSoln = testQPOASES(problem);
  //Vector<double> osqpSoln = testOSQP(problem);

  problem.run(1575, true);
  Vector<double> mySoln = problem.getSolution();



  Vector<double> qf = 0.5 * mySoln.transpose() * problem.P * mySoln;
  double objective = qf(0,0) + problem.q.transpose() * mySoln;
  printf("my objective: %f\n", objective);

  /*
  qf = 0.5 * osqpSoln.transpose() * problem.P * osqpSoln;
  objective = qf(0,0) + problem.q.transpose() * osqpSoln;
  printf("osqp objective: %f\n", objective);

//  qf = 0.5 * qpOasesSoln.transpose() * problem.P * qpOasesSoln;
//  objective = qf(0,0) + problem.q.transpose() * qpOasesSoln;
//  printf("qpoases objective: %f\n", objective);

  Vector<double> soln_diff = mySoln - osqpSoln;
  fprintf(stderr, "diff (me to osqp): %f (rel %f), frac %f\n", soln_diff.norm(), mySoln.norm(), soln_diff.norm()/ mySoln.norm());
  printf("diff max (me to osqp) %g, %g\n", soln_diff.minCoeff(), soln_diff.maxCoeff());

//  soln_diff = mySoln - qpOasesSoln;
//  fprintf(stderr, "diff (me to qpoases): %f (rel %f), frac %f\n", soln_diff.norm(), mySoln.norm(), soln_diff.norm()/ mySoln.norm());
//  printf("diff max (me to qpoases) %g, %g\n", soln_diff.minCoeff(), soln_diff.maxCoeff());

//  Vector<double> Ax = problem.A * osqpSoln;
//  for(int i = 0; i < problem.m; i++) {
//    if(Ax(i,0) < problem.l(i,0) - 2e-3) {
//      printf("lb error %d (%.3f %.3f, %.3f)!\n", i, problem.l[i], problem.u[i], Ax[i]);
//    }
//    if(Ax(i,0) > problem.u(i,0) + 2e-3) {
//      printf("ub error %d (%.3f %.3f, %.3f)!\n", i, problem.l[i], problem.u[i], Ax[i]);
//    }
//  }
   */

}

#define NTPDF double
void test_solver() {
  DenseMatrix<NTPDF> A (100,100);
  A.setRandom();
  A = A * A.transpose();
  Vector<NTPDF> x(100);
  x.setRandom();

  Vector<NTPDF> b = x;

  CholeskyDenseSolver<NTPDF> solver;
  Timer tim;
  for(s64 i = 0; i < 3; i++)
  {
    solver.setup(A);
  }

  solver.solve(x);
  printf("solved in %.3f ms\n", tim.getMs() / 3);

  Vector<NTPDF> diff = A * x - b;
  printf("diff %g, %g\n", diff.minCoeff(), diff.maxCoeff());
}




int main()
{
  //test_solver();
  for(int i = 0; i < 10; i++)
    test_MPC_generator();
  return 0;
}
