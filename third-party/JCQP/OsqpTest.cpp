
#include "OsqpTest.h"
#include <eigen3/Eigen/SparseCore>
#include "OsqpEigen/Solver.hpp"
#include "OsqpEigen/Data.hpp"
#include "OsqpEigen/Settings.hpp"
#include "Timer.h"


Vector<double> testOSQP(QpProblem<double>& problem)
{
  OsqpEigen::Solver solver;
//  for(s64 i = 0; i < problem.n; i++)
//  {
//    for(s64 j = 0; j < problem.n; j++)
//    {
//      problem.P(i,j) += 1e-8;
//    }
//  }
//
//  for(s64 i = 0; i < problem.m; i++)
//  {
//    for(s64 j = 0; j < problem.n; j++)
//    {
//      problem.A(i,j) += 1e-8;
//    }
//  }



//  for(s64 i = 0; i < problem.n; i++)
//  {
//    for(s64 j = 0; j < problem.n; j++)
//    {
//      pSparse.coeffRef(i,j) -= 1e-8;
//    }
//  }
//
//  for(s64 i = 0; i < problem.m; i++)
//  {
//    for(s64 j = 0; j < problem.n; j++)
//    {
//      aSparse.coeffRef(i,j) -= 1e-8;
//    }
//  }

  solver.data()->setNumberOfVariables(problem.n);
  solver.data()->setNumberOfConstraints(problem.m);


  solver.data()->setLowerBound(problem.l);
  solver.data()->setUpperBound(problem.u);
  solver.settings()->setScaling(0);
  solver.settings()->setAdaptiveRho(false);
  solver.settings()->setRho(7);
  solver.settings()->setWarmStart(false);

  Timer t;
  Eigen::SparseMatrix<double> pSparse = problem.P.sparseView();
  Eigen::SparseMatrix<double> aSparse = problem.A.sparseView();
  solver.data()->setHessianMatrix(pSparse);
  solver.data()->setGradient(problem.q);
  solver.data()->setLinearConstraintsMatrix(aSparse);

  solver.initSolver();
  solver.solve();
  fprintf(stderr, "%.3f\n", t.getMs());

  return solver.getSolution();
}
