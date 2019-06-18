#include "qpOASES.h"
#include "Timer.h"

using qprt = qpOASES::real_t;


void eigenToQP(const DenseMatrix<double> &mat, qprt* m, u64 r, u64 c) {
  u64 idx = 0;
  for(int i = 0; i < r; i++) {
    for(int j = 0; j < c; j++) {
      m[idx++] = mat(i,j);
    }
  }
}

Vector<double> testQPOASES(QpProblem<double>& problem)
{
  qprt* H_qpo  = new qprt[problem.n * problem.n];
  qprt* g_qpo  = new qprt[problem.n];
  qprt* A_qpo  = new qprt[problem.m * problem.n];
  qprt* lb_qpo = new qprt[problem.m];
  qprt* ub_qpo = new qprt[problem.m];
  qprt*  x_opt = new qprt[problem.n];

//  memcpy(H_qpo, problem.P.data(), problem.n * problem.n * sizeof(double));

  eigenToQP(problem.P, H_qpo, problem.n, problem.n);
  eigenToQP(problem.q, g_qpo, problem.n, 1);
  eigenToQP(problem.A, A_qpo, problem.m, problem.n);
  eigenToQP(problem.l, lb_qpo, problem.m, 1);
  eigenToQP(problem.u, ub_qpo, problem.m, 1);
  memset(x_opt, 0, problem.n * sizeof(double));

  qpOASES::QProblem qpp(problem.n, problem.m);
  qpOASES::Options op;
  //op.setToMPC();
  op.printLevel = qpOASES::PL_NONE;
  qpp.setOptions(op);
  qpOASES::int_t nWSR = 10000;
  Timer tim;
  qpp.init(H_qpo, g_qpo, A_qpo, NULL, NULL, lb_qpo, ub_qpo, nWSR);
  qpp.getPrimalSolution(x_opt);
  printf("QPOASES solve time: %.3f ms\n", tim.getMs());

  Vector<double> result(problem.n);
  memcpy(result.data(), x_opt, problem.n * sizeof(double));

  delete[] H_qpo;
  delete[] g_qpo;
  delete[] A_qpo;
  delete[] lb_qpo;
  delete[] ub_qpo;
  delete[] x_opt;

  return result;
}