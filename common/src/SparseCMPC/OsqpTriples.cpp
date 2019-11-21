#include "SparseCMPC/SparseCMPC.h"
#include "osqp.h"
#include <Utilities/Timer.h>

struct OsqpCSC {
  u32 nnz, m, n;
  c_float* values = nullptr;
  c_int* colPtrs = nullptr;
  c_int* rowIdx = nullptr;

  void alloc(u32 _n, u32 _nnz) {
    nnz = _nnz;
    n = _n;
    m = _n;
    values = new c_float[nnz];
    colPtrs = new c_int[n + 1];
    rowIdx = new c_int[nnz];
  }

  void freeAll() {
    delete[] values;
    delete[] colPtrs;
    delete[] rowIdx;
  }
};

static OsqpCSC compress(std::vector<SparseTriple<double>>& entries, u32 m, u32 n) {
  OsqpCSC result;
  u32 nnz = entries.size();
  result.alloc(n, nnz);
  result.m = m;

  u32 i = 0;
  result.colPtrs[0] = 0;
  for(u32 c = 0; c < n; c++) {
    u32 cNNZ = 0;
    while(entries[i].c == c && i < nnz) {
      result.values[i] = entries[i].value;
      result.rowIdx[i] = entries[i].r;
      assert(entries[i].r < m);
      assert(entries[i].c < n);
      assert(i < result.nnz);
      i++;
      cNNZ++;
    }
    result.colPtrs[c + 1] = result.colPtrs[c] + cNNZ;
    assert(c+1 < result.n + 1);
    assert(i == nnz || entries[i].c > c);
  }

  assert(i == nnz);

  return result;
}


void SparseCMPC::runSolverOSQP() {
  Timer timer;
  u32 varCount = 12 * _trajectoryLength + 3 * _bBlockCount;
  //printf("[SparseCMPC] Run with OSQP %d, %d\n", varCount, _constraintCount);
  assert(_constraintCount == _ub.size());
  assert(_constraintCount == _lb.size());
  assert(varCount == _linearCost.size());

  // Quadratic term
  sortAndSumTriples(_costTriples);
  OsqpCSC quadraticCostMatrix = compress(_costTriples, varCount, varCount);

  sortAndSumTriples(_constraintTriples);
  OsqpCSC constraintMatrix = compress(_constraintTriples, _constraintCount, varCount);

  OSQPSettings* settings = (OSQPSettings*)malloc(sizeof(OSQPSettings));
  OSQPWorkspace* workspace;
  (void)workspace;
  OSQPData* data = (OSQPData*)malloc(sizeof(OSQPData));

  data->n = varCount;
  data->m = _constraintCount;
  data->P = csc_matrix(varCount, varCount, quadraticCostMatrix.nnz,
      quadraticCostMatrix.values, quadraticCostMatrix.rowIdx, quadraticCostMatrix.colPtrs);
  data->q = _linearCost.data();
  data->A = csc_matrix(_constraintCount, varCount, constraintMatrix.nnz,
      constraintMatrix.values, constraintMatrix.rowIdx, constraintMatrix.colPtrs);
  data->l = _lb.data();
  data->u = _ub.data();

  //printf("t3: %.3f\n", timer.getMs());
  timer.start();

//# define EPS_ABS (1E-3)
//# define EPS_REL (1E-3)
//# define EPS_PRIM_INF (1E-4)
//# define EPS_DUAL_INF (1E-4)
  osqp_set_default_settings(settings);
  settings->eps_abs = 1e-5;
  settings->eps_rel = 1e-5;
  //settings->max_iter = 300;
  //settings->alpha = 1.0; //todo try me
  workspace = osqp_setup(data, settings);

  //printf("t4: %.3f\n", timer.getMs());
  timer.start();

  osqp_solve(workspace);

  //printf("t5: %.3f\n", timer.getMs());

  _result = Eigen::Matrix<float, Eigen::Dynamic, 1>(varCount);
  for(u32 i = 0; i < varCount; i++) {
    _result[i] = workspace->solution->x[i];
  }
//
//
//
//  osqp_cleanup(workspace);

//  QpProblem<double> solver(varCount, _constraintCount, false);
//  solver.A_triples = std::move(_constraintTriples);
//  solver.P_triples = std::move(_costTriples);
//  solver.u = Vector<double>(_ub.size());
//  solver.l = Vector<double>(_ub.size());
//  for(u32 i = 0; i < _constraintCount; i++) {
//    solver.u[i] = _ub[i];
//    solver.l[i] = _lb[i];
//  }
//  solver.q = Vector<double>(varCount);
//  for(u32 i = 0; i < varCount; i++) {
//    solver.q[i] = _linearCost[i];
//  }
//
//  solver.settings.alpha = 1.5;
//  solver.settings.rho = 1e-3;
//  solver.settings.terminate = 1e-6;
//  solver.settings.sigma = 1e-6;
//
//
//  solver.runFromTriples(-1, true);
//  _result = solver.getSolution().cast<float>();

  quadraticCostMatrix.freeAll();
  constraintMatrix.freeAll();
}
