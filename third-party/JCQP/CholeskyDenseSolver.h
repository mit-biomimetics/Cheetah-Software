/*!
 * @file: CholeskyDenseSolver.h
 *
 * Dense Cholesky Solver.
 */

#ifndef QPSOLVER_CHOLESKYDENSESOLVER_H
#define QPSOLVER_CHOLESKYDENSESOLVER_H

#include "types.h"

template<typename T>
class CholeskyDenseSolver
{
public:
  CholeskyDenseSolver(bool print) : _print(print) { }
  ~CholeskyDenseSolver() {
    delete[] pivots;
    delete[] solve1;
    delete[] solve2;
  }
  void setup(const DenseMatrix<T>& kktMat);
  void solve(Vector<T>& in);

  DenseMatrix<T>& getInternal() { return L; }
  DenseMatrix<T> getReconstructedPermuted();

private:
  DenseMatrix<T> L;
  void solveAVX(Vector<T>& in);
  void setupAVX(T* mat, T* result, T* vec, u64 row);
  T* solve1 = nullptr, *solve2 = nullptr;
  Vector<T> D;
  s64* pivots = nullptr;
  s64 n = 0;
  bool _print;
};


#endif //QPSOLVER_CHOLESKYDENSESOLVER_H
