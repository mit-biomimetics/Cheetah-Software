#ifndef QPSOLVER_CHOLESKYSPARSESOLVER_H
#define QPSOLVER_CHOLESKYSPARSESOLVER_H

#include "types.h"

template<typename T>
struct MatCSC {
  u32 nnz, m, n;
  T* values = nullptr; // x
  u32* colPtrs = nullptr; // p
  u32* rowIdx = nullptr; // i

  void alloc(u32 _n, u32 _nnz) {
    nnz = _nnz;
    n = _n;
    m = _n;
    values = new T[nnz];
    colPtrs = new u32[n + 1];
    rowIdx = new u32[nnz];
  }

  void freeAll() {
    delete[] values;
    delete[] colPtrs;
    delete[] rowIdx;
  }
};


template<typename T>
class CholeskySparseSolver
{
public:
  CholeskySparseSolver() = default;
  void preSetup(const DenseMatrix<T>& kktMat, bool b_print = true);
  void setup(bool b_print = true);
  void solve(Vector<T>& out);
  void amdOrder(MatCSC<T>& mat, u32* perm, u32* iperm);

  ~CholeskySparseSolver() {
    A.freeAll();
    L.freeAll();
    delete[] reverseOrder;
    delete[] nnzLCol;
    delete[] P;
    delete[] D;
    delete[] rP;
    delete[] rD;
    delete[] parent;
    delete[] tempSolve;
  }
private:
  void reorder();
  u32 symbolicFactor();
  void factor();
  void sanityCheck();
  void solveOrder();
  MatCSC<T> A;
  MatCSC<T> L;
  u32 n;
  T* tempSolve = nullptr;
  T* reverseOrder = nullptr;
  u32* nnzLCol = nullptr; // # of nonzeros per column in L
  u32* P = nullptr;       // reorder permutation
  u32* rP = nullptr;      // reorder inverse permutation
  T*   D = nullptr;       // Diagonal
  T*   rD = nullptr;      // inverse diagonal
  s32* parent = nullptr;  // tree
};


#endif //QPSOLVER_CHOLESKYSPARSESOLVER_H
