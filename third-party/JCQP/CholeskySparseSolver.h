#ifndef QPSOLVER_CHOLESKYSPARSESOLVER_H
#define QPSOLVER_CHOLESKYSPARSESOLVER_H

#include "types.h"
#include "SparseMatrixMath.h"




template<typename T>
class CholeskySparseSolver
{
public:
  CholeskySparseSolver() = default;
  void preSetup(const DenseMatrix<T>& kktMat, bool b_print = true);
  void preSetup(const std::vector<SparseTriple<T>>& kktMat, u32 n, bool b_print = true);
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
  SparseTriple<T> A_triple;
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
