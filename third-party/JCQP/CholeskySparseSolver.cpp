/*!
 * @file: CholeskySparseSolver.cpp
 *
 * Sparse Cholesky Solver.
 */

#include "CholeskySparseSolver.h"
#include "SparseMatrixMath.h"
#include "Timer.h"
#include "amd.h"
#include <iostream>
#include <immintrin.h>


template<typename T>
void CholeskySparseSolver<T>::preSetup(const DenseMatrix<T> &kktMat, bool b_print)
{
  Timer tim;
  // get sizes
  A.m = kktMat.rows();
  A.n = kktMat.cols();
  n = A.n;

  // count nnz;
  A.nnz = 0;
  for(u32 c = 0; c < A.n; c++) {
    for(u32 r = 0; r < A.m; r++) {
      if(kktMat(r,c) != 0. && r <= c) A.nnz++;
    }
  }

  if(b_print) printf("CHOLSPARSE PASS1 %.3f ms\n", tim.getMs());
  tim.start();

  if(b_print) printf("nnz %d, fill %.3f\n", A.nnz, (double)A.nnz / (double)(A.n*A.m));

  // allocate
  A.values = new T[A.nnz];
  A.colPtrs = new u32[A.n + 1];
  A.rowIdx = new u32[A.nnz];
  nnzLCol = new u32[n];
  parent = new s32[n];
  D = new T[n];
  rD = new T[n];
  tempSolve = new T[n];
  L.colPtrs = new u32[A.n + 1];

  if(b_print) printf("CHOLSPARSE ALLOC %.3f ms\n", tim.getMs());
  tim.start();

  // convert to sparse
  u32 i = 0;
  A.colPtrs[0] = 0;
  for(u32 c = 0; c < A.n; c++) {
    u32 cNNZ = 0;
    for(u32 r = 0; r < A.m; r++) {
      T v = kktMat(r,c);
      if(v != 0. && r <= c) {
        A.values[i] = v;
        A.rowIdx[i] = r;
        i++;
        cNNZ++;
      }
    }
    A.colPtrs[c + 1] = A.colPtrs[c] + cNNZ;
  }

  assert(i == A.nnz);

  if(b_print) printf("CHOLSPARSE SETUP %.3f ms\n", tim.getMs());
  tim.start();

  sanityCheck();
  if(b_print) printf("CHOLSPARSE CHECK %.3f ms\n", tim.getMs());
}

template<typename T>
void CholeskySparseSolver<T>::preSetup(const std::vector<SparseTriple<T>>& kktMat, u32 _n, bool b_print) {
  Timer tim;

  A.m = _n;
  A.n = _n;
  n = _n;
  A.nnz = kktMat.size();

  if(b_print) {
    printf("PRESETUP:\n\tnnz: %d, fill %.3f\n", A.nnz, (double)A.nnz / (double)(A.n*A.m));
  }

  // allocate
  A.values = new T[A.nnz];
  A.colPtrs = new u32[A.n + 1];
  A.rowIdx = new u32[A.nnz];
  nnzLCol = new u32[n];
  parent = new s32[n];
  D = new T[n];
  rD = new T[n];
  tempSolve = new T[n];
  L.colPtrs = new u32[A.n + 1];
  if(b_print) {
    printf("\tallocate time: %.3f ms\n", tim.getMs());
    tim.start();
  }

  // compress!
  u32 i = 0, j = 0;
  A.colPtrs[0] = 0;
  for(u32 c = 0; c < A.n; c++) {
    u32 cNNZ = 0;
    while(kktMat[j].c == c && j < A.nnz) {
      if(kktMat[j].r <= kktMat[j].c) {
        A.values[i] = kktMat[j].value;
        A.rowIdx[i] = kktMat[j].r;
        cNNZ++;
        i++;
      }
      j++;
    }
    A.colPtrs[c + 1] = A.colPtrs[c] + cNNZ;
    assert(j == A.nnz || kktMat[j].c > c);
  }

  assert(j == A.nnz);
  A.nnz = i;
}

template<typename T>
void CholeskySparseSolver<T>::setup(bool b_print)
{
  Timer tim;

  reorder();
  if(b_print) printf("CHOLSPARSE REORDER %.3f ms\n", tim.getMs());
  tim.start();

  L.nnz = symbolicFactor();
  if(b_print) printf("CHOLSPARSE SYM_FACTOR %.3f ms\n", tim.getMs());
  tim.start();
  if(b_print) printf("factor nnz %d, fill %.3f\n", L.nnz, (double)L.nnz / (double)(A.n*A.m));

  L.values = new T[L.nnz];
  L.rowIdx = new u32[L.nnz];
  factor();
  if(b_print) printf("CHOLSPARSE FACTOR %.3f ms\n", tim.getMs());
  tim.start();

#ifndef JCQP_USE_AVX2
  solveOrder();
  if(b_print) printf("CHOLSPARSE SOLVEORDER %.3f ms\n", tim.getMs());
#endif
}

template<typename T>
void CholeskySparseSolver<T>::sanityCheck()
{
  assert(A.colPtrs[0] == 0);
  for(u32 i = 0; i < n; i++) {
    assert(A.colPtrs[i] < A.colPtrs[i + 1]);
  }


  for(u32 i = 0; i < A.colPtrs[n]; i++) {
    assert(A.rowIdx[i] >= 0 && A.rowIdx[i] < n);
  }
}

template<typename T>
void CholeskySparseSolver<T>::reorder()
{
  P = new u32[n];
  rP = new u32[n];
  amdOrder(A,P,rP);

  MatCSC<T> permuted;
  permuted.alloc(n, A.nnz);

  u32* temp = new u32[n];
  for(u32 j = 0; j < n; j++) temp[j] = 0;

  for(u32 j = 0; j < n; j++) {
    u32 jNew = rP[j];

    for(u32 p = A.colPtrs[j]; p < A.colPtrs[j + 1]; p++) {
      u32 i = A.rowIdx[p];
      if(i > j) continue;
      u32 iNew = rP[i];
      temp[std::max(iNew, jNew)]++;
    }
  }

  // C.ptr
  // temp
  // n
  u32 k = 0;
  for(u32 i = 0; i < n; i++) {
    permuted.colPtrs[i] = k;
    k += temp[i];
    temp[i] = permuted.colPtrs[i];
  }
  permuted.colPtrs[n] = k;

  for(u32 j = 0; j < n; j++) {
    u32 jNew = rP[j];
    for(u32 p = A.colPtrs[j]; p < A.colPtrs[j+1]; p++) {
      u32 i = A.rowIdx[p];
      if(i > j) continue;
      u32 iNew = rP[i];
      u32 s = temp[std::max(iNew, jNew)]++;
      permuted.rowIdx[s] = std::min(iNew, jNew);
      permuted.values[s] = A.values[p];
    }
  }


  delete[] temp;
  A.freeAll();
  A = permuted;
}

template<typename T>
u32 CholeskySparseSolver<T>::symbolicFactor()
{
  s32* temp = new s32[n];
  for(u32 i = 0; i < n; i++) {
    temp[i] = 0;
    nnzLCol[i] = 0;
    parent[i] = -1;
  }

  for(u32 j = 0; j < n; j++) {
    temp[j] = j;
    for(u32 p = A.colPtrs[j]; p < A.colPtrs[j + 1]; p++) {
      s32 i = A.rowIdx[p];
      while((u32)temp[i] != j) {
        if(parent[i] == -1) {
          parent[i] = j;
        }
        nnzLCol[i]++;
        temp[i] = j;
        i = parent[i];
      }
    }
  }

//  for(u32 k = 0; k < n; k++) { // loop over cols
//    parent[k] = -1; // no parent assigned
//    temp[k] = k;    // visit
//    nnzLCol[k] = 0; // no nonzeros
//
//    for(u32 p = A.colPtrs[k]; p < A.colPtrs[k + 1]; p++) { // loop over entries in column
//      s32 i = A.rowIdx[p]; // get row
//      if(i < k) {
//        for(;temp[i] != k; i = parent[i]) { // ascend tree
//          if(parent[i] == -1) {
//            parent[i] = k;
//          }
//          nnzLCol[i]++;
//          temp[i] = k;
//        }
//      }
//    }
//  }
//
  L.colPtrs[0] = 0;
  for(u32 i = 0; i < n; i++) {
    L.colPtrs[i + 1] = L.colPtrs[i] + nnzLCol[i];
  }

  delete[] temp;
  return L.colPtrs[n];
}

template<typename T>
void CholeskySparseSolver<T>::factor()
{

  T* y = new T[n];
  u32* next = new u32[n];
  u32* yIdx = new u32[n];
  u8* colUsed = new u8[n];
  u32* elims = new u32[n];

  L.colPtrs[0] = 0;
  for(s32 i = 0; i < (s32)n; i++) {
    L.colPtrs[i + 1] = L.colPtrs[i] + nnzLCol[i]; // col ptrs
    colUsed[i] = 0; // unused
    y[i] = 0.;
    D[i] = 0.;
    next[i] = L.colPtrs[i];
  }

  D[0] = A.values[0];
  rD[0] = 1./D[0];

  for(s32 k = 1; k < (s32)n; k++) {
    u32 nnzRow = 0;
    u32 end = A.colPtrs[k + 1];
    for(s32 i = A.colPtrs[k]; i < (s32)end; i++) {
      u32 bIdx = A.rowIdx[i];
      if(bIdx == (u32)k) {
        D[k] = A.values[i];
        continue;
      }

      y[bIdx] = A.values[i];
      s32 nextIdx = bIdx;
      if(colUsed[nextIdx] == 0) {
        colUsed[nextIdx] = 1;
        elims[0] = nextIdx;
        u32 nnzElim = 1;
        nextIdx = parent[bIdx];

        while(nextIdx != -1 && nextIdx < k) {
          if(colUsed[nextIdx] == 1) break;
          colUsed[nextIdx] = 1;
          elims[nnzElim++] = nextIdx;
          nextIdx = parent[nextIdx];
        }


        while(nnzElim) {
          yIdx[nnzRow++] = elims[--nnzElim];
        }
      }
    }

    for(s32 i = nnzRow - 1; i >= 0; i--) {
      u32 cSel = yIdx[i];

      u32 end2 = next[cSel];
      T yi = y[cSel];
      for(s32 j = L.colPtrs[cSel]; j < (s32)end2; j++) {
        y[L.rowIdx[j]] -= L.values[j] * yi;
      }

      L.rowIdx[end2] = k;
      L.values[end2] = yi * rD[cSel];
      D[k] -= yi * L.values[end2];
      next[cSel]++;

      y[cSel] = 0.;
      colUsed[cSel] = 0;
    }

    rD[k] = 1./D[k];
  }

  delete[] y;
  delete[] next;
  delete[] yIdx;
  delete[] colUsed;
  delete[] elims;

//  s32* flag = new s32[n];
//  s32* pat = new s32[n];
//  for(s32 k = 0; k < n; k++) {
//    y[k] = 0.;
//    s32 top = n;
//    flag[k] = k;
//    nnzLCol[k] = 0;
//    for(s32 p = A.colPtrs[k]; p < A.colPtrs[k + 1]; p++) {
//      s32 i = A.rowIdx[p]; // get row
//      if(i <= k) {
//        y[i] += A.values[p];
//        u32 len;
//        for(len = 0; flag[i] != k; i = parent[i]) {
//          pat[len++] = i;
//          flag[i] = k;
//        }
//
//        while(len > 0) {
//          pat[--top] = pat[--len];
//        }
//      }
//    }
//
//    D[k] = y[k];
//    y[k] = 0.;
//    for(; top < n; top++) {
//      s32 i = pat[top];
//      T yi = y[i];
//      y[i] = 0.;
//      s32 p;
//      for(p = L.colPtrs[i]; p < L.colPtrs[i] + nnzLCol[i]; p++) {
//        y[L.rowIdx[p]] -= L.values[p] * yi;
//      }
//      T l_ki = yi / D[i];
//      D[k] -= l_ki * yi;
//      L.rowIdx[p] = k;
//      L.values[p] = l_ki;
//      nnzLCol[i]++;
//    }
//    assert(D[k] != 0.);
//  }

  //for(int i = 0; i < n; i++) D[i] = 1. / D[i];
}

#ifdef JCQP_USE_AVX2
__attribute__ ((unused))
static void spTriAVX(u32 n, const u32* p, const u32* i, const double* v, double* x)
{
  constexpr int UNROLL = 1;

  for(u32 c = 0; c < n; c++) {

    __m256d xV = _mm256_set1_pd(x[c]);     // xV = {x[c], x[c], x[c], x[c]}

    s32 a = p[c];
    s32 end = p[c+1];
    s32 vEnd = end - 4 * UNROLL;
    const double* matPtr = v + a;
    const u32* iPtr = i + a;
    for(; a < vEnd; a += 4 * UNROLL) {

      for(int u = 0; u < UNROLL; u++) {
        __m256d mat = _mm256_loadu_pd(matPtr);
        __m256d lhs = _mm256_set_pd(x[iPtr[3]], x[iPtr[2]], x[iPtr[1]], x[iPtr[0]]); // todo gather
        __m256d fma = _mm256_fnmadd_pd(xV, mat, lhs);
        //__m256d fma = _mm256_mul_pd(xV, mat);

        x[iPtr[0]] = fma[0];
        x[iPtr[1]] = fma[1];
        x[iPtr[2]] = fma[2];
        x[iPtr[3]] = fma[3];

        matPtr += 4;
        iPtr += 4;
      }
    }

    for(; a < end; a++) {
      x[iPtr[0]] -= matPtr[0] * x[c];
      matPtr++;
      iPtr++;
    }

  }
}

__attribute__ ((unused))
static void spTriTAVX(u32 n, const u32* p, const u32* i, const double* v, double* x) {
  constexpr int UNROLL = 1; // 440 best
  for(s32 c = n - 1; c >= 0; c--) {
    s32 a = p[c];
    s32 end = p[c+1];
    s32 vEnd = end - 4 * UNROLL;
    const u32* iPtr = i + a;


    __m256d acc = _mm256_set1_pd(0.);

    for(; a < vEnd; a += 4 * UNROLL) {
      for(int u = 0; u < UNROLL; u++) {
        __m256d mat = _mm256_loadu_pd(v);
        __m256d rhs = _mm256_set_pd(x[iPtr[3]], x[iPtr[2]], x[iPtr[1]], x[iPtr[0]]);
        acc = _mm256_fmadd_pd(mat, rhs, acc);

        v += 4;
        iPtr += 4;
      }
    }

    x[c] -= acc[0];
    x[c] -= acc[1];
    x[c] -= acc[2];
    x[c] -= acc[3];

    for(;a < end; a++) {
      x[c] -= *(v++) * x[*(iPtr++)];
    }
  }
}

/*!
 * Backsubstitution solve L^{-T}x = b in place.
 * @param n - matrix size
 * @param p - column pointers of sparse matrix
 * @param i - row indices of sparse matrix
 * @param v - values of sparse matrix
 * @param x - input/output (x,b)
 */
__attribute__ ((unused))
static void spTriTAVX2(u32 n, const u32* p, const u32* i, const double* v, double* x) {
  constexpr int UNROLL = 1;
  for(s32 c = n - 1; c >= 0; c--) { // loop backward over columns
    s32 colStart = p[c]; // start idx of column
    s32 end = p[c+1];    // end idx of column
    s32 vEnd = end - 4 * UNROLL; // end idx of vectorized/unrolled
    const u32* iPtr = i + colStart; // indices

    const double* matPtr = v + colStart; // matrix
    __m256d acc = _mm256_set1_pd(0.);

    for(; colStart < vEnd; colStart += 4 * UNROLL) {
      for(int u = 0; u < UNROLL; u++) {
        __m256d mat = _mm256_loadu_pd(matPtr);
        __m256d rhs = _mm256_set_pd(x[iPtr[3]], x[iPtr[2]], x[iPtr[1]], x[iPtr[0]]);
        acc = _mm256_fmadd_pd(mat, rhs, acc);

        matPtr += 4;
        iPtr += 4;
      }
    } // end for

    // horizontal sum and set x
    x[c] -= acc[0];
    x[c] -= acc[1];
    x[c] -= acc[2];
    x[c] -= acc[3];

    // whatever didn't fit in the vectors
    for(;colStart < end; colStart++) {
      x[c] -= *(matPtr++) * x[*(iPtr++)];
    }
  }
}
#endif

template<typename T>
void CholeskySparseSolver<T>::solve(Vector<T>& outv) {
  T* outP = outv.data();
  T* out = tempSolve;

  for(u32 j = 0; j < n; j++) {
    out[j] = outP[P[j]];
  }

#ifdef JCQP_USE_AVX2
  spTriAVX(n, L.colPtrs, L.rowIdx, L.values, out);
#else
  for(u32 j = 0; j < n; j++) {
    u32 end = L.colPtrs[j + 1];
    for(u32 p = L.colPtrs[j]; p < end; p++) {
      out[L.rowIdx[p]] -= L.values[p] * out[j];
    }
  }
#endif



  for(u32 i = 0; i < n; i++) {
    out[i] *= rD[i];
  }

#ifdef JCQP_USE_AVX2
  spTriTAVX2(n, L.colPtrs, L.rowIdx, L.values, out);
#else
  for(s32 i = n - 1; i >= 0; i--) {
    u32 end = L.colPtrs[i + 1];
    for(u32 j = L.colPtrs[i]; j < end; j++) {
      out[i] -= L.values[j] * out[L.rowIdx[j]];
    }
  }
#endif

  for(u32 j = 0; j < n; j++) {
    outP[P[j]] = out[j];
  }
}

/*!
 * order which sparse matrix is accessed on Ltranspose solve.  not used
 */
template<typename T>
void CholeskySparseSolver<T>::solveOrder()
{
  reverseOrder = new T[L.nnz];
  u32 c = 0;
  for(s32 i = n - 1; i >= 0; i--) {
    u32 end = L.colPtrs[i+1];
    for (u32 j = L.colPtrs[i]; j < end; j++) {
      reverseOrder[c++] = L.values[j];
    }
  }
}

/*!
 * Use AMD library to find minimum degree permutation.
 */
template<typename T>
void CholeskySparseSolver<T>::amdOrder(MatCSC<T> &mat, u32 *perm, u32* iperm)
{
  T amdInfo[AMD_INFO];
  int amdRV = amd_order((int)mat.n, (int*)mat.colPtrs, (int*)mat.rowIdx, (int*)perm, nullptr, amdInfo);
  if(amdRV < 0) throw std::runtime_error("reorder failed");

  for(u32 i = 0; i < mat.n; i++) {
    iperm[perm[i]] = i;
  }
}

template class CholeskySparseSolver<double>;
