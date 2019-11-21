/*!
 * @file: CholeskyDenseSolver.cpp
 *
 * Dense Cholesky Solver.
 */

#include "CholeskyDenseSolver.h"
#include "Timer.h"
#include <cassert>
#include <iostream>
#include <immintrin.h>

static constexpr s64 UNROLL_MATVEC = 8; //! Loop unroll for all inner loops, AVX2 only

/*!
 * Perform initial factorization
 */
template<typename T>
void CholeskyDenseSolver<T>::setup(const DenseMatrix<T> &kktMat)
{

  Timer tim;
  assert(kktMat.cols() == kktMat.rows()); // must be a square matrix
  n = kktMat.rows();

  // setup
  L.resize(n,n); // Lower triangular
  D.resize(n);   // Diagonal
  L = kktMat;    // solved in place for better cache usage
  D.setZero();
  pivots = new s64[n]; // permutations for pivoting

  Vector<T> temp(n); // O(n) temporary storage, used in a few places

  if(_print)
    printf("CHOLDENSE SETUP %.3f ms\n", tim.getMs());

  tim.start();

  // in place dense LDLT algorithm
  for(s64 row = 0; row < n; row++)
  {
    // step 1: find the largest row to pivot with
    T largestDiagValue = -std::numeric_limits<T>::infinity();
    s64 pivot = -1;
    for(s64 rowSel = row; rowSel < n; rowSel++)
    {
      T v = std::abs(L(rowSel, rowSel));
      if(v > largestDiagValue)
      {
        largestDiagValue = v;
        pivot = rowSel;
      }
    }
    assert(pivot != -1);

    // step 2: PIVOT
    pivots[row] = pivot;
    if(pivot != row) // but only if we need to (pivot > row)!
    {
      assert(pivot > row);
      // swap row
      for(s64 i = 0; i < row; i++)
      {
        std::swap(L(row, i), L(pivot, i));
      }

      // swap col
      for(s64 i = pivot + 1; i < n; i++)
      {
        std::swap(L(i, row), L(i, pivot));
      }

      // swap elt
      std::swap(L(row, row), L(pivot, pivot));

      for(s64 i = row + 1; i < pivot; i++)
      {
        std::swap(L(i, row), L(pivot, i));
      }
    }

    // step 3, FACTOR:
    // wikipedia equation: A_ij - L_ij * D_j = sum_{k=1}^{j-1} L_ik * L_jk * D_k
    // wikipedia eqn 2   : D_j = A_jj - sum_{k = 1}^{j-1} L_jk^2 D_k
    // wikipedia "i" is my "row"
    if(row > 0)
    {
      for(s64 i = 0; i < row; i++)
      {
        temp[i] = L(i, i) * L(row, i); // temp[k] = D_k * L_ik
      }

      for(s64 i = 0; i < row; i++)
      {
        L(row, row) -= temp[i] * L(row, i); // D_i = A_ii - sum_k D_k * L_ik
      }

      //T factor = L(row, row);

      // compute ith (row-th) column of L
      // all L's are assumed to be in the row-th column
      // reuse the L_ik * D_k from earlier.
      // instead of solving 1-by-1, compute contributions from columns because matrix is col-major order.
      if(n - row - 1 > 0)
      {
        T* mat = L.data() + row + 1;
        T* result = L.data() + row + 1 + n * row;
        T* vec = temp.data();
        setupAVX(mat, result, vec, row); // O(n^2) inner loop
      }
    }

    // step 4, divide.
    T diag = L(row, row);
    if(std::abs(diag) < 1e-6) {
      printf("ERROR: LDLT fail (poorly conditioned pivot %g)\n", diag);
      assert(false);
    }

    if(n - row - 1 > 0) {
      T mult = T(1) / diag;
      for(s64 i = row + 1; i < n; i++) {
        L(i, row) *= mult;
      }
    }
  }

  if(_print)
    printf("CHOLDENSE FACT %.3f ms\n", tim.getMs());

  tim.start();

  // set up memory:
  solve1 = new T[n * (n - 1)]; // elements of L ordered for first triangular solve
  solve2 = new T[n * (n - 1)]; // elements of L ordered for second triangular solve
  s64 c = 0;
  for(s64 j = 0; j < n; j++)
  {
    for(s64 i = j + 1; i < n; i++)
    {
      solve1[c++] = L(i,j); // i > j
    }
  }
  c = 0;
  for(s64 j = n; j --> 0;)
  {
    for(s64 i = 0; i < j; i++)
    {
      solve2[c++] = L(j, i);
    }
  }

  if(_print)
    printf("CHOLDENSE FIN %.3f ms\n", tim.getMs());
}

/*!
 * Inner O(n^2) loop used in factorization, written with AVX intrinsics
 */
template<>
void CholeskyDenseSolver<double>::setupAVX(double *mat, double *result, double *vec, u64 row)
{
  s64 R = n - row - 1;
  s64 C = row;
  double* matPtr;
  double* resultPtr;
  for(s64 c = 0; c < C; c++)
  {

#ifdef JCQP_USE_AVX2
    __m256d vecData = _mm256_set1_pd(vec[c]); // broadcast to all
#endif
    resultPtr = result;
    matPtr = mat + c * n;
    s64 r = 0;
    r = 0;

#ifdef JCQP_USE_AVX2
    s64 end = R - 4 * UNROLL_MATVEC; // gcc fails to make this optimization.
    for(; r < end; r += 4 * UNROLL_MATVEC)
    {
      for(s64 u = 0; u < UNROLL_MATVEC; u++)
      {
        __m256d matData = _mm256_loadu_pd(matPtr);
        __m256d resultData = _mm256_loadu_pd(resultPtr);
        resultData = _mm256_fnmadd_pd(vecData, matData, resultData);
        _mm256_storeu_pd(resultPtr, resultData);

        matPtr += 4;
        resultPtr += 4;
      }
    }
#endif

    // leftovers that don't fit in a vector.
    for(; r < R; r++)
    {
      *resultPtr = *resultPtr - *matPtr * vec[c];
      resultPtr++;
      matPtr++;
    }
  }
}

/*!
 * Specialization for float, can pack 8 numbers/register
 */
template<>
void CholeskyDenseSolver<float>::setupAVX(float *mat, float *result, float *vec, u64 row)
{
  s64 R = n - row - 1;
  s64 C = row;
  float* matPtr;
  float* resultPtr;
  for(s64 c = 0; c < C; c++)
  {
#ifdef JCQP_USE_AVX2
    __m256 vecData = _mm256_set1_ps(vec[c]); // broadcast to all
#endif
    resultPtr = result;
    matPtr = mat + c * n;
    s64 r = 0;

#ifdef JCQP_USE_AVX2
    s64 end = R - 8 * UNROLL_MATVEC; // gcc fails to make this optimization.
    for(; r < end; r += 8 * UNROLL_MATVEC)
    {
      for(s64 u = 0; u < UNROLL_MATVEC; u++)
      {
        __m256 matData = _mm256_loadu_ps(matPtr);
        __m256 resultData = _mm256_loadu_ps(resultPtr);
        resultData = _mm256_fnmadd_ps(vecData, matData, resultData);
        _mm256_storeu_ps(resultPtr, resultData);

        matPtr += 8;
        resultPtr += 8;
      }
    }
#endif
    // leftovers.
    for(; r < R; r++)
    {
      *resultPtr = *resultPtr - *matPtr * vec[c];
      resultPtr++;
      matPtr++;
    }
  }
}

/*!
 * Reconstruct the original matrix (permuted). Only used for debugging, this is slow
 */
template<typename T>
DenseMatrix<T> CholeskyDenseSolver<T>::getReconstructedPermuted()
{
  DenseMatrix<T> LDebug(n,n);
  for(s64 j = 0; j < n; j++)
  {
    for(s64 i = 0; i < n; i++)
    {
      if(i == j)
      {
        LDebug(i,j) = 1;
      }
      else if (i > j)
      {
        LDebug(i,j) = L(i,j);
      }
    }
  }
  DenseMatrix<T> Reconstructed = LDebug * L.diagonal().asDiagonal() * LDebug.transpose();
  return Reconstructed;
}

/*!
 * public solve function
 */
template<typename T>
void CholeskyDenseSolver<T>::solve(Vector<T> &out)
{
  // only use the AVX solver if we have it available.
#ifdef JCQP_USE_AVX2
  solveAVX(out);
  return;
#endif

  s64 c = 0;
  // first permute:
  for(s64 i = 0; i < n; i++)
  {
    std::swap(out[pivots[i]], out[i]);
  }
  // now temp = b


  // next solve L * y = b in place with temp.  L is column major
  // the inner loop is down the result of the solve.
  for(s64 j = 0; j < n; j++)
  {
    for(s64 i = j + 1; i < n; i++)
    {
      out[i] -= solve1[c++] * out[j]; // i > j
    }
  }

  // next do scaling for inv(D) * y
  // now temp contains y;
  for(s64 i = 0; i < n; i++)
  {
    out[i] /= L(i,i);
  }


  c = 0;
  for(s64 j = n; j --> 0;)
  {
    for(s64 i = 0; i < j; i++)
    {
      out[i] -= solve2[c++] * out[j];
    }
  }

  // finally, untranspose
  for(s64 i = n; i-->0;)
  {
    std::swap(out[pivots[i]], out[i]);
  }
}

#ifdef JCQP_USE_AVX2
template<>
void CholeskyDenseSolver<double>::solveAVX(Vector<double> &out)
{
  // first permute:
  for(s64 i = 0; i < n; i++)
  {
    std::swap(out[pivots[i]], out[i]);
  }

  // O(n^2) loop 1
  double* matPtr = solve1;
  double* lhsPtr;
  for(s64 j = 0; j < n; j++)
  {
    lhsPtr = out.data() + j + 1;
    __m256d rhs = _mm256_set1_pd(out[j]);
    s64 i = j + 1;
    s64 vEnd = n - 4 * UNROLL_MATVEC;
    for(; i < vEnd; i += 4 * UNROLL_MATVEC)
    {
      for(s64 u = 0; u < UNROLL_MATVEC; u++)
      {
        __m256d mat = _mm256_loadu_pd(matPtr);
        __m256d lhs = _mm256_loadu_pd(lhsPtr);
        lhs = _mm256_fnmadd_pd(rhs, mat, lhs);
        _mm256_storeu_pd(lhsPtr, lhs);

        matPtr += 4;
        lhsPtr += 4;
      }
    }

    for(; i < n; i++)
    {
      *lhsPtr = *lhsPtr - *matPtr * out[j];
      lhsPtr++;
      matPtr++;
    }
  }


  // next do scaling for inv(D) * y
  // now temp contains y;
  for(s64 i = 0; i < n; i++)
  {
    out[i] /= L(i,i);
  }

  matPtr = solve2;
  for(s64 j = n; j --> 0;)
  {
    lhsPtr = out.data();
    __m256d rhs = _mm256_set1_pd(out[j]);
    s64 i = 0;
    s64 vEnd = j - 4 * UNROLL_MATVEC;
    for(; i < vEnd; i += 4 * UNROLL_MATVEC)
    {
      for(s64 u = 0; u < UNROLL_MATVEC; u++)
      {
        __m256d mat = _mm256_loadu_pd(matPtr);
        __m256d lhs = _mm256_loadu_pd(lhsPtr);
        lhs = _mm256_fnmadd_pd(rhs, mat, lhs);
        _mm256_storeu_pd(lhsPtr, lhs);

        matPtr += 4;
        lhsPtr += 4;
      }
    }

    for(; i < j; i++)
    {
      *lhsPtr = *lhsPtr - *matPtr * out[j];
      lhsPtr++;
      matPtr++;
    }
  }

  // O(n^2) loop 2
  // finally, untranspose
  for(s64 i = n; i-->0;)
  {
    std::swap(out[pivots[i]], out[i]);
  }
}

template<>
void CholeskyDenseSolver<float>::solveAVX(Vector<float> &out)
{
  // first permute:
  for(s64 i = 0; i < n; i++)
  {
    std::swap(out[pivots[i]], out[i]);
  }

  // O(n^2) loop 1
  float* matPtr = solve1;
  float* lhsPtr;
  for(s64 j = 0; j < n; j++)
  {
    lhsPtr = out.data() + j + 1;
    __m256 rhs = _mm256_set1_ps(out[j]);
    s64 i = j + 1;
    s64 vEnd = n - 8 * UNROLL_MATVEC;
    for(; i < vEnd; i += 8 * UNROLL_MATVEC)
    {
      for(s64 u = 0; u < UNROLL_MATVEC; u++)
      {
        __m256 mat = _mm256_loadu_ps(matPtr);
        __m256 lhs = _mm256_loadu_ps(lhsPtr);
        lhs = _mm256_fnmadd_ps(rhs, mat, lhs);
        _mm256_storeu_ps(lhsPtr, lhs);

        matPtr += 8;
        lhsPtr += 8;
      }
    }

    for(; i < n; i++)
    {
      *lhsPtr = *lhsPtr - *matPtr * out[j];
      lhsPtr++;
      matPtr++;
    }
  }


  // next do scaling for inv(D) * y
  // now temp contains y;
  for(s64 i = 0; i < n; i++)
  {
    out[i] /= L(i,i);
  }

  matPtr = solve2;
  for(s64 j = n; j --> 0;)
  {
    lhsPtr = out.data();
    __m256 rhs = _mm256_set1_ps(out[j]);
    s64 i = 0;
    s64 vEnd = j - 8 * UNROLL_MATVEC;
    for(; i < vEnd; i += 8 * UNROLL_MATVEC)
    {
      for(s64 u = 0; u < UNROLL_MATVEC; u++)
      {
        __m256 mat = _mm256_loadu_ps(matPtr);
        __m256 lhs = _mm256_loadu_ps(lhsPtr);
        lhs = _mm256_fnmadd_ps(rhs, mat, lhs);
        _mm256_storeu_ps(lhsPtr, lhs);

        matPtr += 8;
        lhsPtr += 8;
      }
    }

    for(; i < j; i++)
    {
      *lhsPtr = *lhsPtr - *matPtr * out[j];
      lhsPtr++;
      matPtr++;
    }
  }

  // O(n^2) loop 2
  // finally, untranspose
  for(s64 i = n; i-->0;)
  {
    std::swap(out[pivots[i]], out[i]);
  }
}
#endif


template class CholeskyDenseSolver<double>;
template class CholeskyDenseSolver<float>;
