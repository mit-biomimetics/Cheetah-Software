#include "eigenvalues.h"
#include <eigen3/Eigen/Eigenvalues>

/*!
 * Adjust a matrix so that its eigenvalues are less than a given magnitude.
 * Moved to its on CU because this doesn't change often and takes like 2 minutes to recompile on O3 :P
 */
template<typename T>
DenseMatrix<T> constrainEigenvalueMagnitude(const DenseMatrix<T> &m, T magnitude)
{
  Eigen::EigenSolver<DenseMatrix<T>> eigSolver(m);
  using MatCT = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>;
  MatCT D = eigSolver.eigenvalues().asDiagonal();
  MatCT V = eigSolver.eigenvectors();

  for(s64 i = 0; i < D.rows(); i++) {
    T eigMag = std::abs(D(i,i));
    if(eigMag > magnitude)
      D(i,i) = magnitude * D(i,i) / eigMag;
  }

  MatCT restored = V * D * V.inverse(); // slow, but this matrix is small.
  DenseMatrix<T> result(m.rows(), m.rows());
  for(s64 i = 0; i < m.rows(); i++) {
    for(s64 j = 0; j < m.rows(); j++) {
      result(i,j) = std::real(restored(i,j));
    }
  }

  return result;
}

template DenseMatrix<float> constrainEigenvalueMagnitude(const DenseMatrix<float> &m, float magnitude);
template DenseMatrix<double> constrainEigenvalueMagnitude(const DenseMatrix<double> &m, double magnitude);