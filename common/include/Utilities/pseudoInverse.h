#include <stdio.h>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/SVD>
using namespace std;

/*!
 * Compute the pseudo inverse of a matrix
 * @param matrix : input matrix
 * @param sigmaThreshold : threshold for singular values being zero
 * @param invMatrix : output matrix
 */
template <typename T>
void pseudoInverse(DMat<T> const& matrix, double sigmaThreshold,
                   DMat<T>& invMatrix) {
  if ((1 == matrix.rows()) && (1 == matrix.cols())) {
    invMatrix.resize(1, 1);
    if (matrix.coeff(0, 0) > sigmaThreshold) {
      invMatrix.coeffRef(0, 0) = 1.0 / matrix.coeff(0, 0);
    } else {
      invMatrix.coeffRef(0, 0) = 0.0;
    }
    return;
  }

  Eigen::JacobiSVD<DMat<T>> svd(matrix,
                                Eigen::ComputeThinU | Eigen::ComputeThinV);
  // not sure if we need to svd.sort()... probably not
  int const nrows(svd.singularValues().rows());
  DMat<T> invS;
  invS = DMat<T>::Zero(nrows, nrows);
  for (int ii(0); ii < nrows; ++ii) {
    if (svd.singularValues().coeff(ii) > sigmaThreshold) {
      invS.coeffRef(ii, ii) = 1.0 / svd.singularValues().coeff(ii);
    } else {
      // invS.coeffRef(ii, ii) = 1.0/ sigmaThreshold;
      // printf("sigular value is too small: %f\n",
      // svd.singularValues().coeff(ii));
    }
  }
  invMatrix = svd.matrixV() * invS * svd.matrixU().transpose();
}
