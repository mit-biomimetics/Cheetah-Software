#ifndef QPSOLVER_EIGENVALUES_H
#define QPSOLVER_EIGENVALUES_H

#include "types.h"

template<typename T>
DenseMatrix<T> constrainEigenvalueMagnitude(const DenseMatrix<T> &mat, T magnitude);

#endif //QPSOLVER_EIGENVALUES_H
