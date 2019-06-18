#ifndef QPSOLVER_TYPES_H
#define QPSOLVER_TYPES_H

#include <stdint.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

using u64 = uint64_t;
using u32 = uint32_t;
using u16 = uint16_t;
using  u8 =  uint8_t;
using s64 = int64_t;
using s32 = int32_t;
using s16 = int16_t;
using  s8 =  int8_t;

template<typename T>
using DenseMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T>
using SparseMatrix = Eigen::SparseMatrix<T>;

template<typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

#endif //QPSOLVER_TYPES_H
