/*! @file cTypes.h
 *  @brief Common types that are only valid in C++
 *
 *  This file contains types which are only used in C++ code.  This includes
 * Eigen types, template types, aliases, ...
 */

#ifndef PROJECT_CPPTYPES_H
#define PROJECT_CPPTYPES_H

#include <vector>
#include "cTypes.h"
#include <eigen3/Eigen/Dense>

// Rotation Matrix
template <typename T>
using RotMat = typename Eigen::Matrix<T, 3, 3>;

// 2x1 Vector
template <typename T>
using Vec2 = typename Eigen::Matrix<T, 2, 1>;

// 3x1 Vector
template <typename T>
using Vec3 = typename Eigen::Matrix<T, 3, 1>;

// 4x1 Vector
template <typename T>
using Vec4 = typename Eigen::Matrix<T, 4, 1>;

// 6x1 Vector
template <typename T>
using Vec6 = Eigen::Matrix<T, 6, 1>;

// 10x1 Vector
template <typename T>
using Vec10 = Eigen::Matrix<T, 10, 1>;

// 12x1 Vector
template <typename T>
using Vec12 = Eigen::Matrix<T, 12, 1>;

// 18x1 Vector
template <typename T>
using Vec18 = Eigen::Matrix<T, 18, 1>;

// 28x1 vector
template <typename T>
using Vec28 = Eigen::Matrix<T, 28, 1>;

// 3x3 Matrix
template <typename T>
using Mat3 = typename Eigen::Matrix<T, 3, 3>;

// 4x1 Vector
template <typename T>
using Quat = typename Eigen::Matrix<T, 4, 1>;

// Spatial Vector (6x1, all subspaces)
template <typename T>
using SVec = typename Eigen::Matrix<T, 6, 1>;

// Spatial Transform (6x6)
template <typename T>
using SXform = typename Eigen::Matrix<T, 6, 6>;

// 6x6 Matrix
template <typename T>
using Mat6 = typename Eigen::Matrix<T, 6, 6>;

// 12x12 Matrix
template <typename T>
using Mat12 = typename Eigen::Matrix<T, 12, 12>;

// 18x18 Matrix
template <typename T>
using Mat18 = Eigen::Matrix<T, 18, 18>;

// 28x28 Matrix
template <typename T>
using Mat28 = Eigen::Matrix<T, 28, 28>;

// 3x4 Matrix
template <typename T>
using Mat34 = Eigen::Matrix<T, 3, 4>;

// 3x4 Matrix
template <typename T>
using Mat23 = Eigen::Matrix<T, 2, 3>;

// 4x4 Matrix
template <typename T>
using Mat4 = typename Eigen::Matrix<T, 4, 4>;

// 10x1 Vector
template <typename T>
using MassProperties = typename Eigen::Matrix<T, 10, 1>;

// Dynamically sized vector
template <typename T>
using DVec = typename Eigen::Matrix<T, Eigen::Dynamic, 1>;

// Dynamically sized matrix
template <typename T>
using DMat = typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

// Dynamically sized matrix with spatial vector columns
template <typename T>
using D6Mat = typename Eigen::Matrix<T, 6, Eigen::Dynamic>;

// Dynamically sized matrix with cartesian vector columns
template <typename T>
using D3Mat = typename Eigen::Matrix<T, 3, Eigen::Dynamic>;

// std::vector (a list) of Eigen things
template <typename T>
using vectorAligned = typename std::vector<T, Eigen::aligned_allocator<T>>;

enum class RobotType { CHEETAH_3, MINI_CHEETAH };

#endif  // PROJECT_CPPTYPES_H
