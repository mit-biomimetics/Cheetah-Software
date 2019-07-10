/*! @file orientation_tools.h
 *  @brief Utility functions for 3D rotations
 *
 *  This file contains rotation utilities.  We generally use "coordinate
 * transformations" as opposed to the displacement transformations that are
 * commonly found in graphics.  To describe the orientation of a body, we use a
 * rotation matrix which transforms from world to body coordinates. This is the
 * transpose of the matrix which would rotate the body itself into the correct
 * orientation.
 *
 *  This follows the convention of Roy Featherstone's excellent book, Rigid Body
 * Dynamics Algorithms and the spatial_v2 MATLAB library that comes with it.
 * Note that we don't use the spatial_v2 convention for quaternions!
 */

#ifndef LIBBIOMIMETICS_ORIENTATION_TOOLS_H
#define LIBBIOMIMETICS_ORIENTATION_TOOLS_H

#include "Math/MathUtilities.h"
#include "cppTypes.h"

#include <eigen3/Eigen/Dense>

#include <cmath>
#include <iostream>
#include <type_traits>

namespace ori {

static constexpr double quaternionDerviativeStabilization = 0.1;

enum class CoordinateAxis { X, Y, Z };

/*!
 * Convert radians to degrees
 */
template <typename T>
T rad2deg(T rad) {
  static_assert(std::is_floating_point<T>::value,
                "must use floating point value");
  return rad * T(180) / T(M_PI);
}

/*!
 * Convert degrees to radians
 */
template <typename T>
T deg2rad(T deg) {
  static_assert(std::is_floating_point<T>::value,
                "must use floating point value");
  return deg * T(M_PI) / T(180);
}

/*!
 * Compute rotation matrix for coordinate transformation. Note that
 * coordinateRotation(CoordinateAxis:X, .1) * v will rotate v by -.1 radians -
 * this transforms into a frame rotated by .1 radians!.
 */
template <typename T>
Mat3<T> coordinateRotation(CoordinateAxis axis, T theta) {
  static_assert(std::is_floating_point<T>::value,
                "must use floating point value");
  T s = std::sin(theta);
  T c = std::cos(theta);

  Mat3<T> R;

  if (axis == CoordinateAxis::X) {
    R << 1, 0, 0, 0, c, s, 0, -s, c;
  } else if (axis == CoordinateAxis::Y) {
    R << c, 0, -s, 0, 1, 0, s, 0, c;
  } else if (axis == CoordinateAxis::Z) {
    R << c, s, 0, -s, c, 0, 0, 0, 1;
  }

  return R;
}

/*!
 * Go from rpy to rotation matrix.
 */
template <typename T>
Mat3<typename T::Scalar> rpyToRotMat(const Eigen::MatrixBase<T>& v) {
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 3,
                "must have 3x1 vector");
  Mat3<typename T::Scalar> m = coordinateRotation(CoordinateAxis::X, v[0]) *
                               coordinateRotation(CoordinateAxis::Y, v[1]) *
                               coordinateRotation(CoordinateAxis::Z, v[2]);
  return m;
}

/*!
 * Convert a 3x1 vector to a skew-symmetric 3x3 matrix
 */
template <typename T>
Mat3<typename T::Scalar> vectorToSkewMat(const Eigen::MatrixBase<T>& v) {
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 3,
                "Must have 3x1 matrix");
  Mat3<typename T::Scalar> m;
  m << 0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0;
  return m;
}

/*!
 * Put the skew-symmetric component of 3x3 matrix m into a 3x1 vector
 */
template <typename T>
Vec3<typename T::Scalar> matToSkewVec(const Eigen::MatrixBase<T>& m) {
  static_assert(T::ColsAtCompileTime == 3 && T::RowsAtCompileTime == 3,
                "Must have 3x3 matrix");
  return 0.5 * Vec3<typename T::Scalar>(m(2, 1) - m(1, 2), m(0, 2) - m(2, 0),
                                        (m(1, 0) - m(0, 1)));
}

/*!
 * Convert a coordinate transformation matrix to an orientation quaternion.
 */
template <typename T>
Quat<typename T::Scalar> rotationMatrixToQuaternion(
    const Eigen::MatrixBase<T>& r1) {
  static_assert(T::ColsAtCompileTime == 3 && T::RowsAtCompileTime == 3,
                "Must have 3x3 matrix");
  Quat<typename T::Scalar> q;
  Mat3<typename T::Scalar> r = r1.transpose();
  typename T::Scalar tr = r.trace();
  if (tr > 0.0) {
    typename T::Scalar S = sqrt(tr + 1.0) * 2.0;
    q(0) = 0.25 * S;
    q(1) = (r(2, 1) - r(1, 2)) / S;
    q(2) = (r(0, 2) - r(2, 0)) / S;
    q(3) = (r(1, 0) - r(0, 1)) / S;
  } else if ((r(0, 0) > r(1, 1)) && (r(0, 0) > r(2, 2))) {
    typename T::Scalar S = sqrt(1.0 + r(0, 0) - r(1, 1) - r(2, 2)) * 2.0;
    q(0) = (r(2, 1) - r(1, 2)) / S;
    q(1) = 0.25 * S;
    q(2) = (r(0, 1) + r(1, 0)) / S;
    q(3) = (r(0, 2) + r(2, 0)) / S;
  } else if (r(1, 1) > r(2, 2)) {
    typename T::Scalar S = sqrt(1.0 + r(1, 1) - r(0, 0) - r(2, 2)) * 2.0;
    q(0) = (r(0, 2) - r(2, 0)) / S;
    q(1) = (r(0, 1) + r(1, 0)) / S;
    q(2) = 0.25 * S;
    q(3) = (r(1, 2) + r(2, 1)) / S;
  } else {
    typename T::Scalar S = sqrt(1.0 + r(2, 2) - r(0, 0) - r(1, 1)) * 2.0;
    q(0) = (r(1, 0) - r(0, 1)) / S;
    q(1) = (r(0, 2) + r(2, 0)) / S;
    q(2) = (r(1, 2) + r(2, 1)) / S;
    q(3) = 0.25 * S;
  }
  return q;
}

/*!
 * Convert a quaternion to a rotation matrix.  This matrix represents a
 * coordinate transformation into the frame which has the orientation specified
 * by the quaternion
 */
template <typename T>
Mat3<typename T::Scalar> quaternionToRotationMatrix(
    const Eigen::MatrixBase<T>& q) {
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 4,
                "Must have 4x1 quat");
  typename T::Scalar e0 = q(0);
  typename T::Scalar e1 = q(1);
  typename T::Scalar e2 = q(2);
  typename T::Scalar e3 = q(3);

  Mat3<typename T::Scalar> R;

  R << 1 - 2 * (e2 * e2 + e3 * e3), 2 * (e1 * e2 - e0 * e3),
      2 * (e1 * e3 + e0 * e2), 2 * (e1 * e2 + e0 * e3),
      1 - 2 * (e1 * e1 + e3 * e3), 2 * (e2 * e3 - e0 * e1),
      2 * (e1 * e3 - e0 * e2), 2 * (e2 * e3 + e0 * e1),
      1 - 2 * (e1 * e1 + e2 * e2);
  R.transposeInPlace();
  return R;
}

/*!
 * Convert a quaternion to RPY.  Uses ZYX order (yaw-pitch-roll), but returns
 * angles in (roll, pitch, yaw).
 */
template <typename T>
Vec3<typename T::Scalar> quatToRPY(const Eigen::MatrixBase<T>& q) {
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 4,
                "Must have 4x1 quat");
  Vec3<typename T::Scalar> rpy;
  typename T::Scalar as = std::min(-2. * (q[1] * q[3] - q[0] * q[2]), .99999);
  rpy(2) =
      std::atan2(2 * (q[1] * q[2] + q[0] * q[3]),
                 square(q[0]) + square(q[1]) - square(q[2]) - square(q[3]));
  rpy(1) = std::asin(as);
  rpy(0) =
      std::atan2(2 * (q[2] * q[3] + q[0] * q[1]),
                 square(q[0]) - square(q[1]) - square(q[2]) + square(q[3]));
  return rpy;
}

template <typename T>
Quat<typename T::Scalar> rpyToQuat(const Eigen::MatrixBase<T>& rpy) {
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 3,
                "Must have 3x1 vec");
  Mat3<typename T::Scalar> R = rpyToRotMat(rpy);
  Quat<typename T::Scalar> q = rotationMatrixToQuaternion(R);
  return q;
}

/*!
 * Convert a quaternion to so3.
 */
template <typename T>
Vec3<typename T::Scalar> quatToso3(const Eigen::MatrixBase<T>& q) {
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 4,
                "Must have 4x1 quat");
  Vec3<typename T::Scalar> so3;
  typename T::Scalar theta = 2. * std::acos(q[0]);
  so3[0] = theta * q[1] / std::sin(theta / 2.);
  so3[1] = theta * q[2] / std::sin(theta / 2.);
  so3[2] = theta * q[3] / std::sin(theta / 2.);
  return so3;
}

template <typename T>
Vec3<typename T::Scalar> rotationMatrixToRPY(const Eigen::MatrixBase<T>& R) {
  static_assert(T::ColsAtCompileTime == 3 && T::RowsAtCompileTime == 3,
                "Must have 3x3 matrix");
  Quat<typename T::Scalar> q = rotationMatrixToQuaternion(R);
  Vec3<typename T::Scalar> rpy = quatToRPY(q);
  return rpy;
}
/*!
 * Quaternion derivative calculation, like rqd(q, omega) in MATLAB.
 * the omega is expressed in body frame
 * @param q
 * @param omega
 * @return
 */
template <typename T, typename T2>
Quat<typename T::Scalar> quatDerivative(const Eigen::MatrixBase<T>& q,
                                        const Eigen::MatrixBase<T2>& omega) {
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 4,
                "Must have 4x1 quat");
  static_assert(T2::ColsAtCompileTime == 1 && T2::RowsAtCompileTime == 3,
                "Must have 3x1 omega");
  // first case in rqd
  Mat4<typename T::Scalar> Q;
  Q << q[0], -q[1], -q[2], -q[3], q[1], q[0], -q[3], q[2], q[2], q[3], q[0],
      -q[1], q[3], -q[2], q[1], q[0];

  Quat<typename T::Scalar> qq(
      quaternionDerviativeStabilization * omega.norm() * (1 - q.norm()),
      omega[0], omega[1], omega[2]);
  Quat<typename T::Scalar> dq = 0.5 * Q * qq;
  return dq;
}

/*!
 * Take the product of two quaternions
 */
template <typename T>
Quat<typename T::Scalar> quatProduct(const Eigen::MatrixBase<T>& q1,
                                     const Eigen::MatrixBase<T>& q2) {
  typename T::Scalar r1 = q1[0];
  typename T::Scalar r2 = q2[0];
  Vec3<typename T::Scalar> v1(q1[1], q1[2], q1[3]);
  Vec3<typename T::Scalar> v2(q2[1], q2[2], q2[3]);

  typename T::Scalar r = r1 * r2 - v1.dot(v2);
  Vec3<typename T::Scalar> v = r1 * v2 + r2 * v1 + v1.cross(v2);
  Quat<typename T::Scalar> q(r, v[0], v[1], v[2]);
  return q;
}

/*!
 * Compute new quaternion given:
 * @param quat The old quaternion
 * @param omega The angular velocity (IN INERTIAL COORDINATES!)
 * @param dt The timestep
 * @return
 */
template <typename T, typename T2, typename T3>
Quat<typename T::Scalar> integrateQuat(const Eigen::MatrixBase<T>& quat,
                                       const Eigen::MatrixBase<T2>& omega,
                                       T3 dt) {
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 4,
                "Must have 4x1 quat");
  static_assert(T2::ColsAtCompileTime == 1 && T2::RowsAtCompileTime == 3,
                "Must have 3x1 omega");
  Vec3<typename T::Scalar> axis;
  typename T::Scalar ang = omega.norm();
  if (ang > 0) {
    axis = omega / ang;
  } else {
    axis = Vec3<typename T::Scalar>(1, 0, 0);
  }

  ang *= dt;
  Vec3<typename T::Scalar> ee = std::sin(ang / 2) * axis;
  Quat<typename T::Scalar> quatD(std::cos(ang / 2), ee[0], ee[1], ee[2]);

  Quat<typename T::Scalar> quatNew = quatProduct(quatD, quat);
  quatNew = quatNew / quatNew.norm();
  return quatNew;
}

/*!
 * Compute new quaternion given:
 * @param quat The old quaternion
 * @param omega The angular velocity (IN INERTIAL COORDINATES!)
 * @param dt The timestep
 * @return
 */
template <typename T, typename T2, typename T3>
Quat<typename T::Scalar> integrateQuatImplicit(
    const Eigen::MatrixBase<T>& quat, const Eigen::MatrixBase<T2>& omega,
    T3 dt) {
  static_assert(T::ColsAtCompileTime == 1 && T::RowsAtCompileTime == 4,
                "Must have 4x1 quat");
  static_assert(T2::ColsAtCompileTime == 1 && T2::RowsAtCompileTime == 3,
                "Must have 3x1 omega");
  Vec3<typename T::Scalar> axis;
  typename T::Scalar ang = omega.norm();
  if (ang > 0) {
    axis = omega / ang;
  } else {
    axis = Vec3<typename T::Scalar>(1, 0, 0);
  }

  ang *= dt;
  Vec3<typename T::Scalar> ee = std::sin(ang / 2) * axis;
  Quat<typename T::Scalar> quatD(std::cos(ang / 2), ee[0], ee[1], ee[2]);

  Quat<typename T::Scalar> quatNew = quatProduct(quat, quatD);
  quatNew = quatNew / quatNew.norm();
  return quatNew;
}

template <typename T>
void quaternionToso3(const Quat<T> quat, Vec3<T>& so3) {
  so3[0] = quat[1];
  so3[1] = quat[2];
  so3[2] = quat[3];

  T theta =
      2.0 * asin(sqrt(so3[0] * so3[0] + so3[1] * so3[1] + so3[2] * so3[2]));

  if (fabs(theta) < 0.0000001) {
    so3.setZero();
    return;
  }
  so3 /= sin(theta / 2.0);
  so3 *= theta;
}
template <typename T>
Quat<T> so3ToQuat(Vec3<T>& so3) {
  Quat<T> quat;

  T theta = sqrt(so3[0] * so3[0] + so3[1] * so3[1] + so3[2] * so3[2]);

  if (fabs(theta) < 1.e-6) {
    quat.setZero();
    quat[0] = 1.;
    return quat;
  }
  quat[0] = cos(theta / 2.);
  quat[1] = so3[0] / theta * sin(theta / 2.);
  quat[2] = so3[1] / theta * sin(theta / 2.);
  quat[3] = so3[2] / theta * sin(theta / 2.);
  return quat;
}
}  // namespace ori

#endif  // LIBBIOMIMETICS_ORIENTATION_TOOLS_H
