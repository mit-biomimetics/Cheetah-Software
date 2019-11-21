/*! @file test_orientation_tools.cpp
 *  @brief Test orientation functions
 *
 * Test the orientation related functions in orientation_tools.h
 * Does not check any spatial stuff
 */

#include "Math/orientation_tools.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace ori;

/*!
 * Check radians to degrees
 */
TEST(Orientation, rad2deg) {
  EXPECT_EQ(M_PI / 4., deg2rad(45.));
  EXPECT_EQ(45, rad2deg(M_PI / 4.));
}


/*!
 * Check the "almostEqual" function for comparing two eigen quantities
 */
TEST(Orientation, almostEqual) {
  // check matrix "almostEqual" function
  Mat3<double> a, b;
  a << 1, 2, 3, 4, 5, 6, 7, 8, 9.2;
  b << 1, 2, 3, 4, 5, 6, 7, 8, 9.1;
  EXPECT_EQ(true, almostEqual(a, b, .3));
  EXPECT_EQ(false, almostEqual(a, b, .01));

  DVec<double> qdd(4);
  qdd << 1, 2, 3, 4;

  DVec<double> qdd2 = qdd;

  EXPECT_TRUE(almostEqual(qdd, qdd2, 1e-6));
  qdd2(1) += .2;
  EXPECT_FALSE(almostEqual(qdd, qdd2, 1e-6));

  DMat<double> testDynamicMat(3, 3);
  DMat<double> testDynamicMat2(3, 3);
  testDynamicMat << 1, 2, 3, 4, 5, 6, 7, 8, 9.2;
  testDynamicMat2 << 1, 2, 3, 4, 5, 6, 7, 8, 9.1;
  EXPECT_EQ(true, almostEqual(testDynamicMat, testDynamicMat2, .3));
  EXPECT_EQ(false, almostEqual(testDynamicMat, testDynamicMat2, .01));
}

TEST(Orientation, crossMatrix) {
  Vec3<double> a(15,2,6);
  Vec3<double> b(2,-1,3);
  Vec3<double> axb1 = a.cross(b);
  Mat3<double> ax = crossMatrix(a);
  Vec3<double> axb2 = ax * b;
  EXPECT_TRUE(almostEqual(axb1, axb2, .0001));
}

/*!
 * Check the coordinateRotation algorithm which generates rotaiton matrices for
 * axis rotations
 */
TEST(Orientation, coordinateRotation) {
  // check rotation matrices
  double s = std::sin(.2);
  double c = std::cos(.2);
  Mat3<double> R_ref_x, R_ref_y, R_ref_z;
  R_ref_x << 1, 0, 0, 0, c, s, 0, -s, c;
  R_ref_y << c, 0, -s, 0, 1, 0, s, 0, c;
  R_ref_z << c, s, 0, -s, c, 0, 0, 0, 1;
  EXPECT_EQ(R_ref_x, coordinateRotation(CoordinateAxis::X, .2));
  EXPECT_EQ(R_ref_y, coordinateRotation(CoordinateAxis::Y, .2));
  EXPECT_EQ(R_ref_z, coordinateRotation(CoordinateAxis::Z, .2));
}

/*!
 * Check the skew functions, which go between vectors and skew-symmetric
 * matrices
 */
TEST(Orientation, skew) {
  // check skew vec->mat and mat->vec
  Mat3<double> A, B;
  A << 8, 1, 6, 3, 5, 7, 4, 9, 2;
  B << 0, -3, 2, 3, 0, -1, -2, 1, 0;
  Vec3<double> v(1, 1, 1);
  Vec3<double> w(1, 2, 3);

  EXPECT_EQ(matToSkewVec(A), v);
  EXPECT_EQ(B, vectorToSkewMat(w));
  EXPECT_EQ(w, matToSkewVec(vectorToSkewMat(w)));
  EXPECT_EQ(v, matToSkewVec(vectorToSkewMat(v)));
}

/*!
 * Check the quaternion to rotation matrix and rotation matrix to quaternion
 * functions
 */
TEST(Orientation, quatToRotm) {
  // check R -> q and q -> R
  Quat<double> q(.9672, -0.0672, -0.1653, -0.1808);
  Mat3<double> R;
  R << .8799, .3720, -.2955, -.3276, .9256, .1898, .3441, -.0702, .9363;
  Quat<double> q2 = rotationMatrixToQuaternion(R.transpose());
  Mat3<double> R2 = quaternionToRotationMatrix(q).transpose();
  EXPECT_TRUE(almostEqual(q2, q, .0001));
  EXPECT_TRUE(almostEqual(R2, R, .0001));
}

/*!
 * Test the quaternion to roll-pitch-yaw conversion
 */
TEST(Orientation, quatToRpy) {
  // check q -> rpy
  Quat<double> q(0.9672, -0.0672, -0.1653, -0.1808);
  Vec3<double> rpy1(-0.0748, -0.3513, -0.3564);
  EXPECT_TRUE(almostEqual(rpy1, quatToRPY(q), .1));
}

/*!
 * Test the quaternion to roll-pitch-yaw conversion direction
 */
TEST(Orientation, quatToRpy2) {
  // check q -> rpy
  RotMat<double> R = coordinateRotation(CoordinateAxis::Y, .5);
  Quat<double> q = rotationMatrixToQuaternion(R);
  Vec3<double> rpy = quatToRPY(q);
  Vec3<double> rpy_ref(0, .5, 0);
  EXPECT_TRUE(almostEqual(rpy, rpy_ref, .00001));
}

/*!
 * Check that rpy is in roll-pitch-yaw order
 */
TEST(Orientation, quatToRpy3) {
  // check q -> rpy
  RotMat<double> R = coordinateRotation(CoordinateAxis::Z, .5);
  Quat<double> q = rotationMatrixToQuaternion(R);
  Vec3<double> rpy = quatToRPY(q);
  Vec3<double> rpy_ref(0, 0, 0.5);
  EXPECT_TRUE(almostEqual(rpy, rpy_ref, .00001));
}

/*!
 * Check the quaternion derivative function.
 * This isn't used and probably should be removed
 */
TEST(Orienation, quaternionDerivative) {
  Quat<double> ref(-10.8376, -0.6752, -2.5128, -1.3504);
  Quat<double> q =
      quatDerivative(Quat<double>(1, 2, 3, 4), Vec3<double>(1, 2, 3));
  EXPECT_TRUE(almostEqual(q, ref, .0005));
}

/*!
 * Check the quaternion product
 */
TEST(Orientation, quaternionProduct) {
  Quat<double> q1(1, 2, 3, 4);
  Quat<double> q2(5, 6, 7, 8);
  Quat<double> ref(-60, 12, 30, 24);
  EXPECT_TRUE(almostEqual(ref, quatProduct(q1, q2), .0001));
}

TEST(Orientation, quaternionProductDirection) {
  RotMat<double> R1 = coordinateRotation(CoordinateAxis::X, 7.23) *
                      coordinateRotation(CoordinateAxis::Y, -2.343) *
                      coordinateRotation(CoordinateAxis::Z, 1.2324);

  RotMat<double> R2 = coordinateRotation(CoordinateAxis::Z, .3231) *
                      coordinateRotation(CoordinateAxis::X, -4.2332) *
                      coordinateRotation(CoordinateAxis::Y, -3.3213);

  RotMat<double> R12_ref = R2 * R1;

  Quat<double> Q1 = rotationMatrixToQuaternion(R1);
  Quat<double> Q2 = rotationMatrixToQuaternion(R2);
  Quat<double> Q12 = quatProduct(Q1, Q2);

  RotMat<double> R12 = quaternionToRotationMatrix(Q12);

  EXPECT_TRUE(almostEqual(R12, R12_ref, .0001));
}

/*!
 * Check that the quaternion integration goes in the correct direction
 */
TEST(Orientation, quaternionIntegration) {
  RotMat<double> eye = RotMat<double>::Identity();

  Vec3<double> omegaX(1, 0, 0);
  Vec3<double> omegaY(0, -1, 0);
  Vec3<double> omegaZ(0, 0, 1);

  RotMat<double> rot_x_ref = coordinateRotation(CoordinateAxis::X, .1);
  RotMat<double> rot_y_ref = coordinateRotation(CoordinateAxis::Y, -.1);
  RotMat<double> rot_z_ref = coordinateRotation(CoordinateAxis::Z, .1);

  Quat<double> eyeQ = rotationMatrixToQuaternion(eye);

  Quat<double> rot_x_quat = integrateQuat(eyeQ, omegaX, .1);
  Quat<double> rot_y_quat = integrateQuat(eyeQ, omegaY, .1);
  Quat<double> rot_z_quat = integrateQuat(eyeQ, omegaZ, .1);

  EXPECT_TRUE(
      almostEqual(quaternionToRotationMatrix(rot_x_quat), rot_x_ref, .0001));
  EXPECT_TRUE(
      almostEqual(quaternionToRotationMatrix(rot_y_quat), rot_y_ref, .0001));
  EXPECT_TRUE(
      almostEqual(quaternionToRotationMatrix(rot_z_quat), rot_z_ref, .0001));
}

/*!
 * Check rpy to rotation matrix
 */
TEST(Orientation, rpyToRotMat) {
  Vec3<double> rpy(deg2rad(20.), deg2rad(34.), deg2rad(160.));
  Mat3<double> R = rpyToRotMat(rpy);
  Quat<double> q = rotationMatrixToQuaternion(R);
  Vec3<double> rpyAgain = quatToRPY(q);
  EXPECT_TRUE(almostEqual(rpyAgain, rpy, .001));
  std::cout << "1: " << rpy.transpose() << "\n2: " << rpyAgain.transpose()
            << "\n";
}

/*!
 * Check all 6 possible conversions between rotation matrix, quaternion, and rpy
 */
TEST(Orientation, allOrientationConversions) {
  Vec3<double> refRPY(deg2rad(20.), deg2rad(34.), deg2rad(160.));
  // do all rpy -> conversions
  Quat<double> quatFromRPY = rpyToQuat(refRPY);
  Mat3<double> rFromRPY = rpyToRotMat(refRPY);

  // do all quat -> conversions
  Mat3<double> rFromQuat = quaternionToRotationMatrix(quatFromRPY);
  Vec3<double> rpyFromQuat = quatToRPY(quatFromRPY);

  // do all r -> conversions
  Vec3<double> rpyFromR = rotationMatrixToRPY(rFromRPY);
  Quat<double> quatFromR = rotationMatrixToQuaternion(rFromQuat);

  EXPECT_TRUE(almostEqual(refRPY, rpyFromQuat, .0001));
  EXPECT_TRUE(almostEqual(refRPY, rpyFromR, .0001));

  EXPECT_TRUE(almostEqual(rFromQuat, rFromRPY, .0001));

  EXPECT_TRUE(almostEqual(quatFromR, quatFromRPY, .0001));
}