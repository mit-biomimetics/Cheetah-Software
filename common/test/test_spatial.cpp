/*! @file test_spatial.cpp
 *  @brief Test spatial vector/transform maniuplation functions
 *
 * Test the spatial utilities
 * Doesn't test the algorithms.
 */

#include "Dynamics/SpatialInertia.h"
#include "Dynamics/spatial.h"
#include "Math/orientation_tools.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace ori;
using namespace spatial;

/*!
 * Check the spatialRotation function, which generates spatial transforms for
 * coordinate axis rotations
 */
TEST(Spatial, axisRotation) {
  // test X which rotates around an axis
  SXform<double> X1, X2;
  X1 << 0.8384, 0.4580, -0.2955, 0, 0, 0, -0.4183, 0.8882, 0.1898, 0, 0, 0,
      0.3494, -0.0355, 0.9363, 0, 0, 0, 0, 0, 0, 0.8384, 0.4580, -0.2955, 0, 0,
      0, -0.4183, 0.8882, 0.1898, 0, 0, 0, 0.3494, -0.0355, 0.9363;
  X2 = spatialRotation(CoordinateAxis::X, .2) *
       spatialRotation(CoordinateAxis::Y, .3) *
       spatialRotation(CoordinateAxis::Z, .5);

  EXPECT_TRUE(almostEqual(X1, X2, .001));
}

/*!
 * Check the motionCrossMatrix function, which behaves like crm() in MATLAB
 */
TEST(Spatial, crm) {
  // test motion cross product matrix
  SVec<double> v1, v2, v3, v4, v5, v6;
  v1 << 1, 2, 3, 4, 5, 6;
  v2 << 6, 5, 4, 3, 2, 1;
  v3 = motionCrossMatrix(v1) * v2;
  v4 << -7, 14, -7, -14, 28, -14;
  v6 = motionCrossMatrix(v1) * v1;
  v5 << 0, 0, 0, 0, 0, 0;
  EXPECT_TRUE(almostEqual(v3, v4, .001));
  EXPECT_TRUE(almostEqual(v6, v5, .0001));
}

/*!
 * Check the forceCrossMatrix function, which behaves like crf() in MATLAB
 * Uses crm(-v).transpose() = crf(v) identity
 */
TEST(Spatial, crf) {
  // test force cross product matrix
  SVec<double> v1, v2;
  v1 << 1, 2, 3, 4, 5, 6;
  v2 << 6, 5, 4, 3, 2, 1;
  Mat6<double> m1 = motionCrossMatrix(-v1);
  m1.transposeInPlace();
  Mat6<double> m2 = forceCrossMatrix(v1);
  EXPECT_TRUE(almostEqual(m1, m2, .0001));
}

/*!
 * Test motionCrossProduct, a function to compute crm(v1)*v2 with fewer
 * multiplies
 */
TEST(Spatial, crm_prod) {
  // test motion cross product
  SVec<double> v1, v2, v3, v4;
  v1 << 1, 2, 3, 4, 5, 6;
  v2 << 6, 5, 4, 3, 2, 1;
  v3 = motionCrossMatrix(v1) * v2;
  v4 = motionCrossProduct(v1, v2);
  EXPECT_TRUE(almostEqual(v3, v4, .0001));
}

/*!
 * Test forceCrossProduct, a function to compute crf(v1)*v2 with fewer
 * multiplies
 */
TEST(Spatial, crf_prod) {
  // test force cross product
  SVec<double> v1, v2, v3, v4;
  v1 << 1, 2, 3, 4, 5, 6;
  v2 << 6, 5, 4, 3, 2, 1;
  v3 = forceCrossMatrix(v1) * v2;
  v4 = forceCrossProduct(v1, v2);
  EXPECT_TRUE(almostEqual(v3, v4, .0001));
}

/*!
 * Test spatial inertia
 * Checks that it is built correctly from mass, CoM and rotational inertia (like
 * mcI) Also checks the 4x4 Pseudo Inertia
 */
TEST(Spatial, inertia) {
  // test spatial inertia, mcI, and Pat's PseudoInertia
  Mat3<double> I;
  I << 1, 2, 3, 2, 4, 5, 3, 5, 6;
  Vec3<double> com(10, 11, 12);
  SpatialInertia<double> IS(42, com, I);
  Mat6<double> ref;
  Mat4<double> pref;
  pref << 4204.5, 4618, 5037, 420, 4618, 5083.5, 5539, 462, 5037, 5539, 6047.5,
      504, 420, 462, 504, 42;

  ref << 11131, -4618, -5037, 0, -504, 462, -4618, 10252, -5539, 504, 0, -420,
      -5037, -5539, 9288, -462, 420, 0, 0, 504, -462, 42, 0, 0, -504, 0, 420, 0,
      42, 0, 462, -420, 0, 0, 0, 42;

  SpatialInertia<double> IS2(IS.getPseudoInertia());

  EXPECT_TRUE(almostEqual(ref, IS.getMatrix(), .001));
  EXPECT_EQ(42, IS.getMass());
  EXPECT_TRUE(almostEqual(I, IS.getInertiaTensor(), .00001));
  EXPECT_TRUE(almostEqual(com, IS.getCOM(), .00001));
  EXPECT_TRUE(almostEqual(IS.getPseudoInertia(), pref, .1));
  EXPECT_TRUE(almostEqual(IS.getMatrix(), IS2.getMatrix(), .00001));
}

/*!
 * Check the flipAlongAxis method of spatial inertias
 */
TEST(Spatial, inertia_flips) {
  // test flipping inertias around axes
  MassProperties<double> a, aref;
  a << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;
  aref << 1, 2, -3, 4, 5, 6, 7, -8, 9, -10;
  EXPECT_TRUE(almostEqual(
      SpatialInertia<double>(aref).getMatrix(),
      SpatialInertia<double>(a).flipAlongAxis(CoordinateAxis::Y).getMatrix(),
      .0001));
}

/*!
 * Check the homogeneous transformation (pluho)
 * Also check the spatial transformation functions (plux)
 */
TEST(Spatial, pluho_and_plux) {
  // test homogeneous transformations and spatial transformations
  RotMat<double> R = coordinateRotation(CoordinateAxis::X, 1.0) *
                     coordinateRotation(CoordinateAxis::Y, 2.0) *
                     coordinateRotation(CoordinateAxis::Z, 3.0);
  Vec3<double> r(4, 5, 6);
  Mat6<double> X = createSXform(R, r);
  Mat6<double> Xref, X2;
  Mat4<double> H = sxformToHomogeneous(X);
  Mat4<double> Href;
  X2 = homogeneousToSXform(H);
  Xref << 0.4120, -0.0587, -0.9093, 0, 0, 0, -0.8337, -0.4269, -0.3502, 0, 0, 0,
      -0.3676, 0.9024, -0.2248, 0, 0, 0, -4.1941, 6.1091, -2.2948, 0.4120,
      -0.0587, -0.9093, 0.8106, -3.6017, 2.4610, -0.8337, -0.4269, -0.3502,
      -6.5385, -1.3064, 5.4477, -0.3676, 0.9024, -0.2248;
  Href << 0.4120, -0.0587, -0.9093, 4.1015, -0.8337, -0.4269, -0.3502, 7.5706,
      -0.3676, 0.9024, -0.2248, -1.6923, 0, 0, 0, 1.0000;
  EXPECT_TRUE(almostEqual(Xref, X, .001));
  EXPECT_TRUE(almostEqual(R, rotationFromSXform(X), .001));
  EXPECT_TRUE(almostEqual(r, translationFromSXform(X), .001));
  EXPECT_TRUE(almostEqual(H, Href, .001));
  EXPECT_TRUE(almostEqual(X, X2, .00001));
}

/*!
 * Check invertSXform, which computes inverse(X) quickly if X is a plucker
 * coordinate transform
 */
TEST(Spatial, invert_sxform) {
  RotMat<double> R = coordinateRotation(CoordinateAxis::X, 1.0) *
                     coordinateRotation(CoordinateAxis::Y, 2.0) *
                     coordinateRotation(CoordinateAxis::Z, 3.0);
  Vec3<double> r(4, 5, 6);
  Mat6<double> X = createSXform(R, r);
  Mat6<double> Xi_ref = X.inverse();
  Mat6<double> Xi = invertSXform(X);
  EXPECT_TRUE(almostEqual(Xi_ref, Xi, .00001));
}

/*!
 * Test the jointXform and jointMotionSubspace functions, which are similar to
 * jcalc
 */
TEST(Spatial, jcalc) {
  Mat6<double> Xr, Xp, Xr_ref, Xp_ref;
  SVec<double> phi_r, phi_p, phi_r_ref, phi_p_ref;

  Xr = jointXform(JointType::Revolute, CoordinateAxis::Y, std::asin(.707));
  Xp = jointXform(JointType::Prismatic, CoordinateAxis::Z, .2);
  phi_r = jointMotionSubspace<double>(JointType::Revolute, CoordinateAxis::Y);
  phi_p = jointMotionSubspace<double>(JointType::Prismatic, CoordinateAxis::Z);

  Xr_ref << 0.7072, 0, -0.7070, 0, 0, 0, 0, 1.0000, 0, 0, 0, 0, 0.7070, 0,
      0.7072, 0, 0, 0, 0, 0, 0, 0.7072, 0, -0.7070, 0, 0, 0, 0, 1.0000, 0, 0, 0,
      0, 0.7070, 0, 0.7072;

  phi_r_ref << 0, 1, 0, 0, 0, 0;

  Xp_ref << 1.0000, 0, 0, 0, 0, 0, 0, 1.0000, 0, 0, 0, 0, 0, 0, 1.0000, 0, 0, 0,
      0, 0.2000, 0, 1.0000, 0, 0, -0.2000, 0, 0, 0, 1.0000, 0, 0, 0, 0, 0, 0,
      1.0000;

  phi_p_ref << 0, 0, 0, 0, 0, 1;

  EXPECT_TRUE(almostEqual(Xr, Xr_ref, .001));
  EXPECT_TRUE(almostEqual(Xp, Xp_ref, .001));
  EXPECT_TRUE(almostEqual(phi_r, phi_r_ref, .000001));
  EXPECT_TRUE(almostEqual(phi_p, phi_p_ref, .000001));
}

/*!
 * Test functions for converting between spatial inertias and the minimal
 * representation "mass property vector", which contains 10 elements.
 */
TEST(Spatial, mass_properties) {
  MassProperties<double> a;
  a << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;
  Mat6<double> I, I_ref;
  I_ref << 5, 10, 9, 0, -4, 3, 10, 6, 8, 4, 0, -2, 9, 8, 7, -3, 2, 0, 0, 4, -3,
      1, 0, 0, -4, 0, 2, 0, 1, 0, 3, -2, 0, 0, 0, 1;
  SpatialInertia<double> IS(a);
  EXPECT_TRUE(almostEqual(I_ref, IS.getMatrix(), .001));
  EXPECT_TRUE(almostEqual(a, IS.asMassPropertyVector(), .001));
}

/*!
 * Test utility function for computing the inertia of a uniform rectangular box
 */
TEST(Spatial, box_inertia) {
  // test inertia of uniformly distributed box
  Mat3<double> I_ref;
  I_ref << 2.0833, 0, 0, 0, 1.66667, 0, 0, 0, 1.083333;
  Mat3<double> I_calc = rotInertiaOfBox(1., Vec3<double>(2, 3, 4));
  EXPECT_TRUE(almostEqual(I_ref, I_calc, .001));
}

/*!
 * Test utility function for converting spatial velocity to linear velocity at a
 * point
 */
TEST(Spatial, velocityConver) {
  Vec3<double> vRef(-2, 17, 0);
  SVec<double> vs;
  vs << 1, 2, 3, 4, 5, 6;
  Vec3<double> p(7, 8, 9);
  Vec3<double> v = spatialToLinearVelocity(vs, p);
  EXPECT_TRUE(almostEqual(v, vRef, .00001));
}

/*!
 * Test utility function for transforming a point by a spatial transform
 */
TEST(Spatial, pointTransform) {
  Vec3<double> p0(1, 2, 3);
  Vec3<double> pxRef(-3.0026, -1.9791, -3.7507);
  Mat3<double> R = coordinateRotation(CoordinateAxis::X, .2) *
                   coordinateRotation(CoordinateAxis::Y, .3) *
                   coordinateRotation(CoordinateAxis::Z, .5);

  Vec3<double> r(4, 5, 6);
  Mat6<double> X = createSXform(R, r);

  Vec3<double> px = sXFormPoint(X, p0);
  EXPECT_TRUE(almostEqual(px, pxRef, .0005));
}

/*!
 * Test utility function for converting force at a point to a spatial force
 */
TEST(Spatial, forceToSpatialForce) {
  Vec3<double> fLinear(1, 2, 3);
  Vec3<double> point(7, 8, 9);
  SVec<double> fRef;
  fRef << 6, -12, 6, 1, 2, 3;
  SVec<double> fSpatial = forceToSpatialForce(fLinear, point);
  EXPECT_TRUE(almostEqual(fSpatial, fRef, .0005));
}

/*!
 * Test utility function for converting the spatial velocity to the velocity at
 * a point
 */
TEST(Spatial, spatialToLinearVelocity) {
  SVec<double> vspat;
  vspat << 1.93, 2.34, 3.345, -4.23, 5.8383, 6.921;
  Vec3<double> vlin = vspat.tail(3);
  Vec3<double> vang = vspat.head(3);
  Vec3<double> point;
  point << -23.23, 2.638, 9.324;

  Vec3<double> vpoint = vlin + vang.cross(point);
  Vec3<double> vpoint2 = spatialToLinearVelocity(vspat, point);

  EXPECT_TRUE(almostEqual(vpoint2, vpoint, 1e-8));
}

/*!
 * Test utility function for converting the spatial acceleration to the linear
 * acceleration of a point
 */
TEST(Spatial, spatialToLinearAcceleration) {
  SVec<double> vspat, aspat;
  // Top spinning with angular velocity w on a skateboard that is traveling at
  // velocity v and is currently at point p
  double w = 5;
  double v = 1;
  double p = 1;

  vspat << 0, 0, w, v, -w * p, 0;
  aspat << 0, 0, 0, 0, -w * v, 0;

  Vec3<double> p2;
  p2 << p, 0, 0;

  Vec3<double> a1 = spatialToLinearAcceleration(aspat, vspat);

  Vec3<double> a2 = spatialToLinearAcceleration(aspat, vspat, p2);

  Vec3<double> a1_expected, a2_expected;
  a1_expected << w * w * p, 0, 0;
  a2_expected << 0, 0, 0;

  EXPECT_TRUE(almostEqual(a1, a1_expected, 1e-8));
  EXPECT_TRUE(almostEqual(a2, a2_expected, 1e-8));

  vspat << 1.93, 2.34, 3.345, -4.23, 5.8383, 6.921;
  aspat << -5.164, 68.4, 1.56879, -98.44, 8.14, 6.324;
  p2 << 54.797, -6.1654, 3.64587;

  Vec3<double> w1, v1, wd1, vd1;
  Vec3<double> w2, v2, wd2, vd2;

  w1 = vspat.head(3);
  v1 = vspat.tail(3);
  wd1 = aspat.head(3);
  vd1 = aspat.tail(3);

  w2 = w1;
  wd2 = w2;
  v2 = v1 + w1.cross(p2);
  vd2 = vd1 + wd1.cross(p2);

  a1_expected = vd1 + w1.cross(v1);
  a2_expected = vd2 + w2.cross(v2);

  a1 = spatialToLinearAcceleration(aspat, vspat);
  a2 = spatialToLinearAcceleration(aspat, vspat, p2);

  EXPECT_TRUE(almostEqual(a1, a1_expected, 1e-8));
  EXPECT_TRUE(almostEqual(a2, a2_expected, 1e-8));
}
