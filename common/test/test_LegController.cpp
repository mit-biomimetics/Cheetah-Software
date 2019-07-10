/*! @file test_LegController.cpp
 *  Todo: it would be nice to test these against the full dynamics model!
 */

#include "Controllers/LegController.h"
#include "Dynamics/MiniCheetah.h"
#include "Dynamics/Quadruped.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

/*!
 * Test jacobian by verifying that
 * fwd_kin(q + dq) = fwd_kin(q) + J * qd
 */
TEST(LegController, JacobianAndFwdKinematics1) {
  Quadruped<double> quadruped = buildMiniCheetah<double>();
  Vec3<double> q0(1, 2, 3);
  Vec3<double> dq(.001, .002, -.001);
  Mat3<double> J;
  Vec3<double> p0, p1Ref, p1;

  Vec3<double> q1 = q0 + dq;
  computeLegJacobianAndPosition(quadruped, q0, &J, &p0, 0);
  computeLegJacobianAndPosition(quadruped, q1, (Mat3<double>*)nullptr, &p1Ref,
                                0);

  p1 = p0 + J * dq;

  EXPECT_TRUE(almostEqual(p1, p1Ref, .001 * .001));
}

/*!
 * Test jacobian by finding it analytically
 */
TEST(LegController, JacobianAndFwdKinematics2) {
  Quadruped<double> quadruped = buildMiniCheetah<double>();
  Vec3<double> q0(4, 5, 6);
  Vec3<double> p0;
  Mat3<double> Jref;
  Mat3<double> J;

  computeLegJacobianAndPosition(quadruped, q0, &Jref, &p0, 0);
  double d = .001;

  for (int dim = 0; dim < 3; dim++) {
    Vec3<double> dq = Vec3<double>::Zero();
    dq(dim) = d;
    Vec3<double> q1 = q0 + dq;
    Vec3<double> p1;
    computeLegJacobianAndPosition(quadruped, q1, (Mat3<double>*)nullptr, &p1,
                                  0);
    Vec3<double> dp = p1 - p0;
    J.block<3, 1>(0, dim) = dp / d;
  }

  EXPECT_TRUE(almostEqual(J, Jref, .001));
}

/*!
 * Check that the foot is in the right spot when all joints are at zero
 */
TEST(LegController, FwdKinematicsLegSign) {
  Quadruped<double> quadruped = buildMiniCheetah<double>();
  Vec3<double> q(0, 0, 0);
  Vec3<double> p;
  computeLegJacobianAndPosition(quadruped, q, (Mat3<double>*)nullptr, &p, 0);

  Vec3<double> pRef(0, -quadruped._abadLinkLength,
                    -quadruped._hipLinkLength - quadruped._kneeLinkLength);

  EXPECT_TRUE(almostEqual(pRef, p, .00001));

  computeLegJacobianAndPosition(quadruped, q, (Mat3<double>*)nullptr, &p, 1);

  Vec3<double> pRef2(0, quadruped._abadLinkLength,
                     -quadruped._hipLinkLength - quadruped._kneeLinkLength);

  EXPECT_TRUE(almostEqual(pRef2, p, .00001));
}