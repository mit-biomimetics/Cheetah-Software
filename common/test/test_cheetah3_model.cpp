/*! @file test_dynamics.cpp
 *  @brief Test dynamics algorithms
 *
 * Test the dynamics algorithms in DynamicsSimulator and models of Cheetah 3
 */

#include "Dynamics/FloatingBaseModel.h"
#include "Dynamics/Quadruped.h"
#include "Utilities/utilities.h"
//#include "DynamicsSimulator.h"
#include "include/Dynamics/Cheetah3.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace spatial;

#include <stdio.h>
using namespace std;

/*!
 * Creates a cheetah model and runs forward kinematics and ABA
 * Doesn't test anyting - this is just to make sure it doesn't crash
 */
TEST(Cheetah3, simulatorDynamicsDoesntCrashCheetah3) {
  FloatingBaseModel<double> cheetah = buildCheetah3<double>().buildModel();
  DVec<double> tau(12);
  FBModelStateDerivative<double> _dstate;
  _dstate.qdd = DVec<double>::Zero(cheetah._nDof - 6);

  cheetah.forwardKinematics();
  cheetah.runABA(tau, _dstate);
}

/*!
 * Run the contact inertia algorithm for cheetah 3
 * Set a weird body orientation, velocity, q, dq, and tau
 * Checks that it matches the result from J H^(-1) J^T
 */
TEST(Cheetah3, simulatorContactInertiaCheetah3) {
  FloatingBaseModel<double> cheetahModel = buildCheetah3<double>().buildModel();

  RotMat<double> rBody = coordinateRotation(CoordinateAxis::X, .123) *
                         coordinateRotation(CoordinateAxis::Z, .232) *
                         coordinateRotation(CoordinateAxis::Y, .111);
  SVec<double> bodyVel;
  bodyVel << 1, 2, 3, 4, 5, 6;
  FBModelState<double> x;
  DVec<double> q(12);
  DVec<double> dq(12);
  DVec<double> tauref(12);
  for (size_t i = 0; i < 12; i++) {
    q[i] = i + 1;
    dq[i] = (i + 1) * 2;
    tauref[i] = (i + 1) * -30.;
  }

  // set state
  x.bodyOrientation = rotationMatrixToQuaternion(rBody.transpose());
  x.bodyVelocity = bodyVel;
  x.bodyPosition = Vec3<double>(6, 7, 8);
  x.q = q;
  x.qd = dq;

  // do aba
  cheetahModel.setState(x);

  cheetahModel.contactJacobians();
  DMat<double> H = cheetahModel.massMatrix();
  D3Mat<double> J0 = cheetahModel._Jc[15];
  DMat<double> LambdaInv1 = J0 * H.colPivHouseholderQr().solve(J0.transpose());

  D6Mat<double> forceOnly = D6Mat<double>::Zero(6, 3);
  forceOnly.bottomRows<3>() = Mat3<double>::Identity();

  DMat<double> LambdaInv2 = cheetahModel.invContactInertia(15, forceOnly);

  EXPECT_TRUE(almostEqual(LambdaInv1, LambdaInv2, .001));
}

/*!
 * Check a test force for cheetah 3
 * Set a weird body orientation, velocity, q, dq, and tau
 * Checks the result matches qdd = H^(-1) J^T F
 */
TEST(Cheetah3, simulatorTestForceCheetah3) {
  FloatingBaseModel<double> cheetahModel = buildCheetah3<double>().buildModel();

  RotMat<double> rBody = coordinateRotation(CoordinateAxis::X, .123) *
                         coordinateRotation(CoordinateAxis::Z, .232) *
                         coordinateRotation(CoordinateAxis::Y, .111);
  SVec<double> bodyVel;
  bodyVel << 1, 2, 3, 4, 5, 6;
  FBModelState<double> x;
  DVec<double> q(12);
  DVec<double> dq(12);
  DVec<double> tauref(12);
  for (size_t i = 0; i < 12; i++) {
    q[i] = i + 1;
    dq[i] = (i + 1) * 2;
    tauref[i] = (i + 1) * -30.;
  }

  // set state
  x.bodyOrientation = rotationMatrixToQuaternion(rBody.transpose());
  x.bodyVelocity = bodyVel;
  x.bodyPosition = Vec3<double>(6, 7, 8);
  x.q = q;
  x.qd = dq;

  FBModelStateDerivative<double> dx;
  dx.qdd.setZero(12);

  cheetahModel.setState(x);

  cheetahModel.contactJacobians();
  DMat<double> H = cheetahModel.massMatrix();
  DMat<double> J0 = cheetahModel._Jc[15].bottomRows(1);

  Vec3<double> zforce;
  zforce << 0, 0, 1;

  double foot_accel_in_z = cheetahModel.applyTestForce(15, zforce, dx);

  DVec<double> qddFull = H.colPivHouseholderQr().solve(J0.transpose());
  DMat<double> LambdaInv = J0 * qddFull;

  SVec<double> dBodyVelocity = qddFull.head<6>();
  DVec<double> qdd = qddFull.tail(12);

  double foot_accel_2 = LambdaInv(0, 0);

  double LambdaInv2 = cheetahModel.invContactInertia(15, zforce);

  EXPECT_TRUE(almostEqual(dBodyVelocity, dx.dBodyVelocity, .001));
  EXPECT_TRUE(almostEqual(qdd, dx.qdd, 1e-6));
  EXPECT_TRUE(abs(foot_accel_2 - foot_accel_in_z) < .001);
  EXPECT_TRUE(abs(LambdaInv2 - foot_accel_in_z) < .001);

  qdd(0) += 5;
  EXPECT_FALSE(almostEqual(qdd, dx.qdd, 1e-6));
}
