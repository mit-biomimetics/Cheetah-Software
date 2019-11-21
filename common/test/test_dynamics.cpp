/*! @file test_dynamics.cpp
 *  @brief Test dynamics algorithms
 *
 * Test the dynamics algorithms in DynamicsSimulator and models of Cheetah 3
 */

#include "Dynamics/Cheetah3.h"
#include "Dynamics/DynamicsSimulator.h"
#include "Dynamics/FloatingBaseModel.h"
#include "Dynamics/Quadruped.h"
#include "Utilities/utilities.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace spatial;

#include <stdio.h>
using namespace std;

/*!
 * Generate a model of the Cheetah 3 and check its total mass, tree structure
 * and sum of all spatial inertias (including rotors)
 * Compared against MATLAB model
 */
TEST(Dynamics, cheetah3model) {
  FloatingBaseModel<double> cheetah = buildCheetah3<double>().buildModel();
  cheetah.check();

  // check total mass
  EXPECT_TRUE(fpEqual(41.0737, cheetah.totalNonRotorMass(), .0001));
  EXPECT_TRUE(fpEqual(6.3215, cheetah.totalRotorMass(), .0001));

  // check tree structure
  std::vector<int> parentRef{0, 0, 0,  0, 0,  0,  5, 6,  7,
                             5, 9, 10, 5, 12, 13, 5, 15, 16};
  EXPECT_TRUE(vectorEqual(parentRef, cheetah.getParentVector()));

  // this is kind of stupid, but a reasonable sanity check for inertias
  Mat6<double> inertiaSumRef;
  inertiaSumRef << 0.4352, -0.0000, -0.0044, 0, 0.6230, -0.0000, -0.0000,
      1.0730, 0.0000, -0.6230, 0, -0.0822, -0.0044, 0.0000, 1.0071, 0.0000,
      0.0822, 0, 0, -0.6230, 0.0000, 42.6540, 0, 0, 0.6230, 0, 0.0822, 0,
      42.6540, 0, -0.0000, -0.0822, 0, 0, 0, 42.6540;

  Mat6<double> inertiaSum = Mat6<double>::Zero();

  for (size_t i = 0; i < 18; i++) {
    inertiaSum += cheetah.getBodyInertiaVector()[i].getMatrix() +
                  cheetah.getRotorInertiaVector()[i].getMatrix() * 0.25;
  }
  EXPECT_TRUE(almostEqual(inertiaSum, inertiaSumRef, .0001));
}

/*!
 * Test the model transforms
 */
TEST(Dynamics, cheetah3ModelTransforms) {
  FloatingBaseModel<double> cheetah = buildCheetah3<double>().buildModel();
  Mat6<double> XTotal = Mat6<double>::Identity();
  Mat6<double> XRotTotal = Mat6<double>::Identity();
  for (size_t i = 0; i < 18; i++) {
    XTotal = XTotal + cheetah._Xtree[i];
    XRotTotal = XRotTotal + cheetah._Xrot[i];
  }
  Mat6<double> Xtr, Xrtr;
  Xtr << 11.0000, 0.0000, 0, 0, 0, 0, -0.0000, 11.0000, 0, 0, 0, 0, 0, 0,
      19.0000, 0, 0, 0, 0, -1.3680, 0, 11.0000, 0.0000, 0, 1.3680, 0, 0.0000,
      -0.0000, 11.0000, 0, 0.0000, 0, 0, 0, 0, 19.0000;
  Xrtr << 11.0000, 0.0000, 0, 0, 0, 0, -0.0000, 11.0000, 0, 0, 0, 0, 0, 0,
      19.0000, 0, 0, 0, 0, 0, 0.0000, 11.0000, 0.0000, 0, 0, 0, 0, -0.0000,
      11.0000, 0, 0, 0, 0, 0, 0, 19.0000;

  EXPECT_TRUE(almostEqual(Xtr, XTotal, .0005));
  EXPECT_TRUE(almostEqual(Xrtr, XRotTotal, .0005));
}

/*!
 * Run the articulated body algorithm (and forward kinematics) on Cheetah 3
 * Set a weird body orientation, velocity, q, dq, and tau
 * Checks that quatD, pd, vd, and qdd match MATLAB
 */
TEST(Dynamics, simulatorDynamicsABANoExternalForceCheetah3) {
  FloatingBaseModel<double> cheetahModel = buildCheetah3<double>().buildModel();
  DynamicsSimulator<double> sim(cheetahModel, true);

  RotMat<double> rBody = coordinateRotation(CoordinateAxis::X, .123) *
                         coordinateRotation(CoordinateAxis::Z, .232) *
                         coordinateRotation(CoordinateAxis::Y, .111);

  SVec<double> bodyVel;
  bodyVel << 1, 2, 3, 4, 5, 6;
  FBModelState<double> x;
  DVec<double> q(12);
  DVec<double> dq(12);
  DVec<double> tau(12);
  for (size_t i = 0; i < 12; i++) {
    q[i] = i + 1;
    dq[i] = (i + 1) * 2;
    tau[i] = (i + 1) * -30.;
  }

  // set state
  x.bodyOrientation = rotationMatrixToQuaternion(rBody.transpose());
  x.bodyVelocity = bodyVel;
  x.bodyPosition = Vec3<double>(6, 7, 8);
  x.q = q;
  x.qd = dq;

  // do aba
  sim.setState(x);
  sim.forwardKinematics();
  sim.runABA(tau);

  // check:
  Vec3<double> pdRef(4.3717, 4.8598, 5.8541);
  SVec<double> vbdRef;
  vbdRef << 455.5224, 62.1684, -70.56, -11.4505, 10.2354, -15.2394;

  DVec<double> qddRef(12);
  qddRef << -0.7924, -0.2205, -0.9163, -1.6136, -0.4328, -1.6911, -2.9878,
      -0.9358, -2.6194, -3.3773, -1.3235, -3.1598;
  qddRef *= 1000;

  EXPECT_TRUE(almostEqual(pdRef, sim.getDState().dBodyPosition, .001));
  EXPECT_TRUE(almostEqual(vbdRef, sim.getDState().dBodyVelocity, .001));

  for (size_t i = 0; i < 12; i++) {
    // the qdd's are large - see qddRef, so we're only accurate to within ~1.
    EXPECT_TRUE(fpEqual(sim.getDState().qdd[i], qddRef[i], 3.));
  }

  // check the integrator for the body (orientation, position, and spatial
  // velocity). we use a huge timestep here so that any error in the integrator
  // isn't multiplied by something small
  sim.integrate(2.);

  Quat<double> quat2Ref(-0.8962, -0.0994, -0.2610, -0.3446);
  Vec3<double> x2Ref(14.7433, 16.7196, 19.7083);
  SVec<double> v2Ref;
  v2Ref << 912.0449, 126.3367, -138.1201, -18.9011, 25.4709, -24.4788;

  EXPECT_TRUE(almostEqual(quat2Ref, sim.getState().bodyOrientation, .0002));
  EXPECT_TRUE(almostEqual(x2Ref, sim.getState().bodyPosition, .0002));
  EXPECT_TRUE(almostEqual(v2Ref, sim.getState().bodyVelocity, .0002));
}

/*!
 * Run the RNEA and component H/Cqd/G algorithms on Cheetah 3
 * Set a weird body orientation, velocity, q, dq, qdd
 * Checks that genForce matches MATLAB, and CRBA/Cqd/G agree the output as well
 */
TEST(Dynamics, inverseDynamicsNoContacts) {
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
  cheetahModel.setState(x);

  // Construct state derivative
  SVec<double> vbd;
  vbd << 455.5224, 62.1684, -70.56, -11.4505, 10.2354, -15.2394;

  DVec<double> qdd(12);
  qdd << -0.7924, -0.2205, -0.9163, -1.6136, -0.4328, -1.6911, -2.9878, -0.9358,
      -2.6194, -3.3773, -1.3235, -3.1598;
  qdd *= 1000;

  FBModelStateDerivative<double> dx;
  dx.dBodyVelocity = vbd;
  dx.qdd = qdd;

  // Vector form of state derivative
  DVec<double> nu_dot(18);
  nu_dot.head(6) = dx.dBodyVelocity;
  nu_dot.tail(12) = qdd;

  // Components of the equations
  DMat<double> H = cheetahModel.massMatrix();
  DVec<double> Cqd = cheetahModel.generalizedCoriolisForce();
  DVec<double> G = cheetahModel.generalizedGravityForce();

  // Compute ID two different ways
  DVec<double> genForce = cheetahModel.inverseDynamics(dx);
  DVec<double> component_ID_result = H * nu_dot + Cqd + G;

  SVec<double> zero6x1 = SVec<double>::Zero();
  DVec<double> tauReturn = genForce.tail(12);
  SVec<double> fReturn = genForce.head(6);

  EXPECT_TRUE(almostEqual(zero6x1, fReturn, .03));
  EXPECT_TRUE(almostEqual(tauref, tauReturn, .01));
  EXPECT_TRUE(almostEqual(component_ID_result, genForce, 1e-6));
}

/*!
 * Computes Jacobians and bias forces for contacts
 * Checks that finite differenced contact positions equal velocities
 * and finite differenced contact velocities equal J*nu_dot + Jdot*nu
 */
TEST(Dynamics, contactJacobians) {
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
  cheetahModel.setState(x);

  // Construct state derivative
  SVec<double> vbd;
  vbd << 45.5224, 62.1684, -70.56, -11.4505, 10.2354, -15.2394;

  DVec<double> qdd(12);
  qdd << -0.7924, -0.2205, -0.9163, -1.6136, -0.4328, -1.6911, -2.9878, -0.9358,
      -2.6194, -3.3773, -1.3235, -3.1598;
  qdd *= 400;

  FBModelStateDerivative<double> dx;
  dx.dBodyVelocity = vbd;
  dx.qdd = qdd;

  // Vector form of state derivative
  DVec<double> nu(18), nu_dot(18);
  nu.head(6) = x.bodyVelocity;
  nu.tail(12) = x.qd;

  nu_dot.head(6) = dx.dBodyVelocity;
  nu_dot.tail(12) = qdd;

  cheetahModel.contactJacobians();

  // Compute Positions
  Vec3<double> pos0, posf, vel0, vel0_fromJ, velf, acc0, J0dqd;
  pos0 = cheetahModel._pGC[15];
  vel0 = cheetahModel._vGC[15];
  D3Mat<double> J0 = cheetahModel._Jc[15];
  vel0_fromJ = J0 * nu;
  acc0 = J0 * nu_dot + cheetahModel._Jcdqd[15];

  double dt = 1e-9;

  // Integrate the state
  Mat3<double> Rup = quaternionToRotationMatrix(x.bodyOrientation);
  dx.dBodyPosition = Rup.transpose() * x.bodyVelocity.tail(3);

  Vec3<double> omega0 = Rup.transpose() * x.bodyVelocity.head(3);
  Vec3<double> axis;
  double ang = omega0.norm();
  axis = omega0 / ang;

  ang *= dt;
  Vec3<double> ee = std::sin(ang / 2) * axis;
  Quat<double> quatD(std::cos(ang / 2), ee[0], ee[1], ee[2]);
  Quat<double> quatNew = quatProduct(quatD, x.bodyOrientation);
  quatNew = quatNew / quatNew.norm();

  // actual integration
  x.q += x.qd * dt;
  x.qd += dx.qdd * dt;
  x.bodyVelocity += dx.dBodyVelocity * dt;
  x.bodyPosition += dx.dBodyPosition * dt;
  x.bodyOrientation = quatNew;

  cheetahModel.setState(x);
  cheetahModel.contactJacobians();

  // Final
  posf = cheetahModel._pGC[15];
  velf = cheetahModel._vGC[15];

  // Finite Difference
  Vec3<double> vel_fd, acc_fd;
  vel_fd = (posf - pos0) / dt;
  acc_fd = (velf - vel0) / dt;

  EXPECT_TRUE(almostEqual(vel0, vel0_fromJ, 1e-6));
  EXPECT_TRUE(almostEqual(vel_fd, vel0, 1e-5));
  EXPECT_TRUE(almostEqual(acc_fd, acc0, 1e-3));
}

/*!
 * Run the articulated body algorithm (and forward kinematics) on Cheetah 3
 * Set a weird body orientation, velocity, q, dq, and tau
 * Sets external spatial forces on all bodies
 * Checks that quatD, pd, vd, and qdd match MATLAB
 */
TEST(Dynamics, simulatorDynamicsWithExternalForceCheetah3) {
  FloatingBaseModel<double> cheetahModel = buildCheetah3<double>().buildModel();
  DynamicsSimulator<double> sim(cheetahModel, true);

  RotMat<double> rBody = coordinateRotation(CoordinateAxis::X, .123) *
                         coordinateRotation(CoordinateAxis::Z, .232) *
                         coordinateRotation(CoordinateAxis::Y, .111);
  SVec<double> bodyVel;
  bodyVel << 1, 2, 3, 4, 5, 6;
  FBModelState<double> x;
  DVec<double> q(12);
  DVec<double> dq(12);
  DVec<double> tau(12);
  for (size_t i = 0; i < 12; i++) {
    q[i] = i + 1;
    dq[i] = (i + 1) * 2;
    tau[i] = (i + 1) * -30.;
  }

  // set state
  x.bodyOrientation = rotationMatrixToQuaternion(rBody.transpose());
  x.bodyVelocity = bodyVel;
  x.bodyPosition = Vec3<double>(6, 7, 8);
  x.q = q;
  x.qd = dq;

  // generate external forces
  vectorAligned<SVec<double>> forces(18);
  for (size_t i = 0; i < 18; i++) {
    for (size_t j = 0; j < 6; j++) {
      forces[i][j] = i + j + 1;
    }
  }

  // do aba
  sim.setState(x);
  sim.setAllExternalForces(forces);
  sim.step(0.0, tau, 5e5, 5e3);

  // check:
  Vec3<double> pdRef(4.3717, 4.8598, 5.8541);
  SVec<double> vbdRef;
  vbdRef << 806.6664, 44.1266, 33.8287, -10.1360, 2.1066, 12.1677;

  DVec<double> qddRef(12);
  qddRef << -0.5630, -0.1559, -1.0182, -0.9012, -0.3585, -1.6170, -2.0995,
      -0.8959, -2.7325, -2.0808, -1.3134, -3.1206;
  qddRef *= 1000;

  EXPECT_TRUE(almostEqual(pdRef, sim.getDState().dBodyPosition, .001));
  EXPECT_TRUE(almostEqual(vbdRef, sim.getDState().dBodyVelocity, .001));

  for (size_t i = 0; i < 12; i++) {
    // the qdd's are large - see qddRef, so we're only accurate to within ~1.
    EXPECT_TRUE(fpEqual(sim.getDState().qdd[i], qddRef[i], 3.));
  }
}

/*!
 * Run the articulated body algorithm (and forward kinematics) on Cheetah 3
 * Set a weird body orientation, velocity, q, dq, and tau
 * Checks that foot position and velocities match MATLAB
 */
TEST(Dynamics, simulatorFootPosVelCheetah3) {
  FloatingBaseModel<double> cheetahModel = buildCheetah3<double>().buildModel();
  DynamicsSimulator<double> sim(cheetahModel);

  RotMat<double> rBody = coordinateRotation(CoordinateAxis::X, .123) *
                         coordinateRotation(CoordinateAxis::Z, .232) *
                         coordinateRotation(CoordinateAxis::Y, .111);
  SVec<double> bodyVel;
  bodyVel << 1, 2, 3, 4, 5, 6;
  FBModelState<double> x;
  DVec<double> q(12);
  DVec<double> dq(12);
  DVec<double> tau(12);
  for (size_t i = 0; i < 12; i++) {
    q[i] = i + 1;
    dq[i] = (i + 1) * 2;
    tau[i] = (i + 1) * -30.;
  }

  // set state
  x.bodyOrientation = rotationMatrixToQuaternion(rBody.transpose());
  x.bodyVelocity = bodyVel;
  x.bodyPosition = Vec3<double>(6, 7, 8);
  x.q = q;
  x.qd = dq;

  // generate external forces
  vectorAligned<SVec<double>> forces(18);
  for (size_t i = 0; i < 18; i++) {
    for (size_t j = 0; j < 6; j++) {
      forces[i][j] = i + j + 1;
    }
  }

  // fwd kin is included in this
  sim.setState(x);
  sim.setAllExternalForces(forces);
  sim.step(0.0, tau, 5e5, 5e3);

  // Vec3<double> bodypRef1ML(6.25, 6.8271, 8.155);
  Vec3<double> bodypRef1ML(6.260735, 6.812413, 8.056675);
  Vec3<double> footpRefML(5.1594, 7.3559, 7.674);
  // Vec3<double> bodyvRef1ML(5.1989, 5.4008, 5.1234);
  Vec3<double> bodyvRef1ML(5.028498, 5.540016, 5.083884);
  Vec3<double> footvRefML(-9.3258, -0.1926, 26.3323);

  // I add the body points in a different order, so comparing them is kind of
  // annoying. this just tests one body point and one foot point.
  //EXPECT_TRUE(almostEqual(bodypRef1ML, sim.getModel()._pGC.at(2), .0005));
  EXPECT_TRUE(almostEqual(footpRefML, sim.getModel()._pGC.at(15), .0005));
  //EXPECT_TRUE(almostEqual(bodyvRef1ML, sim.getModel()._vGC.at(2), .0005));
  EXPECT_TRUE(almostEqual(footvRefML, sim.getModel()._vGC.at(15), .0005));
}
