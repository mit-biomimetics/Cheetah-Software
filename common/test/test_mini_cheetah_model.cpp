/*! @file test_dynamics.cpp
 *  @brief Test dynamics algorithms
 *
 * Test the model of Mini Cheetah
 */

#include "Dynamics/DynamicsSimulator.h"
#include "Dynamics/FloatingBaseModel.h"
#include "Dynamics/MiniCheetah.h"
#include "Dynamics/Quadruped.h"
#include "Utilities/utilities.h"

#include <iostream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Utilities/Utilities_print.h"

using namespace spatial;

/*!
 * Test the total mass, tree, and spatial inertia sum
 */
TEST(MiniCheetah, miniCheetahModel1) {
  FloatingBaseModel<double> cheetah = buildMiniCheetah<double>().buildModel();

  // masses
  EXPECT_TRUE(fpEqual(8.2520, cheetah.totalNonRotorMass(), .0001));
  EXPECT_TRUE(fpEqual(0.66, cheetah.totalRotorMass(), .0001));

  // check tree structure
  std::vector<int> parentRef{0, 0, 0,  0, 0,  0,  5, 6,  7,
                             5, 9, 10, 5, 12, 13, 5, 15, 16};
  EXPECT_TRUE(vectorEqual(parentRef, cheetah.getParentVector()));

  // this is kind of stupid, but a reasonable sanity check for inertias
  Mat6<double> inertiaSumRef;
  inertiaSumRef <<
      // 0.0272,      0, 0.0001,       0, 0.0663, 0,
      // 0, 0.3758,    0.0, -0.0663,      0, 0,
      // 0.0001, 0.0000, 0.0497, 0, 0, 0,
      // 0, -0.0663, 0, 8.4170, 0, 0,
      // 0.0663, 0, 0, 0, 8.4170, 0,
      // 0, 0, 0, 0, 0, 8.4170;

      0.0272,
      0, 0, 0, 0.0663, 0, 0, 0.05, 0, -0.066336, 0, 0, 0, 0, 0.0497, 0, 0, 0, 0,
      -0.066336, 0, 8.417, 0, 0, 0.066336, 0, 0, 0, 8.417, 0, 0, 0, 0, 0, 0,
      8.417;

  Mat6<double> inertiaSum = Mat6<double>::Zero();

  for (size_t i = 0; i < 18; i++) {
    inertiaSum += cheetah.getBodyInertiaVector()[i].getMatrix() +
                  cheetah.getRotorInertiaVector()[i].getMatrix() * 0.25;
  }

  EXPECT_TRUE(almostEqual(inertiaSum, inertiaSumRef, .0003));
}

/*!
 * Test the model transforms
 */
TEST(MiniCheetah, miniCheetahModel2) {
  FloatingBaseModel<double> cheetah = buildMiniCheetah<double>().buildModel();
  Mat6<double> XTotal = Mat6<double>::Identity();
  Mat6<double> XRotTotal = Mat6<double>::Identity();
  for (size_t i = 0; i < 18; i++) {
    XTotal = XTotal + cheetah._Xtree[i];
    XRotTotal = XRotTotal + cheetah._Xrot[i];
  }
  Mat6<double> Xtr, Xrtr;
  Xtr << 11.0000, 0.0000, 0, 0, 0, 0, -0.0000, 11.0000, 0, 0, 0, 0, 0, 0,
      19.0000, 0, 0, 0, 0, -0.8360, 0, 11.0000, 0.0000, 0, 0.8360, 0, 0.0000,
      -0.0000, 11.0000, 0, 0, 0, 0, 0, 0, 19.0000;

  Xrtr << 11.0000, 0.0000, 0, 0, 0, 0, -0.0000, 11.0000, 0, 0, 0, 0, 0, 0,
      19.0000, 0, 0, 0, 0, 0, 0, 11.0000, 0.0000, 0, 0, 0, 0.0000, -0.0000,
      11.0000, 0, 0.0000, 0, 0, 0, 0, 19.0000;

  EXPECT_TRUE(almostEqual(Xtr, XTotal, .0005));
  EXPECT_TRUE(almostEqual(Xrtr, XRotTotal, .0005));
}

/*!
 * Creates a cheetah model and runs forward kinematics and ABA
 * Doesn't test anyting - this is just to make sure it doesn't crash
 */
TEST(MiniCheetah, simulatorDynamicsDoesntCrashMiniCheetah) {
  FloatingBaseModel<double> cheetah = buildMiniCheetah<double>().buildModel();
  DynamicsSimulator<double> sim(cheetah);
  DVec<double> tau(12);
  sim.forwardKinematics();
  sim.runABA(tau);
}

/*!
 * Run the articulated body algorithm (and forward kinematics) on Cheetah 3
 * Set a weird body orientation, velocity, q, dq, and tau
 * Checks that quatD, pd, vd, and qdd match MATLAB
 */
TEST(MiniCheetah, simulatorDynamicsABANoExternalForceMiniCheetah) {
  FloatingBaseModel<double> cheetahModel =
      buildMiniCheetah<double>().buildModel();
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
    tau[i] = (i + 1) * -3.;
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
  vbdRef << 1217.31, -182.793, -171.524, -4.09393, 33.6798, -59.975;

  DVec<double> qddRef(12);
  qddRef << -1946.36, -1056.62, -1470.57, -3133.03, -1263.89, -2626.9, -6637.38,
      -4071.88, -4219.14, -7554.08, -3392.44, -5112.85;

  EXPECT_TRUE(almostEqual(pdRef, sim.getDState().dBodyPosition, .001));
  EXPECT_TRUE(almostEqual(vbdRef, sim.getDState().dBodyVelocity, 1));

  for (size_t i = 0; i < 12; i++) {
    // the qdd's are large - see qddRef, so we're only accurate to within ~1.
    EXPECT_TRUE(fpEqual(sim.getDState().qdd[i], qddRef[i], 3.));
  }
}

/*!
 * Run the articulated body algorithm (and forward kinematics) on Cheetah 3
 * Set a weird body orientation, velocity, q, dq, and tau
 * Sets external spatial forces on all bodies
 * Checks that quatD, pd, vd, and qdd match MATLAB
 */
TEST(MiniCheetah, simulatorDynamicsWithExternalForceMiniCheetah) {
  FloatingBaseModel<double> cheetahModel =
      buildMiniCheetah<double>().buildModel();
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
    tau[i] = (i + 1) * -3.;
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
      forces[i][j] = .3 * (i + j + 1);
    }
  }

  // do aba
  sim.setState(x);
  sim.setAllExternalForces(forces);
  sim.step(0.0, tau, 5e5, 5e3);

  // check:
  Vec3<double> pdRef(4.3717, 4.8598, 5.8541);
  SVec<double> vbdRef;
  vbdRef << 3153.46, 42.6931, 264.584, -3.80573, -53.3519, 68.5713;

  DVec<double> qddRef(12);
  qddRef << -1211.85, -1167.39, -2024.94, 394.414, -297.854, -2292.29, -1765.92,
      -4040.78, -4917.41, 149.368, -2508.2, -4969.02;

  EXPECT_TRUE(almostEqual(pdRef, sim.getDState().dBodyPosition, .001));
  EXPECT_TRUE(almostEqual(vbdRef, sim.getDState().dBodyVelocity, .01));

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
TEST(MiniCheetah, simulatorFootPosVelMiniCheetah) {
  FloatingBaseModel<double> cheetahModel =
      buildMiniCheetah<double>().buildModel();
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

  Vec3<double> footpRefML(5.4937, 7.1459, 7.8096);
  Vec3<double> footvRefML(-2.5284, 2.0944, 16.3732);

  // I add the body points in a different order, so comparing them is kind of
  // annoying. this just one foot point.
  EXPECT_TRUE(almostEqual(footpRefML, sim.getModel()._pGC.at(15), .0005));
  EXPECT_TRUE(almostEqual(footvRefML, sim.getModel()._vGC.at(15), .0005));
}

/*!
 * Check that the hip location convention is correct
 */
TEST(MiniCheetah, hipLocationConvention) {
  Vec3<double> hipLocationRef[4];
  hipLocationRef[0] = Vec3<double>(0.19, -0.049, 0.0);
  hipLocationRef[1] = Vec3<double>(0.19, 0.049, 0.0);
  hipLocationRef[2] = Vec3<double>(-0.19, -0.049, 0.0);
  hipLocationRef[3] = Vec3<double>(-0.19, 0.049, 0.0);

  Vec3<double> hipLocations[4];

  auto quadruped = buildMiniCheetah<double>();

  for (int i = 0; i < 4; i++) {
    hipLocations[i] = quadruped.getHipLocation(i);
  }

  for (int i = 0; i < 4; i++) {
    // std::cout << "ref: " << hipLocationRef[i].transpose() << "\nact: " <<
    // hipLocations[i].transpose() << "\n\n";
    EXPECT_TRUE(almostEqual(hipLocations[i], hipLocationRef[i], .0001));
  }
}

TEST(MiniCheetah, ContactPositionVelocity) {
  FloatingBaseModel<double> cheetahModel =
      buildMiniCheetah<double>().buildModel();
  DynamicsSimulator<double> sim(cheetahModel);

  RotMat<double> rBody = coordinateRotation(CoordinateAxis::X, 0.0) *
                         coordinateRotation(CoordinateAxis::Z, 0.0) *
                         coordinateRotation(CoordinateAxis::Y, 0.0);

  SVec<double> bodyVel;
  bodyVel << 0., 0., 0., 1.0, 0., 0.;
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
  cheetahModel.setState(x);
  cheetahModel.forwardKinematics();  // compute forward kinematics

  double length(0.19);
  double width(0.049);

  Vec3<double> FR_abd_ref(x.bodyPosition[0] + length, x.bodyPosition[1] - width,
                          x.bodyPosition[2]);

  Vec3<double> FL_abd_ref(x.bodyPosition[0] + length, x.bodyPosition[1] + width,
                          x.bodyPosition[2]);

  Vec3<double> HR_abd_ref(x.bodyPosition[0] - length, x.bodyPosition[1] - width,
                          x.bodyPosition[2]);

  Vec3<double> HL_abd_ref(x.bodyPosition[0] - length, x.bodyPosition[1] + width,
                          x.bodyPosition[2]);

  // for(size_t i(0); i<8; ++i){
  // printf("%lu th\n", i);
  // pretty_print(cheetahModel._pGC.at(i), std::cout, "cp pos");
  //}

  // I add the body points in a different order, so comparing them is kind of
  // annoying. this just one foot point.
  EXPECT_TRUE(
      almostEqual(FR_abd_ref, cheetahModel._pGC.at(linkID::FR_abd), .0005));
  EXPECT_TRUE(
      almostEqual(FL_abd_ref, cheetahModel._pGC.at(linkID::FL_abd), .0005));
  EXPECT_TRUE(
      almostEqual(HR_abd_ref, cheetahModel._pGC.at(linkID::HR_abd), .0005));
  EXPECT_TRUE(
      almostEqual(HL_abd_ref, cheetahModel._pGC.at(linkID::HL_abd), .0005));
}

TEST(MiniCheetah, InertiaProperty) {
  FloatingBaseModel<double> cheetahModel =
      buildMiniCheetah<double>().buildModel();

  SVec<double> bodyVel;
  FBModelState<double> x;
  DVec<double> q(12);
  DVec<double> dq(12);
  DVec<double> tau(12);

  for (size_t i = 0; i < 12; i++) {
    q[i] = 0.;
    dq[i] = 0.;
    tau[i] = 0.;
  }
  q[1] = -M_PI / 2.;
  q[4] = -M_PI / 2.;
  q[7] = M_PI / 2.;
  q[10] = M_PI / 2.;

  Vec3<double> rpy_ori;
  rpy_ori.setZero();
  // set state
  x.bodyOrientation = ori::rpyToQuat(rpy_ori);
  x.bodyVelocity.setZero();
  x.bodyPosition.setZero();
  x.q = q;
  x.qd = dq;

  cheetahModel.setState(x);
  cheetahModel.massMatrix();

  DMat<double> H = cheetahModel.getMassMatrix();

  // pretty_print(H, std::cout, "H");
  // exit(0);
}
