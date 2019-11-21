/*! @file test_utilities.cpp
 *  @brief Test Utilities functions
 *
 * Test the various utilities
 */

#include "SimUtilities/ImuSimulator.h"
#include "cppTypes.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

// check the imu data when stationary and at the origin
TEST(ImuSimulator, passThrough) {
  // no noise
  SimulatorControlParameters simParams;
  simParams.vectornav_imu_quat_noise = 0.f;
  simParams.vectornav_imu_accelerometer_noise = 0.f;
  simParams.vectornav_imu_gyro_noise = 0.f;


  ImuSimulator<double> imuSim(simParams, 0);

  FBModelState<double> xfb;
  FBModelStateDerivative<double> xfbd;

  Quat<double> eye;
  eye << 1, 0, 0, 0;

  SVec<double> v0;
  v0 << 0, 0, 0, 0, 0, 0;

  // stationary, upright
  xfb.bodyOrientation = eye;
  xfb.bodyVelocity = v0;
  xfbd.dBodyVelocity = v0;

  VectorNavData vn;

  // update data
  imuSim.updateVectornav(xfb, xfbd, &vn);

  Vec3<float> imuZero(0, 0, 0);
  Vec3<float> imuGravity(0, 0, 9.81);
  Quat<float> imuLevel;
  // imuLevel << 1,0,0,0;
  imuLevel << 0, 0, 0,
      1;  // vectornave quaternion description sequence is (x, y, z, w)

  EXPECT_TRUE(almostEqual(imuZero, vn.gyro, .00001f));
  EXPECT_TRUE(almostEqual(imuLevel, vn.quat, .00001f));
  EXPECT_TRUE(almostEqual(imuGravity, vn.accelerometer, .00001f));
}

// test with
// no omega or body acceleration, but linear, velocity, orientation, and
// position
TEST(ImuSimulator, orientation) {
  // no noise
  SimulatorControlParameters simParams;
  simParams.vectornav_imu_quat_noise = 0.f;
  simParams.vectornav_imu_accelerometer_noise = 0.f;
  simParams.vectornav_imu_gyro_noise = 0.f;

  ImuSimulator<double> imuSim(simParams, 0);

  FBModelState<double> xfb;
  FBModelStateDerivative<double> xfbd;

  // pitch the robot down:
  Quat<double> quat =
      rotationMatrixToQuaternion(coordinateRotation(CoordinateAxis::Y, .2));
  Vec3<double> position(2, 4, 3);
  SVec<double> velocity;
  velocity << 0, 0, 0, 2, 3, 6;
  SVec<double> acceleration = SVec<double>::Zero();

  xfb.bodyVelocity = velocity;
  xfb.bodyPosition = position;
  xfb.bodyOrientation = quat;
  xfbd.dBodyPosition = velocity.tail<3>();
  xfbd.dBodyVelocity = acceleration;

  VectorNavData vn;

  // update data
  imuSim.updateVectornav(xfb, xfbd, &vn);

  Vec3<float> imuZero(0, 0, 0);
  Vec3<float> imuGravity(0, 0, 9.81);
  // Quat<float> imuLevel = quat.cast<float>();
  Quat<float> imuLevel;
  Quat<float> sim_quat = quat.cast<float>();
  // vectornave quaternion description sequence is (x, y, z, w)
  imuLevel << sim_quat[1], sim_quat[2], sim_quat[3], sim_quat[0];

  EXPECT_TRUE(almostEqual(imuZero, vn.gyro, .00001f));
  EXPECT_TRUE(almostEqual(imuLevel, vn.quat, .00001f));
  EXPECT_TRUE(vn.accelerometer[0] < .0001f);
}

// test with omega
TEST(ImuSimulator, omega) {
  // no noise
  SimulatorControlParameters simParams;
  simParams.vectornav_imu_quat_noise = 0.f;
  simParams.vectornav_imu_accelerometer_noise = 0.f;
  simParams.vectornav_imu_gyro_noise = 0.f;

  ImuSimulator<double> imuSim(simParams, 0);

  FBModelState<double> xfb;
  FBModelStateDerivative<double> xfbd;

  // pitch the robot down:
  Quat<double> quat =
      rotationMatrixToQuaternion(coordinateRotation(CoordinateAxis::Y, .2));
  Vec3<double> position(2, 4, 3);
  SVec<double> velocity;
  velocity << 1, 2, 3, 0, 0, 0;
  SVec<double> acceleration = SVec<double>::Zero();

  xfb.bodyVelocity = velocity;
  xfb.bodyPosition = position;
  xfb.bodyOrientation = quat;
  xfbd.dBodyPosition = velocity.tail<3>();
  xfbd.dBodyVelocity = acceleration;

  VectorNavData vn;

  // update data
  imuSim.updateVectornav(xfb, xfbd, &vn);

  Vec3<float> imuZero(0, 0, 0);
  Vec3<float> imuGravity(0, 0, 9.81);
  Quat<float> imuLevel;
  Quat<float> sim_quat = quat.cast<float>();
  imuLevel << sim_quat[1], sim_quat[2], sim_quat[3], sim_quat[0];
  Vec3<float> omegaB = xfb.bodyVelocity.head<3>().cast<float>();

  EXPECT_TRUE(almostEqual(omegaB, vn.gyro, .00001f));
  EXPECT_TRUE(almostEqual(imuLevel, vn.quat, .00001f));
}

// test running in a circle omega cross v terms
TEST(ImuSimulator, omegaCrossV) {
  // no noise
  SimulatorControlParameters simParams;
  simParams.vectornav_imu_quat_noise = 0.f;
  simParams.vectornav_imu_accelerometer_noise = 0.f;
  simParams.vectornav_imu_gyro_noise = 0.f;

  ImuSimulator<double> imuSim(simParams, 0);

  FBModelState<double> xfb;
  FBModelStateDerivative<double> xfbd;

  // robot is upright, running in clockwise circle when viewed from above, at 12
  // o'clock:
  Quat<double> quat = rotationMatrixToQuaternion(Mat3<double>::Identity());
  Vec3<double> position(2, 4, 3);
  SVec<double> velocity;
  velocity << 0, 0, -2, 3, 0, 0;  // clockwise
  SVec<double> acceleration =
      SVec<double>::Zero();  // zero spatial acceleration

  xfb.bodyVelocity = velocity;
  xfb.bodyPosition = position;
  xfb.bodyOrientation = quat;
  xfbd.dBodyPosition = velocity.tail<3>();
  xfbd.dBodyVelocity = acceleration;

  VectorNavData vn;

  // update data
  imuSim.updateVectornav(xfb, xfbd, &vn);

  Vec3<float> imuZero(0, 0, 0);
  Vec3<float> imuGravity(0, 0, 9.81);
  Quat<float> imuLevel;
  Quat<float> sim_quat = quat.cast<float>();
  imuLevel << sim_quat[1], sim_quat[2], sim_quat[3], sim_quat[0];
  Vec3<float> omegaB = xfb.bodyVelocity.head<3>().cast<float>();

  // we should see -y acceleration (man sitting inside robot will feel pulled
  // toward the outside)
  Vec3<float> acc_ref(0, -6, 9.81);

  EXPECT_TRUE(almostEqual(omegaB, vn.gyro, .00001f));
  EXPECT_TRUE(almostEqual(imuLevel, vn.quat, .00001f));
  EXPECT_TRUE(almostEqual(vn.accelerometer, acc_ref, .0001f));
}

// test noise:
TEST(ImuSimulator, noise) {
  SimulatorControlParameters simParamsNoise;
  simParamsNoise.vectornav_imu_quat_noise = 0.005f;
  simParamsNoise.vectornav_imu_accelerometer_noise = 0.002f;
  simParamsNoise.vectornav_imu_gyro_noise = 0.005f;

  SimulatorControlParameters simParams;
  simParams.vectornav_imu_quat_noise = 0.f;
  simParams.vectornav_imu_accelerometer_noise = 0.f;
  simParams.vectornav_imu_gyro_noise = 0.f;


  ImuSimulator<double> imuSimNoise(simParamsNoise, 0);

  FBModelState<double> xfb;
  FBModelStateDerivative<double> xfbd;

  Quat<double> quat = rotationMatrixToQuaternion(
      ori::coordinateRotation(CoordinateAxis::X, 2.34) *
      ori::coordinateRotation(CoordinateAxis::Y, -5.3) *
      ori::coordinateRotation(CoordinateAxis::Z, -.3));

  Vec3<double> position(2, 4, 3);
  SVec<double> velocity;
  velocity << 2, 1, -2, 3, 3, 2;  // clockwise
  SVec<double> acceleration =
      SVec<double>::Zero();  // zero spatial acceleration

  xfb.bodyVelocity = velocity;
  xfb.bodyPosition = position;
  xfb.bodyOrientation = quat;
  xfbd.dBodyPosition = velocity.tail<3>();
  xfbd.dBodyVelocity = acceleration;

  VectorNavData vn;

  Vec3<float> rpy_ref = ori::quatToRPY(quat.cast<float>());
  Vec3<float> rpy_avg(0, 0, 0);

  constexpr size_t nIterations = 1000;
  for (size_t i = 0; i < nIterations; i++) {
    imuSimNoise.updateVectornav(xfb, xfbd, &vn);
    Quat<float> sim_quat;
    sim_quat << vn.quat[3], vn.quat[0], vn.quat[1], vn.quat[2];
    Vec3<float> rpy = ori::quatToRPY(sim_quat);
    rpy_avg += rpy / float(nIterations);
  }

  EXPECT_TRUE(almostEqual(rpy_avg, rpy_ref, .01));
}
