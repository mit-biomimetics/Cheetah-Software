/*! @file test_actuatorModel.cpp
 *  @brief Test the actuator model of the mini cheetah and cheetah 3 robots
 *
 *
 */

#include "Dynamics/ActuatorModel.h"
#include "Dynamics/Cheetah3.h"
#include "Dynamics/MiniCheetah.h"
#include "Dynamics/Quadruped.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

TEST(ActuatorModel, miniCheetah) {
  Quadruped<double> quad = buildMiniCheetah<double>();
  auto actuatorModels = quad.buildActuatorModels();
  auto& hipModel = actuatorModels[1];

  // first, disable friction
  hipModel.setFriction(false);

  double tauMaxPositive = 0, tauMaxNegative = 0;
  double maxQdMaxTorque = 0, qdMax = 0;

  // check our max torque in the positive direction:
  for (double tau = 0; tau < 200; tau += .1) {
    double tauAct = hipModel.getTorque(tau, 0);
    if (!fpEqual(tau, tauAct, .0001)) {
      tauMaxPositive = tauAct;
      break;
    }
  }

  for (double tau = 0; tau > -200; tau -= .1) {
    double tauAct = hipModel.getTorque(tau, 0);
    if (!fpEqual(tau, tauAct, .0001)) {
      tauMaxNegative = tauAct;
      break;
    }
  }

  for (double qd = 0; qd < 40; qd += .01) {
    double tauAct = hipModel.getTorque(18, qd);
    if (!fpEqual(18., tauAct, .0001)) {
      maxQdMaxTorque = qd;
      break;
    }
  }

  for (double qd = 0; qd < 60; qd += .01) {
    double tauAct = hipModel.getTorque(18, qd);
    if (tauAct <= 0) {
      qdMax = qd;
      break;
    }
  }

  EXPECT_TRUE(fpEqual(tauMaxPositive, 18., .0001));
  EXPECT_TRUE(fpEqual(tauMaxNegative, -18., .0001));
  EXPECT_TRUE(fpEqual(maxQdMaxTorque, 28.47, .02));
  EXPECT_TRUE(fpEqual(qdMax, 40., .0001));
}

TEST(ActuatorModel, cheetah3) {
  Quadruped<double> quad = buildCheetah3<double>();
  auto actuatorModels = quad.buildActuatorModels();
  auto& hipModel = actuatorModels[1];

  // first, disable friction
  hipModel.setFriction(false);

  double tauMaxPositive = 0, tauMaxNegative = 0;
  double maxQdMaxTorque = 0, qdMax = 0;

  // check our max torque in the positive direction:
  for (double tau = 0; tau < 300; tau += .1) {
    double tauAct = hipModel.getTorque(tau, 0);
    if (!fpEqual(tau, tauAct, .0001)) {
      tauMaxPositive = tauAct;
      break;
    }
  }

  for (double tau = 0; tau > -300; tau -= .1) {
    double tauAct = hipModel.getTorque(tau, 0);
    if (!fpEqual(tau, tauAct, .0001)) {
      tauMaxNegative = tauAct;
      break;
    }
  }

  for (double qd = 0; qd < 40; qd += .01) {
    double tauAct = hipModel.getTorque(208, qd);
    if (!fpEqual(208., tauAct, .0001)) {
      maxQdMaxTorque = qd;
      break;
    }
  }

  for (double qd = 0; qd < 60; qd += .01) {
    double tauAct = hipModel.getTorque(208, qd);
    if (tauAct <= 0) {
      qdMax = qd;
      break;
    }
  }

  EXPECT_TRUE(fpEqual(tauMaxPositive, 208.6, .1));
  EXPECT_TRUE(fpEqual(tauMaxNegative, -208.6, .1));
  EXPECT_TRUE(fpEqual(maxQdMaxTorque, 8.44, .02));
  EXPECT_TRUE(fpEqual(qdMax, 15.94, .1));
}