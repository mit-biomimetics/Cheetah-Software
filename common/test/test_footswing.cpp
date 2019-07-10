#include <include/Controllers/FootSwingTrajectory.h>
#include "cppTypes.h"
#include "Math/MathUtilities.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"


TEST(FootSwing, fstest) {
  FootSwingTrajectory<double> traj;
  Vec3<double> p0(1,1,1);
  traj.setInitialPosition(p0);
  Vec3<double> pf(2,2,1.2);
  traj.setFinalPosition(pf);
  traj.setHeight(1);

  // initial
  traj.computeSwingTrajectoryBezier(0, 0.5);
  EXPECT_TRUE(almostEqual(p0, traj.getPosition(), 0.001));

  // midpoint
  traj.computeSwingTrajectoryBezier(0.5, 0.5);
  Vec3<double> pMid = (p0 + pf)/2;
  pMid[2] = 2;
  EXPECT_TRUE(almostEqual(pMid, traj.getPosition(), 0.0001));

  // final
  traj.computeSwingTrajectoryBezier(1, 0.5);
  EXPECT_TRUE(almostEqual(pf, traj.getPosition(), 0.00001));




}

TEST(FootSwing, fsderiv) {
  FootSwingTrajectory<double> traj;
  Vec3<double> _p0(1,1,1);
  traj.setInitialPosition(_p0);
  Vec3<double> pf(2,2,1.2);
  traj.setFinalPosition(pf);
  traj.setHeight(1);


  // velocity/acceleration:
  double phasePerSecond = 2.; // 0.5 second swing.
  double t0 = 0.27;
  double dt = 0.001;
  double t1 = t0 + dt;
  double ph0 = t0 * phasePerSecond;
  double ph1 = t1 * phasePerSecond;

  traj.computeSwingTrajectoryBezier(ph0, 0.5);
  auto p0 = traj.getPosition();
  auto v0 = traj.getVelocity();
  auto a0 = traj.getAcceleration();


  traj.computeSwingTrajectoryBezier(ph1, 0.5);
  auto p1 = traj.getPosition();
  auto v1 = traj.getVelocity();
  auto a1 = traj.getAcceleration();

  Vec3<double> vdiff = (p1 - p0) / dt;
  Vec3<double> vref = (v0 + v1) / 2;
  Vec3<double> adiff = (v1 - v0) / dt;
  Vec3<double> aref = (a0 + a1) / 2;

  EXPECT_TRUE(almostEqual(adiff, aref, 0.001));
  EXPECT_TRUE(almostEqual(vdiff, vref, 0.001));
}