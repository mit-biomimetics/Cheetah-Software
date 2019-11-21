/*! @file Cheetah3.h
 *  @brief Utility function to build a Cheetah 3 Quadruped object
 *
 *  This file is based on Cheetah3FullRotorModel_mex.m (originally written by
 * Pat) and builds a model of the Cheetah 3 robot.  The inertia parameters of
 * legs and rotors were determined through an experimental procedure described
 * in "Linear Matrix Inequalities for Physically-Consistent Inertial Parameter
 *      Identification: A Statistical Perspective on the Mass Distribution" by
 * Wensing, Kim, Slotine. (see https://arxiv.org/abs/1701.04395)
 *
 *  It turns out that the parameters are not fully observable when the base is
 * fixed, which is described in "Observability in Inertial Parameter
 * Identification" by Wensing, Niemeyer, Slotine (see
 * https://arxiv.org/abs/1711.03896).
 *
 *  However, these estimates are still very good and were confirmed against a
 * CAD model of the robot.
 */

#ifndef LIBBIOMIMETICS_CHEETAH3_H
#define LIBBIOMIMETICS_CHEETAH3_H

#include "Dynamics/spatial.h"
#include "FloatingBaseModel.h"
#include "Quadruped.h"
#include "cppTypes.h"

using namespace spatial;

/*!
 * Generate a Quadruped model of Cheetah 3
 */
template <typename T>
Quadruped<T> buildCheetah3() {
  Quadruped<T> cheetah;
  cheetah._robotType = RobotType::CHEETAH_3;

  cheetah._bodyMass = 26.60;
  cheetah._bodyLength = .6;
  cheetah._bodyWidth = .256;
  cheetah._bodyHeight = .2;
  cheetah._abadGearRatio = 7.66666666667;
  cheetah._hipGearRatio = 7.66666666667;
  cheetah._kneeGearRatio = 8.846;
  cheetah._abadLinkLength = 0.045;
  cheetah._hipLinkLength = 0.342;
  cheetah._kneeLinkY_offset = 0.0;
  cheetah._kneeLinkLength = 0.345;
  cheetah._maxLegLength = 0.687;

  cheetah._batteryV = 65;
  cheetah._motorKT = 0.266;
  cheetah._motorR = 0.45;
  cheetah._jointDamping = .1;
  cheetah._jointDryFriction = 1;
  cheetah._motorTauMax = 27.2;

  MassProperties<T> abadMassProperties, hipMassProperties, kneeMassProperties,
      abadRotorMassProperties, hipRotorMassProperties, kneeRotorMassProperties;

  abadMassProperties << 1.645278752937798e+00, 3.987427689943171e-24,
      1.050694793045323e-01, -7.248117018543008e-03, 6.871761154771191e-03,
      2.715028710162546e-04, 6.820782113449669e-03, 4.545198314447007e-04,
      1.000360164583147e-25, -2.201422454848303e-25;

  hipMassProperties << 1.071200401412615e+00, -7.564034441661159e-04,
      -3.180527291037499e-02, -1.047496808654130e-01, 2.020976253260794e-02,
      1.856851411013688e-02, 3.115448656660994e-03, 1.294666102054298e-04,
      1.217274116646696e-03, -4.116533457993085e-04;

  kneeMassProperties << 9.019384095786234e-01, 2.215928611851693e-02,
      3.623331361414656e-03, -4.411934868569409e-02, 2.269554734974269e-02,
      2.714252268535285e-02, 5.566745484613744e-03, -3.238036214456428e-04,
      -2.194717054345885e-03, 1.274627163554897e-03;

  abadRotorMassProperties << 4.999939189070923e-01, 1.184093324655296e-24,
      1.976769699734076e-03, -3.681440077402916e-04, 9.908339941420835e-04,
      4.955917167197880e-04, 4.956200895280653e-04, 3.533262919985919e-08,
      5.751153538955154e-27, 1.342701882012297e-27;

  hipRotorMassProperties << 5.405045213157799e-01, -3.164547557293410e-03,
      -6.776800604201722e-03, 7.158586831736347e-04, 4.352341568156850e-04,
      1.015373256571613e-03, 9.126776037286181e-04, -3.097121515937353e-05,
      -2.154046949222089e-04, 1.306997274233568e-04;

  kneeRotorMassProperties << 5.398676189202086e-01, -2.571556872957194e-04,
      -7.550002224983261e-03, 1.125959428120149e-03, 7.035267590972641e-04,
      8.693209041255046e-04, 4.122075637011218e-04, -3.986925582339967e-05,
      -2.295396403157565e-04, 2.324567952553503e-05;

  cheetah._abadInertia = SpatialInertia<T>(abadMassProperties);
  cheetah._hipInertia = SpatialInertia<T>(hipMassProperties);
  cheetah._kneeInertia = SpatialInertia<T>(kneeMassProperties);
  cheetah._abadRotorInertia = SpatialInertia<T>(abadRotorMassProperties);
  cheetah._hipRotorInertia = SpatialInertia<T>(hipRotorMassProperties);
  cheetah._kneeRotorInertia = SpatialInertia<T>(kneeRotorMassProperties);
  Vec3<T> bodyCOM(0, 0, 0);
  Vec3<T> bodyDims(cheetah._bodyLength, cheetah._bodyWidth,
                   cheetah._bodyHeight);
  cheetah._bodyInertia = SpatialInertia<T>(
      cheetah._bodyMass, bodyCOM, rotInertiaOfBox(cheetah._bodyMass, bodyDims));

  // this doesn't generalize to the mini cheetah?
  cheetah._abadLocation =
      Vec3<T>(cheetah._bodyLength, cheetah._bodyWidth, 0) * 0.5;

  // note that this is wrong to match a small bug in the actual simulator, just
  // to test that I get the right answer. TODO fix!
  cheetah._abadRotorLocation = Vec3<T>(0, -cheetah._bodyWidth / 2., 0);

  cheetah._hipLocation = Vec3<T>(0, cheetah._abadLinkLength, 0);
  cheetah._hipRotorLocation = Vec3<T>(0, cheetah._abadLinkLength, 0);
  cheetah._kneeLocation = Vec3<T>(0, 0, -cheetah._hipLinkLength);
  cheetah._kneeRotorLocation = Vec3<T>(0, 0, 0);
  return cheetah;
}

#endif  // LIBBIOMIMETICS_CHEETAH3_H
