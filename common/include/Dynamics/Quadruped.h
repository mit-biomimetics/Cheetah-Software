/*! @file Quadruped.h
 *  @brief Data structure containing parameters for quadruped robot
 *
 *  This file contains the Quadruped class.  This stores all the parameters for
 * a quadruped robot.  There are utility functions to generate Quadruped objects
 * for Cheetah 3 (and eventually mini-cheetah). There is a buildModel() method
 * which can be used to create a floating-base dynamics model of the quadruped.
 */

#ifndef LIBBIOMIMETICS_QUADRUPED_H
#define LIBBIOMIMETICS_QUADRUPED_H

#include "Dynamics/ActuatorModel.h"
#include "Dynamics/FloatingBaseModel.h"
#include "Dynamics/SpatialInertia.h"

#include <eigen3/Eigen/StdVector>

#include <vector>

namespace cheetah {
constexpr size_t num_act_joint = 12;
constexpr size_t num_q = 19;
constexpr size_t dim_config = 18;
constexpr size_t num_leg = 4;
constexpr size_t num_leg_joint = 3;
constexpr float servo_rate = 0.001;
}  // namespace cheetah

namespace linkID {
constexpr size_t FR = 9;   // Front Right Foot
constexpr size_t FL = 11;  // Front Left Foot
constexpr size_t HR = 13;  // Hind Right Foot
constexpr size_t HL = 15;  // Hind Left Foot

constexpr size_t FR_abd = 2;  // Front Right Abduction
constexpr size_t FL_abd = 0;  // Front Left Abduction
constexpr size_t HR_abd = 3;  // Hind Right Abduction
constexpr size_t HL_abd = 1;  // Hind Left Abduction
}  // namespace linkID

using std::vector;

/*!
 * Representation of a quadruped robot's physical properties.
 *
 * When viewed from the top, the quadruped's legs are:
 *
 * FRONT
 * 2 1   RIGHT
 * 4 3
 * BACK
 *
 */
template <typename T>
class Quadruped {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  RobotType _robotType;
  T _bodyLength, _bodyWidth, _bodyHeight, _bodyMass;
  T _abadGearRatio, _hipGearRatio, _kneeGearRatio;
  T _abadLinkLength, _hipLinkLength, _kneeLinkLength, _maxLegLength;
  T _motorKT, _motorR, _batteryV;
  T _motorTauMax;
  T _jointDamping, _jointDryFriction;
  SpatialInertia<T> _abadInertia, _hipInertia, _kneeInertia, _abadRotorInertia,
      _hipRotorInertia, _kneeRotorInertia, _bodyInertia;
  Vec3<T> _abadLocation, _abadRotorLocation, _hipLocation, _hipRotorLocation,
      _kneeLocation, _kneeRotorLocation;
  FloatingBaseModel<T> buildModel();
  bool buildModel(FloatingBaseModel<T>& model);
  std::vector<ActuatorModel<T>> buildActuatorModels();

  static T getSideSign(int leg) {
    const T sideSigns[4] = {-1, 1, -1, 1};
    assert(leg >= 0 && leg < 4);
    return sideSigns[leg];
  }

  /*!
   * Get location of the hip for the given leg in robot frame
   */
  Vec3<T> getHipLocation(int leg) {
    assert(leg >= 0 && leg < 4);
    Vec3<T> pHip((leg == 0 || leg == 1) ? _abadLocation(0) : -_abadLocation(0),
                 (leg == 1 || leg == 3) ? _abadLocation(1) : -_abadLocation(1),
                 _abadLocation(2));
    return pHip;
  }
};

template <typename T, typename T2>
Vec3<T> withLegSigns(const Eigen::MatrixBase<T2>& v, int legID);

#endif  // LIBBIOMIMETICS_QUADRUPED_H
