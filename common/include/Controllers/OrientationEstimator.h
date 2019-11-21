/*! @file OrientationEstimator.h
 *  @brief All Orientation Estimation Algorithms
 *
 *  This file will contain all orientation algorithms.
 *  Orientation estimators should compute:
 *  - orientation: a quaternion representing orientation
 *  - rBody: coordinate transformation matrix (satisfies vBody = Rbody * vWorld)
 *  - omegaBody: angular velocity in body frame
 *  - omegaWorld: angular velocity in world frame
 *  - rpy: roll pitch yaw
 */
#ifndef PROJECT_ORIENTATIONESTIMATOR_H
#define PROJECT_ORIENTATIONESTIMATOR_H

#include "Controllers/StateEstimatorContainer.h"

/*!
 * "Cheater" estimator for orientation which always returns the correct value in simulation
 */
template <typename T>
class CheaterOrientationEstimator : public GenericEstimator<T> {
 public:
  virtual void run();
  virtual void setup() {}
};

/*!
 * Estimator for the VectorNav IMU.  The VectorNav provides an orientation already and
 * we just return that.
 */
template <typename T>
class VectorNavOrientationEstimator : public GenericEstimator<T> {
 public:
  virtual void run();
  virtual void setup() {}
  
 protected:
  bool _b_first_visit = true;
  Quat<T> _ori_ini_inv;
};


#endif  // PROJECT_ORIENTATIONESTIMATOR_H
