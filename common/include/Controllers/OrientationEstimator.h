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

template <typename T>
class CheaterOrientationEstimator : public GenericEstimator<T> {
 public:
  virtual void run();
  virtual void setup() {}
};

template <typename T>
class VectorNavOrientationEstimator : public GenericEstimator<T> {
 public:
  virtual void run();
  virtual void setup() {}
};


#endif  // PROJECT_ORIENTATIONESTIMATOR_H
