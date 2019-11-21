/*! @file OrientationEstimator.cpp
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

#include "Controllers/OrientationEstimator.h"

/*!
 * Get quaternion, rotation matrix, angular velocity (body and world),
 * rpy, acceleration (world, body) by copying from cheater state data
 */
template <typename T>
void CheaterOrientationEstimator<T>::run() {
  this->_stateEstimatorData.result->orientation =
      this->_stateEstimatorData.cheaterState->orientation.template cast<T>();
  this->_stateEstimatorData.result->rBody = ori::quaternionToRotationMatrix(
      this->_stateEstimatorData.result->orientation);
  this->_stateEstimatorData.result->omegaBody =
      this->_stateEstimatorData.cheaterState->omegaBody.template cast<T>();
  this->_stateEstimatorData.result->omegaWorld =
      this->_stateEstimatorData.result->rBody.transpose() *
      this->_stateEstimatorData.result->omegaBody;
  this->_stateEstimatorData.result->rpy =
      ori::quatToRPY(this->_stateEstimatorData.result->orientation);
  this->_stateEstimatorData.result->aBody =
      this->_stateEstimatorData.cheaterState->acceleration.template cast<T>();
  this->_stateEstimatorData.result->aWorld =
      this->_stateEstimatorData.result->rBody.transpose() *
      this->_stateEstimatorData.result->aBody;
}

/*!
 * Get quaternion, rotation matrix, angular velocity (body and world),
 * rpy, acceleration (world, body) from vector nav IMU
 */
template <typename T>
void VectorNavOrientationEstimator<T>::run() {
  this->_stateEstimatorData.result->orientation[0] =
      this->_stateEstimatorData.vectorNavData->quat[3];
  this->_stateEstimatorData.result->orientation[1] =
      this->_stateEstimatorData.vectorNavData->quat[0];
  this->_stateEstimatorData.result->orientation[2] =
      this->_stateEstimatorData.vectorNavData->quat[1];
  this->_stateEstimatorData.result->orientation[3] =
      this->_stateEstimatorData.vectorNavData->quat[2];

  if(_b_first_visit){
    Vec3<T> rpy_ini = ori::quatToRPY(this->_stateEstimatorData.result->orientation);
    rpy_ini[0] = 0;
    rpy_ini[1] = 0;
    _ori_ini_inv = rpyToQuat(-rpy_ini);
    _b_first_visit = false;
  }
  this->_stateEstimatorData.result->orientation = 
    ori::quatProduct(_ori_ini_inv, this->_stateEstimatorData.result->orientation);

  this->_stateEstimatorData.result->rpy =
      ori::quatToRPY(this->_stateEstimatorData.result->orientation);


  this->_stateEstimatorData.result->rBody = ori::quaternionToRotationMatrix(
      this->_stateEstimatorData.result->orientation);

  this->_stateEstimatorData.result->omegaBody =
      this->_stateEstimatorData.vectorNavData->gyro.template cast<T>();

  this->_stateEstimatorData.result->omegaWorld =
      this->_stateEstimatorData.result->rBody.transpose() *
      this->_stateEstimatorData.result->omegaBody;

  this->_stateEstimatorData.result->aBody =
      this->_stateEstimatorData.vectorNavData->accelerometer.template cast<T>();
  this->_stateEstimatorData.result->aWorld =
      this->_stateEstimatorData.result->rBody.transpose() *
      this->_stateEstimatorData.result->aBody;
}


template class CheaterOrientationEstimator<float>;
template class CheaterOrientationEstimator<double>;

template class VectorNavOrientationEstimator<float>;
template class VectorNavOrientationEstimator<double>;
