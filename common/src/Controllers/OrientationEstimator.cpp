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

#include "Controllers/OrientationEstimator.h"

template<typename T>
void CheaterOrientationEstimator<T>::run()  {
  this->_stateEstimatorData.result->orientation = this->_stateEstimatorData.cheaterState->orientation.template cast<T>();
  this->_stateEstimatorData.result->rBody = ori::quaternionToRotationMatrix(this->_stateEstimatorData.result->orientation);
  this->_stateEstimatorData.result->omegaBody = this->_stateEstimatorData.cheaterState->omegaBody.template cast<T>();
  this->_stateEstimatorData.result->omegaWorld = this->_stateEstimatorData.result->rBody.transpose() * this->_stateEstimatorData.result->omegaBody;
  this->_stateEstimatorData.result->rpy = ori::quatToRPY(this->_stateEstimatorData.result->orientation);
  this->_stateEstimatorData.result->aBody = this->_stateEstimatorData.cheaterState->acceleration.template cast<T>();
  this->_stateEstimatorData.result->aWorld = this->_stateEstimatorData.result->rBody.transpose() * this->_stateEstimatorData.result->aBody;
}

template<typename T>
void VectorNavOrientationEstimator<T>::run() {
  this->_stateEstimatorData.result->orientation[0] = 
      this->_stateEstimatorData.vectorNavData->quat[3];
  this->_stateEstimatorData.result->orientation[1] = 
      this->_stateEstimatorData.vectorNavData->quat[0];
  this->_stateEstimatorData.result->orientation[2] = 
      this->_stateEstimatorData.vectorNavData->quat[1];
  this->_stateEstimatorData.result->orientation[3] = 
      this->_stateEstimatorData.vectorNavData->quat[2];
  this->_stateEstimatorData.result->rBody = ori::quaternionToRotationMatrix(this->_stateEstimatorData.result->orientation);
  this->_stateEstimatorData.result->omegaBody = this->_stateEstimatorData.vectorNavData->gyro.template cast <T>();
  this->_stateEstimatorData.result->omegaWorld = this->_stateEstimatorData.result->rBody.transpose() * this->_stateEstimatorData.result->omegaBody;
  this->_stateEstimatorData.result->rpy = ori::quatToRPY(this->_stateEstimatorData.result->orientation);
  this->_stateEstimatorData.result->aBody =
      this->_stateEstimatorData.vectorNavData->accelerometer.template cast<T>();
  this->_stateEstimatorData.result->aWorld =
      this->_stateEstimatorData.result->rBody.transpose() * this->_stateEstimatorData.result->aBody;
}

template<typename T>
void KvhOrientationEstimator<T>::run() {
  // Need to get the quaternion estimate
  //this->_stateEstimatorData.result->orientation = this->_stateEstimatorData.vectorNavData->quat.template cast<T>();
  

  //this->_stateEstimatorData.result->rBody = ori::quaternionToRotationMatrix(this->_stateEstimatorData.result->orientation);
  //this->_stateEstimatorData.result->omegaBody = this->_stateEstimatorData.kvhImuData->gyro.template cast <T>();
  //this->_stateEstimatorData.result->omegaWorld = this->_stateEstimatorData.result->rBody.transpose() * this->_stateEstimatorData.result->omegaBody;
  //this->_stateEstimatorData.result->rpy = ori::quatToRPY(this->_stateEstimatorData.result->orientation);
}

template class CheaterOrientationEstimator<float>;
template class CheaterOrientationEstimator<double>;

template class VectorNavOrientationEstimator<float>;
template class VectorNavOrientationEstimator<double>;

template class KvhOrientationEstimator<float>;
template class KvhOrientationEstimator<double>;
