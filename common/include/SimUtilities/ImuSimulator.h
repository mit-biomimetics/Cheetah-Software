/*! @file ImuSimulator.h
 *  @brief Simulated IMU
 */

#ifndef PROJECT_IMUSIMULATOR_H
#define PROJECT_IMUSIMULATOR_H

#include "cppTypes.h"
#include "Dynamics/FloatingBaseModel.h"
#include "ControlParameters/SimulatorParameters.h"
#include "SimUtilities/IMUTypes.h"

#include <random>



template <typename T>
class ImuSimulator {
public:

  explicit ImuSimulator(SimulatorControlParameters& simSettings, u64 seed = 0) :
    _simSettings(simSettings),
    _mt(seed),
    _kvhGyroDistribution(-simSettings.kvh_imu_gyro_noise, simSettings.kvh_imu_gyro_noise),
    _kvhAccelerometerDistribution(-simSettings.kvh_imu_accelerometer_noise, simSettings.kvh_imu_gyro_noise),
    _vectornavGyroDistribution(
            -simSettings.vectornav_imu_gyro_noise, 
            simSettings.vectornav_imu_gyro_noise),
    _vectornavAccelerometerDistribution(-simSettings.vectornav_imu_accelerometer_noise, simSettings.vectornav_imu_accelerometer_noise),
    _vectornavQuatDistribution(-simSettings.vectornav_imu_quat_noise, simSettings.vectornav_imu_quat_noise) {
    if(simSettings.vectornav_imu_quat_noise != 0) {
      _vectorNavOrientationNoise = true;
    }
  }

  void updateVectornav(
          const FBModelState<T>& robotState, 
          const FBModelStateDerivative<T>& robotStateD, 
          VectorNavData* data);

  void updateKVH(const FBModelState<T>& robotState, const FBModelStateDerivative<T>& robotStateD, KvhImuData* data);

  void computeAcceleration(const FBModelState<T>& robotState, const FBModelStateDerivative<T>& robotStateD,
          Vec3<float>& acc, std::uniform_real_distribution<float>& dist, const RotMat<float>& R_body);

  void updateCheaterState(const FBModelState<T> &robotState, const FBModelStateDerivative<T>& robotStateD, CheaterState<T>& state);

private:
  SimulatorControlParameters& _simSettings;
  std::mt19937 _mt;
  std::uniform_real_distribution<float> _kvhGyroDistribution;
  std::uniform_real_distribution<float> _kvhAccelerometerDistribution;
  std::uniform_real_distribution<float> _vectornavGyroDistribution;
  std::uniform_real_distribution<float> _vectornavAccelerometerDistribution;
  std::uniform_real_distribution<float> _vectornavQuatDistribution;
  bool _vectorNavOrientationNoise = false;
};
#endif //PROJECT_IMUSIMULATOR_H