/*! @file RobotParameters.cpp
 *  @brief Declaration of various robot parameters
 *
 *  This class contains all the ControlParameters for the robot.
 */

#ifndef PROJECT_ROBOTPARAMETERS_H
#define PROJECT_ROBOTPARAMETERS_H

#include "ControlParameters/ControlParameters.h"

class RobotControlParameters : public ControlParameters {
 public:
  RobotControlParameters()
      : ControlParameters("robot-parameters"),
        INIT_PARAMETER(myValue),
        INIT_PARAMETER(control_mode),
        INIT_PARAMETER(testValue),
        INIT_PARAMETER(controller_dt),
        INIT_PARAMETER(stand_kp_cartesian),
        INIT_PARAMETER(stand_kd_cartesian),
        INIT_PARAMETER(kpCOM),
        INIT_PARAMETER(kdCOM),
        INIT_PARAMETER(kpBase),
        INIT_PARAMETER(kdBase),
        INIT_PARAMETER(cheater_mode),
        INIT_PARAMETER(imu_process_noise_position),
        INIT_PARAMETER(imu_process_noise_velocity),
        INIT_PARAMETER(foot_process_noise_position),
        INIT_PARAMETER(foot_sensor_noise_position),
        INIT_PARAMETER(foot_sensor_noise_velocity),
        INIT_PARAMETER(foot_height_sensor_noise) {}

  DECLARE_PARAMETER(double, myValue)
  DECLARE_PARAMETER(double, control_mode)
  DECLARE_PARAMETER(double, testValue)
  DECLARE_PARAMETER(double, controller_dt)
  DECLARE_PARAMETER(Vec3<double>, stand_kp_cartesian)
  DECLARE_PARAMETER(Vec3<double>, stand_kd_cartesian)
  DECLARE_PARAMETER(Vec3<double>, kpCOM)
  DECLARE_PARAMETER(Vec3<double>, kdCOM)
  DECLARE_PARAMETER(Vec3<double>, kpBase)
  DECLARE_PARAMETER(Vec3<double>, kdBase)

  // state estimator
  DECLARE_PARAMETER(s64, cheater_mode)
  DECLARE_PARAMETER(double, imu_process_noise_position)
  DECLARE_PARAMETER(double, imu_process_noise_velocity)
  DECLARE_PARAMETER(double, foot_process_noise_position)
  DECLARE_PARAMETER(double, foot_sensor_noise_position)
  DECLARE_PARAMETER(double, foot_sensor_noise_velocity)
  DECLARE_PARAMETER(double, foot_height_sensor_noise)
};

#endif  // PROJECT_ROBOTPARAMETERS_H
