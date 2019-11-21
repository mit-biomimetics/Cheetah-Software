/*! @file SimulatorParameters.cpp
 *  @brief Declaration of various simulator parameters
 *
 *  This class contains all the ControlParameters for the simulator.
 *  In most cases, the simulator just loads the control parameters in
 * simulator-defaults.ini and this is okay
 */

#ifndef PROJECT_SIMULATORPARAMETERS_H
#define PROJECT_SIMULATORPARAMETERS_H

#include "ControlParameters/ControlParameters.h"

#define SIMULATOR_DEFAULT_PARAMETERS "/simulator-defaults.yaml"
#define MINI_CHEETAH_DEFAULT_PARAMETERS "/mini-cheetah-defaults.yaml"
#define CHEETAH_3_DEFAULT_PARAMETERS "/cheetah-3-defaults.yaml"

/*!
 * Simulator specific control parameters
 */
class SimulatorControlParameters : public ControlParameters {
 public:
  SimulatorControlParameters()
      : ControlParameters("simulator-parameters"),
        INIT_PARAMETER(vectornav_imu_accelerometer_noise),
        INIT_PARAMETER(vectornav_imu_gyro_noise),
        INIT_PARAMETER(vectornav_imu_quat_noise),
        INIT_PARAMETER(game_controller_deadband),
        INIT_PARAMETER(simulation_speed),
        INIT_PARAMETER(simulation_paused),
        INIT_PARAMETER(high_level_dt),
        INIT_PARAMETER(low_level_dt),
        INIT_PARAMETER(dynamics_dt),
        INIT_PARAMETER(floor_kp),
        INIT_PARAMETER(floor_kd),
        INIT_PARAMETER(use_spring_damper),
        INIT_PARAMETER(sim_state_lcm),
        INIT_PARAMETER(sim_lcm_ttl),
        INIT_PARAMETER(go_home),
        INIT_PARAMETER(home_pos),
        INIT_PARAMETER(home_rpy),
        INIT_PARAMETER(home_kp_lin),
        INIT_PARAMETER(home_kd_lin),
        INIT_PARAMETER(home_kp_ang),
        INIT_PARAMETER(home_kd_ang) {}

  DECLARE_PARAMETER(float, vectornav_imu_accelerometer_noise)
  DECLARE_PARAMETER(float, vectornav_imu_gyro_noise)
  DECLARE_PARAMETER(float, vectornav_imu_quat_noise)

  DECLARE_PARAMETER(float, game_controller_deadband)

  DECLARE_PARAMETER(double, simulation_speed)
  DECLARE_PARAMETER(s64, simulation_paused)
  DECLARE_PARAMETER(double, high_level_dt)
  DECLARE_PARAMETER(double, low_level_dt)
  DECLARE_PARAMETER(double, dynamics_dt)

  DECLARE_PARAMETER(double, floor_kp)
  DECLARE_PARAMETER(double, floor_kd)
  DECLARE_PARAMETER(s64, use_spring_damper)
  DECLARE_PARAMETER(s64, sim_state_lcm)
  DECLARE_PARAMETER(s64, sim_lcm_ttl)

  DECLARE_PARAMETER(s64, go_home)
  DECLARE_PARAMETER(Vec3<double>,home_pos)
  DECLARE_PARAMETER(Vec3<double>,home_rpy)
  DECLARE_PARAMETER(double, home_kp_lin)
  DECLARE_PARAMETER(double, home_kd_lin)
  DECLARE_PARAMETER(double, home_kp_ang)
  DECLARE_PARAMETER(double, home_kd_ang)

};

#endif  // PROJECT_SIMULATORPARAMETERS_H
