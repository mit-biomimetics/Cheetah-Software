/*!
 * @file DesiredStateCommand.h
 * @brief Logic to convert a joystick command into a desired trajectory for the robot
 *
 * This will generate a state trajectory which can easily be used for model predictive controllers
 */

/*========================= Gamepad Control ==========================*/
/**
 *
 */
#ifndef DESIRED_STATE_COMMAND_H
#define DESIRED_STATE_COMMAND_H

#include <iostream>

#include "Controllers/StateEstimatorContainer.h"
#include "cppTypes.h"

#include "SimUtilities/GamepadCommand.h"
#include "robot/include/rt/rt_rc_interface.h"

/**
 *
 */
template <typename T>
struct DesiredStateData {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  DesiredStateData() { zero(); }

  // Zero out all of the data
  void zero();

  // Instantaneous desired state command
  Vec12<T> stateDes;

  Vec12<T> pre_stateDes;

  // Desired future state trajectory (for up to 10 timestep MPC)
  Eigen::Matrix<T, 12, 10> stateTrajDes;
};

/**
 *
 */
template <typename T>
class DesiredStateCommand {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  // Initialize with the GamepadCommand struct
  DesiredStateCommand(GamepadCommand* command, rc_control_settings* rc_command,
                      RobotControlParameters* _parameters,
                      StateEstimate<T>* sEstimate, float _dt) {
    gamepadCommand = command;
    rcCommand = rc_command;
    stateEstimate = sEstimate;
    parameters = _parameters;

    data.stateDes.setZero();
    data.pre_stateDes.setZero();
    leftAnalogStick.setZero();
    rightAnalogStick.setZero();

    dt = _dt;
  }

  void convertToStateCommands();
  void setCommandLimits(T minVelX_in, T maxVelX_in,
                        T minVelY_in, T maxVelY_in, T minTurnRate_in, T maxTurnRate_in);
  void desiredStateTrajectory(int N, Vec10<T> dtVec);
  void printRawInfo();
  void printStateCommandInfo();
  float deadband(float command, T minVal, T maxVal);

  // These should come from the inferface
  T maxRoll = 0.4;
  T minRoll = -0.4;
  T maxPitch = 0.4;
  T minPitch = -0.4;
  T maxVelX = 3.0;
  T minVelX = -3.0;
  //T maxVelX = 5.0;
  //T minVelX = -5.0;
  T maxVelY = 2.0;
  T minVelY = -2.0;
  //T maxVelY = 0.5;
  //T minVelY = -0.5;
  T maxTurnRate = 2.5;
  T minTurnRate = -2.5;
  //T maxTurnRate = 0.5;
  //T minTurnRate = -0.5;


  Vec2<float> leftAnalogStick;
  Vec2<float> rightAnalogStick;

  // Holds the instantaneous desired state and future desired state trajectory
  DesiredStateData<T> data;

  const rc_control_settings* rcCommand;
  const GamepadCommand* gamepadCommand;

  bool trigger_pressed = false;

private:
  StateEstimate<T>* stateEstimate;
  RobotControlParameters* parameters;


  // Dynamics matrix for discrete time approximation
  Mat12<T> A;

  // Control loop timestep change
  T dt;

  // Value cutoff for the analog stick deadband
  T deadbandRegion = 0.075;
  //const T filter = 0.01;
  const T filter = 0.1;

  // Choose how often to print info, every N iterations
  int printNum = 5;  // N*(0.001s) in simulation time

  // Track the number of iterations since last info print
  int printIter = 0;
};

#endif
