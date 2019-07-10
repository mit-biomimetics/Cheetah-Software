/*========================= Gamepad Control ==========================*/
/**
 *
 */
#ifndef DESIRED_STATE_COMMAND_H
#define DESIRED_STATE_COMMAND_H

#include <Controllers/StateEstimatorContainer.h>
#include <cppTypes.h>
#include <iostream>
#include "SimUtilities/GamepadCommand.h"

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
  DesiredStateCommand(GamepadCommand* command, StateEstimate<T>* sEstimate) {
    gamepadCommand = command;
    stateEstimate = sEstimate;
  }

  void convertToStateCommands();
  void desiredStateTrajectory(int N, Vec10<T> dtVec);
  void printRawInfo();
  void printStateCommandInfo();
  float deadband(float command, T minVal, T maxVal);

  // These should come from the inferface
  T maxRoll = 0.4;
  T minRoll = -0.4;
  T maxPitch = 0.4;
  T minPitch = -0.4;
  T maxVelX = 1.0;
  T minVelX = -1.0;
  T maxVelY = 0.5;
  T minVelY = -0.5;
  T maxTurnRate = 2.0;
  T minTurnRate = -2.0;

  // Holds the instantaneous desired state and future desired state trajectory
  DesiredStateData<T> data;

 private:
  GamepadCommand* gamepadCommand;
  StateEstimate<T>* stateEstimate;

  // Dynamics matrix for discrete time approximation
  Mat12<T> A;

  // Control loop timestep change
  T dt = 0.001;

  // Value cutoff for the analog stick deadband
  T deadbandRegion = 0.075;

  // Choose how often to print info, every N iterations
  int printNum = 5;  // N*(0.001s) in simulation time

  // Track the number of iterations since last info print
  int printIter = 0;
};

#endif