/*! @file SimulatorMessage.h
 *  @brief Messages sent to/from the development simulator
 *
 *  These messsages contain all data that is exchanged between the robot program
 * and the simulator using shared memory.   This is basically everything except
 * for debugging logs, which are handled by LCM instead
 */

#ifndef PROJECT_SIMULATORTOROBOTMESSAGE_H
#define PROJECT_SIMULATORTOROBOTMESSAGE_H

#include "ControlParameters/ControlParameterInterface.h"
#include "SimUtilities/GamepadCommand.h"
#include "SimUtilities/IMUTypes.h"
#include "SimUtilities/SpineBoard.h"
#include "SimUtilities/VisualizationData.h"
#include "SimUtilities/ti_boardcontrol.h"
#include "Utilities/SharedMemory.h"

/*!
 * The mode for the simulator
 */
enum class SimulatorMode {
  RUN_CONTROL_PARAMETERS,  // don't run the robot controller, just process
                           // Control Parameters
  RUN_CONTROLLER,          // run the robot controller
  DO_NOTHING,              // just to check connection
  EXIT                     // quit!
};

/*!
 * A plain message from the simulator to the robot
 */
struct SimulatorToRobotMessage {
  GamepadCommand gamepadCommand;  // joystick
  RobotType robotType;  // which robot the simulator thinks we are simulating

  // imu data
  VectorNavData vectorNav;
  CheaterState<double> cheaterState;

  // leg data
  SpiData spiData;
  TiBoardData tiBoardData[4];
  // todo cheetah 3
  ControlParameterRequest controlParameterRequest;

  SimulatorMode mode;
};

/*!
 * A plain message from the robot to the simulator
 */
struct RobotToSimulatorMessage {
  RobotType robotType;
  SpiCommand spiCommand;
  TiBoardCommand tiBoardCommand[4];

  VisualizationData visualizationData;
  CheetahVisualization mainCheetahVisualization;
  ControlParameterResponse controlParameterResponse;
};

/*!
 * All the data shared between the robot and the simulator
 */
struct SimulatorMessage {
  RobotToSimulatorMessage robotToSim;
  SimulatorToRobotMessage simToRobot;
};

/*!
 * A SimulatorSyncronizedMessage is stored in shared memory and is accessed by
 * both the simulator and the robot The simulator and robot take turns have
 * exclusive access to the entire message. The intended sequence is:
 *  - robot: waitForSimulator()
 *  - simulator: *simulates robot* (simulator can read/write, robot cannot do
 * anything)
 *  - simulator: simDone()
 *  - simulator: waitForRobot()
 *  - robot: *runs controller*    (robot can read/write, simulator cannot do
 * anything)
 *  - robot: robotDone();
 *  - robot: waitForSimulator()
 *  ...
 */
struct SimulatorSyncronizedMessage : public SimulatorMessage {
  /*!
   * The init() method should only be called *after* shared memory is connected!
   */
  void init() {
    robotToSimSemaphore.init(0);
    simToRobotSemaphore.init(0);
  }

  void waitForSimulator() { simToRobotSemaphore.decrement(); }

  void simulatorIsDone() { simToRobotSemaphore.increment(); }

  void waitForRobot() { robotToSimSemaphore.decrement(); }

  bool tryWaitForRobot() { return robotToSimSemaphore.tryDecrement(); }

  bool waitForRobotWithTimeout() {
    return robotToSimSemaphore.decrementTimeout(1, 0);
  }

  void robotIsDone() { robotToSimSemaphore.increment(); }

 private:
  SharedMemorySemaphore robotToSimSemaphore, simToRobotSemaphore;
};

#endif  // PROJECT_SIMULATORTOROBOTMESSAGE_H
