/*! @file SimulationBridge.h
 *  @brief  The SimulationBridge runs a RobotController and connects it to a
 * Simulator, using shared memory. It is the simulation version of the
 * HardwareBridge.
 */

#ifndef PROJECT_SIMULATIONDRIVER_H
#define PROJECT_SIMULATIONDRIVER_H

#include <thread>

#include "ControlParameters/RobotParameters.h"
#include "RobotRunner.h"
#include "SimUtilities/SimulatorMessage.h"
#include "Types.h"
#include "Utilities/PeriodicTask.h"
#include "Utilities/SharedMemory.h"

class SimulationBridge {
 public:
  explicit SimulationBridge(RobotType robot, RobotController* robot_ctrl) : 
    _robot(robot) {
     _fakeTaskManager = new PeriodicTaskManager;
    _robotRunner = new RobotRunner(robot_ctrl, _fakeTaskManager, 0, "robot-task");
    _userParams = robot_ctrl->getUserControlParameters();

 }
  void run();
  void handleControlParameters();
  void runRobotControl();
  ~SimulationBridge() {
    delete _fakeTaskManager;
    delete _robotRunner;
  }
  void run_sbus();

 private:
  PeriodicTaskManager taskManager;
  bool _firstControllerRun = true;
  PeriodicTaskManager* _fakeTaskManager = nullptr;
  RobotType _robot;
  RobotRunner* _robotRunner = nullptr;
  SimulatorMode _simMode;
  SharedMemoryObject<SimulatorSyncronizedMessage> _sharedMemory;
  RobotControlParameters _robotParams;
  ControlParameters* _userParams = nullptr;
  u64 _iterations = 0;

  std::thread* sbus_thread;
};

#endif  // PROJECT_SIMULATIONDRIVER_H
