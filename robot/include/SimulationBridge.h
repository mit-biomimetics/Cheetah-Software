/*! @file SimulatorDriver.h
 *  @brief  The SimulatorDriver runs a RobotController and connects it to a Simulator, using shared memory.
 *  It is the simulation version of the HardwareDriver.
 */

#ifndef PROJECT_SIMULATIONDRIVER_H
#define PROJECT_SIMULATIONDRIVER_H

#include "Types.h"
#include "Utilities/SharedMemory.h"
#include "SimUtilities/SimulatorMessage.h"
#include "ControlParameters/RobotParameters.h"
#include "RobotController.h"
#include "Utilities/PeriodicTask.h"
#include <thread>

class SimulationBridge {
  public:
    explicit SimulationBridge(RobotType robot) : _robot(robot) { }
    void run();
    void handleControlParameters();
    void runRobotControl();
    ~SimulationBridge() {
      delete _fakeTaskManager;
      delete _robotController;
    }
    void run_sbus();

  private:
    PeriodicTaskManager taskManager;
    bool _firstControllerRun = true;
    PeriodicTaskManager* _fakeTaskManager = nullptr;
    RobotType _robot;
    RobotController* _robotController = nullptr;
    SimulatorMode _simMode;
    SharedMemoryObject<SimulatorSyncronizedMessage> _sharedMemory;
    RobotControlParameters _robotParams;
    u64 _iterations = 0;

    std::thread* sbus_thread;
};


#endif //PROJECT_SIMULATIONDRIVER_H
