/*!
 * @file HardwareBridge.h
 * @brief Interface between robot code and robot hardware
 *
 * This class initializes the hardware of both robots and allows the robot controller to access it
 */

#ifndef PROJECT_HARDWAREBRIDGE_H
#define PROJECT_HARDWAREBRIDGE_H


#define MAX_STACK_SIZE 16384 // 16KB  of stack
#define TASK_PRIORITY 49     // linux priority, this is not the nice value

#include <string>
#include <lcm-cpp.hpp>
#include "Utilities/PeriodicTask.h"
#include "RobotController.h"
#include "gamepad_lcmt.hpp"
#include "control_parameter_respones_lcmt.hpp"
#include "control_parameter_request_lcmt.hpp"


class HardwareBridge {
public:
  HardwareBridge() : statusTask(&taskManager, 0.5f), _interfaceLCM(getLcmUrl(255)),
  _visualizationLCM(getLcmUrl(255)) { }
  void prefaultStack();
  void setupScheduler();
  void initError(const char* reason, bool printErrno = false);
  void initCommon();
  ~HardwareBridge() {
    delete _robotController;
  }
  void handleGamepadLCM(const lcm::ReceiveBuffer* rbuf,
                        const std::string& chan,
                        const gamepad_lcmt* msg);

  void handleInterfaceLCM();
  void handleControlParameter(const lcm::ReceiveBuffer* rbuf,
                              const std::string& chan,
                              const control_parameter_request_lcmt* msg);

  void publishVisualizationLCM();
  void run_sbus();

protected:
  PeriodicTaskManager taskManager;
  PrintTaskStatus     statusTask;
  GamepadCommand      _gamepadCommand;
  VisualizationData   _visualizationData;
  CheetahVisualization _mainCheetahVisualization;
  lcm::LCM            _interfaceLCM;
  lcm::LCM            _visualizationLCM;
  control_parameter_respones_lcmt _parameter_response_lcmt;
  SpiData _spiData;
  SpiCommand _spiCommand;


  bool _firstRun = true;
  RobotController* _robotController = nullptr;
  RobotControlParameters _robotParams;
  u64 _iterations = 0;
  std::thread _interfaceLcmThread;
  volatile bool _interfaceLcmQuit = false;

  int _port;
};

class MiniCheetahHardwareBridge : public HardwareBridge {
public:
  MiniCheetahHardwareBridge();
  void runSpi();
  void initHardware();
  void run();
  void abort(const std::string& reason);
  void abort(const char* reason);

private:
  VectorNavData _vectorNavData;
  lcm::LCM _spiLcm;
};

#endif //PROJECT_HARDWAREBRIDGE_H
