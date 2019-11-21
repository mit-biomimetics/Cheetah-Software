/*!
 * @file HardwareBridge.cpp
 * @brief Interface between robot code and robot hardware
 *
 * This class initializes the hardware of both robots and allows the robot
 * controller to access it
 */
#ifdef linux 

#include <sys/mman.h>
#include <unistd.h>
#include <cstring>
#include <thread>
#include "Configuration.h"

#include "HardwareBridge.h"
#include "rt/rt_interface_lcm.h"
#include "rt/rt_sbus.h"
#include "rt/rt_spi.h"
#include "rt/rt_vectornav.h"

#define USE_MICROSTRAIN

/*!
 * If an error occurs during initialization, before motors are enabled, print
 * error and exit.
 * @param reason Error message string
 * @param printErrno If true, also print C errno
 */
void HardwareBridge::initError(const char* reason, bool printErrno) {
  printf("FAILED TO INITIALIZE HARDWARE: %s\n", reason);

  if (printErrno) {
    printf("Error: %s\n", strerror(errno));
  }

  exit(-1);
}

/*!
 * All hardware initialization steps that are common between Cheetah 3 and Mini Cheetah
 */
void HardwareBridge::initCommon() {
  printf("[HardwareBridge] Init stack\n");
  prefaultStack();
  printf("[HardwareBridge] Init scheduler\n");
  setupScheduler();
  if (!_interfaceLCM.good()) {
    initError("_interfaceLCM failed to initialize\n", false);
  }

  printf("[HardwareBridge] Subscribe LCM\n");
  _interfaceLCM.subscribe("interface", &HardwareBridge::handleGamepadLCM, this);
  _interfaceLCM.subscribe("interface_request",
                          &HardwareBridge::handleControlParameter, this);

  printf("[HardwareBridge] Start interface LCM handler\n");
  _interfaceLcmThread = std::thread(&HardwareBridge::handleInterfaceLCM, this);
}

/*!
 * Run interface LCM
 */
void HardwareBridge::handleInterfaceLCM() {
  while (!_interfaceLcmQuit) _interfaceLCM.handle();
}

/*!
 * Writes to a 16 KB buffer on the stack. If we are using 4K pages for our
 * stack, this will make sure that we won't have a page fault when the stack
 * grows.  Also mlock's all pages associated with the current process, which
 * prevents the cheetah software from being swapped out.  If we do run out of
 * memory, the robot program will be killed by the OOM process killer (and
 * leaves a log) instead of just becoming unresponsive.
 */
void HardwareBridge::prefaultStack() {
  printf("[Init] Prefault stack...\n");
  volatile char stack[MAX_STACK_SIZE];
  memset(const_cast<char*>(stack), 0, MAX_STACK_SIZE);
  if (mlockall(MCL_CURRENT | MCL_FUTURE) == -1) {
    initError(
        "mlockall failed.  This is likely because you didn't run robot as "
        "root.\n",
        true);
  }
}

/*!
 * Configures the scheduler for real time priority
 */
void HardwareBridge::setupScheduler() {
  printf("[Init] Setup RT Scheduler...\n");
  struct sched_param params;
  params.sched_priority = TASK_PRIORITY;
  if (sched_setscheduler(0, SCHED_FIFO, &params) == -1) {
    initError("sched_setscheduler failed.\n", true);
  }
}

/*!
 * LCM Handler for gamepad message
 */
void HardwareBridge::handleGamepadLCM(const lcm::ReceiveBuffer* rbuf,
                                      const std::string& chan,
                                      const gamepad_lcmt* msg) {
  (void)rbuf;
  (void)chan;
  _gamepadCommand.set(msg);
}

/*!
 * LCM Handler for control parameters
 */
void HardwareBridge::handleControlParameter(
    const lcm::ReceiveBuffer* rbuf, const std::string& chan,
    const control_parameter_request_lcmt* msg) {
  (void)rbuf;
  (void)chan;
  if (msg->requestNumber <= _parameter_response_lcmt.requestNumber) {
    // nothing to do!
    printf(
        "[HardwareBridge] Warning: the interface has run a ControlParameter "
        "iteration, but there is no new request!\n");
    // return;
  }

  // sanity check
  s64 nRequests = msg->requestNumber - _parameter_response_lcmt.requestNumber;
  if (nRequests != 1) {
    printf("[ERROR] Hardware bridge: we've missed %ld requests\n",
           nRequests - 1);
  }

  switch (msg->requestKind) {
    case (s8)ControlParameterRequestKind::SET_USER_PARAM_BY_NAME: {
      if(!_userControlParameters) {
        printf("[Warning] Got user param %s, but not using user parameters!\n",
               (char*)msg->name);
      } else {
        std::string name((char*)msg->name);
        ControlParameter& param = _userControlParameters->collection.lookup(name);

        // type check
        if ((s8)param._kind != msg->parameterKind) {
          throw std::runtime_error(
              "type mismatch for parameter " + name + ", robot thinks it is " +
              controlParameterValueKindToString(param._kind) +
              " but received a command to set it to " +
              controlParameterValueKindToString(
                  (ControlParameterValueKind)msg->parameterKind));
        }

        // do the actual set
        ControlParameterValue v;
        memcpy(&v, msg->value, sizeof(v));
        param.set(v, (ControlParameterValueKind)msg->parameterKind);

        // respond:
        _parameter_response_lcmt.requestNumber =
            msg->requestNumber;  // acknowledge that the set has happened
        _parameter_response_lcmt.parameterKind =
            msg->parameterKind;  // just for debugging print statements
        memcpy(_parameter_response_lcmt.value, msg->value, 64);
        //_parameter_response_lcmt.value = _parameter_request_lcmt.value; // just
        //for debugging print statements
        strcpy((char*)_parameter_response_lcmt.name,
               name.c_str());  // just for debugging print statements
        _parameter_response_lcmt.requestKind = msg->requestKind;

        printf("[User Control Parameter] set %s to %s\n", name.c_str(),
               controlParameterValueToString(
                   v, (ControlParameterValueKind)msg->parameterKind)
                   .c_str());
      }
    } break;

    case (s8)ControlParameterRequestKind::SET_ROBOT_PARAM_BY_NAME: {
      std::string name((char*)msg->name);
      ControlParameter& param = _robotParams.collection.lookup(name);

      // type check
      if ((s8)param._kind != msg->parameterKind) {
        throw std::runtime_error(
            "type mismatch for parameter " + name + ", robot thinks it is " +
            controlParameterValueKindToString(param._kind) +
            " but received a command to set it to " +
            controlParameterValueKindToString(
                (ControlParameterValueKind)msg->parameterKind));
      }

      // do the actual set
      ControlParameterValue v;
      memcpy(&v, msg->value, sizeof(v));
      param.set(v, (ControlParameterValueKind)msg->parameterKind);

      // respond:
      _parameter_response_lcmt.requestNumber =
          msg->requestNumber;  // acknowledge that the set has happened
      _parameter_response_lcmt.parameterKind =
          msg->parameterKind;  // just for debugging print statements
      memcpy(_parameter_response_lcmt.value, msg->value, 64);
      //_parameter_response_lcmt.value = _parameter_request_lcmt.value; // just
      //for debugging print statements
      strcpy((char*)_parameter_response_lcmt.name,
             name.c_str());  // just for debugging print statements
      _parameter_response_lcmt.requestKind = msg->requestKind;

      printf("[Robot Control Parameter] set %s to %s\n", name.c_str(),
             controlParameterValueToString(
                 v, (ControlParameterValueKind)msg->parameterKind)
                 .c_str());

    } break;

    default: {
      throw std::runtime_error("parameter type unsupported");
    }
    break;
  }
  _interfaceLCM.publish("interface_response", &_parameter_response_lcmt);
}


MiniCheetahHardwareBridge::MiniCheetahHardwareBridge(RobotController* robot_ctrl, bool load_parameters_from_file)
    : HardwareBridge(robot_ctrl), _spiLcm(getLcmUrl(255)), _microstrainLcm(getLcmUrl(255)) {
  _load_parameters_from_file = load_parameters_from_file;
}

/*!
 * Main method for Mini Cheetah hardware
 */
void MiniCheetahHardwareBridge::run() {
  initCommon();
  initHardware();

  if(_load_parameters_from_file) {
    printf("[Hardware Bridge] Loading parameters from file...\n");

    try {
      _robotParams.initializeFromYamlFile(THIS_COM "config/mini-cheetah-defaults.yaml");
    } catch(std::exception& e) {
      printf("Failed to initialize robot parameters from yaml file: %s\n", e.what());
      exit(1);
    }

    if(!_robotParams.isFullyInitialized()) {
      printf("Failed to initialize all robot parameters\n");
      exit(1);
    }

    printf("Loaded robot parameters\n");

    if(_userControlParameters) {
      try {
        _userControlParameters->initializeFromYamlFile(THIS_COM "config/mc-mit-ctrl-user-parameters.yaml");
      } catch(std::exception& e) {
        printf("Failed to initialize user parameters from yaml file: %s\n", e.what());
        exit(1);
      }

      if(!_userControlParameters->isFullyInitialized()) {
        printf("Failed to initialize all user parameters\n");
        exit(1);
      }

      printf("Loaded user parameters\n");
    } else {
      printf("Did not load user parameters because there aren't any\n");
    }
  } else {
    printf("[Hardware Bridge] Loading parameters over LCM...\n");
    while (!_robotParams.isFullyInitialized()) {
      printf("[Hardware Bridge] Waiting for robot parameters...\n");
      usleep(1000000);
    }

    if(_userControlParameters) {
      while (!_userControlParameters->isFullyInitialized()) {
        printf("[Hardware Bridge] Waiting for user parameters...\n");
        usleep(1000000);
      }
    }
  }



  printf("[Hardware Bridge] Got all parameters, starting up!\n");

  _robotRunner =
      new RobotRunner(_controller, &taskManager, _robotParams.controller_dt, "robot-control");

  _robotRunner->driverCommand = &_gamepadCommand;
  _robotRunner->spiData = &_spiData;
  _robotRunner->spiCommand = &_spiCommand;
  _robotRunner->robotType = RobotType::MINI_CHEETAH;
  _robotRunner->vectorNavData = &_vectorNavData;
  _robotRunner->controlParameters = &_robotParams;
  _robotRunner->visualizationData = &_visualizationData;
  _robotRunner->cheetahMainVisualization = &_mainCheetahVisualization;

  _robotRunner->init();
  _firstRun = false;

  // init control thread

  statusTask.start();

  // spi Task start
  PeriodicMemberFunction<MiniCheetahHardwareBridge> spiTask(
      &taskManager, .002, "spi", &MiniCheetahHardwareBridge::runSpi, this);
  spiTask.start();

  // microstrain
  if(_microstrainInit)
    _microstrainThread = std::thread(&MiniCheetahHardwareBridge::runMicrostrain, this);

  // robot controller start
  _robotRunner->start();

  // visualization start
  PeriodicMemberFunction<MiniCheetahHardwareBridge> visualizationLCMTask(
      &taskManager, .0167, "lcm-vis",
      &MiniCheetahHardwareBridge::publishVisualizationLCM, this);
  visualizationLCMTask.start();

  // rc controller
  _port = init_sbus(false);  // Not Simulation
  PeriodicMemberFunction<HardwareBridge> sbusTask(
      &taskManager, .005, "rc_controller", &HardwareBridge::run_sbus, this);
  sbusTask.start();

  // temporary hack: microstrain logger
  PeriodicMemberFunction<MiniCheetahHardwareBridge> microstrainLogger(
      &taskManager, .001, "microstrain-logger", &MiniCheetahHardwareBridge::logMicrostrain, this);
  microstrainLogger.start();

  for (;;) {
    usleep(1000000);
    // printf("joy %f\n", _robotRunner->driverCommand->leftStickAnalog[0]);
  }
}

/*!
 * Receive RC with SBUS
 */
void HardwareBridge::run_sbus() {
  if (_port > 0) {
    int x = receive_sbus(_port);
    if (x) {
      sbus_packet_complete();
    }
  }
}

void MiniCheetahHardwareBridge::runMicrostrain() {
  while(true) {
    _microstrainImu.run();

#ifdef USE_MICROSTRAIN
    _vectorNavData.accelerometer = _microstrainImu.acc;
    _vectorNavData.quat[0] = _microstrainImu.quat[1];
    _vectorNavData.quat[1] = _microstrainImu.quat[2];
    _vectorNavData.quat[2] = _microstrainImu.quat[3];
    _vectorNavData.quat[3] = _microstrainImu.quat[0];
    _vectorNavData.gyro = _microstrainImu.gyro;
#endif
  }


}

void MiniCheetahHardwareBridge::logMicrostrain() {
  _microstrainImu.updateLCM(&_microstrainData);
  _microstrainLcm.publish("microstrain", &_microstrainData);
}

/*!
 * Initialize Mini Cheetah specific hardware
 */
void MiniCheetahHardwareBridge::initHardware() {
  _vectorNavData.quat << 1, 0, 0, 0;
#ifndef USE_MICROSTRAIN
  printf("[MiniCheetahHardware] Init vectornav\n");
  if (!init_vectornav(&_vectorNavData)) {
    printf("Vectornav failed to initialize\n");
    //initError("failed to initialize vectornav!\n", false);
  }
#endif

  init_spi();
  _microstrainInit = _microstrainImu.tryInit(0, 921600);
}

/*!
 * Run Mini Cheetah SPI
 */
void MiniCheetahHardwareBridge::runSpi() {
  spi_command_t* cmd = get_spi_command();
  spi_data_t* data = get_spi_data();

  memcpy(cmd, &_spiCommand, sizeof(spi_command_t));
  spi_driver_run();
  memcpy(&_spiData, data, sizeof(spi_data_t));

  _spiLcm.publish("spi_data", data);
  _spiLcm.publish("spi_command", cmd);
}

/*!
 * Send LCM visualization data
 */
void HardwareBridge::publishVisualizationLCM() {
  cheetah_visualization_lcmt visualization_data;
  for (int i = 0; i < 3; i++) {
    visualization_data.x[i] = _mainCheetahVisualization.p[i];
  }

  for (int i = 0; i < 4; i++) {
    visualization_data.quat[i] = _mainCheetahVisualization.quat[i];
    visualization_data.rgba[i] = _mainCheetahVisualization.color[i];
  }

  for (int i = 0; i < 12; i++) {
    visualization_data.q[i] = _mainCheetahVisualization.q[i];
  }

  _visualizationLCM.publish("main_cheetah_visualization", &visualization_data);
}

Cheetah3HardwareBridge::Cheetah3HardwareBridge(RobotController *rc) : HardwareBridge(rc) {

}

void Cheetah3HardwareBridge::run() {
  initCommon();
  // initHardware(); todo

  printf("[Hardware Bridge] Loading parameters over LCM...\n");
  while (!_robotParams.isFullyInitialized()) {
    printf("[Hardware Bridge] Waiting for robot parameters...\n");
    usleep(1000000);
  }

  if(_userControlParameters) {
    while (!_userControlParameters->isFullyInitialized()) {
      printf("[Hardware Bridge] Waiting for user parameters...\n");
      usleep(1000000);
    }
  }

  printf("[Hardware Bridge] Got all parameters, starting up!\n");

  _robotRunner =
      new RobotRunner(_controller, &taskManager, _robotParams.controller_dt, "robot-control");

  _robotRunner->driverCommand = &_gamepadCommand;
  _robotRunner->tiBoardData = _tiBoardData;
  _robotRunner->tiBoardCommand = _tiBoardCommand;
  _robotRunner->robotType = RobotType::CHEETAH_3;
  _robotRunner->controlParameters = &_robotParams;
  _robotRunner->visualizationData = &_visualizationData;
  _robotRunner->cheetahMainVisualization = &_mainCheetahVisualization;

  _robotRunner->init();
  _firstRun = false;

  // init control thread

  statusTask.start();

  // Ecat Task start todo
//  PeriodicMemberFunction<Cheetah3HardwareBridge> ecatTask(
//      &taskManager, .002, "ecat", &Cheetah3HardwareBridge::runEcat, this);
//  ecatTask.start();

  // microstrain
//  if(_microstrainInit)
//    _microstrainThread = std::thread(&MiniCheetahHardwareBridge::runMicrostrain, this);

  // robot controller start
  _robotRunner->start();

  // visualization start
  PeriodicMemberFunction<Cheetah3HardwareBridge> visualizationLCMTask(
      &taskManager, .0167, "lcm-vis",
      &MiniCheetahHardwareBridge::publishVisualizationLCM, this);
  visualizationLCMTask.start();

  // rc controller
  _port = init_sbus(false);  // Not Simulation
  PeriodicMemberFunction<HardwareBridge> sbusTask(
      &taskManager, .005, "rc_controller", &HardwareBridge::run_sbus, this);
  sbusTask.start();

  // temporary hack: microstrain logger
//  PeriodicMemberFunction<MiniCheetahHardwareBridge> microstrainLogger(
//      &taskManager, .001, "microstrain-logger", &MiniCheetahHardwareBridge::logMicrostrain, this);
//  microstrainLogger.start();

  for (;;) {
    usleep(1000000);
    // printf("joy %f\n", _robotRunner->driverCommand->leftStickAnalog[0]);
  }
}

#endif
