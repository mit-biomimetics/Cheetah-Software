#include "Simulation.h"
#include "Dynamics/Quadruped.h"
#include "ParamHandler.hpp"

#include <Configuration.h>
#include <include/GameController.h>
#include <unistd.h>
#include <fstream>

// if DISABLE_HIGH_LEVEL_CONTROL is defined, the simulator will run freely,
// without trying to connect to a robot
//#define DISABLE_HIGH_LEVEL_CONTROL

/*!
 * Initialize the simulator here.  It is _not_ okay to block here waiting for
 * the robot to connect. Use firstRun() instead!
 */
Simulation::Simulation(RobotType robot, Graphics3D* window,
                       SimulatorControlParameters& params, ControlParameters& userParams, std::function<void(void)> uiUpdate)
    : _simParams(params), _userParams(userParams), _tau(12) {
  _uiUpdate = uiUpdate;
  // init parameters
  printf("[Simulation] Load parameters...\n");
  _simParams
      .lockMutex();  // we want exclusive access to the simparams at this point
  if (!_simParams.isFullyInitialized()) {
    printf(
        "[ERROR] Simulator parameters are not fully initialized.  You forgot: "
        "\n%s\n",
        _simParams.generateUnitializedList().c_str());
    throw std::runtime_error("simulator not initialized");
  }

  // init LCM
  if (_simParams.sim_state_lcm) {
    printf("[Simulation] Setup LCM...\n");
    _lcm = new lcm::LCM(getLcmUrl(_simParams.sim_lcm_ttl));
    if (!_lcm->good()) {
      printf("[ERROR] Failed to set up LCM\n");
      throw std::runtime_error("lcm bad");
    }
  }

  // init quadruped info
  printf("[Simulation] Build quadruped...\n");
  _robot = robot;
  _quadruped = _robot == RobotType::MINI_CHEETAH ? buildMiniCheetah<double>()
                                                 : buildCheetah3<double>();
  printf("[Simulation] Build actuator model...\n");
  _actuatorModels = _quadruped.buildActuatorModels();
  _window = window;

  // init graphics
  if (_window) {
    printf("[Simulation] Setup Cheetah graphics...\n");
    Vec4<float> truthColor, seColor;
    truthColor << 0.2, 0.4, 0.2, 0.6;
    seColor << .75,.75,.75, 1.0;
    _simRobotID = _robot == RobotType::MINI_CHEETAH ? window->setupMiniCheetah(truthColor, true, true)
                                                    : window->setupCheetah3(truthColor, true, true);
    _controllerRobotID = _robot == RobotType::MINI_CHEETAH
                             ? window->setupMiniCheetah(seColor, false, false)
                             : window->setupCheetah3(seColor, false, false);
  }

  // init rigid body dynamics
  printf("[Simulation] Build rigid body model...\n");
  _model = _quadruped.buildModel();
  _robotDataModel = _quadruped.buildModel();
  _simulator =
      new DynamicsSimulator<double>(_model, (bool)_simParams.use_spring_damper);
  _robotDataSimulator = new DynamicsSimulator<double>(_robotDataModel, false);

  DVec<double> zero12(12);
  for (u32 i = 0; i < 12; i++) {
    zero12[i] = 0.;
  }

  // set some sane defaults:
  _tau = zero12;
  _robotControllerState.q = zero12;
  _robotControllerState.qd = zero12;
  FBModelState<double> x0;
  x0.bodyOrientation = rotationMatrixToQuaternion(
      ori::coordinateRotation(CoordinateAxis::Z, 0.));
  // Mini Cheetah
  x0.bodyPosition.setZero();
  x0.bodyVelocity.setZero();
  x0.q = zero12;
  x0.qd = zero12;

  // Mini Cheetah Initial Posture
  // x0.bodyPosition[2] = -0.49;
  // Cheetah 3
  // x0.bodyPosition[2] = -0.34;
  // x0.q[0] = -0.807;
  // x0.q[1] = -1.2;
  // x0.q[2] = 2.4;

  // x0.q[3] = 0.807;
  // x0.q[4] = -1.2;
  // x0.q[5] = 2.4;

  // x0.q[6] = -0.807;
  // x0.q[7] = -1.2;
  // x0.q[8] = 2.4;

  // x0.q[9] = 0.807;
  // x0.q[10] = -1.2;
  // x0.q[11] = 2.4;

  // Initial (Mini Cheetah stand)
  // x0.bodyPosition[2] = -0.185;
  // Cheetah 3
  // x0.bodyPosition[2] = -0.075;

  // x0.q[0] = -0.03;
  // x0.q[1] = -0.79;
  // x0.q[2] = 1.715;

  // x0.q[3] = 0.03;
  // x0.q[4] = -0.79;
  // x0.q[5] = 1.715;

  // x0.q[6] = -0.03;
  // x0.q[7] = -0.72;
  // x0.q[8] = 1.715;

  // x0.q[9] = 0.03;
  // x0.q[10] = -0.72;
  // x0.q[11] = 1.715;

  // Cheetah lies on the ground
  //x0.bodyPosition[2] = -0.45;
  x0.bodyPosition[2] = 0.05;
  x0.q[0] = -0.7;
  x0.q[1] = 1.;
  x0.q[2] = 2.715;

  x0.q[3] = 0.7;
  x0.q[4] = 1.;
  x0.q[5] = 2.715;

  x0.q[6] = -0.7;
  x0.q[7] = -1.0;
  x0.q[8] = -2.715;

  x0.q[9] = 0.7;
  x0.q[10] = -1.0;
  x0.q[11] = -2.715;


  setRobotState(x0);
  _robotDataSimulator->setState(x0);

  printf("[Simulation] Setup low-level control...\n");
  // init spine:
  if (_robot == RobotType::MINI_CHEETAH) {
    for (int leg = 0; leg < 4; leg++) {
      _spineBoards[leg].init(Quadruped<float>::getSideSign(leg), leg);
      _spineBoards[leg].data = &_spiData;
      _spineBoards[leg].cmd = &_spiCommand;
      _spineBoards[leg].resetData();
      _spineBoards[leg].resetCommand();
    }
  } else if (_robot == RobotType::CHEETAH_3) {
    // init ti board
    for (int leg = 0; leg < 4; leg++) {
      _tiBoards[leg].init(Quadruped<float>::getSideSign(leg));
      _tiBoards[leg].set_link_lengths(_quadruped._abadLinkLength,
                                      _quadruped._hipLinkLength,
                                      _quadruped._kneeLinkLength);
      _tiBoards[leg].reset_ti_board_command();
      _tiBoards[leg].reset_ti_board_data();
      _tiBoards[leg].run_ti_board_iteration();
    }
  } else {
    assert(false);
  }

  // init shared memory
  printf("[Simulation] Setup shared memory...\n");
  _sharedMemory.createNew(DEVELOPMENT_SIMULATOR_SHARED_MEMORY_NAME, true);
  _sharedMemory().init();

  // shared memory fields:
  _sharedMemory().simToRobot.robotType = _robot;
  _window->_drawList._visualizationData =
      &_sharedMemory().robotToSim.visualizationData;

  // load robot control parameters
  printf("[Simulation] Load control parameters...\n");
  if (_robot == RobotType::MINI_CHEETAH) {
    _robotParams.initializeFromYamlFile(getConfigDirectoryPath() +
                                        MINI_CHEETAH_DEFAULT_PARAMETERS);
  } else if (_robot == RobotType::CHEETAH_3) {
    _robotParams.initializeFromYamlFile(getConfigDirectoryPath() +
                                        CHEETAH_3_DEFAULT_PARAMETERS);
  } else {
    assert(false);
  }

  if (!_robotParams.isFullyInitialized()) {
    printf("Not all robot control parameters were initialized. Missing:\n%s\n",
           _robotParams.generateUnitializedList().c_str());
    throw std::runtime_error("not all parameters initialized from ini file");
  }
  // init IMU simulator
  printf("[Simulation] Setup IMU simulator...\n");
  _imuSimulator = new ImuSimulator<double>(_simParams);

  _simParams.unlockMutex();
  printf("[Simulation] Ready!\n");
}

void Simulation::sendControlParameter(const std::string& name,
                                      ControlParameterValue value,
                                      ControlParameterValueKind kind, bool isUser) {
  ControlParameterRequest& request =
      _sharedMemory().simToRobot.controlParameterRequest;
  ControlParameterResponse& response =
      _sharedMemory().robotToSim.controlParameterResponse;

  // first check no pending message
  assert(request.requestNumber == response.requestNumber);

  // new message
  request.requestNumber++;

  // message data
  request.requestKind = isUser ? ControlParameterRequestKind::SET_USER_PARAM_BY_NAME : ControlParameterRequestKind::SET_ROBOT_PARAM_BY_NAME;
  strcpy(request.name, name.c_str());
  request.value = value;
  request.parameterKind = kind;
  printf("%s\n", request.toString().c_str());

  // run robot:
  _robotMutex.lock();
  _sharedMemory().simToRobot.mode = SimulatorMode::RUN_CONTROL_PARAMETERS;
  _sharedMemory().simulatorIsDone();

  // wait for robot code to finish
  if (_sharedMemory().waitForRobotWithTimeout()) {
  } else {
    handleControlError();
    request.requestNumber = response.requestNumber; // so if we come back we won't be off by 1
    _robotMutex.unlock();
    return;
  }

  //_sharedMemory().waitForRobot();
  _robotMutex.unlock();

  // verify response is good
  assert(response.requestNumber == request.requestNumber);
  assert(response.parameterKind == request.parameterKind);
  assert(std::string(response.name) == request.name);
}

/*!
 * Report a control error.  This doesn't throw and exception and will return so you can clean up
 */
void Simulation::handleControlError() {
  _wantStop = true;
  _running = false;
  _connected = false;
  _uiUpdate();
  if(!_sharedMemory().robotToSim.errorMessage[0]) {
    printf(
      "[ERROR] Control code timed-out!\n");
    _errorCallback("Control code has stopped responding without giving an error message.\nIt has likely crashed - "
                   "check the output of the control code for more information");

  } else {
    printf("[ERROR] Control code has an error!\n");
    _errorCallback("Control code has an error:\n" + std::string(_sharedMemory().robotToSim.errorMessage));
  }

}

/*!
 * Called before the simulator is run the first time.  It's okay to put stuff in
 * here that blocks on having the robot connected.
 */
void Simulation::firstRun() {
  // connect to robot
  _robotMutex.lock();
  _sharedMemory().simToRobot.mode = SimulatorMode::DO_NOTHING;
  _sharedMemory().simulatorIsDone();

  printf("[Simulation] Waiting for robot...\n");

  // this loop will check to see if the robot is connected at 10 Hz
  // doing this in a loop allows us to click the "stop" button in the GUI
  // and escape from here before the robot code connects, if needed
  while (!_sharedMemory().tryWaitForRobot()) {
    if (_wantStop) {
      return;
    }
    usleep(100000);
  }
  printf("Success! the robot is alive\n");
  _connected = true;
  _uiUpdate();
  _robotMutex.unlock();

  // send all control parameters
  printf("[Simulation] Send robot control parameters to robot...\n");
  for (auto& kv : _robotParams.collection._map) {
    sendControlParameter(kv.first, kv.second->get(kv.second->_kind),
                         kv.second->_kind, false);
  }

  for (auto& kv : _userParams.collection._map) {
    sendControlParameter(kv.first, kv.second->get(kv.second->_kind),
                         kv.second->_kind, true);
  }
}

/*!
 * Take a single timestep of dt seconds
 */
void Simulation::step(double dt, double dtLowLevelControl,
                      double dtHighLevelControl) {
  // Low level control (if needed)
  if (_currentSimTime >= _timeOfNextLowLevelControl) {
    lowLevelControl();
    _timeOfNextLowLevelControl = _timeOfNextLowLevelControl + dtLowLevelControl;
  }

  // High level control
  if (_currentSimTime >= _timeOfNextHighLevelControl) {
#ifndef DISABLE_HIGH_LEVEL_CONTROL
    highLevelControl();
#endif
    _timeOfNextHighLevelControl =
        _timeOfNextHighLevelControl + dtHighLevelControl;
  }

  // actuator model:
  if (_robot == RobotType::MINI_CHEETAH) {
    for (int leg = 0; leg < 4; leg++) {
      for (int joint = 0; joint < 3; joint++) {
        _tau[leg * 3 + joint] = _actuatorModels[joint].getTorque(
            _spineBoards[leg].torque_out[joint],
            _simulator->getState().qd[leg * 3 + joint]);
      }
    }
  } else if (_robot == RobotType::CHEETAH_3) {
    for (int leg = 0; leg < 4; leg++) {
      for (int joint = 0; joint < 3; joint++) {
        _tau[leg * 3 + joint] = _actuatorModels[joint].getTorque(
            _tiBoards[leg].data->tau_des[joint],
            _simulator->getState().qd[leg * 3 + joint]);
      }
    }
  } else {
    assert(false);
  }

  // dynamics
  _currentSimTime += dt;

  // Set Homing Information
  RobotHomingInfo<double> homing;
  homing.active_flag = _simParams.go_home;
  homing.position = _simParams.home_pos;
  homing.rpy = _simParams.home_rpy;
  homing.kp_lin = _simParams.home_kp_lin;
  homing.kd_lin = _simParams.home_kd_lin;
  homing.kp_ang = _simParams.home_kp_ang;
  homing.kd_ang = _simParams.home_kd_ang;
  _simulator->setHoming(homing);

  _simulator->step(dt, _tau, _simParams.floor_kp, _simParams.floor_kd);
}

void Simulation::lowLevelControl() {
  if (_robot == RobotType::MINI_CHEETAH) {
    // update spine board data:
    for (int leg = 0; leg < 4; leg++) {
      _spiData.q_abad[leg] = _simulator->getState().q[leg * 3 + 0];
      _spiData.q_hip[leg] = _simulator->getState().q[leg * 3 + 1];
      _spiData.q_knee[leg] = _simulator->getState().q[leg * 3 + 2];

      _spiData.qd_abad[leg] = _simulator->getState().qd[leg * 3 + 0];
      _spiData.qd_hip[leg] = _simulator->getState().qd[leg * 3 + 1];
      _spiData.qd_knee[leg] = _simulator->getState().qd[leg * 3 + 2];
    }

    // run spine board control:
    for (auto& spineBoard : _spineBoards) {
      spineBoard.run();
    }

  } else if (_robot == RobotType::CHEETAH_3) {
    // update data
    for (int leg = 0; leg < 4; leg++) {
      for (int joint = 0; joint < 3; joint++) {
        _tiBoards[leg].data->q[joint] =
            _simulator->getState().q[leg * 3 + joint];
        _tiBoards[leg].data->dq[joint] =
            _simulator->getState().qd[leg * 3 + joint];
      }
    }

    // run control
    for (auto& tiBoard : _tiBoards) {
      tiBoard.run_ti_board_iteration();
    }
  } else {
    assert(false);
  }
}



void Simulation::highLevelControl() {
  // send joystick data to robot:
  _sharedMemory().simToRobot.gamepadCommand = _window->getDriverCommand();
  _sharedMemory().simToRobot.gamepadCommand.applyDeadband(
      _simParams.game_controller_deadband);

  // send IMU data to robot:
  _imuSimulator->updateCheaterState(_simulator->getState(),
                                    _simulator->getDState(),
                                    _sharedMemory().simToRobot.cheaterState);

  _imuSimulator->updateVectornav(_simulator->getState(),
                                   _simulator->getDState(),
                                   &_sharedMemory().simToRobot.vectorNav);


  // send leg data to robot
  if (_robot == RobotType::MINI_CHEETAH) {
    _sharedMemory().simToRobot.spiData = _spiData;
  } else if (_robot == RobotType::CHEETAH_3) {
    for (int i = 0; i < 4; i++) {
      _sharedMemory().simToRobot.tiBoardData[i] = *_tiBoards[i].data;
    }
  } else {
    assert(false);
  }

  // signal to the robot that it can start running
  // the _robotMutex is used to prevent qt (which runs in its own thread) from
  // sending a control parameter while the robot code is already running.
  _robotMutex.lock();
  _sharedMemory().simToRobot.mode = SimulatorMode::RUN_CONTROLLER;
  _sharedMemory().simulatorIsDone();

  // wait for robot code to finish (and send LCM while waiting)
  if (_lcm) {
    buildLcmMessage();
    _lcm->publish(SIM_LCM_NAME, &_simLCM);
  }

  // first make sure we haven't killed the robot code
  if (_wantStop) return;

  // next try waiting at most 1 second:
  if (_sharedMemory().waitForRobotWithTimeout()) {
  } else {
    handleControlError();
    _robotMutex.unlock();
    return;
  }
  _robotMutex.unlock();

  // update
  if (_robot == RobotType::MINI_CHEETAH) {
    _spiCommand = _sharedMemory().robotToSim.spiCommand;

    // pretty_print(_spiCommand.q_des_abad, "q des abad", 4);
    // pretty_print(_spiCommand.q_des_hip, "q des hip", 4);
    // pretty_print(_spiCommand.q_des_knee, "q des knee", 4);
  } else if (_robot == RobotType::CHEETAH_3) {
    for (int i = 0; i < 4; i++) {
      _tiBoards[i].command = _sharedMemory().robotToSim.tiBoardCommand[i];
    }
  } else {
    assert(false);
  }

  _highLevelIterations++;
}

void Simulation::buildLcmMessage() {
  _simLCM.time = _currentSimTime;
  _simLCM.timesteps = _highLevelIterations;
  auto& state = _simulator->getState();
  auto& dstate = _simulator->getDState();

  Vec3<double> rpy = ori::quatToRPY(state.bodyOrientation);
  RotMat<double> Rbody = ori::quaternionToRotationMatrix(state.bodyOrientation);
  Vec3<double> omega = Rbody.transpose() * state.bodyVelocity.head<3>();
  Vec3<double> v = Rbody.transpose() * state.bodyVelocity.tail<3>();

  for (size_t i = 0; i < 4; i++) {
    _simLCM.quat[i] = state.bodyOrientation[i];
  }

  for (size_t i = 0; i < 3; i++) {
    _simLCM.vb[i] = state.bodyVelocity[i + 3];  // linear velocity in body frame
    _simLCM.rpy[i] = rpy[i];
    for (size_t j = 0; j < 3; j++) {
      _simLCM.R[i][j] = Rbody(i, j);
    }
    _simLCM.omegab[i] = state.bodyVelocity[i];
    _simLCM.omega[i] = omega[i];
    _simLCM.p[i] = state.bodyPosition[i];
    _simLCM.v[i] = v[i];
    _simLCM.vbd[i] = dstate.dBodyVelocity[i + 3];
  }

  for (size_t leg = 0; leg < 4; leg++) {
    for (size_t joint = 0; joint < 3; joint++) {
      _simLCM.q[leg][joint] = state.q[leg * 3 + joint];
      _simLCM.qd[leg][joint] = state.qd[leg * 3 + joint];
      _simLCM.qdd[leg][joint] = dstate.qdd[leg * 3 + joint];
      _simLCM.tau[leg][joint] = _tau[leg * 3 + joint];
      size_t gcID = _simulator->getModel()._footIndicesGC.at(leg);
      _simLCM.p_foot[leg][joint] = _simulator->getModel()._pGC.at(gcID)[joint];
      _simLCM.f_foot[leg][joint] = _simulator->getContactForce(gcID)[joint];
    }
  }
}

/*!
 * Add an infinite collision plane to the simulator
 * @param mu          : friction of the plane
 * @param resti       : restitution coefficient
 * @param height      : height of plane
 * @param addToWindow : if true, also adds graphics for the plane
 */
void Simulation::addCollisionPlane(double mu, double resti, double height,
                                   double sizeX, double sizeY, double checkerX,
                                   double checkerY, bool addToWindow) {
  _simulator->addCollisionPlane(mu, resti, height);
  if (addToWindow && _window) {
    _window->lockGfxMutex();
    Checkerboard checker(sizeX, sizeY, checkerX, checkerY);

    size_t graphicsID = _window->_drawList.addCheckerboard(checker, true);
    _window->_drawList.buildDrawList();
    _window->_drawList.updateCheckerboard(height, graphicsID);
    _window->unlockGfxMutex();
  }
}

/*!
 * Add an box collision to the simulator
 * @param mu          : location of the box
 * @param resti       : restitution coefficient
 * @param depth       : depth (x) of box
 * @param width       : width (y) of box
 * @param height      : height (z) of box
 * @param pos         : position of box
 * @param ori         : orientation of box
 * @param addToWindow : if true, also adds graphics for the plane
 */
void Simulation::addCollisionBox(double mu, double resti, double depth,
                                 double width, double height,
                                 const Vec3<double>& pos,
                                 const Mat3<double>& ori, bool addToWindow,
                                 bool transparent) {
  _simulator->addCollisionBox(mu, resti, depth, width, height, pos, ori);
  if (addToWindow && _window) {
    _window->lockGfxMutex();
    _window->_drawList.addBox(depth, width, height, pos, ori, transparent);
    _window->unlockGfxMutex();
  }
}

void Simulation::addCollisionMesh(double mu, double resti, double grid_size,
                                  const Vec3<double>& left_corner_loc,
                                  const DMat<double>& height_map,
                                  bool addToWindow, bool transparent) {
  _simulator->addCollisionMesh(mu, resti, grid_size, left_corner_loc,
                               height_map);
  if (addToWindow && _window) {
    _window->lockGfxMutex();
    _window->_drawList.addMesh(grid_size, left_corner_loc, height_map,
                               transparent);
    _window->unlockGfxMutex();
  }
}

/*!
 * Runs the simulator in the current thread until the _running variable is set
 * to false. Updates graphics at 60 fps if desired. Runs simulation at the
 * desired speed
 * @param dt
 */
void Simulation::runAtSpeed(std::function<void(std::string)> errorCallback, bool graphics) {
  _errorCallback = errorCallback;
  firstRun();  // load the control parameters

  // if we requested to stop, stop.
  if (_wantStop) return;
  assert(!_running);
  _running = true;
  Timer frameTimer;
  Timer freeRunTimer;
  u64 desiredSteps = 0;
  u64 steps = 0;

  double frameTime = 1. / 60.;
  double lastSimTime = 0;

  printf(
      "[Simulator] Starting run loop (dt %f, dt-low-level %f, dt-high-level %f "
      "speed %f graphics %d)...\n",
      _simParams.dynamics_dt, _simParams.low_level_dt, _simParams.high_level_dt,
      _simParams.simulation_speed, graphics);

  while (_running) {
    double dt = _simParams.dynamics_dt;
    double dtLowLevelControl = _simParams.low_level_dt;
    double dtHighLevelControl = _simParams.high_level_dt;
    _desiredSimSpeed = (_window && _window->wantTurbo()) ? 100.f : _simParams.simulation_speed;
    if(_window && _window->wantSloMo()) {
      _desiredSimSpeed /= 10.;
    }
    u64 nStepsPerFrame = (u64)(((1. / 60.) / dt) * _desiredSimSpeed);
    if (!_window->IsPaused() && steps < desiredSteps) {
      _simParams.lockMutex();   
      step(dt, dtLowLevelControl, dtHighLevelControl);
      _simParams.unlockMutex();
      steps++;
    } else {
      double timeRemaining = frameTime - frameTimer.getSeconds();
      if (timeRemaining > 0) {
        usleep((u32)(timeRemaining * 1e6));
      }
    }
    if (frameTimer.getSeconds() > frameTime) {
      double realElapsedTime = frameTimer.getSeconds();
      frameTimer.start();
      if (graphics && _window) {
        double simRate = (_currentSimTime - lastSimTime) / realElapsedTime;
        lastSimTime = _currentSimTime;
        sprintf(_window->infoString,
                "[Simulation Run %5.2fx]\n"
                "real-time:  %8.3f\n"
                "sim-time:   %8.3f\n"
                "rate:       %8.3f\n",
                _desiredSimSpeed, freeRunTimer.getSeconds(), _currentSimTime,
                simRate);
        updateGraphics();
      }
      if (!_window->IsPaused() && (desiredSteps - steps) < nStepsPerFrame)
        desiredSteps += nStepsPerFrame;
    }
  }
}

void Simulation::loadTerrainFile(const std::string& terrainFileName,
                                 bool addGraphics) {
  printf("load terrain %s\n", terrainFileName.c_str());
  ParamHandler paramHandler(terrainFileName);

  if (!paramHandler.fileOpenedSuccessfully()) {
    printf("[ERROR] could not open yaml file for terrain\n");
    throw std::runtime_error("yaml bad");
  }

  std::vector<std::string> keys = paramHandler.getKeys();

  for (auto& key : keys) {
    auto load = [&](double& val, const std::string& name) {
      if (!paramHandler.getValue<double>(key, name, val))
        throw std::runtime_error("terrain read bad: " + key + " " + name);
    };

    auto loadVec = [&](double& val, const std::string& name, size_t idx) {
      std::vector<double> v;
      if (!paramHandler.getVector<double>(key, name, v))
        throw std::runtime_error("terrain read bad: " + key + " " + name);
      val = v.at(idx);
    };

    auto loadArray = [&](double* val, const std::string& name, size_t idx) {
      std::vector<double> v;
      if (!paramHandler.getVector<double>(key, name, v))
        throw std::runtime_error("terrain read bad: " + key + " " + name);
      assert(v.size() == idx);
      for (size_t i = 0; i < idx; i++) val[i] = v[i];
    };

    printf("terrain element %s\n", key.c_str());
    std::string typeName;
    paramHandler.getString(key, "type", typeName);
    if (typeName == "infinite-plane") {
      double mu, resti, height, gfxX, gfxY, checkerX, checkerY;
      load(mu, "mu");
      load(resti, "restitution");
      load(height, "height");
      loadVec(gfxX, "graphicsSize", 0);
      loadVec(gfxY, "graphicsSize", 1);
      loadVec(checkerX, "checkers", 0);
      loadVec(checkerY, "checkers", 1);
      addCollisionPlane(mu, resti, height, gfxX, gfxY, checkerX, checkerY,
                        addGraphics);
    } else if (typeName == "box") {
      double mu, resti, depth, width, height, transparent;
      double pos[3];
      double ori[3];
      load(mu, "mu");
      load(resti, "restitution");
      load(depth, "depth");
      load(width, "width");
      load(height, "height");
      loadArray(pos, "position", 3);
      loadArray(ori, "orientation", 3);
      load(transparent, "transparent");

      Mat3<double> R_box = ori::rpyToRotMat(Vec3<double>(ori));
      R_box.transposeInPlace();  // collisionBox uses "rotation" matrix instead
                                 // of "transformation"
      addCollisionBox(mu, resti, depth, width, height, Vec3<double>(pos), R_box,
                      addGraphics, transparent != 0.);
    } else if (typeName == "stairs") {
      double mu, resti, rise, run, stepsDouble, width, transparent;
      double pos[3];
      double ori[3];
      load(mu, "mu");
      load(resti, "restitution");
      load(rise, "rise");
      load(width, "width");
      load(run, "run");
      load(stepsDouble, "steps");
      loadArray(pos, "position", 3);
      loadArray(ori, "orientation", 3);
      load(transparent, "transparent");

      Mat3<double> R = ori::rpyToRotMat(Vec3<double>(ori));
      Vec3<double> pOff(pos);
      R.transposeInPlace();  // "graphics" rotation matrix

      size_t steps = (size_t)stepsDouble;

      double heightOffset = rise / 2;
      double runOffset = run / 2;
      for (size_t step = 0; step < steps; step++) {
        Vec3<double> p(runOffset, 0, heightOffset);
        p = R * p + pOff;

        addCollisionBox(mu, resti, run, width, heightOffset * 2, p, R,
                        addGraphics, transparent != 0.);

        heightOffset += rise / 2;
        runOffset += run;
      }
    } else if (typeName == "mesh") {
      double mu, resti, transparent, grid;
      Vec3<double> left_corner;
      std::vector<std::vector<double> > height_map_2d;
      load(mu, "mu");
      load(resti, "restitution");
      load(transparent, "transparent");
      load(grid, "grid");
      loadVec(left_corner[0], "left_corner_loc", 0);
      loadVec(left_corner[1], "left_corner_loc", 1);
      loadVec(left_corner[2], "left_corner_loc", 2);

      int x_len(0);
      int y_len(0);
      bool file_input(false);
      paramHandler.getBoolean(key, "heightmap_file", file_input);
      if (file_input) {
        // Read from text file
        std::string file_name;
        paramHandler.getString(key, "heightmap_file_name", file_name);
        std::ifstream f_height;
        f_height.open(THIS_COM "/config/" + file_name);
        if (!f_height.good()) {
          std::cout << "file reading error: "
                    << THIS_COM "../config/" + file_name << std::endl;
        }
        int i(0);
        int j(0);
        double tmp;

        std::string line;
        std::vector<double> height_map_vec;
        while (getline(f_height, line)) {
          std::istringstream iss(line);
          j = 0;
          while (iss >> tmp) {
            height_map_vec.push_back(tmp);
            ++j;
          }
          y_len = j;
          height_map_2d.push_back(height_map_vec);
          height_map_vec.clear();
          // printf("y len: %d\n", y_len);
          ++i;
        }
        x_len = i;

      } else {
        paramHandler.get2DArray(key, "height_map", height_map_2d);
        x_len = height_map_2d.size();
        y_len = height_map_2d[0].size();
        // printf("x, y len: %d, %d\n", x_len, y_len);
      }

      DMat<double> height_map(x_len, y_len);
      for (int i(0); i < x_len; ++i) {
        for (int j(0); j < y_len; ++j) {
          height_map(i, j) = height_map_2d[i][j];
          // printf("height (%d, %d) : %f\n", i, j, height_map(i,j) );
        }
      }
      addCollisionMesh(mu, resti, grid, left_corner, height_map, addGraphics,
                       transparent != 0.);

    } else {
      throw std::runtime_error("unknown terrain " + typeName);
    }
  }
}

void Simulation::updateGraphics() {
  _robotControllerState.bodyOrientation =
      _sharedMemory().robotToSim.mainCheetahVisualization.quat.cast<double>();
  _robotControllerState.bodyPosition =
      _sharedMemory().robotToSim.mainCheetahVisualization.p.cast<double>();
  for (int i = 0; i < 12; i++)
    _robotControllerState.q[i] =
        _sharedMemory().robotToSim.mainCheetahVisualization.q[i];
  _robotDataSimulator->setState(_robotControllerState);
  _robotDataSimulator->forwardKinematics();  // calc all body positions
  _window->_drawList.updateRobotFromModel(*_simulator, _simRobotID, true);
  _window->_drawList.updateRobotFromModel(*_robotDataSimulator,
                                          _controllerRobotID, false);
  _window->_drawList.updateAdditionalInfo(*_simulator);
  _window->update();
}


