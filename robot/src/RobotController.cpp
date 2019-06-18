#include "Controllers/OrientationEstimator.h"
#include "RobotController.h"
#include "Dynamics/Cheetah3.h"
#include "Dynamics/MiniCheetah.h"
#include "Controllers/ContactEstimator.h"

#include <Utilities/Utilities_print.h>
#include <WBC_States/BodyCtrl/BodyCtrlTest.hpp>
#include <WBC_States/JPosCtrl/JPosCtrlTest.hpp>
#include <WBC_States/WBICTrot/WBICTrotTest.hpp>
#include <WBC_States/BackFlip/BackFlipTest.hpp>
#include <WBC_States/Bounding/BoundingTest.hpp>
#include <Configuration.h>
#include <ParamHandler.hpp>

#include <Utilities/Timer.h>
#include <unistd.h>

#include "rt/rt_interface_lcm.h"
#include <Controllers/PositionVelocityEstimator.h>

#define REMOTE_CTRL false
//#define REMOTE_CTRL true

/**
 * Initializes the robot model, state estimator, leg controller,
 * robot data, and any control logic specific data.
 */
void RobotController::init() {
  printf("[RobotController] initialize\n");

  // Build the appropriate Quadruped object
  if (robotType == RobotType::MINI_CHEETAH) {
    _quadruped = buildMiniCheetah<float>();
  } else {
    _quadruped = buildCheetah3<float>();
  }

  // Initialize the model and robot data
  _model = _quadruped.buildModel();
  _jpos_initializer = new JPosInitializer<float>(3.);
  _data = new Cheetah_Data<float>();
  _extra_data = new Cheetah_Extra_Data<float>();

  // Always initialize the leg controller and state entimator
  _legController = new LegController<float>(_quadruped);
  _stateEstimator = new StateEstimatorContainer<float>(
    cheaterState, kvhImuData, vectorNavData,
    _legController->datas, &_stateEstimate, controlParameters);
  initializeStateEstimator(false);

  // Choose which control logic option to use
  //CONTROL_OPTION = ControlLogicName::CONTROL_FSM;
  CONTROL_OPTION = ControlLogicName::WBC;

  // Control logic specific initializations
  switch (CONTROL_OPTION) {
  case ControlLogicName::CONTROL_FSM:
    // Initialize Control FSM logic
    initializeControlOptionControlFSM();
    break;

  case ControlLogicName::WBC:
    // Initialize WBC logic
    initializeControlOptionWBC();
    break;
  }
}


/**
 * Runs the overall robot control system by calling each of the major components
 * to run each of their respective steps.
 */
void RobotController::run() {
  // Increment the iteration counter
  iter++;

  // Update the data from the robot
  setupStep();

  // Run the chosen logic during normal control operation
  switch (CONTROL_OPTION) {
  case ControlLogicName::CONTROL_FSM:
    // Run Control FSM logic
    runControlOptionControlFSM();
    break;

  case ControlLogicName::WBC:
    // Run WBC logic
    runControlOptionWBC();
    break;
  }

  // Visualization (will make this into a separate function later)
  for (int leg = 0; leg < 4; leg++) {
    for (int joint = 0; joint < 3; joint++) {
      cheetahMainVisualization->q[leg * 3 + joint] = _legController->datas[leg].q[joint];
    }
  }
  cheetahMainVisualization->p.setZero();
  cheetahMainVisualization->p = _stateEstimate.position;

  // Sets the leg controller commands for the robot appropriate commands
  finalizeStep();
}


/**
 * Initializes the Control FSM.
 */
void RobotController::initializeControlOptionControlFSM() {
  // Initialize a new GaitScheduler object
  _gaitScheduler = new GaitScheduler<float>();

  // Initialize the DesiredStateCommand object
  _desiredStateCommand = new DesiredStateCommand<float>(driverCommand, &_stateEstimate);

  // Initialize a new ContactEstimator object
  //_contactEstimator = new ContactEstimator<double>(); //_contactEstimator->initialize();

  // Initializes the Control FSM with all the required data
  _controlFSM = new ControlFSM<float>(&_quadruped, _stateEstimator, _legController, _gaitScheduler, _desiredStateCommand, controlParameters);

}


/**
 * Calculate the commands for the leg controllers using the ControlFSM logic.
 */
void RobotController::runControlOptionControlFSM() {
  // Run the state estimator step
  _stateEstimator->run(cheetahMainVisualization);

  // Find the current gait schedule
  _gaitScheduler->step();

  // Find the desired state trajectory
  _desiredStateCommand->convertToStateCommands();

  // Run the Control FSM code
  _controlFSM->runFSM();
}


/**
 * Initializes the WBC.
 */
void RobotController::initializeControlOptionWBC() {
  // Get the requested test name from the parameters
  ParamHandler handler(THIS_COM"robot/WBC_States/config/ROBOT_test_setup.yaml");
  std::string test_name;
  if (!handler.getString("test_name", test_name)) {
    printf("cannot find test name\n");
    exit(0);
  }

  // Run the requested test from the string
  if (test_name == "body_ctrl") {
    _wbc_state = new BodyCtrlTest<float>(&_model, robotType);
  } else if (test_name == "wbic_trot") {
    _wbc_state = new WBICTrotTest<float>(&_model, robotType);
  } else if (test_name == "bounding") {
    _wbc_state = new BoundingTest<float>(&_model, robotType);
  } else if (test_name == "back_flip") {
    _wbc_state = new BackFlipTest<float>(&_model, robotType);
  }
}


/**
 * Calculate the commands for the leg controllers using the WBC logic.
 */
void RobotController::runControlOptionWBC() {
  // Run the state estimator step
  _stateEstimator->run(cheetahMainVisualization);

  cheetahMainVisualization->p = _stateEstimate.position;
  //printf("p robot %.3f %.3f %.3f\n", _stateEstimate.position[0], _stateEstimate.position[1],
      //_stateEstimate.position[2]);

  // for debugging the visualizations from robot code
  //testDebugVisualization();
  //StepLocationVisualization();
  //BodyPathVisualization();
  //BodyPathArrowVisualization();

  // RC controller;
#if (REMOTE_CTRL)
  get_main_control_settings(&main_control_settings);
  //printf("%f\n", main_control_settings.mode);
#else
  main_control_settings.mode = 11;
#endif

  if (main_control_settings.mode == 0 ) { // E-Stop
    for (int leg = 0; leg < 4; leg++) {
      for (size_t jidx(0); jidx < cheetah::num_leg_joint; ++jidx) {
        _legController->commands[leg].tauFeedForward[jidx] = 0.;
        _legController->commands[leg].qDes[jidx] = _legController->datas[leg].q[jidx];
        _legController->commands[leg].qdDes[jidx] = 0.;
      }
      _legController->commands[leg].kpJoint.setZero();
      _legController->commands[leg].kdJoint.setZero();
    }
  } else {
    // ======= WBC state command computation  =============== //
    // Commenting out WBC for now to test Locomotion control
    Mat3<float> kpMat;
    Mat3<float> kdMat;
    if (!_jpos_initializer->IsInitialized(_legController)) {
      //if (false) {  // TEST
      kpMat << 5, 0, 0,
            0, 5, 0,
            0, 0, 5;

      kdMat << 0.1, 0, 0,
            0, 0.1, 0,
            0, 0, 0.1;
      _ini_yaw = _stateEstimator->getResult().rpy[2];
      for (int leg = 0; leg < 4; leg++) {
        _legController->commands[leg].kpJoint = kpMat;
        _legController->commands[leg].kdJoint = kdMat;
      }
    } else {
      Vec3<float> rpy = _stateEstimator->getResult().rpy;
      rpy[2] -= _ini_yaw;
      Quat<float> quat_ori = ori::rpyToQuat(rpy);

      // Quaternion (Orientation setting
      for (size_t i(0); i < 4; ++i) { _data->body_ori[i] = quat_ori[i]; }
      // Angular velocity
      for (int i(0); i < 3; ++i) {
        _data->ang_vel[i] = _stateEstimator->getResult().omegaBody[i];
      }

      if (_cheaterModeEnabled) {
        _data->global_body_pos[0] = cheaterState->position[0];
        _data->global_body_pos[1] = cheaterState->position[1];
        _data->global_body_pos[2] = cheaterState->position[2] + 0.5; // ground is -0.5
        _data->cheater_mode = true;
      }

      for (int leg(0); leg < 4; ++leg) {
        for (int jidx(0); jidx < 3; ++jidx) {
          _data->jpos[3 * leg + jidx] = _legController->datas[leg].q[jidx];
          _data->jvel[3 * leg + jidx] = _legController->datas[leg].qd[jidx];
        }
      }
      _data->mode = main_control_settings.mode;

      // Mode Setting
      _data->dir_command[0] = driverCommand->leftStickAnalog[1];
      _data->dir_command[1] = driverCommand->leftStickAnalog[0];

      // Orientation
      _data->ori_command[0] = driverCommand->rightTriggerAnalog;
      _data->ori_command[0] -= driverCommand->leftTriggerAnalog;

      _data->ori_command[1] = driverCommand->rightStickAnalog[1];
      _data->ori_command[2] = driverCommand->rightStickAnalog[0];

#if (REMOTE_CTRL)
      // Remote control
      if (main_control_settings.variable[0] == 4) { // Body Posture Control
        _data->ori_command[0] = main_control_settings.rpy_des[0];
        _data->ori_command[1] = main_control_settings.rpy_des[1];
        _data->ori_command[2] = main_control_settings.rpy_des[2];
      } else {
        _data->dir_command[0] = main_control_settings.v_des[0];
        _data->dir_command[1] = main_control_settings.v_des[1];
        _data->ori_command[2] = main_control_settings.omega_des[2];
      }
#endif
      //Timer timer;
      _wbc_state->GetCommand(_data, _legController->commands, _extra_data);
      //std::cout<< "wbc computation: " << timer.getMs()<<std::endl;
    }
    // === End of WBC state command computation  =========== //
  }
}


void RobotController::StepLocationVisualization() {
  // Cones
  //visualizationData->num_cones = 20*2;
  int num_step = _extra_data->num_step;
  visualizationData->num_cones = num_step;
  for (int j(0) ; j < num_step  ; j++)
  {
    ConeVisualization cone;
    cone.radius = 0.03;
    cone.direction << 0, 0 , 0.05;
    //cone.direction += .2f * j *Vec3<float>(sinf(.4f * j * t + j * .2f), cosf(.4f * j * t + j * .2f), 0);
    cone.point_position <<
                        _extra_data->loc_x[j],
                                    _extra_data->loc_y[j],
                                    (_extra_data->loc_z[j] - 0.5);  // Ground is -0.5
    cone.color << .6 , .2 ,  .4, .6;
    visualizationData->cones[j] = cone;
  }
}
void RobotController::BodyPathVisualization() {
  // Path visualization
  PathVisualization path;
  path.num_points = _extra_data->num_path_pt;
  for (size_t j = 0 ; j < path.num_points ; j++)
  {
    path.position[j] <<
                     _extra_data->path_x[j],
                                 _extra_data->path_y[j],
                                 _extra_data->path_z[j] - 0.5; //Ground is -0.5
  }
  //pretty_print(_extra_data->path_x, "path x", path.num_points);
  //pretty_print(_extra_data->path_y, "path y", path.num_points);
  path.color << 0.6 ,  0.2, 0.05 ,  1;

  visualizationData->num_paths = 1;
  visualizationData->paths[0] = path;
}

void RobotController::BodyPathArrowVisualization() {
  visualizationData->num_arrows = _extra_data->num_middle_pt;

  double ar_len(0.07);
  double yaw;
  for (int i(0); i < _extra_data->num_middle_pt; ++i) {
    visualizationData->arrows[i].base_position <<
        _extra_data->mid_x[i], _extra_data->mid_y[i], _extra_data->mid_z[i] - 0.5 ;
    //visualizationData->arrows[i].direction <<
    //_extra_data->mid_ori_roll[i],
    //_extra_data->mid_ori_pitch[i],
    //_extra_data->mid_ori_yaw[i];

    yaw = _extra_data->mid_ori_yaw[i];

    visualizationData->arrows[i].direction << ar_len*cos(yaw), ar_len*(sin(yaw)), 0.;

    visualizationData->arrows[i].head_width = 0.02;
    visualizationData->arrows[i].head_length = 0.03;
    visualizationData->arrows[i].shaft_width = 0.01;

    visualizationData->arrows[i].color <<
                                       0.8, 0.3, 0.1, 1;
  }
}

void RobotController::testDebugVisualization() {
  // Todo, move some of this into the robot controller.

  float t = (float)_iterations / 1000;

  cheetahMainVisualization->q.setZero();
  cheetahMainVisualization->p.setZero();
  cheetahMainVisualization->quat.setZero();
  cheetahMainVisualization->color << 1, 1, 1, 0.5f;
  cheetahMainVisualization->p[2] = sin(t * 10);

  return;


  // Test sphere visualization
  SphereVisualization sphere;
  sphere.position[0] = 0; sphere.position[1] = 5 * std::sin(t); sphere.position[2] = 5 * std::cos(t);
  sphere.radius = 1;
  sphere.color[0] = 1; sphere.color[1] = 1; sphere.color[2] = 0; sphere.color[3] = 1;
  visualizationData->num_spheres = 1;
  visualizationData->spheres[0] = sphere;


  // Cones
  visualizationData->num_cones = 5;
  for (size_t j = 0 ; j < 5  ; j++)
  {
    ConeVisualization cone;
    cone.radius = (j + 1) / 10.;
    cone.direction << 0, 0 , .3 * j + .2;
    cone.direction += .2f * j * Vec3<float>(sinf(.4f * j * t + j * .2f), cosf(.4f * j * t + j * .2f), 0);
    cone.point_position <<  3 + j, 3,  0;
    cone.color << .4 , .1 ,  .2, .6;
    visualizationData->cones[j] = cone;
  }

  // boxes
  visualizationData->num_blocks = 1;
  BlockVisualization block;
  block.corner_position << -5, -5,  0;
  block.dimension <<  1, .2, .3;
  block.color <<  1,  0,  0, std::fmod((float)_iterations / 1000, 1.f) ;
  visualizationData->blocks[0] = block;


  // Test Arrow Visualization
  visualizationData->num_arrows = 1;

  visualizationData->arrows[0].base_position << 1, 1, 1;

  visualizationData->arrows[0].direction << 1, 1, 1;

  visualizationData->arrows[0].head_width = 0.1;
  visualizationData->arrows[0].head_length = 0.2;
  visualizationData->arrows[0].shaft_width = 0.05;

  visualizationData->arrows[0].color << std::fmod((float)_iterations / 1000 + 0.5f, 1.f), std::fmod((float)_iterations / 1000, 1.f) , 1, .6;

  // Test Path visualization
  PathVisualization path;
  path.num_points = 150;
  for (size_t j = 0 ; j < path.num_points ; j++)
  {
    path.position[j] << j * 2.0 / path.num_points, sin(j * 10. / path.num_points + ((float) _iterations) / 1000) ,  .5;
  }
  path.color << 0 ,  0, 1 ,  1;

  visualizationData->num_paths = 1;
  visualizationData->paths[0] = path;
}


void RobotController::setupStep() {
  // Update the leg data
  if (robotType == RobotType::MINI_CHEETAH) {
    _legController->updateData(spiData);
  } else if (robotType == RobotType::CHEETAH_3) {
    _legController->updateData(tiBoardData);
  } else {
    assert(false);
  }

  // Setup the leg controller for a new iteration
  _legController->zeroCommand();
  _legController->setEnabled(true);
  _legController->setMaxTorqueCheetah3(208.5);

  // state estimator
  // check transition to cheater mode:
  if (!_cheaterModeEnabled && controlParameters->cheater_mode) {
    printf("[RobotController] Transitioning to Cheater Mode...\n");
    initializeStateEstimator(true);
    // todo any configuration
    _cheaterModeEnabled = true;
  }

  // check transition from cheater mode:
  if (_cheaterModeEnabled && !controlParameters->cheater_mode) {
    printf("[RobotController] Transitioning from Cheater Mode...\n");
    initializeStateEstimator(false);
    // todo any configuration
    _cheaterModeEnabled = false;
  }


  // todo safety checks, sanity checks, etc...
}


void RobotController::finalizeStep() {
  if (robotType == RobotType::MINI_CHEETAH) {
    _legController->updateCommand(spiCommand);
  } else if (robotType == RobotType::CHEETAH_3) {
    _legController->updateCommand(tiBoardCommand);
  } else {
    assert(false);
  }
  _iterations++;
}

void RobotController::initializeStateEstimator(bool cheaterMode) {
  _stateEstimator->removeAllEstimators();
  _stateEstimator->addEstimator<ContactEstimator<float>>();
  if (cheaterMode) {
    _stateEstimator->addEstimator<CheaterOrientationEstimator<float>>();
  } else if (robotType == RobotType::MINI_CHEETAH) {
    _stateEstimator->addEstimator<VectorNavOrientationEstimator<float>>();
  } else if (robotType == RobotType::CHEETAH_3) {
    _stateEstimator->addEstimator<KvhOrientationEstimator<float>>();
  } else {
    assert(false);
  }
}

RobotController::~RobotController() {
  delete _legController;
  delete _stateEstimator;
  delete _wbc_state;
  delete _data;
  delete _extra_data;
  delete _jpos_initializer;
}

void RobotController::cleanup() {

}
