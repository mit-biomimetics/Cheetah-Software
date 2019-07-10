#include "MIT_Controller.hpp"
#include "rt/rt_interface_lcm.h"
MIT_Controller::MIT_Controller():RobotController(){  }

//#define RC_ESTOP
/**
 * Initializes the Control FSM.
 */
void MIT_Controller::initializeController() {
  // Initialize a new GaitScheduler object
  _gaitScheduler = new GaitScheduler<float>();

  // Initialize the DesiredStateCommand object
  _desiredStateCommand =
      new DesiredStateCommand<float>(_driverCommand, _stateEstimate);

  // Initialize a new ContactEstimator object
  //_contactEstimator = new ContactEstimator<double>();
  ////_contactEstimator->initialize();

  // Initializes the Control FSM with all the required data
  _controlFSM = new ControlFSM<float>(_quadruped, _stateEstimator,
                                      _legController, _gaitScheduler,
                                      _desiredStateCommand, _controlParameters, _visualizationData, &userParameters);
}

/**
 * Calculate the commands for the leg controllers using the ControlFSM logic.
 */
void MIT_Controller::runController() {

  bool rcEstop = false;

#ifdef RC_ESTOP
  get_main_control_settings(&main_control_settings);
  if(main_control_settings.mode == 0) {
    rcEstop = true;
    printf("ESTOP!\n");
  }
#endif

  if(!rcEstop) {
    // Find the current gait schedule
    _gaitScheduler->step();

    // Find the desired state trajectory
    _desiredStateCommand->convertToStateCommands();

    // Run the Control FSM code
    _controlFSM->runFSM();
  } else {
    printf("ESTOP!\n");
  }

}


