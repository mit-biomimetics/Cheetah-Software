/*============================ Control FSM ============================*/
/**
 * The Finite State Machine that manages the robot's controls. Handles
 * calls to the FSM State functions and manages transitions between all
 * of the states.
 */

#include "ControlFSM.h"
#include <rt/rt_rc_interface.h>

/**
 * Constructor for the Control FSM. Passes in all of the necessary
 * data and stores it in a struct. Initializes the FSM with a starting
 * state and operating mode.
 *
 * @param _quadruped the quadruped information
 * @param _stateEstimator contains the estimated states
 * @param _legController interface to the leg controllers
 * @param _gaitScheduler controls scheduled foot contact modes
 * @param _desiredStateCommand gets the desired COM state trajectories
 * @param controlParameters passes in the control parameters from the GUI
 */
template <typename T>
ControlFSM<T>::ControlFSM(Quadruped<T>* _quadruped,
                          StateEstimatorContainer<T>* _stateEstimator,
                          LegController<T>* _legController,
                          GaitScheduler<T>* _gaitScheduler,
                          DesiredStateCommand<T>* _desiredStateCommand,
                          RobotControlParameters* controlParameters,
                          VisualizationData* visualizationData,
                          MIT_UserParameters* userParameters)
{
  // Add the pointers to the ControlFSMData struct
  data._quadruped = _quadruped;
  data._stateEstimator = _stateEstimator;
  data._legController = _legController;
  data._gaitScheduler = _gaitScheduler;
  data._desiredStateCommand = _desiredStateCommand;
  data.controlParameters = controlParameters;
  data.visualizationData = visualizationData;
  data.userParameters = userParameters;

  // Initialize and add all of the FSM States to the state list
  statesList.invalid = nullptr;
  statesList.passive = new FSM_State_Passive<T>(&data);
  statesList.jointPD = new FSM_State_JointPD<T>(&data);
  statesList.impedanceControl = new FSM_State_ImpedanceControl<T>(&data);
  statesList.standUp = new FSM_State_StandUp<T>(&data);
  statesList.balanceStand = new FSM_State_BalanceStand<T>(&data);
  statesList.locomotion = new FSM_State_Locomotion<T>(&data);
  statesList.recoveryStand = new FSM_State_RecoveryStand<T>(&data);
  statesList.vision = new FSM_State_Vision<T>(&data);
  statesList.backflip = new FSM_State_BackFlip<T>(&data);
  statesList.frontJump = new FSM_State_FrontJump<T>(&data);

  safetyChecker = new SafetyChecker<T>(&data);

  // Initialize the FSM with the Passive FSM State
  initialize();
}

/**
 * Initialize the Control FSM with the default settings. SHould be set to
 * Passive state and Normal operation mode.
 */
template <typename T>
void ControlFSM<T>::initialize() {
  // Initialize a new FSM State with the control data
  currentState = statesList.passive;

  // Enter the new current state cleanly
  currentState->onEnter();

  // Initialize to not be in transition
  nextState = currentState;

  // Initialize FSM mode to normal operation
  operatingMode = FSM_OperatingMode::NORMAL;
}

/**
 * Called each control loop iteration. Decides if the robot is safe to
 * run controls and checks the current state for any transitions. Runs
 * the regular state behavior if all is normal.
 */
template <typename T>
void ControlFSM<T>::runFSM() {
  // Publish state estimator data to other computer
  //for(size_t i(0); i<3; ++i){
    //_state_estimator.p[i] = data._stateEstimator->getResult().position[i];
    //_state_estimator.quat[i] = data._stateEstimator->getResult().orientation[i];
  //}
    //_state_estimator.quat[3] = data._stateEstimator->getResult().orientation[3];
  //state_estimator_lcm.publish("state_estimator_ctrl_pc", &_state_estimator);

  // Check the robot state for safe operation
  operatingMode = safetyPreCheck();

  if(data.controlParameters->use_rc){
    int rc_mode = data._desiredStateCommand->rcCommand->mode;
    if(rc_mode == RC_mode::RECOVERY_STAND){
      data.controlParameters->control_mode = K_RECOVERY_STAND;

    } else if(rc_mode == RC_mode::LOCOMOTION){
      data.controlParameters->control_mode = K_LOCOMOTION;

    } else if(rc_mode == RC_mode::QP_STAND){
      data.controlParameters->control_mode = K_BALANCE_STAND;

    } else if(rc_mode == RC_mode::VISION){
      data.controlParameters->control_mode = K_VISION;

    } else if(rc_mode == RC_mode::BACKFLIP || rc_mode == RC_mode::BACKFLIP_PRE){
      data.controlParameters->control_mode = K_BACKFLIP;
   }
      //data.controlParameters->control_mode = K_FRONTJUMP;
    //std::cout<< "control mode: "<<data.controlParameters->control_mode<<std::endl;
  }

  // Run the robot control code if operating mode is not unsafe
  if (operatingMode != FSM_OperatingMode::ESTOP) {
    // Run normal controls if no transition is detected
    if (operatingMode == FSM_OperatingMode::NORMAL) {
      // Check the current state for any transition
      nextStateName = currentState->checkTransition();

      // Detect a commanded transition
      if (nextStateName != currentState->stateName) {
        // Set the FSM operating mode to transitioning
        operatingMode = FSM_OperatingMode::TRANSITIONING;

        // Get the next FSM State by name
        nextState = getNextState(nextStateName);

        // Print transition initialized info
        //printInfo(1);

      } else {
        // Run the iteration for the current state normally
        currentState->run();
      }
    }

    // Run the transition code while transition is occuring
    if (operatingMode == FSM_OperatingMode::TRANSITIONING) {
      transitionData = currentState->transition();

      // Check the robot state for safe operation
      safetyPostCheck();

      // Run the state transition
      if (transitionData.done) {
        // Exit the current state cleanly
        currentState->onExit();

        // Print finalizing transition info
        //printInfo(2);

        // Complete the transition
        currentState = nextState;

        // Enter the new current state cleanly
        currentState->onEnter();

        // Return the FSM to normal operation mode
        operatingMode = FSM_OperatingMode::NORMAL;
      }
    } else {
      // Check the robot state for safe operation
      safetyPostCheck();
    }

  } else { // if ESTOP
    currentState = statesList.passive;
    currentState->onEnter();
    nextStateName = currentState->stateName;
  }

  // Print the current state of the FSM
  printInfo(0);

  // Increase the iteration counter
  iter++;
}

/**
 * Checks the robot state for safe operation conditions. If it is in
 * an unsafe state, it will not run the normal control code until it
 * is safe to operate again.
 *
 * @return the appropriate operating mode
 */
template <typename T>
FSM_OperatingMode ControlFSM<T>::safetyPreCheck() {
  // Check for safe orientation if the current state requires it
  if (currentState->checkSafeOrientation && data.controlParameters->control_mode != K_RECOVERY_STAND) {
    if (!safetyChecker->checkSafeOrientation()) {
      operatingMode = FSM_OperatingMode::ESTOP;
      std::cout << "broken: Orientation Safety Ceck FAIL" << std::endl;
    }
  }

  // Default is to return the current operating mode
  return operatingMode;
}

/**
 * Checks the robot state for safe operation commands after calculating the
 * control iteration. Prints out which command is unsafe. Each state has
 * the option to enable checks for commands that it cares about.
 *
 * Should this EDamp / EStop or just continue?
 * Should break each separate check into its own function for clarity
 *
 * @return the appropriate operating mode
 */
template <typename T>
FSM_OperatingMode ControlFSM<T>::safetyPostCheck() {
  // Check for safe desired foot positions
  if (currentState->checkPDesFoot) {
    safetyChecker->checkPDesFoot();
  }

  // Check for safe desired feedforward forces
  if (currentState->checkForceFeedForward) {
    safetyChecker->checkForceFeedForward();
  }

  // Default is to return the current operating mode
  return operatingMode;
}

/**
 * Returns the approptiate next FSM State when commanded.
 *
 * @param  next commanded enumerated state name
 * @return next FSM state
 */
template <typename T>
FSM_State<T>* ControlFSM<T>::getNextState(FSM_StateName stateName) {
  // Choose the correct FSM State by enumerated state name
  switch (stateName) {
    case FSM_StateName::INVALID:
      return statesList.invalid;

    case FSM_StateName::PASSIVE:
      return statesList.passive;

    case FSM_StateName::JOINT_PD:
      return statesList.jointPD;

    case FSM_StateName::IMPEDANCE_CONTROL:
      return statesList.impedanceControl;

    case FSM_StateName::STAND_UP:
      return statesList.standUp;

    case FSM_StateName::BALANCE_STAND:
      return statesList.balanceStand;

    case FSM_StateName::LOCOMOTION:
      return statesList.locomotion;

    case FSM_StateName::RECOVERY_STAND:
      return statesList.recoveryStand;

    case FSM_StateName::VISION:
      return statesList.vision;

    case FSM_StateName::BACKFLIP:
      return statesList.backflip;

    case FSM_StateName::FRONTJUMP:
      return statesList.frontJump;

    default:
      return statesList.invalid;
  }
}

/**
 * Prints Control FSM info at regular intervals and on important events
 * such as transition initializations and finalizations. Separate function
 * to not clutter the actual code.
 *
 * @param printing mode option for regular or an event
 */
template <typename T>
void ControlFSM<T>::printInfo(int opt) {
  switch (opt) {
    case 0:  // Normal printing case at regular intervals
      // Increment printing iteration
      printIter++;

      // Print at commanded frequency
      if (printIter == printNum) {
        std::cout << "[CONTROL FSM] Printing FSM Info...\n";
        std::cout
            << "---------------------------------------------------------\n";
        std::cout << "Iteration: " << iter << "\n";
        if (operatingMode == FSM_OperatingMode::NORMAL) {
          std::cout << "Operating Mode: NORMAL in " << currentState->stateString
                    << "\n";

        } else if (operatingMode == FSM_OperatingMode::TRANSITIONING) {
          std::cout << "Operating Mode: TRANSITIONING from "
                    << currentState->stateString << " to "
                    << nextState->stateString << "\n";

        } else if (operatingMode == FSM_OperatingMode::ESTOP) {
          std::cout << "Operating Mode: ESTOP\n";
        }
        std::cout << "Gait Type: " << data._gaitScheduler->gaitData.gaitName
                  << "\n";
        std::cout << std::endl;

        // Reset iteration counter
        printIter = 0;
      }

      // Print robot info about the robot's status
      // data._gaitScheduler->printGaitInfo();
      // data._desiredStateCommand->printStateCommandInfo();

      break;

    case 1:  // Initializing FSM State transition
      std::cout << "[CONTROL FSM] Transition initialized from "
                << currentState->stateString << " to " << nextState->stateString
                << "\n"
                << std::endl;

      break;

    case 2:  // Finalizing FSM State transition
      std::cout << "[CONTROL FSM] Transition finalizing from "
                << currentState->stateString << " to " << nextState->stateString
                << "\n"
                << std::endl;

      break;
  }
}

// template class ControlFSM<double>; This should be fixed... need to make
// RobotRunner a template
template class ControlFSM<float>;
