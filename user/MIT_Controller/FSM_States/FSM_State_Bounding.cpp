/*============================ Bounding =============================*/
/**
 * FSM State for robot bounding. Manages the contact specific logic
 * and handles calling the interfaces to the controllers. This state
 * should be independent of controller, gait, and desired trajectory.
 */

#include "FSM_State_Bounding.h"
#include <Utilities/Timer.h>
#include <Controllers/WBC_Ctrl/BoundingCtrl/BoundingCtrl.hpp>

/**
 * Constructor for the FSM State that passes in state specific info to
 * the generic FSM State constructor.
 *
 * @param _controlFSMData holds all of the relevant control data
 */
template <typename T>
FSM_State_Bounding<T>::FSM_State_Bounding(
    ControlFSMData<T>* _controlFSMData)
    : FSM_State<T>(_controlFSMData, 
        FSM_StateName::BOUNDING, "BOUNDING") {
  // Set the safety checks
  this->turnOnAllSafetyChecks();
  this->checkPDesFoot = false;

  // Initialize GRF and footstep locations to 0s
  this->footFeedForwardForces = Mat34<T>::Zero();
  this->footstepLocations = Mat34<T>::Zero();
  bounding_ctrl = new BoundingCtrl<T>(_controlFSMData->_quadruped->buildModel());
}

template <typename T>
void FSM_State_Bounding<T>::onEnter() {
  // Default is to not transition
  this->nextStateName = this->stateName;

  // Reset the transition data
  this->transitionData.zero();
}

/**
 * Calls the functions to be executed on each control loop iteration.
 */
template <typename T>
void FSM_State_Bounding<T>::run() {
  // Call the bounding control logic for this iteration
  BoundingControlStep();
}

/**
 * Manages which states can be transitioned into either by the user
 * commands or state event triggers.
 *
 * @return the enumerated FSM state name to transition into
 */
template <typename T>
FSM_StateName FSM_State_Bounding<T>::checkTransition() {
  // Get the next state
  iter++;

  // Switch FSM control mode
  switch ((int)this->_data->controlParameters->control_mode) {
    case K_BOUNDING:
       break;

    case K_BALANCE_STAND:
      // Requested change to BALANCE_STAND
      this->nextStateName = FSM_StateName::BALANCE_STAND;

      // Transition time is immediate
      this->transitionDuration = 0.0;

      break;

    case K_PASSIVE:
      // Requested change to BALANCE_STAND
      this->nextStateName = FSM_StateName::PASSIVE;

      // Transition time is immediate
      this->transitionDuration = 0.0;

      break;

    case K_STAND_UP:
      this->nextStateName = FSM_StateName::STAND_UP;
      this->transitionDuration = 0.;
      break;


    default:
      std::cout << "[CONTROL FSM] Bad Request: Cannot transition from "
                << K_BOUNDING << " to "
                << this->_data->controlParameters->control_mode << std::endl;
  }

  // Return the next state name to the FSM
  return this->nextStateName;
}

/**
 * Handles the actual transition for the robot between states.
 * Returns true when the transition is completed.
 *
 * @return true if transition is complete
 */
template <typename T>
TransitionData<T> FSM_State_Bounding<T>::transition() {
  // Switch FSM control mode
  switch (this->nextStateName) {
    case FSM_StateName::BALANCE_STAND:
      BoundingControlStep();

      iter++;
      if (iter >= this->transitionDuration * 1000) {
        this->transitionData.done = true;
      } else {
        this->transitionData.done = false;
      }

      break;

    case FSM_StateName::PASSIVE:
      this->turnOffAllSafetyChecks();

      this->transitionData.done = true;

      break;

    case FSM_StateName::STAND_UP:
      this->transitionData.done = true;
      break;


    default:
      std::cout << "[CONTROL FSM] Something went wrong in transition"
                << std::endl;
  }

  // Return the transition data to the FSM
  return this->transitionData;
}

/**
 * Cleans up the state information on exiting the state.
 */
template <typename T>
void FSM_State_Bounding<T>::onExit() {
  // Nothing to clean up when exiting
  iter = 0;
}

/**
 * Calculate the commands for the leg controllers for each of the feet by
 * calling the appropriate balance controller and parsing the results for
 * each stance or swing leg.
 */
template <typename T>
void FSM_State_Bounding<T>::BoundingControlStep() {
  // StateEstimate<T> stateEstimate = this->_data->_stateEstimator->getResult();

  // Contact state logic
  // estimateContact();

  // Run the balancing controllers to get GRF and next step locations

  void* dummy = nullptr;
  bounding_ctrl->run(dummy, *this->_data);
}

// template class FSM_State_Bounding<double>;
template class FSM_State_Bounding<float>;
