/*============================= Stand Up ==============================*/
/**
 * Transitionary state that is called for the robot to stand up into
 * balance control mode.
 */

#include "FSM_State_StandUp.h"

/**
 * Constructor for the FSM State that passes in state specific info to
 * the generic FSM State constructor.
 *
 * @param _controlFSMData holds all of the relevant control data
 */
template <typename T>
FSM_State_StandUp<T>::FSM_State_StandUp(ControlFSMData<T>* _controlFSMData):
  FSM_State<T>(_controlFSMData, FSM_StateName::STAND_UP, "STAND_UP") {
  // Do nothing
  // Set the pre controls safety checks
  this->checkSafeOrientation = false;

  // Post control safety checks
  this->checkPDesFoot = false;
  this->checkForceFeedForward = false;
}


template <typename T>
void FSM_State_StandUp<T>::onEnter() {
  // Default is to not transition
  this->nextStateName = this->stateName;

  // Reset the transition data
  this->transitionData.zero();

  // Reset iteration counter
  iter = 0;

}


/**
 * Calls the functions to be executed on each control loop iteration.
 */
template <typename T>
void FSM_State_StandUp<T>::run() {
  //float h0 = 0.2;
  //float heightDes = h0 + (0.45 - h0) * ((iter / 1000.0) / 0.5);
}


/**
 * Manages which states can be transitioned into either by the user
 * commands or state event triggers.
 *
 * @return the enumerated FSM state name to transition into
 */
template <typename T>
FSM_StateName FSM_State_StandUp<T>::checkTransition() {
  this->nextStateName = this->stateName;
  iter++;

  // Switch FSM control mode
  switch ((int)this->_data->controlParameters->control_mode) {
  case K_STAND_UP:
    // Normal operation for state based transitions

    // After 0.5s, the robot should be standing and can go to balance
    if (iter / 1000.0 > 0.5) {
      this->nextStateName = FSM_StateName::BALANCE_STAND;

      // Notify the control parameters of the mode switch
      this->_data->controlParameters->control_mode = K_BALANCE_STAND;
    }
    break;

  case K_PASSIVE:  // normal c
    this->nextStateName = FSM_StateName::PASSIVE;
    break;

  default:
    std::cout << "[CONTROL FSM] Bad Request: Cannot transition from " << K_PASSIVE << " to " << this->_data->controlParameters->control_mode << std::endl;

  }

  // Get the next state
  return this->nextStateName;
}


/**
 * Handles the actual transition for the robot between states.
 * Returns true when the transition is completed.
 *
 * @return true if transition is complete
 */
template <typename T>
TransitionData<T> FSM_State_StandUp<T>::transition() {
  // Finish Transition
  switch (this->nextStateName) {

  case FSM_StateName::PASSIVE:  // normal
    this->transitionData.done = true;
    break;

  case FSM_StateName::BALANCE_STAND:
    this->transitionData.done = true;
    break;

  default:
    std::cout << "[CONTROL FSM] Something went wrong in transition" << std::endl;

  }

  // Return the transition data to the FSM
  return this->transitionData;
}


/**
 * Cleans up the state information on exiting the state.
 */
template <typename T>
void FSM_State_StandUp<T>::onExit() {
  // Nothing to clean up when exiting
}



//template class FSM_State_StandUp<double>;
template class FSM_State_StandUp<float>;