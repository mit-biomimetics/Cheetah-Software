/*========================= Impedance Control =========================*/
/**
 * FSM State that allows PD Impedance control in cartesian space for
 * each of the legs.
 */

#include "FSM_State_ImpedanceControl.h"

/**
 * Constructor for the FSM State that passes in state specific info to
 * the generif FSM State constructor.
 *
 * @param _controlFSMData holds all of the relevant control data
 */
template <typename T>
FSM_State_ImpedanceControl<T>::FSM_State_ImpedanceControl(
    ControlFSMData<T>* _controlFSMData)
    : FSM_State<T>(_controlFSMData, FSM_StateName::IMPEDANCE_CONTROL,
                   "IMPEDANCE_CONTROL") {
  // Do nothing here yet
}

template <typename T>
void FSM_State_ImpedanceControl<T>::onEnter() {
  // Default is to not transition
  this->nextStateName = this->stateName;

  // Reset the transition data
  this->transitionData.zero();
}

/**
 * Calls the functions to be executed on each control loop iteration.
 */
template <typename T>
void FSM_State_ImpedanceControl<T>::run() {
  // Do nothing, all commands should begin as zeros
}

/**
 * Manages which states can be transitioned into either by the user
 * commands or state event triggers.
 *
 * @return the enumerated FSM state name to transition into
 */
template <typename T>
FSM_StateName FSM_State_ImpedanceControl<T>::checkTransition() {
  // Get the next state
  // Switch FSM control mode
  switch ((int)this->_data->controlParameters->control_mode) {
    case K_IMPEDANCE_CONTROL:
      // Normal operation for state based transitions
      break;

    case K_BALANCE_STAND:
      // Requested change to balance stand
      this->nextStateName = FSM_StateName::BALANCE_STAND;

      // Transition instantaneously to locomotion state on request
      this->transitionData.tDuration = 0.0;

      // Set the next gait in the scheduler to
      this->_data->_gaitScheduler->gaitData._nextGait = GaitType::STAND;
      break;

    default:
      std::cout << "[CONTROL FSM] Bad Request: Cannot transition from " << 0
                << " to " << this->_data->controlParameters->control_mode
                << std::endl;
  }

  return this->nextStateName;
}

/**
 * Handles the actual transition for the robot between states.
 * Returns true when the transition is completed.
 *
 * @return true if transition is complete
 */
template <typename T>
TransitionData<T> FSM_State_ImpedanceControl<T>::transition() {
  this->transitionData.done = true;

  // Return the transition data to the FSM
  return this->transitionData;
}

/**
 * Cleans up the state information on exiting the state.
 */
template <typename T>
void FSM_State_ImpedanceControl<T>::onExit() {
  // Nothing to clean up when exiting
}

// template class FSM_State_ImpedanceControl<double>;
template class FSM_State_ImpedanceControl<float>;