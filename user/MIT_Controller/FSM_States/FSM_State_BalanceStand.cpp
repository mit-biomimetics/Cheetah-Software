/*=========================== Balance Stand ===========================*/
/**
 * FSM State that forces all legs to be on the ground and uses the QP
 * Balance controller for instantaneous balance control.
 */

#include "FSM_State_BalanceStand.h"
#include <Controllers/WBC_Ctrl/LocomotionCtrl/LocomotionCtrl.hpp>

/**
 * Constructor for the FSM State that passes in state specific info to
 * the generic FSM State constructor.
 *
 * @param _controlFSMData holds all of the relevant control data
 */
template <typename T>
FSM_State_BalanceStand<T>::FSM_State_BalanceStand(
    ControlFSMData<T>* _controlFSMData)
    : FSM_State<T>(_controlFSMData, FSM_StateName::BALANCE_STAND,"BALANCE_STAND") {
  // Set the pre controls safety checks
  this->turnOnAllSafetyChecks();
  // Turn off Foot pos command since it is set in WBC as operational task
  this->checkPDesFoot = false;


  // Initialize GRF to 0s
  this->footFeedForwardForces = Mat34<T>::Zero();

  _wbc_ctrl = new LocomotionCtrl<T>(_controlFSMData->_quadruped->buildModel());
  _wbc_data = new LocomotionCtrlData<T>();

  _wbc_ctrl->setFloatingBaseWeight(1000.);
}

template <typename T>
void FSM_State_BalanceStand<T>::onEnter() {
  // Default is to not transition
  this->nextStateName = this->stateName;

  // Reset the transition data
  this->transitionData.zero();

  // Always set the gait to be standing in this state
  this->_data->_gaitScheduler->gaitData._nextGait = GaitType::STAND;
  
  _ini_body_pos = (this->_data->_stateEstimator->getResult()).position;

  if(_ini_body_pos[2] < 0.2) {
    _ini_body_pos[2] = 0.3;
  }

  last_height_command = _ini_body_pos[2];

  _ini_body_ori_rpy = (this->_data->_stateEstimator->getResult()).rpy;
  _body_weight = this->_data->_quadruped->_bodyMass * 9.81;
}

/**
 * Calls the functions to be executed on each control loop iteration.
 */
template <typename T>
void FSM_State_BalanceStand<T>::run() {
  Vec4<T> contactState;
  contactState<< 0.5, 0.5, 0.5, 0.5;
  this->_data->_stateEstimator->setContactPhase(contactState);
  BalanceStandStep();
}

/**
 * Manages which states can be transitioned into either by the user
 * commands or state event triggers.
 *
 * @return the enumerated FSM state name to transition into
 */
template <typename T>
FSM_StateName FSM_State_BalanceStand<T>::checkTransition() {
  // Get the next state
  _iter++;

  // Switch FSM control mode
  switch ((int)this->_data->controlParameters->control_mode) {
    case K_BALANCE_STAND:
      // Normal operation for state based transitions

      // Need a working state estimator for this
      /*if (velocity > v_max) {
        // Notify the State of the upcoming next state
        this->nextStateName = FSM_StateName::LOCOMOTION;

        // Transition instantaneously to locomotion state on request
        this->transitionDuration = 0.0;

        // Set the next gait in the scheduler to
        this->_data->_gaitScheduler->gaitData._nextGait = GaitType::TROT;

      }*/

      // TEST: in place to show automatic non user requested transitions
      /*if (_iter >= 5458) {
        this->nextStateName = FSM_StateName::LOCOMOTION;
        this->_data->controlParameters->control_mode = K_LOCOMOTION;
        this->transitionDuration = 0.0;
        this->_data->_gaitScheduler->gaitData._nextGait =
            GaitType::AMBLE;  // TROT; // Or get whatever is in
                              // main_control_settings
        _iter = 0;
      }*/
      break;

    case K_LOCOMOTION:
      // Requested change to balance stand
      this->nextStateName = FSM_StateName::LOCOMOTION;

      // Transition instantaneously to locomotion state on request
      this->transitionDuration = 0.0;

      // Set the next gait in the scheduler to
      this->_data->_gaitScheduler->gaitData._nextGait = GaitType::TROT;
      break;

    case K_PASSIVE:
      this->nextStateName = FSM_StateName::PASSIVE;
      // Transition time is immediate
      this->transitionDuration = 0.0;

      break;

    case K_VISION:
      this->nextStateName = FSM_StateName::VISION;
      // Transition time is immediate
      this->transitionDuration = 0.0;
      break;

    case K_RECOVERY_STAND:
      this->nextStateName = FSM_StateName::RECOVERY_STAND;
      // Transition time is immediate
      this->transitionDuration = 0.0;
      break;

    case K_BACKFLIP:
      this->nextStateName = FSM_StateName::BACKFLIP;
      this->transitionDuration = 0.;
      break;

    default:
      std::cout << "[CONTROL FSM] Bad Request: Cannot transition from "
                << K_BALANCE_STAND << " to "
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
TransitionData<T> FSM_State_BalanceStand<T>::transition() {
  // Switch FSM control mode
  switch (this->nextStateName) {
    case FSM_StateName::LOCOMOTION:
      BalanceStandStep();

      _iter++;
      if (_iter >= this->transitionDuration * 1000) {
        this->transitionData.done = true;
      } else {
        this->transitionData.done = false;
      }

      break;

    case FSM_StateName::PASSIVE:
      this->turnOffAllSafetyChecks();
      this->transitionData.done = true;
      break;

    case FSM_StateName::RECOVERY_STAND:
      this->transitionData.done = true;
      break;

    case FSM_StateName::BACKFLIP:
      this->transitionData.done = true;
      break;

    case FSM_StateName::VISION:
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
void FSM_State_BalanceStand<T>::onExit() {
  _iter = 0;
}

/**
 * Calculate the commands for the leg controllers for each of the feet.
 */
template <typename T>
void FSM_State_BalanceStand<T>::BalanceStandStep() {

  _wbc_data->pBody_des = _ini_body_pos;
  _wbc_data->vBody_des.setZero();
  _wbc_data->aBody_des.setZero();

  _wbc_data->pBody_RPY_des = _ini_body_ori_rpy;
  if(this->_data->controlParameters->use_rc){
    const rc_control_settings* rc_cmd = this->_data->_desiredStateCommand->rcCommand;
    // Orientation
    _wbc_data->pBody_RPY_des[0] = rc_cmd->rpy_des[0]*1.4;
    _wbc_data->pBody_RPY_des[1] = rc_cmd->rpy_des[1]*0.46;
    _wbc_data->pBody_RPY_des[2] -= rc_cmd->rpy_des[2];

    // Height
    _wbc_data->pBody_des[2] += 0.12 * rc_cmd->height_variation;
  }else{
    // Orientation
    _wbc_data->pBody_RPY_des[0] = 
     0.6* this->_data->_desiredStateCommand->gamepadCommand->leftStickAnalog[0];
     _wbc_data->pBody_RPY_des[1] = 
      0.6*this->_data->_desiredStateCommand->gamepadCommand->rightStickAnalog[0];
    _wbc_data->pBody_RPY_des[2] -= 
      this->_data->_desiredStateCommand->gamepadCommand->rightStickAnalog[1];
    
    // Height
    _wbc_data->pBody_des[2] += 
      0.12 * this->_data->_desiredStateCommand->gamepadCommand->rightStickAnalog[0];
  }
  _wbc_data->vBody_Ori_des.setZero();

  for(size_t i(0); i<4; ++i){
    _wbc_data->pFoot_des[i].setZero();
    _wbc_data->vFoot_des[i].setZero();
    _wbc_data->aFoot_des[i].setZero();
    _wbc_data->Fr_des[i].setZero();
    _wbc_data->Fr_des[i][2] = _body_weight/4.;
    _wbc_data->contact_state[i] = true;
  }
  
  if(this->_data->_desiredStateCommand->trigger_pressed) {
    _wbc_data->pBody_des[2] = 0.05;

    if(last_height_command - _wbc_data->pBody_des[2] > 0.001) {
      _wbc_data->pBody_des[2] = last_height_command - 0.001;
    }
  }
  last_height_command = _wbc_data->pBody_des[2];

  _wbc_ctrl->run(_wbc_data, *this->_data);
}

// template class FSM_State_BalanceStand<double>;
template class FSM_State_BalanceStand<float>;
