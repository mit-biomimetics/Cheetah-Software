/*============================ Locomotion =============================*/
/**
 * FSM State for robot locomotion. Manages the contact specific logic
 * and handles calling the interfaces to the controllers. This state
 * should be independent of controller, gait, and desired trajectory.
 */

#include "FSM_State_Locomotion.h"
#include <Utilities/Timer.h>

/**
 * Constructor for the FSM State that passes in state specific info to
 * the generic FSM State constructor.
 *
 * @param _controlFSMData holds all of the relevant control data
 */
template <typename T>
FSM_State_Locomotion<T>::FSM_State_Locomotion(
    ControlFSMData<T>* _controlFSMData)
    : FSM_State<T>(_controlFSMData, FSM_StateName::LOCOMOTION, "LOCOMOTION") {
  // Set the safety checks
  this->turnOnAllSafetyChecks();

  // Initialize GRF and footstep locations to 0s
  this->footFeedForwardForces = Mat34<T>::Zero();
  this->footstepLocations = Mat34<T>::Zero();
}

template <typename T>
void FSM_State_Locomotion<T>::onEnter() {
  // Default is to not transition
  this->nextStateName = this->stateName;

  // Reset the transition data
  this->transitionData.zero();
}

/**
 * Calls the functions to be executed on each control loop iteration.
 */
template <typename T>
void FSM_State_Locomotion<T>::run() {
  // Call the locomotion control logic for this iteration
  LocomotionControlStep();
}

/**
 * Manages which states can be transitioned into either by the user
 * commands or state event triggers.
 *
 * @return the enumerated FSM state name to transition into
 */
template <typename T>
FSM_StateName FSM_State_Locomotion<T>::checkTransition() {
  // Get the next state
  iter++;

  // Switch FSM control mode
  switch ((int)this->_data->controlParameters->control_mode) {
    case K_LOCOMOTION:
      // Normal operation for state based transitions

      // Need a working state estimator for automatic state based transition
      /*if (velocity < v_min) {
        this->nextStateName = FSM_StateName::BALANCE_STAND;

        // Transition over the duration of one period
        this->transitionDuration =
      this->_data->_gaitScheduler->gaitData.periodTimeNominal;

        // Notify the gait scheduler that the robot is transitioning to stand
        this->_data->_gaitScheduler->gaitData._nextGait =
      GaitType::TRANSITION_TO_STAND;

      }*/

      // in place to show automatic non user requested transitions
      /*if (iter >= 1387) {  // 2058) {
        // Set the next state to be BALANCE_STAND
        this->nextStateName = FSM_StateName::BALANCE_STAND;

        // Transition time is 1 gait period
        this->transitionDuration =
            this->_data->_gaitScheduler->gaitData.periodTimeNominal;

        // Signal the gait scheduler that we are transitioning to standing
        this->_data->_gaitScheduler->gaitData._nextGait =
            GaitType::TRANSITION_TO_STAND;

        // Notify the control parameters of the mode switch
        this->_data->controlParameters->control_mode = K_BALANCE_STAND;

        // Reset iteration counter
        iter = 0;
      }
       */
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

    default:
      std::cout << "[CONTROL FSM] Bad Request: Cannot transition from "
                << K_LOCOMOTION << " to "
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
TransitionData<T> FSM_State_Locomotion<T>::transition() {
  // Switch FSM control mode
  switch (this->nextStateName) {
    case FSM_StateName::BALANCE_STAND:
      LocomotionControlStep();

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
void FSM_State_Locomotion<T>::onExit() {
  // Nothing to clean up when exiting
  iter = 0;
}

/**
 * Calculate the commands for the leg controllers for each of the feet by
 * calling the appropriate balance controller and parsing the results for
 * each stance or swing leg.
 */
template <typename T>
void FSM_State_Locomotion<T>::LocomotionControlStep() {
  // StateEstimate<T> stateEstimate = this->_data->_stateEstimator->getResult();

  // Contact state logic
  // estimateContact();

  // Run the balancing controllers to get GRF and next step locations
  //Timer t;
  cMPCOld.run<T>(*this->_data);
  //printf("entire run time %.3f\n", t.getMs());
//  this->runControls();
//
//  // Calculate appropriate control actions for each leg to be sent out
//  for (int leg = 0; leg < 4; leg++) {
//    // The actual contact logic should come from the contact estimator later
//    // rather than schedule
//    if (this->_data->_gaitScheduler->gaitData.contactStateScheduled(leg)) {
//      // Leg is in contact
//
//      // Impedance control for the stance leg
//      StanceLegImpedanceControl(leg);
//
//      // Stance leg Ground Reaction Force command
//      this->_data->_legController->commands[leg].forceFeedForward =
//          this->footFeedForwardForces.col(leg);
//
//    } else if (!this->_data->_gaitScheduler->gaitData.contactStateScheduled(
//                   leg)) {
//      // Leg is not in contact
//
//      // Swing leg trajectory
//      // TODO
//      this->footstepLocations.col(leg)
//          << 0., 0.,
//          -this->_data->_quadruped->_maxLegLength / 2. +
//            sin(2.*M_PI * this->_data->_gaitScheduler->gaitData.phaseSwing(leg) );
//
//      // Swing leg impedance control
//      this->cartesianImpedanceControl(
//          leg, this->footstepLocations.col(leg), Vec3<T>::Zero(),
//          this->_data->controlParameters->stand_kp_cartesian,
//          this->_data->controlParameters->stand_kd_cartesian);
//
//      // Feedforward torques for swing leg tracking
//      // TODO
//
//    } else {
//      std::cout << "[CONTROL ERROR] Undefined scheduled contact state\n"
//                << std::endl;
//    }
//
//    // Singularity barrier calculation (maybe an overall safety checks
//    // function?)
//    // TODO
//  }
}

/**
 * Stance leg logic for impedance control. Prevent leg slipping and
 * bouncing, as well as tracking the foot velocity during high speeds.
 */
template <typename T>
void FSM_State_Locomotion<T>::StanceLegImpedanceControl(int leg) {
  // Impedance control for the stance leg
  this->cartesianImpedanceControl(
      leg, this->footstepLocations.col(leg), Vec3<T>::Zero(),
      this->_data->controlParameters->stand_kp_cartesian,
      this->_data->controlParameters->stand_kd_cartesian);
}

// template class FSM_State_Locomotion<double>;
template class FSM_State_Locomotion<float>;
