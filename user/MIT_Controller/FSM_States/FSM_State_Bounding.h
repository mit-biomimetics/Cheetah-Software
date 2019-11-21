#ifndef FSM_STATE_BOUNDING_H
#define FSM_STATE_BOUNDING_H

#include "FSM_State.h"

template<typename T> class WBC_Ctrl;
/**
 *
 */
template <typename T>
class FSM_State_Bounding : public FSM_State<T> {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  FSM_State_Bounding(ControlFSMData<T>* _controlFSMData);

  // Behavior to be carried out when entering a state
  void onEnter();

  // Run the normal behavior for the state
  void run();

  // Checks for any transition triggers
  FSM_StateName checkTransition();

  // Manages state specific transitions
  TransitionData<T> transition();

  // Behavior to be carried out when exiting a state
  void onExit();

 private:
  // Keep track of the control iterations
  int iter = 0;
  WBC_Ctrl<T> * bounding_ctrl;

  // Parses contact specific controls to the leg controller
  void BoundingControlStep();

};

#endif  // FSM_STATE_BOUNDING_H
