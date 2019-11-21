#ifndef FSM_STATE_BALANCESTAND_H
#define FSM_STATE_BALANCESTAND_H

#include "FSM_State.h"

template<typename T> class WBC_Ctrl;
template<typename T> class LocomotionCtrlData;

/**
 *
 */
template <typename T>
class FSM_State_BalanceStand : public FSM_State<T> {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  FSM_State_BalanceStand(ControlFSMData<T>* _controlFSMData);

  // Behavior to be carried out when entering a state
  void onEnter() override;

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
  int _iter = 0;

  // Parses contact specific controls to the leg controller
  void BalanceStandStep();

  WBC_Ctrl<T> * _wbc_ctrl;
  LocomotionCtrlData<T> * _wbc_data;

  T last_height_command = 0;

  Vec3<T> _ini_body_pos;
  Vec3<T> _ini_body_ori_rpy;
  T _body_weight;
};

#endif  // FSM_STATE_BALANCESTAND_H
