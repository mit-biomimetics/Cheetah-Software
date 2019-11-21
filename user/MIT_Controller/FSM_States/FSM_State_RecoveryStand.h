#ifndef FSM_STATE_RECOVERY_STANDUP_H
#define FSM_STATE_RECOVERY_STANDUP_H

#include "FSM_State.h"

/**
 *
 */
template <typename T>
class FSM_State_RecoveryStand : public FSM_State<T> {
 public:
  FSM_State_RecoveryStand(ControlFSMData<T>* _controlFSMData);

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

  TransitionData<T> testTransition();

 private:
  // Keep track of the control iterations
  int iter = 0;
  int _motion_start_iter = 0;

  static constexpr int StandUp = 0;
  static constexpr int FoldLegs = 1;
  static constexpr int RollOver = 2;

  unsigned long long _state_iter;
  int _flag = FoldLegs;

  // JPos
  Vec3<T> fold_jpos[4];
  Vec3<T> stand_jpos[4];
  Vec3<T> rolling_jpos[4];
  Vec3<T> initial_jpos[4];
  Vec3<T> zero_vec3;

  Vec3<T> f_ff;

  // iteration setup
  //const int rollover_ramp_iter = 300;
  //const int rollover_settle_iter = 300;

  //const int fold_ramp_iter = 1000;
  //const int fold_settle_iter = 1000;

  //const int standup_ramp_iter = 500;
  //const int standup_settle_iter = 500;

  // 0.5 kHz
  const int rollover_ramp_iter = 150;
  const int rollover_settle_iter = 150;

  //const int fold_ramp_iter = 500;
  //const int fold_settle_iter = 500;
  const int fold_ramp_iter = 400;
  const int fold_settle_iter = 700;

  const int standup_ramp_iter = 250;
  const int standup_settle_iter = 250;

  void _RollOver(const int & iter);
  void _StandUp(const int & iter);
  void _FoldLegs(const int & iter);

  bool _UpsideDown();
  void _SetJPosInterPts(
      const size_t & curr_iter, size_t max_iter, int leg, 
      const Vec3<T> & ini, const Vec3<T> & fin);

};

#endif  // FSM_STATE_RECOVERY_STANDUP_H
