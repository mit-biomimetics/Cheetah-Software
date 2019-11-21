#ifndef FSM_STATE_TWO_CONTACT_STAND_H
#define FSM_STATE_TWO_CONTACT_STAND_H

#include "FSM_State.h"
#include <Controllers/BalanceController/BalanceControllerVBL.hpp>
#include <Controllers/BalanceController/ReferenceGRF.hpp>
#include <Controllers/WBC_Ctrl/LocomotionCtrl/LocomotionCtrl.hpp>
#include <Controllers/WBC_Ctrl/WBC_Ctrl.hpp>

/**
 *
 */
template <typename T>
class FSM_State_TwoContactStand : public FSM_State<T> {
 public:
  FSM_State_TwoContactStand(ControlFSMData<T>* _controlFSMData);

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
  //
  void get_desired_state();
  void get_foot_locations();
  void set_ctrl_params();
  void get_current_state();
  void compute_force_ff();
  void debug_visualization();
  void liftLeg(int leg, Vec3<T> q, Vec3<T> qd);
  void rpyToR(Mat3<float> &R, double* rpy_in);
  void eulerToQuat(double* rpy_in, double* quat_in);
  void quatToEuler(double* quat_in, double* rpy_in);

  // Keep track of the control iterations
  int iter = 0;

  // Higher Level Robot Body Controllers
  BalanceControllerVBL balanceControllerVBL; 
  ReferenceGRF refGRF;

  // Whole body control
  WBC_Ctrl<T>* _wbc_ctrl;
  LocomotionCtrlData<T>* _wbc_data;

  // LQR Weights
  double x_weights[3], xdot_weights[3], R_weights[3], omega_weights[3];
  double alpha_control, beta_control;

  // Actual state of the body using estimator or cheater state - convert to Vec3 probably best for some of these
  double p_act[3], v_act[3], quat_act[4], se_xfb[13];
  double rpy_yaw_offset[3], quat_yaw_offset[4];

  // Foot positions accounting for initial yaw
  double pFeet_yaw[12], pFeet_yaw_des[12], pFeet_yaw_world[12];

  // Foot positions in world frame
  double pFeet[12], pFeet_world[12], pFeet_world_des[12];

  // Quadruped Model
  FloatingBaseModel<T> model;
  FBModelState<T> state;
  double Ig_in[3], mass_in, p_COM[3];
  Vec3<T> c_body, c_world;
  Mat18<float> H;
  DVec<T> G_ff;
  Mat3<T> Jleg;

  // Feet relative to COM
  Vec3<T> pFeetVec;
  Vec3<T> pFeetVecBody;

  // Desired state of the body
  double pFeet_des[12], pweight, convert = 3.14159/180;
  double p_des[3], v_des[3], a_des[3], rpy[3], omegaDes[3];

  // Whole Body controller
  Vec3<float> pBody_des;
  Vec3<float> vBody_des;
  Vec3<float> aBody_des;
  Vec3<float> pBody_RPY_des;
  Vec3<float> vBody_Ori_des;
  Vec3<float> pFoot_des[4];
  Vec3<float> vFoot_des[4];
  Vec3<float> aFoot_des[4];
  Vec3<float> Fr_des[4];
  Vec4<float> contact_state;

  // Joint positions for legs not in contact
  Vec3<T> q_lift_leg, qd_lift_leg;
  Mat3<T> kpMat, kdMat;

  // Contact Data
  double minForce, maxForce, mu_ctrl = 0.45;
  double minForces[4], maxForces[4];
  double contactStateScheduled[4] = {1.0, 1.0, 1.0, 1.0};
  Vec4<T> conPhase;
  double threshold = 0.2;
  double ramp_end = 0.5;

  // Control Input
  double f_ref_z[4], f_ref_world[12], fOpt[12];

  // Leg Impedance Control
  Vec3<double> impedance_kp;
  Vec3<double> impedance_kd;
  bool ESTOP = false;

  // Check to see if desired position has changed
  double p_des_prev[2] = {99.0,99.0};

  // Account for non-zero initial yaw
  double ini_yaw;
  Mat3<float> rBody_yaw; // rBody adjusted for initial yaw


};

#endif  // FSM_STATE_TWO_CONTACT_STAND_H
