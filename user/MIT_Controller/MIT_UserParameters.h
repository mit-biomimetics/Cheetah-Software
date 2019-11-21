#ifndef PROJECT_MITUSERPARAMETERS_H
#define PROJECT_MITUSERPARAMETERS_H

#include "ControlParameters/ControlParameters.h"

class MIT_UserParameters : public ControlParameters {
public:
  MIT_UserParameters()
      : ControlParameters("user-parameters"),
        INIT_PARAMETER(cmpc_gait),
        INIT_PARAMETER(cmpc_x_drag),
        INIT_PARAMETER(cmpc_use_sparse),
        INIT_PARAMETER(use_wbc),
        INIT_PARAMETER(cmpc_bonus_swing),
        INIT_PARAMETER(Kp_body),
        INIT_PARAMETER(Kd_body),
        INIT_PARAMETER(Kp_ori),
        INIT_PARAMETER(Kd_ori),
        INIT_PARAMETER(Kp_foot),
        INIT_PARAMETER(Kd_foot),
        INIT_PARAMETER(Kp_joint),
        INIT_PARAMETER(Kd_joint),
        //INIT_PARAMETER(Kp_joint_swing),
        //INIT_PARAMETER(Kd_joint_swing),
        INIT_PARAMETER(Q_pos),
        INIT_PARAMETER(Q_vel),
        INIT_PARAMETER(Q_ori),
        INIT_PARAMETER(Q_ang),
        INIT_PARAMETER(R_control),
        INIT_PARAMETER(R_prev),
        INIT_PARAMETER(two_leg_orient),
        INIT_PARAMETER(stance_legs),
        INIT_PARAMETER(use_jcqp),
        INIT_PARAMETER(jcqp_max_iter),
        INIT_PARAMETER(jcqp_rho),
        INIT_PARAMETER(jcqp_sigma),
        INIT_PARAMETER(jcqp_alpha),
        INIT_PARAMETER(jcqp_terminate),
        INIT_PARAMETER(Swing_Kp_cartesian),
        INIT_PARAMETER(Swing_Kd_cartesian),
        INIT_PARAMETER(Swing_Kp_joint),
        INIT_PARAMETER(Swing_Kd_joint),
        INIT_PARAMETER(Swing_step_offset),
        INIT_PARAMETER(Swing_traj_height),
        INIT_PARAMETER(Swing_use_tau_ff),
        INIT_PARAMETER(RPC_Q_p),
        INIT_PARAMETER(RPC_Q_theta),
        INIT_PARAMETER(RPC_Q_dp),
        INIT_PARAMETER(RPC_Q_dtheta),
        INIT_PARAMETER(RPC_R_r),
        INIT_PARAMETER(RPC_R_f),
        INIT_PARAMETER(RPC_H_r_trans),
        INIT_PARAMETER(RPC_H_r_rot),
        INIT_PARAMETER(RPC_H_theta0),
        INIT_PARAMETER(RPC_H_phi0),
        INIT_PARAMETER(RPC_mass),
        INIT_PARAMETER(RPC_inertia),
        INIT_PARAMETER(RPC_gravity),
        INIT_PARAMETER(RPC_mu),
        INIT_PARAMETER(RPC_filter),
        INIT_PARAMETER(RPC_use_pred_comp),
        INIT_PARAMETER(RPC_use_async_filt),
        INIT_PARAMETER(RPC_visualize_pred),
        INIT_PARAMETER(RPC_use_separate),
        INIT_PARAMETER(des_p),
        INIT_PARAMETER(des_theta),
        INIT_PARAMETER(des_dp),
        INIT_PARAMETER(des_dtheta),
        INIT_PARAMETER(des_theta_max),
        INIT_PARAMETER(des_dp_max),
        INIT_PARAMETER(des_dtheta_max),
        INIT_PARAMETER(gait_type),
        INIT_PARAMETER(gait_period_time),
        INIT_PARAMETER(gait_switching_phase),
        INIT_PARAMETER(gait_override),
        INIT_PARAMETER(gait_max_leg_angle),
        INIT_PARAMETER(gait_max_stance_time),
        INIT_PARAMETER(gait_min_stance_time)

  {}

  DECLARE_PARAMETER(double, cmpc_gait);
  DECLARE_PARAMETER(double, cmpc_x_drag);
  DECLARE_PARAMETER(double, cmpc_use_sparse);
  DECLARE_PARAMETER(double, use_wbc);
  DECLARE_PARAMETER(double, cmpc_bonus_swing);

  DECLARE_PARAMETER(Vec3<double>, Kp_body);
  DECLARE_PARAMETER(Vec3<double>, Kd_body);

  DECLARE_PARAMETER(Vec3<double>, Kp_ori);
  DECLARE_PARAMETER(Vec3<double>, Kd_ori);

  DECLARE_PARAMETER(Vec3<double>, Kp_foot);
  DECLARE_PARAMETER(Vec3<double>, Kd_foot);

  DECLARE_PARAMETER(Vec3<double>, Kp_joint);
  DECLARE_PARAMETER(Vec3<double>, Kd_joint);

  DECLARE_PARAMETER(Vec3<double>, Q_pos);
  DECLARE_PARAMETER(Vec3<double>, Q_vel);
  DECLARE_PARAMETER(Vec3<double>, Q_ori);
  DECLARE_PARAMETER(Vec3<double>, Q_ang);
  DECLARE_PARAMETER(double, R_control);
  DECLARE_PARAMETER(double, R_prev);
  DECLARE_PARAMETER(Vec3<double>, two_leg_orient);
  DECLARE_PARAMETER(double, stance_legs);

  //DECLARE_PARAMETER(Vec3<double>, Kp_joint_swing);
  //DECLARE_PARAMETER(Vec3<double>, Kd_joint_swing);

  DECLARE_PARAMETER(double, use_jcqp);
  DECLARE_PARAMETER(double, jcqp_max_iter);
  DECLARE_PARAMETER(double, jcqp_rho);
  DECLARE_PARAMETER(double, jcqp_sigma);
  DECLARE_PARAMETER(double, jcqp_alpha);
  DECLARE_PARAMETER(double, jcqp_terminate);

  // Swing leg parameters
  DECLARE_PARAMETER(Vec3<double>, Swing_Kp_cartesian);
  DECLARE_PARAMETER(Vec3<double>, Swing_Kd_cartesian);
  DECLARE_PARAMETER(Vec3<double>, Swing_Kp_joint);
  DECLARE_PARAMETER(Vec3<double>, Swing_Kd_joint);
  DECLARE_PARAMETER(Vec3<double>, Swing_step_offset);
  DECLARE_PARAMETER(double, Swing_traj_height);
  DECLARE_PARAMETER(double, Swing_use_tau_ff);


  // Parameters used for RPC
  DECLARE_PARAMETER(Vec3<double>, RPC_Q_p);
  DECLARE_PARAMETER(Vec3<double>, RPC_Q_theta);
  DECLARE_PARAMETER(Vec3<double>, RPC_Q_dp);
  DECLARE_PARAMETER(Vec3<double>, RPC_Q_dtheta);
  DECLARE_PARAMETER(Vec3<double>, RPC_R_r);
  DECLARE_PARAMETER(Vec3<double>, RPC_R_f);
  DECLARE_PARAMETER(Vec3<double>, RPC_H_r_trans);
  DECLARE_PARAMETER(Vec3<double>, RPC_H_r_rot);
  DECLARE_PARAMETER(Vec3<double>, RPC_H_theta0);
  DECLARE_PARAMETER(Vec3<double>, RPC_H_phi0);
  DECLARE_PARAMETER(double, RPC_mass);
  DECLARE_PARAMETER(Vec3<double>, RPC_inertia);
  DECLARE_PARAMETER(Vec3<double>, RPC_gravity);
  DECLARE_PARAMETER(double, RPC_mu);
  DECLARE_PARAMETER(Vec3<double>, RPC_filter);
  DECLARE_PARAMETER(double, RPC_use_pred_comp);
  DECLARE_PARAMETER(double, RPC_use_async_filt);
  DECLARE_PARAMETER(double, RPC_visualize_pred);
  DECLARE_PARAMETER(double, RPC_use_separate);

  // Desired states
  DECLARE_PARAMETER(Vec3<double>, des_p);
  DECLARE_PARAMETER(Vec3<double>, des_theta);
  DECLARE_PARAMETER(Vec3<double>, des_dp);
  DECLARE_PARAMETER(Vec3<double>, des_dtheta);
  DECLARE_PARAMETER(Vec3<double>, des_theta_max);
  DECLARE_PARAMETER(Vec3<double>, des_dp_max);
  DECLARE_PARAMETER(Vec3<double>, des_dtheta_max);

  // Gait Scheduler
  DECLARE_PARAMETER(double, gait_type);
  DECLARE_PARAMETER(double, gait_period_time);
  DECLARE_PARAMETER(double, gait_switching_phase);
  DECLARE_PARAMETER(double, gait_override);
  DECLARE_PARAMETER(double, gait_max_leg_angle);
  DECLARE_PARAMETER(double, gait_max_stance_time);
  DECLARE_PARAMETER(double, gait_min_stance_time);

};

#endif //PROJECT_MITUSERPARAMETERS_H
