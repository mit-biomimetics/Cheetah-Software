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
        INIT_PARAMETER(two_leg_threshold),
        INIT_PARAMETER(two_leg_ramp),
        INIT_PARAMETER(stance_legs),
        INIT_PARAMETER(use_jcqp),
        INIT_PARAMETER(jcqp_max_iter),
        INIT_PARAMETER(jcqp_rho),
        INIT_PARAMETER(jcqp_sigma),
        INIT_PARAMETER(jcqp_alpha),
        INIT_PARAMETER(jcqp_terminate)

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
  DECLARE_PARAMETER(double, two_leg_threshold);
  DECLARE_PARAMETER(double, two_leg_ramp);
  DECLARE_PARAMETER(double, stance_legs);

  //DECLARE_PARAMETER(Vec3<double>, Kp_joint_swing);
  //DECLARE_PARAMETER(Vec3<double>, Kd_joint_swing);

  DECLARE_PARAMETER(double, use_jcqp);
  DECLARE_PARAMETER(double, jcqp_max_iter);
  DECLARE_PARAMETER(double, jcqp_rho);
  DECLARE_PARAMETER(double, jcqp_sigma);
  DECLARE_PARAMETER(double, jcqp_alpha);
  DECLARE_PARAMETER(double, jcqp_terminate);
};

#endif //PROJECT_MITUSERPARAMETERS_H
