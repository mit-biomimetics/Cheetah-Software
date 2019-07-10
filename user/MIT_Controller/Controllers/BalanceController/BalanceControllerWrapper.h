#ifndef BALANCECONTROLLERWRAPPER_H
#define BALANCECONTROLLERWRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

void balanceControl_init();
void balanceControl_set_desiredTrajectoryData(
    double* rpy_des_in, double* p_des_in, double* omegab_des_in,
    double* v_des_in);  //, double *vdot_des_in);
void balanceControl_set_desired_swing_pos(double* pFeet_des_in);
void balanceControl_set_actual_swing_pos(double* pFeet_act_in);
void balanceControl_set_PDgains(double* Kp_COM_in, double* Kd_COM_in,
                                double* Kp_Base_in, double* Kd_Base_in);
void balanceControl_set_wrench_weights(double* COM_weights_in,
                                       double* Base_weights_in);
void balanceControl_set_QP_options(double use_hard_constraint_pitch_in);
void balanceControl_updateProblemData(double* xfb_in, double* pFeet_in,
                                      double yaw_in);
void balanceControl_SetContactData(double* contact_state_in,
                                   double* min_forces_in,
                                   double* max_forces_in);
void balanceControl_solveQP_nonThreaded();
void balanceControl_set_friction(double mu_in);
void balanceControl_set_alpha_control(double alpha_control_in);
void balanceControl_set_mass(double mass_in);
void balanceControl_publish_data_lcm();
double balanceControl_get_fOpt_matlab(int index);

void bcPrint();

#ifdef __cplusplus
}
#endif

#endif
