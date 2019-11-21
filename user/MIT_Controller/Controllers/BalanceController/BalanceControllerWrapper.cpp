#include "BalanceControllerWrapper.h"
#include <cstdlib>
#include <iostream>
#include "BalanceController.hpp"

#ifdef __cplusplus
extern "C" {
#endif

static BalanceController* BalanceControllerObj = NULL;
double xOpt_force[12];

void balanceControl_init() {
  if (BalanceControllerObj == NULL) {
    BalanceControllerObj = new BalanceController();
    std::cout << "debug initialize balance controller\n";
  }
}

void balanceControl_set_desiredTrajectoryData(
    double* rpy_des_in, double* p_des_in, double* omegab_des_in,
    double* v_des_in)  //, double* vdot_des_in)
{
  BalanceControllerObj->set_desiredTrajectoryData(
      rpy_des_in, p_des_in, omegab_des_in, v_des_in);  //, vdot_des_in );
}

void balanceControl_set_desired_swing_pos(double* pFeet_des_in) {
  BalanceControllerObj->set_desired_swing_pos(pFeet_des_in);
}

void balanceControl_set_actual_swing_pos(double* pFeet_act_in) {
  BalanceControllerObj->set_actual_swing_pos(pFeet_act_in);
}

void balanceControl_set_PDgains(double* Kp_COM_in, double* Kd_COM_in,
                                double* Kp_Base_in, double* Kd_Base_in) {
  BalanceControllerObj->set_PDgains(Kp_COM_in, Kd_COM_in, Kp_Base_in,
                                    Kd_Base_in);
}

void balanceControl_set_wrench_weights(double* COM_weights_in,
                                       double* Base_weights_in) {
  BalanceControllerObj->set_wrench_weights(COM_weights_in, Base_weights_in);
}

void balanceControl_set_QP_options(double use_hard_constraint_pitch_in) {
  BalanceControllerObj->set_QP_options(use_hard_constraint_pitch_in);
}

void balanceControl_set_friction(double mu_in) {
  BalanceControllerObj->set_friction(mu_in);
}

void balanceControl_set_alpha_control(double alpha_control_in) {
  BalanceControllerObj->set_alpha_control(alpha_control_in);
}

void balanceControl_set_mass(double mass_in) {
  BalanceControllerObj->set_mass(mass_in);
}

void balanceControl_updateProblemData(double* xfb_in, double* pFeet_in,
                                      double yaw_in) {
  double p_des[3], p_act[3], v_des[3], v_act[3], O_err[3];

  BalanceControllerObj->updateProblemData(xfb_in, pFeet_in, p_des, p_act, v_des,
                                          v_act, O_err, yaw_in);
}

void balanceControl_SetContactData(double* contact_state_in,
                                   double* min_forces_in,
                                   double* max_forces_in) {
  BalanceControllerObj->SetContactData(contact_state_in, min_forces_in,
                                       max_forces_in);
}

void balanceControl_solveQP_nonThreaded() {
  BalanceControllerObj->solveQP_nonThreaded(xOpt_force);
}

void balanceControl_publish_data_lcm() {
  if (BalanceControllerObj != NULL) {
    BalanceControllerObj->publish_data_lcm();
  }
}

double balanceControl_get_fOpt_matlab(int index) {
  return xOpt_force[index - 1];
}

void bcPrint() { printf("balance controller here\n\n"); }

#ifdef __cplusplus
}
#endif
