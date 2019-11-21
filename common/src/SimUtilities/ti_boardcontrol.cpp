/*! @file ti_boardcontrol.cpp
 *  @brief TI Board Code, used to simulate the TI board
 *
 *  This is mostly a copy of the exact code that runs on the TI Board
 */

#include <algorithm>
#include <cmath>
#include <iostream>

#include "SimUtilities/ti_boardcontrol.h"

/*!
 * Initialize TI board
 */
void TI_BoardControl::init(float side_sign) {
  this->_side_sign = side_sign;
  data = &data_structure;
}

/*!
 * Set the leg link lengths (abad, hip, knee)
 */
void TI_BoardControl::set_link_lengths(float l1, float l2, float l3) {
  link_lengths_set = true;
  this->_l1 = l1;
  this->_l2 = l2;
  this->_l3 = l3;
}

/*!
 * Reset data for this board
 */
void TI_BoardControl::reset_ti_board_data() {
  std::fill_n(data_structure.tau_des, 3, 0);
  data_structure.loop_count_ti = 0;
  data_structure.ethercat_count_ti = 0;
  data_structure.microtime_ti = 0;
}

/*!
 * Reset command for this board
 */
void TI_BoardControl::reset_ti_board_command() {
  command.position_des[0] = 0;
  command.position_des[1] = 0;
  command.position_des[2] = 0;

  std::fill_n(command.velocity_des, 3, 0);
  std::fill_n(command.kp, 3, 0);
  std::fill_n(command.kd, 3, 0);
  std::fill_n(data_structure.tau_des, 3, 0);
}

/*!
 * Run TI board control
 */
void TI_BoardControl::run_ti_board_iteration() {
  if (!link_lengths_set) {
    std::cout << "[TI Board] Error.  Link lengths haven't been set.\n";
    return;
  }

  if (command.enable == 1 && command.max_torque > 1e-9) {
    impedanceControl(_side_sign, data->q, data->dq, command.position_des,
                     command.velocity_des, command.kp, command.kd,
                     command.force_ff, command.tau_ff, data->position,
                     data->velocity, data->force, data->tau_des);

    // joint position feedback.
    for(int i = 0; i < 3; i++) {
      data->tau_des[i] += (command.q_des[i] - data->q[i]) * command.kp_joint[i];
      data->tau_des[i] += (command.qd_des[i] - data->dq[i]) * command.kd_joint[i];
    }



    for (int i = 0; i < 3; ++i) {
      if (std::abs(data->tau_des[i]) > command.max_torque) {
        data->tau_des[i] =
            data->tau_des[i] / std::abs(data->tau_des[i]) * command.max_torque;
      }
    }
  } else {
    float Jacobian[3][3];
    kinematics(_side_sign, data->q, data->dq, data->position, data->velocity,
               Jacobian);
    data->force[0] = 0;
    data->force[1] = 0;
    data->force[2] = 0;
    data->tau_des[0] = 0;
    data->tau_des[1] = 0;
    data->tau_des[2] = 0;
  }
  data->loop_count_ti += 1;
}

/*!
 * Cartesian impedance control
 * @param side_sign : which side of robot
 * @param q  : joint position
 * @param dq : joint velocity
 * @param position_des : cartesian desired position
 * @param velocity_des : cartesian desired velocity
 * @param kp : cartesian stiffness
 * @param kd : cartesian damping
 * @param force_bias : cartesian force
 * @param torque_bias : torque
 * @param position : current position of foot
 * @param velocity : current velocity of foot
 * @param force : output force
 * @param torque : output torque
 */
void TI_BoardControl::impedanceControl(
    const float side_sign, const float q[3], const float dq[3],
    const float position_des[3], const float velocity_des[3], const float kp[3],
    const float kd[3], const float force_bias[3], const float torque_bias[3],
    float *position, float *velocity, float *force, float *torque) {
  float Jacobian[3][3];
  kinematics(side_sign, q, dq, position, velocity, Jacobian);

  force[0] = kp[0] * (position_des[0] - position[0]) +
             kd[0] * (velocity_des[0] - velocity[0]) + force_bias[0];
  force[1] = kp[1] * (position_des[1] - position[1]) +
             kd[1] * (velocity_des[1] - velocity[1]) + force_bias[1];
  force[2] = kp[2] * (position_des[2] - position[2]) +
             kd[2] * (velocity_des[2] - velocity[2]) + force_bias[2];

  torque[0] = (force[0] * Jacobian[0][0] + force[1] * Jacobian[1][0] +
               force[2] * Jacobian[2][0]) +
              torque_bias[0];
  torque[1] = (force[0] * Jacobian[0][1] + force[1] * Jacobian[1][1] +
               force[2] * Jacobian[2][1]) +
              torque_bias[1];
  torque[2] = (force[0] * Jacobian[0][2] + force[1] * Jacobian[1][2] +
               force[2] * Jacobian[2][2]) +
              torque_bias[2];
}

/*!
 * Kinematics and jacobian
 * @param side_sign : side of leg
 * @param q : joint positions
 * @param dq : joint velocities
 * @param p : position of foot
 * @param v : velocity of foot
 * @param J : jacobian output
 */
void TI_BoardControl::kinematics(const float side_sign, const float q[3],
                                 const float dq[3], float *p, float *v,
                                 float J[][3]) {

  const float s1 = sin(q[0]), s2 = sin(q[1]), s3 = sin(q[2]);
  const float c1 = cos(q[0]), c2 = cos(q[1]), c3 = cos(q[2]);
  const float c23 = c2 * c3 - s2 * s3;
  const float s23 = s2 * c3 + c2 * s3;

  p[0] = _l3 * s23 + _l2 * s2;
  p[1] = _l1 * side_sign * c1 + _l3 * (s1 * c23) + _l2 * c2 * s1;
  p[2] = _l1 * side_sign * s1 - _l3 * (c1 * c23) - _l2 * c1 * c2;

  J[0][0] = 0;
  J[0][1] = _l3 * c23 + _l2 * c2;
  J[0][2] = _l3 * c23;
  J[1][0] = _l3 * c1 * c23 + _l2 * c1 * c2 - _l1 * side_sign * s1;
  J[1][1] = -_l3 * s1 * s23 - _l2 * s1 * s2;
  J[1][2] = -_l3 * s1 * s23;
  J[2][0] = _l3 * s1 * c23 + _l2 * c2 * s1 + _l1 * side_sign * c1;
  J[2][1] = _l3 * c1 * s23 + _l2 * c1 * s2;
  J[2][2] = _l3 * c1 * s23;

  v[0] = 0;
  v[1] = 0;
  v[2] = 0;
  int i, j;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      v[i] += J[i][j] * dq[j];
    }
  }
}
