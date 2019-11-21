/*! @file ti_boardcontrol.h
 *  @brief TI Board Code, used to simulate the TI board
 *
 *  This is mostly a copy of the exact code that runs on the TI Board
 */

#ifndef TI_BOARDCONTROL_H
#define TI_BOARDCONTROL_H

#include "cTypes.h"

/*!
 * Command sent to TI board
 */
struct TiBoardCommand {
  float position_des[3];
  float velocity_des[3];
  float kp[3];
  float kd[3];

  u32 enable; // uint32t
  u32 zero;   // uint32t

  float force_ff[3];
  float tau_ff[3];

  float zero_offset[3];

  float q_des[3];
  float qd_des[3];
  float kp_joint[3];
  float kd_joint[3];

  float max_torque; // todo this isn't used.

};

/*!
 * Data received from TI board
 */
struct TiBoardData {
  float position[3];
  float velocity[3];
  float force[3];
  float q[3];
  float dq[3];
  float tau[3];
  float tau_des[3]; // todo this is encoder ticks?
  u32 loop_count_ti;
  u32 ethercat_count_ti;
  u32 microtime_ti;
};

/*!
 * TI Board control code
 */
class TI_BoardControl {
 public:
  TI_BoardControl() = default;
  void init(float side_sign);
  void run_ti_board_iteration();
  void reset_ti_board_data();
  void reset_ti_board_command();
  void set_link_lengths(float l1, float l2, float l3);
  TiBoardCommand command;
  TiBoardData data_structure;

  // for Ben's TI board code that uses pointers
  TiBoardData* data;

 private:
  void kinematics(const float side_sign, const float q[3], const float dq[3],
                  float* p, float* v, float J[][3]);
  void impedanceControl(const float side_sign, const float q[3],
                        const float dq[3], const float position_des[3],
                        const float velocity_des[3], const float kp[3],
                        const float kd[3], const float force_bias[3],
                        const float torque_bias[3], float* position,
                        float* velocity, float* force, float* torque);

  float _side_sign;
  float _l1, _l2, _l3;
  bool link_lengths_set = false;
};

#endif  // TI_BOARDCONTROL_H
