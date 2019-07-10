/*! @file SpineBoard.cpp
 *  @brief Spine Board Code, used to simulate the SpineBoard.
 */

#include "SimUtilities/SpineBoard.h"
#include <stdio.h>

void SpineBoard::init(float sideSign, s32 board) {
  printf(
      "[Cheetah Control] Hello! I am the SPINE board for leg %d (0-indexed) "
      "with side sign %.3f\n",
      board, sideSign);
  this->board_num = board;
  this->side_sign = sideSign;
}

void SpineBoard::resetData() {
  if (data == nullptr) {
    printf(
        "[ERROR: SPINE board] reset_spine_board_data called when "
        "cheetahlcm_spi_data_t* was null\n");
    return;
  }

  data->flags[board_num] = 0;
  data->qd_abad[board_num] = 0.f;
  data->qd_hip[board_num] = 0.f;
  data->qd_knee[board_num] = 0.f;
  data->q_abad[board_num] = 0.f;
  data->q_hip[board_num] = 0.f;
  data->q_knee[board_num] = 0.f;
  data->spi_driver_status = 0;
}

void SpineBoard::resetCommand() {
  if (cmd == nullptr) {
    printf(
        "[ERROR: SPINE board] reset_spine_board_command called when "
        "cheetahlcm_spi_command_t* was null\n");
    return;
  }

  cmd->flags[board_num] = 0;
  cmd->kd_abad[board_num] = 0.f;
  cmd->kd_hip[board_num] = 0.f;
  cmd->kd_knee[board_num] = 0.f;
  cmd->kp_abad[board_num] = 0.f;
  cmd->kp_hip[board_num] = 0.f;
  cmd->kp_knee[board_num] = 0.f;
  cmd->qd_des_abad[board_num] = 0.f;
  cmd->qd_des_hip[board_num] = 0.f;
  cmd->qd_des_knee[board_num] = 0.f;
  cmd->q_des_abad[board_num] = 0.f;
  cmd->q_des_hip[board_num] = 0.f;
  cmd->q_des_knee[board_num] = 0.f;
  cmd->tau_abad_ff[board_num] = 0.f;
  cmd->tau_hip_ff[board_num] = 0.f;
  cmd->tau_hip_ff[board_num] = 0.f;
  cmd->tau_knee_ff[board_num] = 0.f;
}

void SpineBoard::run() {
  iter_counter++;
  if (cmd == nullptr || data == nullptr) {
    printf(
        "[ERROR: SPINE board] run_spine_board_iteration called with null "
        "command or data!\n");
    torque_out[0] = 0.f;
    torque_out[1] = 0.f;
    torque_out[2] = 0.f;
    return;
  }

  /// Check abad softstop ///
  if (data->q_abad[board_num] > q_limit_p[0]) {
    torque_out[0] = kp_softstop * (q_limit_p[0] - data->q_abad[board_num]) -
                    kd_softstop * (data->qd_abad[board_num]) +
                    cmd->tau_abad_ff[board_num];
  } else if (data->q_abad[board_num] < q_limit_n[0]) {
    torque_out[0] = kp_softstop * (q_limit_n[0] - data->q_abad[board_num]) -
                    kd_softstop * (data->qd_abad[board_num]) +
                    cmd->tau_abad_ff[board_num];
  } else {
    torque_out[0] = cmd->kp_abad[board_num] *
                        (cmd->q_des_abad[board_num] - data->q_abad[board_num]) +
                    cmd->kd_abad[board_num] * (cmd->qd_des_abad[board_num] -
                                               data->qd_abad[board_num]) +
                    cmd->tau_abad_ff[board_num];
  }

  /// Check hip softstop ///
  if (data->q_hip[board_num] > q_limit_p[1]) {
    torque_out[1] = kp_softstop * (q_limit_p[1] - data->q_hip[board_num]) -
                    kd_softstop * (data->qd_hip[board_num]) +
                    cmd->tau_hip_ff[board_num];
  } else if (data->q_hip[board_num] < q_limit_n[1]) {
    torque_out[1] = kp_softstop * (q_limit_n[1] - data->q_hip[board_num]) -
                    kd_softstop * (data->qd_hip[board_num]) +
                    cmd->tau_hip_ff[board_num];
  } else {
    torque_out[1] = cmd->kp_hip[board_num] *
                        (cmd->q_des_hip[board_num] - data->q_hip[board_num]) +
                    cmd->kd_hip[board_num] *
                        (cmd->qd_des_hip[board_num] - data->qd_hip[board_num]) +
                    cmd->tau_hip_ff[board_num];
  }

  /// No knee softstop right now ///
  torque_out[2] = cmd->kp_knee[board_num] *
                      (cmd->q_des_knee[board_num] - data->q_knee[board_num]) +
                  cmd->kd_knee[board_num] *
                      (cmd->qd_des_knee[board_num] - data->qd_knee[board_num]) +
                  cmd->tau_knee_ff[board_num];

  const float* torque_limits = disabled_torque;

  if (cmd->flags[board_num] & 0b1) {
    if (cmd->flags[board_num] & 0b10)
      torque_limits = wimp_torque;
    else
      torque_limits = max_torque;
  }

  for (int i = 0; i < 3; i++) {
    if (torque_out[i] > torque_limits[i]) torque_out[i] = torque_limits[i];
    if (torque_out[i] < -torque_limits[i]) torque_out[i] = -torque_limits[i];
  }

  //    if(! (iter_counter % 100))
  //    {
  //        printf("SPINE BOARD %d: flags: 0x%hx, torque out: %.3f %.3f %.3f\n",
  //        board_num, cmd->flags[board_num],
  //               torque_out[0], torque_out[1], torque_out[2]);
  //    }
}
