/*! @file SpineBoard.h
 *  @brief Spine Board Code, used to simulate the SpineBoard.
 *
 *  This is mostly a copy of the exact code that runs on the SpineBoard
 */

#ifndef PROJECT_SPINEBOARD_H
#define PROJECT_SPINEBOARD_H

#include "cTypes.h"

/*!
 * Command to spine board
 */
struct SpiCommand {
  float q_des_abad[4];
  float q_des_hip[4];
  float q_des_knee[4];

  float qd_des_abad[4];
  float qd_des_hip[4];
  float qd_des_knee[4];

  float kp_abad[4];
  float kp_hip[4];
  float kp_knee[4];

  float kd_abad[4];
  float kd_hip[4];
  float kd_knee[4];

  float tau_abad_ff[4];
  float tau_hip_ff[4];
  float tau_knee_ff[4];

  int32_t flags[4];
};

/*!
 * Data from spine board
 */
struct SpiData {
  float q_abad[4];
  float q_hip[4];
  float q_knee[4];
  float qd_abad[4];
  float qd_hip[4];
  float qd_knee[4];
  int32_t flags[4];
  int32_t spi_driver_status;
};

/*!
 * Spine board control logic
 */
class SpineBoard {
 public:
  SpineBoard() {}
  void init(float side_sign, s32 board);
  void run();
  void resetData();
  void resetCommand();
  SpiCommand* cmd = nullptr;
  SpiData* data = nullptr;
  float torque_out[3];

 private:
  float side_sign;
  s32 board_num;
  const float max_torque[3] = {17.f, 17.f, 26.f};  // TODO CHECK WITH BEN
  const float wimp_torque[3] = {6.f, 6.f, 6.f};    // TODO CHECK WITH BEN
  const float disabled_torque[3] = {0.f, 0.f, 0.f};
  const float q_limit_p[3] = {1.5f, 5.0f, 0.f};
  const float q_limit_n[3] = {-1.5f, -5.0f, 0.f};
  const float kp_softstop = 100.f;
  const float kd_softstop = 0.4f;
  s32 iter_counter = 0;
};

#endif  // PROJECT_SPINEBOARD_H
