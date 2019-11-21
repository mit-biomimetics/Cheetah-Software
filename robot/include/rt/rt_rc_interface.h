/**
 * @file rt_rc_interface.h
 *
 */
#ifndef _RT_RC_INTERFACE
#define _RT_RC_INTERFACE

class rc_control_settings {
  public:
    double     mode;
    double     p_des[2]; // (x, y) -1 ~ 1
    double     height_variation; // -1 ~ 1
    double     v_des[3]; // -1 ~ 1 * (scale 0.5 ~ 1.5)
    double     rpy_des[3]; // -1 ~ 1
    double     omega_des[3]; // -1 ~ 1
    double     variable[3];
};


namespace RC_mode{
  constexpr int OFF = 0;
  constexpr int QP_STAND = 3;
  constexpr int BACKFLIP_PRE = 4;
  constexpr int BACKFLIP = 5;
  constexpr int VISION = 6;
  constexpr int LOCOMOTION = 11;
  constexpr int RECOVERY_STAND = 12;

  // Experiment Mode
  constexpr int TWO_LEG_STANCE_PRE = 20;
  constexpr int TWO_LEG_STANCE = 21;
};

void sbus_packet_complete();

void get_rc_control_settings(void* settings);
//void get_rc_channels(void* settings);

void* v_memcpy(void* dest, volatile void* src, size_t n);

float deadband(float command, float deadbandRegion, float minVal, float maxVal);

#endif
