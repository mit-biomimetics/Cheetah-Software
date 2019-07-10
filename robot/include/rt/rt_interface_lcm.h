/**
 * @file rt_interface_lcm.h
 *
 */

#ifndef _RT_INTERFACE_LCM
#define _RT_INTERFACE_LCM

#include <lcm/lcm-cpp.hpp>
class gui_main_control_settings_t;
class rc_channels_t;

class Handler {
 public:
  /**
   * @brief      Handler for main control settings from MATLAB interface. Never
   * call this yourself, leave it to LCM
   *
   * @param[in]  rbuf     The receive buffer
   * @param[in]  channel  The name of the LCM Channel
   * @param[in]  msg      Pointer to the LCM message
   */
  void main_control_settings_handler(const lcm::ReceiveBuffer* rbuf,
                                     const std::string& channel,
                                     const gui_main_control_settings_t* msg);

  void rc_channels_handler(const lcm::ReceiveBuffer* rbuf,
                           const std::string& channel,
                           const rc_channels_t* msg);
};

void sbus_packet_complete();

void get_main_control_settings(void* settings);
int get_iterations_since_last_lcm();
void get_rc_channels(void* settings);

void init_interface_lcm(lcm::LCM* main_lcm);
void control_iteration_lcm();
void* v_memcpy(void* dest, volatile void* src, size_t n);

#endif
