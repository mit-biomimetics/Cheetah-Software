/*!
 * @file rt_sbus.h
 * @brief Communication with RC controller receiver
 */

#ifndef _rt_sbus
#define _rt_sbus

#include <stdint.h>

void unpack_sbus_data(uint8_t sbus_data[], uint16_t *channels);

int read_sbus_data(int port, uint8_t *sbus_data);

int read_sbus_channel(int channel);

int receive_sbus(int port);

int init_sbus(int is_simulator);


enum TaranisSwitchState {
  SWITCH_UP = 0,
  SWITCH_MIDDLE = 1,
  SWITCH_DOWN = 2,
};

struct Taranis_X7_data {
  TaranisSwitchState left_upper_switch, left_lower_left_switch, left_lower_right_switch,
  right_upper_switch, right_lower_left_switch, right_lower_right_switch;

  float left_stick[2];
  float right_stick[2];
  float knobs[2];
};

void update_taranis_x7(Taranis_X7_data* data);


#endif
