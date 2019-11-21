/*!
 * @file rt_sbus.cpp
 * @brief Communication with RC controller receiver
 */

#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <string>

#ifdef linux
#define termios asmtermios

#include <asm/termios.h>

#undef termios
#endif

#include <termios.h>

//#include "rt/rt_interface_lcm.h"
#include "rt/rt_sbus.h"
#include "rt/rt_serial.h"

pthread_mutex_t sbus_data_m;

uint16_t channels[18];
uint16_t channel_data[18];

/**@brief Name of SBUS serial port in simulator*/
#define K_SBUS_PORT_SIM "/dev/ttyUSB0"
/**@brief Name of SBUS serial port on the mini cheetah*/
#define K_SBUS_PORT_MC "/dev/ttyS4"

/*!
 * Unpack sbus message into channels
 */
void unpack_sbus_data(uint8_t sbus_data[], uint16_t *channels_) {
  if ((sbus_data[0] == 0xF) && (sbus_data[24] == 0x0)) {
    channels_[0] = ((sbus_data[1]) | ((sbus_data[2] & 0x7) << 8));
    channels_[1] = (sbus_data[2] >> 3) | ((sbus_data[3] & 0x3F) << 5);
    channels_[2] = ((sbus_data[3] & 0xC0) >> 6) | (sbus_data[4] << 2) |
                   ((sbus_data[5] & 0x1) << 10);
    channels_[3] = ((sbus_data[5] & 0xFE) >> 1) | ((sbus_data[6] & 0xF) << 7);
    channels_[4] = ((sbus_data[6] & 0xF0) >> 4) | ((sbus_data[7] & 0x7F) << 4);
    channels_[5] = ((sbus_data[7] & 0x80) >> 7) | (sbus_data[8] << 1) |
                   ((sbus_data[9] & 0x3) << 9);
    channels_[6] = ((sbus_data[9] & 0xFC) >> 2) | ((sbus_data[10] & 0x1F) << 6);
    channels_[7] = ((sbus_data[10] & 0xE0) >> 5) | (sbus_data[11] << 3);

    channels_[8] = ((sbus_data[12]) | ((sbus_data[13] & 0x7) << 8));
    channels_[9] = (sbus_data[13] >> 3) | ((sbus_data[14] & 0x3F) << 5);
    channels_[10] = ((sbus_data[14] & 0xC0) >> 6) | (sbus_data[15] << 2) |
                    ((sbus_data[16] & 0x1) << 10);
    channels_[11] =
        ((sbus_data[16] & 0xFE) >> 1) | ((sbus_data[17] & 0xF) << 7);
    channels_[12] =
        ((sbus_data[17] & 0xF0) >> 4) | ((sbus_data[18] & 0x7F) << 4);
    channels_[13] = ((sbus_data[18] & 0x80) >> 7) | (sbus_data[19] << 1) |
                    ((sbus_data[20] & 0x3) << 9);
    channels_[14] =
        ((sbus_data[20] & 0xFC) >> 2) | ((sbus_data[21] & 0x1F) << 6);
    channels_[15] = ((sbus_data[21] & 0xE0) >> 5) | (sbus_data[22] << 3);

    channels_[16] = (sbus_data[23] & 0x80) >> 7;
    channels_[17] = (sbus_data[23] & 0x40) >> 6;

    pthread_mutex_lock(&sbus_data_m);
    for (int i = 0; i < 18; i++) {
       // printf("[%d] %d ", i, channels_[i]);
      channel_data[i] = channels_[i];
    }
    //printf("\n\n");
    pthread_mutex_unlock(&sbus_data_m);

    // for(int i = 0; i < 18; i++) {
    //   printf("[%2d] %04d ", i, channel_data[i]);
    // }
    // printf("\n");
    // for(int i = 0; i < 24; i++) {
    //   printf("[%2d] %04d ", i, sbus_data[i]);
    // }
    // printf("\n\n");

  } else {
    // printf("Bad Packet\n");
  }
}

/*!
 * Read data from serial port
 */
int read_sbus_data(int port, uint8_t *sbus_data) {
  uint8_t packet_full = 0;
  uint8_t read_byte[1] = {0};
  int timeout_counter = 0;
  // int n = read(fd1, read_buff, sizeof(read_buff));
  while ((!packet_full) && (timeout_counter < 50)) {
    timeout_counter++;
    // Read a byte //
    while(read(port, read_byte, sizeof(read_byte)) != 1) {
      
    }

    // Shift the buffer //
    for (int i = 0; i < 24; i++) {
      sbus_data[i] = sbus_data[i + 1];
    }
    sbus_data[24] = read_byte[0];

    // Check for the correct start and stop bytes ///
    if ((sbus_data[0] == 15) && (sbus_data[24] == 0)) {
      // unpack_sbus_data(sbus_data_buff, channels);
      packet_full = 1;
    }
  }
  return packet_full;
}

/*!
 * Get sbus channel
 */
int read_sbus_channel(int channel) {
  pthread_mutex_lock(&sbus_data_m);
  int value = channel_data[channel];
  pthread_mutex_unlock(&sbus_data_m);
  return value;
}

/*!
 * Receive serial and find packets
 */
int receive_sbus(int port) {
  uint16_t read_buff[25] = {0};
  int x = read_sbus_data(port, (uint8_t *)read_buff);
  if (x) {
    unpack_sbus_data((uint8_t *)read_buff, channels);
  } else {
    printf("SBUS tried read 50 bytes without seeing a packet\n");
  }
  return x;
}

/*!
 * Initialize SBUS serial port
 */
int init_sbus(int is_simulator) {
  // char *port1;
  std::string port1;
  if (is_simulator) {
    port1 = K_SBUS_PORT_SIM;
  } else {
    port1 = K_SBUS_PORT_MC;
  }

  if (pthread_mutex_init(&sbus_data_m, NULL) != 0) {
    printf("Failed to initialize sbus data mutex.\n");
  }

  int fd1 = open(port1.c_str(), O_RDWR | O_NOCTTY | O_SYNC);
  if (fd1 < 0) {
    printf("Error opening %s: %s\n", port1.c_str(), strerror(errno));
  } else {
    init_serial_for_sbus(fd1, 100000);
#ifdef linux
    //set_interface_attribs_custom_baud(fd1, 100000, 0, 0);
#endif
  }
  return fd1;
}

static float scale_joystick(uint16_t in) {
  return (in - 172) * 2.f / (1811.f - 172.f) - 1.f;
}

static TaranisSwitchState map_switch(uint16_t in) {
  switch(in) {
    case 1811:
      return TaranisSwitchState::SWITCH_DOWN;
    case 992:
      return TaranisSwitchState::SWITCH_MIDDLE;
    case 172:
      return TaranisSwitchState::SWITCH_UP;
    default:
      printf("[SBUS] switch returned bad value %d\n", in);
      return TaranisSwitchState::SWITCH_UP;
  }
}


void update_taranis_x7(Taranis_X7_data* data) {
  pthread_mutex_lock(&sbus_data_m);
  data->left_stick[0] = scale_joystick(channel_data[3]);
  data->left_stick[1] = scale_joystick(channel_data[0]);
  data->right_stick[0] = scale_joystick(channel_data[1]);
  data->right_stick[1] = scale_joystick(channel_data[2]);
  data->left_lower_left_switch = map_switch(channel_data[4]);
  data->left_lower_right_switch = map_switch(channel_data[5]);
  data->left_upper_switch = map_switch(channel_data[6]);
  data->right_lower_left_switch = map_switch(channel_data[7]);
  data->right_lower_right_switch = map_switch(channel_data[8]);
  data->right_upper_switch = map_switch(channel_data[9]);
  data->knobs[0] = scale_joystick(channel_data[10]);
  data->knobs[1] = scale_joystick(channel_data[11]);

  pthread_mutex_unlock(&sbus_data_m);
}
