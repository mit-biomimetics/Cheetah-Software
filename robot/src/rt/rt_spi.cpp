/*!
 * @file rt_spi.h
 * @brief SPI communication to spine board
 */
#ifdef linux

#include <byteswap.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <string.h>

#include <linux/spi/spidev.h>
#include "rt/rt_spi.h"
#include <lcm/lcm-cpp.hpp>

unsigned char spi_mode = SPI_MODE_0;
unsigned char spi_bits_per_word = 8;
unsigned int spi_speed = 6000000;
uint8_t lsb = 0x01;

int spi_1_fd = -1;
int spi_2_fd = -1;

int spi_open();

static spine_cmd_t g_spine_cmd;
static spine_data_t g_spine_data;

spi_command_t spi_command_drv;
spi_data_t spi_data_drv;
spi_torque_t spi_torque;

pthread_mutex_t spi_mutex;

const float max_torque[3] = {17.f, 17.f, 26.f};  // TODO CHECK WITH BEN
const float wimp_torque[3] = {6.f, 6.f, 6.f};    // TODO CHECK WITH BEN
const float disabled_torque[3] = {0.f, 0.f, 0.f};

// only used for actual robot
const float abad_side_sign[4] = {-1.f, -1.f, 1.f, 1.f};
const float hip_side_sign[4] = {-1.f, 1.f, -1.f, 1.f};
const float knee_side_sign[4] = {-.6429f, .6429f, -.6429f, .6429f};

// only used for actual robot
const float abad_offset[4] = {0.f, 0.f, 0.f, 0.f};
const float hip_offset[4] = {M_PI / 2.f, -M_PI / 2.f, -M_PI / 2.f, M_PI / 2.f};
const float knee_offset[4] = {K_KNEE_OFFSET_POS, -K_KNEE_OFFSET_POS,
                              -K_KNEE_OFFSET_POS, K_KNEE_OFFSET_POS};

/*!
 * Compute SPI message checksum
 * @param data : input
 * @param len : length (in 32-bit words)
 * @return
 */
uint32_t xor_checksum(uint32_t *data, size_t len) {
  uint32_t t = 0;
  for (size_t i = 0; i < len; i++) t = t ^ data[i];
  return t;
}

/*!
 * Emulate the spi board to estimate the torque.
 */
void fake_spine_control(spi_command_t *cmd, spi_data_t *data,
                        spi_torque_t *torque_out, int board_num) {
  torque_out->tau_abad[board_num] =
      cmd->kp_abad[board_num] *
          (cmd->q_des_abad[board_num] - data->q_abad[board_num]) +
      cmd->kd_abad[board_num] *
          (cmd->qd_des_abad[board_num] - data->qd_abad[board_num]) +
      cmd->tau_abad_ff[board_num];

  torque_out->tau_hip[board_num] =
      cmd->kp_hip[board_num] *
          (cmd->q_des_hip[board_num] - data->q_hip[board_num]) +
      cmd->kd_hip[board_num] *
          (cmd->qd_des_hip[board_num] - data->qd_hip[board_num]) +
      cmd->tau_hip_ff[board_num];

  torque_out->tau_knee[board_num] =
      cmd->kp_knee[board_num] *
          (cmd->q_des_knee[board_num] - data->q_knee[board_num]) +
      cmd->kd_knee[board_num] *
          (cmd->qd_des_knee[board_num] - data->qd_knee[board_num]) +
      cmd->tau_knee_ff[board_num];

  const float *torque_limits = disabled_torque;

  if (cmd->flags[board_num] & 0b1) {
    if (cmd->flags[board_num] & 0b10)
      torque_limits = wimp_torque;
    else
      torque_limits = max_torque;
  }

  if (torque_out->tau_abad[board_num] > torque_limits[0])
    torque_out->tau_abad[board_num] = torque_limits[0];
  if (torque_out->tau_abad[board_num] < -torque_limits[0])
    torque_out->tau_abad[board_num] = -torque_limits[0];

  if (torque_out->tau_hip[board_num] > torque_limits[1])
    torque_out->tau_hip[board_num] = torque_limits[1];
  if (torque_out->tau_hip[board_num] < -torque_limits[1])
    torque_out->tau_hip[board_num] = -torque_limits[1];

  if (torque_out->tau_knee[board_num] > torque_limits[2])
    torque_out->tau_knee[board_num] = torque_limits[2];
  if (torque_out->tau_knee[board_num] < -torque_limits[2])
    torque_out->tau_knee[board_num] = -torque_limits[2];
}

/*!
 * Initialize SPI
 */
void init_spi() {
  // check sizes:
  size_t command_size = sizeof(spi_command_t);
  size_t data_size = sizeof(spi_data_t);

  memset(&spi_command_drv, 0, sizeof(spi_command_drv));
  memset(&spi_data_drv, 0, sizeof(spi_data_drv));

  if (pthread_mutex_init(&spi_mutex, NULL) != 0)
    printf("[ERROR: RT SPI] Failed to create spi data mutex\n");

  if (command_size != K_EXPECTED_COMMAND_SIZE) {
    printf("[RT SPI] Error command size is %ld, expected %d\n", command_size,
           K_EXPECTED_COMMAND_SIZE);
  } else
    printf("[RT SPI] command size good\n");

  if (data_size != K_EXPECTED_DATA_SIZE) {
    printf("[RT SPI] Error data size is %ld, expected %d\n", data_size,
           K_EXPECTED_DATA_SIZE);
  } else
    printf("[RT SPI] data size good\n");

  printf("[RT SPI] Open\n");
  spi_open();
}

/*!
 * Open SPI device
 */
int spi_open() {
  int rv = 0;
  spi_1_fd = open("/dev/spidev2.0", O_RDWR);
  if (spi_1_fd < 0) perror("[ERROR] Couldn't open spidev 2.0");
  spi_2_fd = open("/dev/spidev2.1", O_RDWR);
  if (spi_2_fd < 0) perror("[ERROR] Couldn't open spidev 2.1");

  rv = ioctl(spi_1_fd, SPI_IOC_WR_MODE, &spi_mode);
  if (rv < 0) perror("[ERROR] ioctl spi_ioc_wr_mode (1)");

  rv = ioctl(spi_2_fd, SPI_IOC_WR_MODE, &spi_mode);
  if (rv < 0) perror("[ERROR] ioctl spi_ioc_wr_mode (2)");

  rv = ioctl(spi_1_fd, SPI_IOC_RD_MODE, &spi_mode);
  if (rv < 0) perror("[ERROR] ioctl spi_ioc_rd_mode (1)");

  rv = ioctl(spi_2_fd, SPI_IOC_RD_MODE, &spi_mode);
  if (rv < 0) perror("[ERROR] ioctl spi_ioc_rd_mode (2)");

  rv = ioctl(spi_1_fd, SPI_IOC_WR_BITS_PER_WORD, &spi_bits_per_word);
  if (rv < 0) perror("[ERROR] ioctl spi_ioc_wr_bits_per_word (1)");

  rv = ioctl(spi_2_fd, SPI_IOC_WR_BITS_PER_WORD, &spi_bits_per_word);
  if (rv < 0) perror("[ERROR] ioctl spi_ioc_wr_bits_per_word (2)");

  rv = ioctl(spi_1_fd, SPI_IOC_RD_BITS_PER_WORD, &spi_bits_per_word);
  if (rv < 0) perror("[ERROR] ioctl spi_ioc_rd_bits_per_word (1)");

  rv = ioctl(spi_2_fd, SPI_IOC_RD_BITS_PER_WORD, &spi_bits_per_word);
  if (rv < 0) perror("[ERROR] ioctl spi_ioc_rd_bits_per_word (2)");

  rv = ioctl(spi_1_fd, SPI_IOC_WR_MAX_SPEED_HZ, &spi_speed);
  if (rv < 0) perror("[ERROR] ioctl spi_ioc_wr_max_speed_hz (1)");
  rv = ioctl(spi_2_fd, SPI_IOC_WR_MAX_SPEED_HZ, &spi_speed);
  if (rv < 0) perror("[ERROR] ioctl spi_ioc_wr_max_speed_hz (2)");

  rv = ioctl(spi_1_fd, SPI_IOC_RD_MAX_SPEED_HZ, &spi_speed);
  if (rv < 0) perror("[ERROR] ioctl spi_ioc_rd_max_speed_hz (1)");
  rv = ioctl(spi_2_fd, SPI_IOC_RD_MAX_SPEED_HZ, &spi_speed);
  if (rv < 0) perror("[ERROR] ioctl spi_ioc_rd_max_speed_hz (2)");

  rv = ioctl(spi_1_fd, SPI_IOC_RD_LSB_FIRST, &lsb);
  if (rv < 0) perror("[ERROR] ioctl spi_ioc_rd_lsb_first (1)");

  rv = ioctl(spi_2_fd, SPI_IOC_RD_LSB_FIRST, &lsb);
  if (rv < 0) perror("[ERROR] ioctl spi_ioc_rd_lsb_first (2)");
  return rv;
}

int spi_driver_iterations = 0;

/*!
 * convert spi command to spine_cmd_t
 */
void spi_to_spine(spi_command_t *cmd, spine_cmd_t *spine_cmd, int leg_0) {
  for (int i = 0; i < 2; i++) {
    // spine_cmd->q_des_abad[i] = (cmd->q_des_abad[i+leg_0] +
    // abad_offset[i+leg_0]) * abad_side_sign[i+leg_0]; spine_cmd->q_des_hip[i]
    // = (cmd->q_des_hip[i+leg_0] + hip_offset[i+leg_0]) *
    // hip_side_sign[i+leg_0]; spine_cmd->q_des_knee[i] =
    // (cmd->q_des_knee[i+leg_0] + knee_offset[i+leg_0]) /
    // knee_side_sign[i+leg_0];
    spine_cmd->q_des_abad[i] =
        (cmd->q_des_abad[i + leg_0] * abad_side_sign[i + leg_0]) +
        abad_offset[i + leg_0];
    spine_cmd->q_des_hip[i] =
        (cmd->q_des_hip[i + leg_0] * hip_side_sign[i + leg_0]) +
        hip_offset[i + leg_0];
    spine_cmd->q_des_knee[i] =
        (cmd->q_des_knee[i + leg_0] / knee_side_sign[i + leg_0]) +
        knee_offset[i + leg_0];

    spine_cmd->qd_des_abad[i] =
        cmd->qd_des_abad[i + leg_0] * abad_side_sign[i + leg_0];
    spine_cmd->qd_des_hip[i] =
        cmd->qd_des_hip[i + leg_0] * hip_side_sign[i + leg_0];
    spine_cmd->qd_des_knee[i] =
        cmd->qd_des_knee[i + leg_0] / knee_side_sign[i + leg_0];

    spine_cmd->kp_abad[i] = cmd->kp_abad[i + leg_0];
    spine_cmd->kp_hip[i] = cmd->kp_hip[i + leg_0];
    spine_cmd->kp_knee[i] = cmd->kp_knee[i + leg_0];

    spine_cmd->kd_abad[i] = cmd->kd_abad[i + leg_0];
    spine_cmd->kd_hip[i] = cmd->kd_hip[i + leg_0];
    spine_cmd->kd_knee[i] = cmd->kd_knee[i + leg_0];

    spine_cmd->tau_abad_ff[i] =
        cmd->tau_abad_ff[i + leg_0] * abad_side_sign[i + leg_0];
    spine_cmd->tau_hip_ff[i] =
        cmd->tau_hip_ff[i + leg_0] * hip_side_sign[i + leg_0];
    spine_cmd->tau_knee_ff[i] =
        cmd->tau_knee_ff[i + leg_0] * knee_side_sign[i + leg_0];

    spine_cmd->flags[i] = cmd->flags[i + leg_0];
  }
  spine_cmd->checksum = xor_checksum((uint32_t *)spine_cmd, 32);
}

/*!
 * convert spine_data_t to spi data
 */
void spine_to_spi(spi_data_t *data, spine_data_t *spine_data, int leg_0) {
  for (int i = 0; i < 2; i++) {
    data->q_abad[i + leg_0] = (spine_data->q_abad[i] - abad_offset[i + leg_0]) *
                              abad_side_sign[i + leg_0];
    data->q_hip[i + leg_0] = (spine_data->q_hip[i] - hip_offset[i + leg_0]) *
                             hip_side_sign[i + leg_0];
    data->q_knee[i + leg_0] = (spine_data->q_knee[i] - knee_offset[i + leg_0]) *
                              knee_side_sign[i + leg_0];

    data->qd_abad[i + leg_0] =
        spine_data->qd_abad[i] * abad_side_sign[i + leg_0];
    data->qd_hip[i + leg_0] = spine_data->qd_hip[i] * hip_side_sign[i + leg_0];
    data->qd_knee[i + leg_0] =
        spine_data->qd_knee[i] * knee_side_sign[i + leg_0];

    data->flags[i + leg_0] = spine_data->flags[i];
  }

  uint32_t calc_checksum = xor_checksum((uint32_t *)spine_data, 14);
  if (calc_checksum != (uint32_t)spine_data->checksum)
    printf("SPI ERROR BAD CHECKSUM GOT 0x%hx EXPECTED 0x%hx\n", calc_checksum,
           spine_data->checksum);
}

/*!
 * send receive data and command from spine
 */
void spi_send_receive(spi_command_t *command, spi_data_t *data) {
  // update driver status flag
  spi_driver_iterations++;
  data->spi_driver_status = spi_driver_iterations << 16;

  // transmit and receive buffers
  uint16_t tx_buf[K_WORDS_PER_MESSAGE];
  uint16_t rx_buf[K_WORDS_PER_MESSAGE];

  for (int spi_board = 0; spi_board < 2; spi_board++) {
    // copy command into spine type:
    spi_to_spine(command, &g_spine_cmd, spi_board * 2);

    // pointers to command/data spine array
    uint16_t *cmd_d = (uint16_t *)&g_spine_cmd;
    uint16_t *data_d = (uint16_t *)&g_spine_data;

    // zero rx buffer
    memset(rx_buf, 0, K_WORDS_PER_MESSAGE * sizeof(uint16_t));

    // copy into tx buffer flipping bytes
    for (int i = 0; i < K_WORDS_PER_MESSAGE; i++)
      tx_buf[i] = (cmd_d[i] >> 8) + ((cmd_d[i] & 0xff) << 8);
    // tx_buf[i] = __bswap_16(cmd_d[i]);

    // each word is two bytes long
    size_t word_len = 2;  // 16 bit word

    // spi message struct
    struct spi_ioc_transfer spi_message[1];

    // zero message struct.
    memset(spi_message, 0, 1 * sizeof(struct spi_ioc_transfer));

    // set up message struct
    for (int i = 0; i < 1; i++) {
      spi_message[i].bits_per_word = spi_bits_per_word;
      spi_message[i].cs_change = 1;
      spi_message[i].delay_usecs = 0;
      spi_message[i].len = word_len * 66;
      spi_message[i].rx_buf = (uint64_t)rx_buf;
      spi_message[i].tx_buf = (uint64_t)tx_buf;
    }

    // do spi communication
    int rv = ioctl(spi_board == 0 ? spi_1_fd : spi_2_fd, SPI_IOC_MESSAGE(1),
                   &spi_message);
    (void)rv;

    // flip bytes the other way
    for (int i = 0; i < 30; i++)
      data_d[i] = (rx_buf[i] >> 8) + ((rx_buf[i] & 0xff) << 8);
    // data_d[i] = __bswap_16(rx_buf[i]);

    // copy back to data
    spine_to_spi(data, &g_spine_data, spi_board * 2);
  }
}

/*!
 * Run SPI
 */
void spi_driver_run() {
  // do spi board calculations
  for (int i = 0; i < 4; i++) {
    fake_spine_control(&spi_command_drv, &spi_data_drv, &spi_torque, i);
  }

  // in here, the driver is good
  pthread_mutex_lock(&spi_mutex);
  spi_send_receive(&spi_command_drv, &spi_data_drv);
  pthread_mutex_unlock(&spi_mutex);
}

/*!
 * Get the spi command
 */
spi_command_t *get_spi_command() {
  return &spi_command_drv;
}

/*!
 * Get the spi data
 */
spi_data_t *get_spi_data() { return &spi_data_drv; }

#endif
