#include <pthread.h>
#include <rt/rt_spi_lcm.h>

lcm::LCM *min_lcm_spi = NULL;

// just adds the spi command handler to the main lcm
void init_spi_lcm(lcm::LCM *main_lcm) {
  printf("[RT SPI LCM] Initializing...\n");
  min_lcm_spi = main_lcm;
}

void publish_spi_data(spi_data_t *data) {
  min_lcm_spi->publish("CHEETAH_spi_data", data);
}

void publish_spi_torque(spi_torque_t *data) {
  if (min_lcm_spi == NULL) return;
  min_lcm_spi->publish("CHEETAH_spi_torque", data);
}
