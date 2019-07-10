#ifndef _rt_spi_lcm
#define _rt_spi_lcm

#include <spi_data_t.hpp>
#include <spi_torque_t.hpp>

#include <lcm/lcm-cpp.hpp>

void init_spi_lcm(lcm::LCM *main_lcm);

void publish_spi_data(spi_data_t *data);
void publish_spi_torque(spi_torque_t *data);

#endif
