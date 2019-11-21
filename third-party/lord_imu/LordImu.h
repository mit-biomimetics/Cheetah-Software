#ifndef PROJECT_LORDIMU_H
#define PROJECT_LORDIMU_H

#include <string>
#include <mutex>

#include "cTypes.h"
#include "cppTypes.h"
#include "cppTypes.h"
#include "mip_sdk.h"
#include "mip_gx4_25.h"

#include "../../lcm-types/cpp/microstrain_lcmt.hpp"

struct LordImuDeviceInfo {
  std::string modelName;
  std::string modelNumber;
  std::string serialNumber;
  std::string lotNumber;
  std::string deviceOptions;
  u16 firmwareVersion;
  u32 ahrs_base_rate;
  u32 filter_rate;
};

struct LordImuDescriptor {
  u8 descriptor;
  u8 descriptorSet;
};

struct LordImuBasicStatus {
  u16 model;
  u8 status_type;
  u32 flags;
  u32 system_time_ms;
};

struct LordImuDiagnosticStatus {
  u16 model;
  u8 status_type;
  u32 flags;
  u32 system_time_ms;
  u8 has_magnetometer;
  u8 has_pressure;
  u8 power_state;
  u16 gyro_range;
  u16 accel_range;
  u32 junk;

  float temp;
  u32 temp_read_ms;
  u8 temp_error;

  u32 gnss_count;
  u32 gnss_ms;

  u8 imu_enable;
  u32 imu_outgoing_dropped;

  u32 bytes_written;
  u32 bytes_read;
  u32 write_overruns;
  u32 read_overruns;
};

class LordImu {
public:
  void init(u32 port, u32 baud_rate);
  bool tryInit(u32 port, u32 baud_rate);
  void mode_setup();
  void get_device_info();
  void self_test();
  void basic_report();
  void setup_streaming();
  void enable();
  void print_packet_stats();
  void zero_gyro();
  void run();
  void updateLCM(microstrain_lcmt* message);

  u32 invalid_packets = 0;
  u32 timeout_packets = 0;
  u32 unknown_packets = 0;
  u32 good_packets = 0;

  Vec3<float> gyro;
  Vec3<float> acc;
  Vec4<float> quat;


private:
  mip_interface device_interface;
  LordImuDeviceInfo deviceInfo;
  LordImuBasicStatus basicStatus;

};


#endif //PROJECT_LORDIMU_H
