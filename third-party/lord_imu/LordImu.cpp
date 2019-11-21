#include "LordImu.h"
#include <cstdio>
#include <stdexcept>
#include <thread>
#include "../../common/include/cppTypes.h"
#include "../../common/include/Dynamics/spatial.h"
#include "../../common/include/Math/orientation_tools.h"

constexpr u32 IMU_PACKET_TIMEOUT_MS = 1000;
//constexpr u32 MIP_SDK_GX4_25_IMU_DIRECT_MODE = 0x02;
constexpr u32 MIP_SDK_STANDARD_MODE = 0x01;

static std::mutex dataMutex;

bool LordImu::tryInit(u32 port, u32 baud_rate) {
  try {
    init(port, baud_rate);
  } catch(std::exception& e) {
    printf("[LordIMU] failed to initialize: %s\n", e.what());
    return false;
  }

  return true;
}

void LordImu::init(u32 port, u32 baud_rate) {
  printf("[Lord IMU] Open port %d, baud rate %d\n", port, baud_rate);

  if(mip_interface_init(port, baud_rate, &device_interface, IMU_PACKET_TIMEOUT_MS) != MIP_INTERFACE_OK) {
    throw std::runtime_error("Failed to initialize MIP interface for IMU\n");
  }

  printf("[Lord IMU] Port open. Mode setup...\n");
  mode_setup();
  printf("[Lord IMU] Get info...\n");
  get_device_info();
  printf("[Lord IMU] Self test...\n");
  //self_test();
  printf("[Lord IMU] Basic report...\n");
  basic_report();
//  printf("[Lord IMU] Zero Gyro...\n");
//  zero_gyro();
  printf("[Lord IMU] Setup IMU...\n");
  setup_streaming();
  printf("[Lord IMU] Enable Data...\n");
  enable();
}

void LordImu::mode_setup() {
  u8 com_mode = MIP_SDK_STANDARD_MODE;
  printf("[Lord IMU] Set direct mode\n");

  while(mip_system_com_mode(&device_interface, MIP_FUNCTION_SELECTOR_WRITE, &com_mode) != MIP_INTERFACE_OK) {
    printf("failed to set com mode\n");
  }

  while(mip_system_com_mode(&device_interface, MIP_FUNCTION_SELECTOR_READ, &com_mode) != MIP_INTERFACE_OK) {
    printf("failed to read com mode\n");
  }

  if(com_mode != MIP_SDK_STANDARD_MODE) {
    printf("failed to set mode, wanted %d, got %d\n", MIP_SDK_STANDARD_MODE, com_mode);
  } else {
    printf("[Lord IMU] Done\n");
  }

  usleep(100000);


  while(mip_base_cmd_idle(&device_interface) != MIP_INTERFACE_OK){
    printf("idle fail\n");
  }

  usleep(100000);

  printf("[Lord IMU] Ping...\n");

  while(mip_base_cmd_ping(&device_interface) != MIP_INTERFACE_OK){
    printf("fail\n");
  }
}

static std::string dev_string(u16* str) {
  char temp[32];
  char* in = (char*)str;

  u32 i = 0;
  while((!in[i] || in[i] == ' ') && i < 16) {
    i++;
  }

  char* charp = temp;
  while(i < 16) {
    *charp = in[i];
    charp++;
    i++;
  }
  *charp = 0;
  return std::string(temp);
}

void LordImu::get_device_info() {
  base_device_info_field device_info;
  while(mip_base_cmd_get_device_info(&device_interface, &device_info) != MIP_INTERFACE_OK){
    printf("fail\n");
  }

  u16 base_rate = 0;
  while(mip_3dm_cmd_get_ahrs_base_rate(&device_interface, &base_rate) != MIP_INTERFACE_OK){
    printf("fail\n");
  }

  u16 filter_rate = 0;
  while(mip_3dm_cmd_get_filter_base_rate(&device_interface, &filter_rate) != MIP_INTERFACE_OK){
    printf("fail\n");
  }



  deviceInfo.modelName = dev_string(device_info.model_name);
  deviceInfo.modelNumber = dev_string(device_info.model_number);
  deviceInfo.serialNumber = dev_string(device_info.serial_number);
  deviceInfo.lotNumber = dev_string(device_info.lotnumber);
  deviceInfo.deviceOptions = dev_string(device_info.device_options);
  deviceInfo.firmwareVersion = device_info.firmware_version;
  deviceInfo.ahrs_base_rate = base_rate;
  deviceInfo.filter_rate = filter_rate;

  printf("[Lord IMU] Got device info:\n");
  printf("  name: %s\n"
         "  model: %s\n"
         "  serial: %s\n"
         "  lot: %s\n"
         "  options: %s\n"
         "  firmware: %d\n"
         "  ahrs rate: %d Hz\n"
         "  filter rate: %d Hz\n",
         deviceInfo.modelName.c_str(),
         deviceInfo.modelNumber.c_str(),
         deviceInfo.serialNumber.c_str(),
         deviceInfo.lotNumber.c_str(),
         deviceInfo.deviceOptions.c_str(),
         deviceInfo.firmwareVersion,
         deviceInfo.ahrs_base_rate,
         deviceInfo.filter_rate);
}


void LordImu::self_test() {
  u32 bit_result;
  while(mip_base_cmd_built_in_test(&device_interface, &bit_result) != MIP_INTERFACE_OK){
    printf("fail\n");
  }

  if(bit_result) {
    printf("[Lord IMU] Self test failed: 0x%x\n", bit_result);
    throw std::runtime_error("self test fail\n");
  } else {
    printf("[Lord IMU] Self test passed\n");
  }
}


void LordImu::basic_report() {
//  while(mip_3dm_cmd_hw_specific_imu_device_status(&device_interface, 6237, 1, &basicStatus) != MIP_INTERFACE_OK){
//    printf("fail\n");
//  }
}
static LordImu* gLordImu;
static void filter_callback(void* user_ptr, u8* packet, u16 packet_size, u8 callback_type) {
  (void)user_ptr;
  (void)packet_size;

  mip_field_header* field_header;
  u8* field_data;
  u16 field_offset = 0;
  mip_filter_attitude_quaternion quat;

//  Mat3<float> R = Mat3<float>::Identity();
  //R << 0, -1, 0, -1, 0, 0, 0, 0, 1;

  switch(callback_type) {
    case MIP_INTERFACE_CALLBACK_VALID_PACKET:
      while(mip_get_next_field(packet, &field_header,
          &field_data, &field_offset) == MIP_OK) {
        switch(field_header->descriptor) {
          case MIP_FILTER_DATA_ATT_QUATERNION:
          {
            memcpy(&quat, field_data, sizeof(mip_filter_attitude_quaternion));
            mip_filter_attitude_quaternion_byteswap(&quat);
            dataMutex.lock();
            gLordImu->quat = Vec4<float>(quat.q);
            Mat3<float> g_R_imu, r_R_imup;

            g_R_imu << 0, 1, 0, 1, 0, 0, 0, 0, -1;
            r_R_imup << 0, 0, 1, 0, -1, 0, 1, 0, 0;

            Vec4<float> ql = ori::rotationMatrixToQuaternion(g_R_imu.transpose());
            Vec4<float> qr = ori::rotationMatrixToQuaternion(r_R_imup);
            gLordImu->quat = ori::quatProduct(ql, ori::quatProduct(gLordImu->quat, qr));

            gLordImu->good_packets++;
            dataMutex.unlock();
          }
            break;
          default:
            printf("[Lord IMU] Unknown FILTER packet %d\n", field_header->descriptor);
            break;
        }
      }
      break;
    case MIP_INTERFACE_CALLBACK_CHECKSUM_ERROR:
      gLordImu->invalid_packets++;
      break;
    case MIP_INTERFACE_CALLBACK_TIMEOUT:
      gLordImu->timeout_packets++;
      break;
    default:
      gLordImu->unknown_packets++;
      break;
  }

}



static void ahrs_callback(void* user_ptr, u8* packet, u16 packet_size, u8 callback_type) {
  (void)user_ptr;
  (void)packet_size;

  mip_field_header* field_header;
  u8* field_data;
  u16 field_offset = 0;
  mip_ahrs_scaled_accel accel;
  mip_ahrs_scaled_gyro gyro;



  switch(callback_type) {
    case MIP_INTERFACE_CALLBACK_VALID_PACKET:
      //gLordImu->good_packets++;
      while(mip_get_next_field(packet, &field_header,
          &field_data, &field_offset) == MIP_OK) {
        switch(field_header->descriptor) {
          case MIP_AHRS_DATA_ACCEL_SCALED:
            memcpy(&accel, field_data, sizeof(mip_ahrs_scaled_accel));
            mip_ahrs_scaled_accel_byteswap(&accel);
            //gLordImu->acc = Vec3<float>(accel.scaled_accel);
            dataMutex.lock();
            gLordImu->acc[0] = 9.81f * accel.scaled_accel[2];
            gLordImu->acc[1] = -9.81f * accel.scaled_accel[1];
            gLordImu->acc[2] = 9.81f * accel.scaled_accel[0];
            gLordImu->good_packets++;
            dataMutex.unlock();
            break;
          case MIP_AHRS_DATA_GYRO_SCALED:
            dataMutex.lock();
            memcpy(&gyro, field_data, sizeof(mip_ahrs_scaled_gyro));
            mip_ahrs_scaled_gyro_byteswap(&gyro);
            //gLordImu->gyro = Vec3<float>(gyro.scaled_gyro);
            gLordImu->gyro[0] = gyro.scaled_gyro[2];
            gLordImu->gyro[1] = -gyro.scaled_gyro[1];
            gLordImu->gyro[2] = gyro.scaled_gyro[0];
            gLordImu->good_packets++;
            dataMutex.unlock();
            break;
          default:
            printf("[Lord IMU] Unknown AHRS packet %d\n", field_header->descriptor);
            break;
        }
      }
      break;
    case MIP_INTERFACE_CALLBACK_CHECKSUM_ERROR:
      gLordImu->invalid_packets++;
      break;
    case MIP_INTERFACE_CALLBACK_TIMEOUT:
      gLordImu->timeout_packets++;
      break;
    default:
      gLordImu->unknown_packets++;
      break;
  }
}


void LordImu::setup_streaming() {
  gLordImu = this;

  u8 enable = MIP_3DM_CONING_AND_SCULLING_DISABLE;

  while(mip_3dm_cmd_coning_sculling_compensation(&device_interface, MIP_FUNCTION_SELECTOR_WRITE, &enable) != MIP_INTERFACE_OK){
    printf("fail\n");
  }

  // callbacks
  if(mip_interface_add_descriptor_set_callback(&device_interface,
      MIP_FILTER_DATA_SET, nullptr, filter_callback) != MIP_INTERFACE_OK) {
    throw std::runtime_error("failed to set IMU filter callback");
  }

  if(mip_interface_add_descriptor_set_callback(&device_interface,
      MIP_AHRS_DATA_SET, nullptr, ahrs_callback) != MIP_INTERFACE_OK) {
    throw std::runtime_error("failed to set IMU ahrs callback");
  }

  printf("[Lord IMU] Setup message types...\n");
  u8 data_types[3] = {MIP_AHRS_DATA_GYRO_SCALED, MIP_AHRS_DATA_ACCEL_SCALED};
  u16 data_downsampling[3] = {1, 1};

  u8 num_entries = 2;
  while(mip_3dm_cmd_ahrs_message_format(&device_interface,
      MIP_FUNCTION_SELECTOR_WRITE, &num_entries, data_types, data_downsampling) != MIP_INTERFACE_OK) {
    printf("fail\n");
  }

  data_types[0] = {MIP_FILTER_DATA_ATT_QUATERNION};
  num_entries = 1;

  while(mip_3dm_cmd_filter_message_format(&device_interface,
       MIP_FUNCTION_SELECTOR_WRITE, &num_entries, data_types, data_downsampling) != MIP_INTERFACE_OK) {
    printf("fail\n");
  }
}

void LordImu::enable() {

//  u8 hs = MIP_FILTER_HEADING_SOURCE_MAGNETOMETER;
//  while(mip_filter_heading_source(&device_interface, MIP_FUNCTION_SELECTOR_WRITE, &hs) != MIP_INTERFACE_OK) {
//    printf("hs fail\n");
//  }

  while(mip_filter_reset_filter(&device_interface) != MIP_INTERFACE_OK){
    printf("reset fail\n");
  }

  float angles[3];
  angles[0] = angles[1] = angles[2] = 0;

  while(mip_filter_set_init_attitude(&device_interface, angles) != MIP_INTERFACE_OK){
    printf("init fail\n");
  }

  u8 option = 0; // no magnetometer
  while(mip_filter_heading_source(&device_interface, 1, &option) != MIP_INTERFACE_OK) {
    printf("setup heading fail\n");
  }

  u8 enable_value = 1;
  while(mip_3dm_cmd_continuous_data_stream(&device_interface,
      MIP_FUNCTION_SELECTOR_WRITE,
      MIP_3DM_AHRS_DATASTREAM,
      &enable_value) != MIP_INTERFACE_OK) {
    printf("fail\n");
  }

  while(mip_3dm_cmd_continuous_data_stream(&device_interface,
                                           MIP_FUNCTION_SELECTOR_WRITE,
                                           MIP_3DM_INS_DATASTREAM,
                                           &enable_value) != MIP_INTERFACE_OK) {
    printf("fail\n");
  }

  printf("[Lord IMU] Ready to run!\n");
}

void LordImu::run() {
  for(u32 i = 0; i < 32; i++)
    mip_interface_update(&device_interface);
  usleep(100);
}

void LordImu::updateLCM(microstrain_lcmt *message) {
  dataMutex.lock();
  for(u32 i = 0; i < 4; i++) {
    message->quat[i] = quat[i];
  }

  Vec3<float> rpy = ori::quatToRPY(quat);
  for(u32 i = 0; i < 3; i++) {
    message->rpy[i] = rpy[i];
    message->acc[i] = acc[i];
    message->omega[i] = gyro[i];
  }

  message->good_packets = good_packets;
  message->bad_packets = invalid_packets + unknown_packets + timeout_packets;
  dataMutex.unlock();
}

void LordImu::print_packet_stats() {
  printf("-- IMU Packets --\n"
         "good: %d\n"
         "invalid: %d\n"
         "timeout: %d\n"
         "unknown: %d\n",
         good_packets, invalid_packets,
         timeout_packets, unknown_packets);
}

void LordImu::zero_gyro() {
  u16 duration = 5000; //milliseconds
  float bias_vector[3];

  printf("please wait...\n");
  while(mip_3dm_cmd_capture_gyro_bias(&device_interface, duration, bias_vector) != MIP_INTERFACE_OK){
    printf("fail\n");
  }

  printf("Gyro Bias Captured:\nbias_vector[0] = %f\nbias_vector[1] = %f\nbias_vector[2] = %f\n\n", bias_vector[0], bias_vector[1], bias_vector[2]);
}
