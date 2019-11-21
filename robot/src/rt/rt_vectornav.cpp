/*!
 * @file rt_vectornav.cpp
 * @brief VectorNav IMU communication
 */

#ifdef linux

#include <inttypes.h>
#include <pthread.h>
#include <stdio.h>
#include <string>
#include <stdexcept>

#include <lcm/lcm-cpp.hpp>

#include "SimUtilities/IMUTypes.h"
#include "Utilities/utilities.h"
#include "rt/rt_vectornav.h"
#include "vectornav_lcmt.hpp"

#define K_MINI_CHEETAH_VECTOR_NAV_SERIAL "/dev/ttyS0"

//#define PRINT_VECTORNAV_DEBUG

int processErrorReceived(const std::string& errorMessage, VnError errorCode);
void vectornav_handler(void* userData, VnUartPacket* packet,
                       size_t running_index);
/*!
 * VectorNav Driver data
 */
typedef struct {
  VnSensor vs;
  BinaryOutputRegister bor;
} vn_sensor;

vn_sensor vn;

static lcm::LCM* vectornav_lcm;
vectornav_lcmt vectornav_lcm_data;
static VectorNavData* g_vn_data = nullptr;

/*!
 * Initialize Vectornav communication and set up sensor
 */
bool init_vectornav(VectorNavData* vn_data) {
  g_vn_data = vn_data;
  printf("[Simulation] Setup LCM...\n");
  vectornav_lcm = new lcm::LCM(getLcmUrl(255));
  if (!vectornav_lcm->good()) {
    printf("[ERROR] Failed to set up LCM\n");
    throw std::runtime_error("lcm bad");
  }

  VnError error;
  VpeBasicControlRegister vpeReg;
  ImuFilteringConfigurationRegister filtReg;
  const char SENSOR_PORT[] = K_MINI_CHEETAH_VECTOR_NAV_SERIAL;
  const uint32_t SENSOR_BAUDRATE = 115200;
  char modelNumber[30];
  char strConversions[50];
  uint32_t newHz, oldHz;
  // uint32_t hz_desired = 200;

  printf("[rt_vectornav] init_vectornav()\n");

  // initialize vectornav library
  VnSensor_initialize(&(vn.vs));

  // connect to sensor
  if ((error = VnSensor_connect(&(vn.vs), SENSOR_PORT, SENSOR_BAUDRATE)) !=
      E_NONE) {
    printf("[rt_vectornav] VnSensor_connect failed.\n");
    processErrorReceived("Error connecting to sensor.", error);
    return false;
  }

  // read the sensor's model number
  if ((error = VnSensor_readModelNumber(&(vn.vs), modelNumber,
                                        sizeof(modelNumber))) != E_NONE) {
    printf("[rt_vectornav] VnSensor_readModelNumber failed.\n");
    processErrorReceived("Error reading model number.", error);
    return false;
  }
  printf("Model Number: %s\n", modelNumber);

  // switch the sensor to 1 kHz mode
  if ((error = VnSensor_readAsyncDataOutputFrequency(&(vn.vs), &oldHz)) !=
      E_NONE) {
    printf("[rt_vectornav] VnSensor_readAsyncDataOutputFrequency failed.\n");
    processErrorReceived("Error reading async data output frequency.", error);
    return false;
  }

  // non-zero frequency causes the IMU to output ascii packets at the set
  // frequency, as well as binary
  if ((error = VnSensor_writeAsyncDataOutputFrequency(&(vn.vs), 0, true)) !=
      E_NONE) {
    printf("[rt_vectornav] VnSensor_wrtieAsyncDataOutputFrequency failed.\n");
    processErrorReceived("Error writing async data output frequency.", error);
    return false;
  }
  if ((error = VnSensor_readAsyncDataOutputFrequency(&(vn.vs), &newHz)) !=
      E_NONE) {
    printf("[rt_vectornav] VnSensor_readAsyncDataOutputFrequency failed.\n");
    processErrorReceived("Error reading async data output frequency.", error);
    return false;
  }
  printf("[rt_vectornav] Changed frequency from %d to %d Hz.\n", oldHz, newHz);

  // change to relative heading mode to avoid compass weirdness
  if ((error = VnSensor_readVpeBasicControl(&(vn.vs), &vpeReg)) != E_NONE) {
    printf("[rt_vectornav] VnSensor_ReadVpeBasicControl failed.\n");
    processErrorReceived("Error reading VPE basic control.", error);
    return false;
  }
  strFromHeadingMode(strConversions, (VnHeadingMode)vpeReg.headingMode);
  printf("[rt_vectornav] Sensor was in mode: %s\n", strConversions);
  vpeReg.headingMode = VNHEADINGMODE_RELATIVE;
  if ((error = VnSensor_writeVpeBasicControl(&(vn.vs), vpeReg, true)) !=
      E_NONE) {
    printf("[rt_vectornav] VnSensor_writeVpeBasicControl failed.\n");
    processErrorReceived("Error writing VPE basic control.", error);
    return false;
  }
  if ((error = VnSensor_readVpeBasicControl(&(vn.vs), &vpeReg)) != E_NONE) {
    processErrorReceived("Error reading VPE basic control.", error);
    printf("[rt_vectornav] VnSensor_ReadVpeBasicControl failed.\n");
    return false;
  }
  strFromHeadingMode(strConversions, (VnHeadingMode)vpeReg.headingMode);
  printf("[rt_vectornav] Sensor now id mode: %s\n", strConversions);

  if ((error = VnSensor_readImuFilteringConfiguration(&(vn.vs), &filtReg)) !=
      E_NONE) {
    printf("[rt_vectornav] VnSensor_readGyroCompensation failed.\n");
  }
  printf("[rt_vectornav] AccelWindow: %d\n", filtReg.accelWindowSize);
  //        filtReg.accelWindowSize = 4;    // We're sampling at 200 hz, but the
  //        imu samples at 800 hz. filtReg.accelFilterMode = 3; if((error =
  //        VnSensor_writeImuFilteringConfiguration(&(vn.vs), filtReg, true)) !=
  //        E_NONE)
  //        {
  //            printf("[rt_vectornav] VnSensor_writeGyroCompensation
  //            failed.\n");
  //        }

  // setup binary output message type
  BinaryOutputRegister_initialize(
      &(vn.bor), ASYNCMODE_PORT2,
      4,  // divisor:  output frequency = 800/divisor
      (CommonGroup)(COMMONGROUP_QUATERNION | COMMONGROUP_ANGULARRATE |
                    COMMONGROUP_ACCEL),
      TIMEGROUP_NONE, IMUGROUP_NONE, GPSGROUP_NONE, ATTITUDEGROUP_NONE,
      INSGROUP_NONE, GPSGROUP_NONE);

  if ((error = VnSensor_writeBinaryOutput1(&(vn.vs), &(vn.bor), true)) !=
      E_NONE) {
    printf("[rt_vectornav] VnSensor_writeBinaryOutput1 failed.\n");
    processErrorReceived("Error writing binary output 1.", error);
    return false;
  }

  // setup handler
  VnSensor_registerAsyncPacketReceivedHandler(&(vn.vs), vectornav_handler,
                                              NULL);
  printf("[rt_vectornav] IMU is set up!\n");
  return true;
}

int got_first_vectornav_message = 0;

/*!
 * Got new packet handler for vectornav
 */
void vectornav_handler(void* userData, VnUartPacket* packet,
                       size_t running_index) {
  (void)userData;
  (void)running_index;
  vec4f quat;
  vec3f omega;
  vec3f a;

  if (VnUartPacket_type(packet) != PACKETTYPE_BINARY) {
    printf("[vectornav_handler] got a packet that wasn't binary.\n");
    return;
  }

  if (!VnUartPacket_isCompatible(
          packet,
          (CommonGroup)(COMMONGROUP_QUATERNION | COMMONGROUP_ANGULARRATE |
                        COMMONGROUP_ACCEL),
          TIMEGROUP_NONE, IMUGROUP_NONE, GPSGROUP_NONE, ATTITUDEGROUP_NONE,
          INSGROUP_NONE, GPSGROUP_NONE)) {
    printf("[vectornav_handler] got a packet with the wrong type of data.\n");
    return;
  }

  quat = VnUartPacket_extractVec4f(packet);
  omega = VnUartPacket_extractVec3f(packet);
  a = VnUartPacket_extractVec3f(packet);

  for (int i = 0; i < 4; i++) {
    vectornav_lcm_data.q[i] = quat.c[i];
    g_vn_data->quat[i] = quat.c[i];
  }

  for (int i = 0; i < 3; i++) {
    vectornav_lcm_data.w[i] = omega.c[i];
    vectornav_lcm_data.a[i] = a.c[i];
    g_vn_data->gyro[i] = omega.c[i];
    g_vn_data->accelerometer[i] = a.c[i];
  }

  vectornav_lcm->publish("hw_vectornav", &vectornav_lcm_data);

#ifdef PRINT_VECTORNAV_DEBUG
  char strConversions[50];
  str_vec4f(strConversions, quat);
  printf("[QUAT] %s\n", strConversions);

  str_vec3f(strConversions, omega);
  printf("[OMEGA] %s\n", strConversions);

  str_vec3f(strConversions, a);
  printf("[ACC] %s\n", strConversions);
#endif
}

/*!
 * Error callback for vectornav
 */
int processErrorReceived(const std::string& errorMessage, VnError errorCode) {
  char errorCodeStr[100];
  strFromVnError(errorCodeStr, errorCode);
  printf("%s\nVECTORNAV ERROR: %s\n", errorMessage.c_str(), errorCodeStr);
  return -1;
}
#endif
