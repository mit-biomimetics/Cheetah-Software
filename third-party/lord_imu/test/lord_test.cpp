#include <cstdio>

#include "mip_sdk.h"
#include "mip_gx4_25.h"
#include "../LordImu.h"
#include "Math/orientation_tools.h"


#include <stdio.h>
#include <unistd.h>


int main(int argc, char** argv) {
  u32 com_port, baudrate;
  if(argc != 3) {
    printf("usage: imu-test com-port baudrate\n");
    return 1;
  }

  com_port = std::atoi(argv[1]);
  baudrate = std::atoi(argv[2]);

  LordImu imu;
  if(imu.tryInit(com_port, baudrate)) {
    while(true) {
      imu.run();
      //std::cout << "Q: " << imu.quat.transpose() << "\n";
      Vec3<float> rpy = ori::quatToRPY(imu.quat); (void)rpy;
      Vec3<float> acc = imu.acc; (void)acc;
      Vec3<float> ang = imu.gyro; (void)ang;
      //printf("rpy: %.3f %.3f %.3f\n", rpy[0], rpy[1], rpy[2]);
      //printf("acc: %.3f %.3f %.3f\n", acc[0], acc[1], acc[2]);
      printf("ang: %.3f %.3f %.3f\n", ang[0], ang[1], ang[2]);
    }
  }


  return 0;
}
