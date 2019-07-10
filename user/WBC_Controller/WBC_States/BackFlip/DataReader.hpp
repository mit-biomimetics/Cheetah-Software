#ifndef BACKFLIP_DATA_READER_H
#define BACKFLIP_DATA_READER_H
#include <cppTypes.h>

enum plan_offsets {
  q0_offset = 0,     // x, z, yaw, front hip, front knee, rear hip, rear knee
  qd0_offset = 7,    // x, z, yaw, front hip, front knee, rear hip, rear knee
  tau_offset = 14,   // front hip, front knee, rear hip, rear knee
  force_offset = 18  // front x, front z, rear x, rear z
};

typedef Eigen::Matrix<float, 7, 1> Vector7f;

class DataReader {
 public:
  static const int plan_cols = 22;

  DataReader(const RobotType &);
  void load_control_plan(const char *filename);
  void unload_control_plan();
  float *get_initial_configuration();
  float *get_plan_at_time(int timestep);
  int plan_timesteps = -1;

 private:
  RobotType _type;
  float *plan_buffer;
  bool plan_loaded = false;
};

#endif  // BACKFLIP_DATA_READER_H
