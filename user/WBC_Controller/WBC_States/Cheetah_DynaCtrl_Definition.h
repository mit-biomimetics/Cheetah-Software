#ifndef Cheetah_DYNACORE_CONTROL_DEFINITION
#define Cheetah_DYNACORE_CONTROL_DEFINITION

#include <Configuration.h>
#include <Controllers/LegController.h>
#include <Dynamics/Quadruped.h>

#define CheetahConfigPath THIS_COM"user/WBC_Controller/WBC_States/config/"

template <typename T>
class Cheetah_Data {
 public:
  T ang_vel[3];
  T body_ori[4];
  T jpos[cheetah::num_act_joint];
  T jvel[cheetah::num_act_joint];
  bool foot_contact[4];
  T dir_command[2];
  T ori_command[3];

  T global_body_pos[3];

  int mode;
  bool cheater_mode;
};

template <typename T>
class Cheetah_Extra_Data {
 public:
  int num_step;
  T loc_x[100];  // Large enough number
  T loc_y[100];
  T loc_z[100];

  int num_path_pt;
  T path_x[100];
  T path_y[100];
  T path_z[100];

  int num_middle_pt;
  T mid_ori_roll[20];
  T mid_ori_pitch[20];
  T mid_ori_yaw[20];

  T mid_x[20];
  T mid_y[20];
  T mid_z[20];
};

#endif
