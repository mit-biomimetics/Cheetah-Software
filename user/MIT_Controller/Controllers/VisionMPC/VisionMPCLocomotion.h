#ifndef CHEETAH_SOFTWARE_VISION_MPCLOCOMOTION_H
#define CHEETAH_SOFTWARE_VISION_MPCLOCOMOTION_H

#include <Controllers/FootSwingTrajectory.h>
#include <FSM_States/ControlFSMData.h>
#include "cppTypes.h"

using Eigen::Array4f;
using Eigen::Array4i;


class VisionGait
{
public:
  VisionGait(int nMPC_segments, Vec4<int> offsets, Vec4<int>  durations, const std::string& name="");
  ~VisionGait();
  Vec4<float> getContactState();
  Vec4<float> getSwingState();
  int* mpc_gait();
  void setIterations(int iterationsPerMPC, int currentIteration);
  int _stance;
  int _swing;


private:
  int _nMPC_segments;
  int* _mpc_table;
  Array4i _offsets; // offset in mpc segments
  Array4i _durations; // duration of step in mpc segments
  Array4f _offsetsFloat; // offsets in phase (0 to 1)
  Array4f _durationsFloat; // durations in phase (0 to 1)
  int _iteration;
  int _nIterations;
  float _phase;


};


class VisionMPCLocomotion {
public:
  VisionMPCLocomotion(float _dt, int _iterations_between_mpc, MIT_UserParameters* parameters);
  void initialize();

  template<typename T>
  void run(ControlFSMData<T>& data, 
      const Vec3<T> & vel_cmd, const DMat<T> & height_map, const DMat<int> & idx_map);

  Vec3<float> pBody_des;
  Vec3<float> vBody_des;
  Vec3<float> aBody_des;

  Vec3<float> pBody_RPY_des;
  Vec3<float> vBody_Ori_des;

  Vec3<float> pFoot_des[4];
  Vec3<float> vFoot_des[4];
  Vec3<float> aFoot_des[4];

  Vec3<float> Fr_des[4];

  Vec4<float> contact_state;

private:
  void _UpdateFoothold(Vec3<float> & foot, const Vec3<float> & body_pos,
      const DMat<float> & height_map, const DMat<int> & idx_map);
  void _IdxMapChecking(int x_idx, int y_idx, int & x_idx_selected, int & y_idx_selected, 
      const DMat<int> & idx_map);

  Vec3<float> _fin_foot_loc[4];
  float grid_size = 0.015;

  Vec3<float> v_des_world;
  Vec3<float> rpy_des;
  Vec3<float> v_rpy_des;

  float _body_height = 0.31;
  void updateMPCIfNeeded(int* mpcTable, ControlFSMData<float>& data);
  void solveDenseMPC(int *mpcTable, ControlFSMData<float> &data);
  int iterationsBetweenMPC;
  int horizonLength;
  float dt;
  float dtMPC;
  int iterationCounter = 0;
  Vec3<float> f_ff[4];
  Vec4<float> swingTimes;
  FootSwingTrajectory<float> footSwingTrajectories[4];
  VisionGait trotting, bounding, pronking, galloping, standing, trotRunning;
  Mat3<float> Kp, Kd, Kp_stance, Kd_stance;
  bool firstRun = true;
  bool firstSwing[4];
  float swingTimeRemaining[4];
  float stand_traj[6];
  int current_gait;
  int gaitNumber;

  Vec3<float> world_position_desired;
  Vec3<float> rpy_int;
  Vec3<float> rpy_comp;
  Vec3<float> pFoot[4];
  float trajAll[12*36];

  MIT_UserParameters* _parameters = nullptr;
};


#endif //CHEETAH_SOFTWARE_VISION_MPCLOCOMOTION_H
