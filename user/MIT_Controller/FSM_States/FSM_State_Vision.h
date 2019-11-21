#ifndef FSM_STATE_VISION_H
#define FSM_STATE_VISION_H

#include <Controllers/convexMPC/ConvexMPCLocomotion.h>
#include <Controllers/VisionMPC/VisionMPCLocomotion.h>
#include "FSM_State.h"
#include <thread>
#include <lcm-cpp.hpp>
#include "heightmap_t.hpp"
#include "traversability_map_t.hpp"
#include "velocity_visual_t.hpp"
#include "obstacle_visual_t.hpp"
#include "localization_lcmt.hpp"

template<typename T> class WBC_Ctrl;
template<typename T> class LocomotionCtrlData;
/**
 *
 */
template <typename T>
class FSM_State_Vision : public FSM_State<T> {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  FSM_State_Vision(ControlFSMData<T>* _controlFSMData);

  // Behavior to be carried out when entering a state
  void onEnter();

  // Run the normal behavior for the state
  void run();

  // Checks for any transition triggers
  FSM_StateName checkTransition();

  // Manages state specific transitions
  TransitionData<T> transition();

  // Behavior to be carried out when exiting a state
  void onExit();

 private:
  // Keep track of the control iterations
  int iter = 0;
  VisionMPCLocomotion vision_MPC;
  ConvexMPCLocomotion cMPCOld;

  WBC_Ctrl<T> * _wbc_ctrl;
  LocomotionCtrlData<T> * _wbc_data;

  Vec3<T> _ini_body_pos;
  Vec3<T> _ini_body_ori_rpy;
  Vec3<T> zero_vec3;
  Vec3<T> _global_robot_loc;
  Vec3<T> _robot_rpy;

  size_t x_size = 100;
  size_t y_size = 100;
  double grid_size = 0.015;

  DMat<T> _height_map;
  DMat<int> _idx_map;

  void handleHeightmapLCM(const lcm::ReceiveBuffer* rbuf, const std::string& chan, const heightmap_t* msg);
  void handleIndexmapLCM(const lcm::ReceiveBuffer* rbuf, const std::string& chan, const traversability_map_t* msg);
  void handleLocalization(const lcm::ReceiveBuffer* rbuf, const std::string& chan, const localization_lcmt* msg);
  bool _b_localization_data = false;
  void visionLCMThread() { while (true) { _visionLCM.handle(); } }

  lcm::LCM _visionLCM;
  std::thread _visionLCMThread;

  vectorAligned< Vec3<T> > _obs_list; // loc, height
  obstacle_visual_t _obs_visual_lcm;

  void _updateStateEstimator();
  void _JPosStand();
  void _UpdateObstacle();
  void _LocomotionControlStep(const Vec3<T> & vel_cmd);
  void _UpdateVelCommand(Vec3<T> & vel_cmd);
  void _RCLocomotionControl();
  void _Visualization(const Vec3<T> & des_vel);
  void _print_obstacle_list();

  Vec3<T> _target_pos;

};

#endif  // FSM_STATE_LOCOMOTION_H
