#ifndef WBIC_TROT_TEST_eheetah
#define WBIC_TROT_TEST_Cheetah

#include <Dynamics/Quadruped.h>
#include <Utilities/filters.h>
#include <WBC_States/Test.hpp>

#include <lcm-cpp.hpp>
#include <thread>
#include "vision_data_t.hpp"
#include "wbc_test_data_t.hpp"

template <typename T>
class StateProvider;
template <typename T>
class Path;

namespace WBICTrotPhase {
constexpr int lift_up = 0;
constexpr int full_contact_1 = 1;
constexpr int frhl_swing_start_trans = 2;
constexpr int frhl_swing = 3;
constexpr int frhl_swing_end_trans = 4;
constexpr int full_contact_2 = 5;
constexpr int flhr_swing_start_trans = 6;
constexpr int flhr_swing = 7;
constexpr int flhr_swing_end_trans = 8;
constexpr int NUM_TROT_PHASE = 9;
};  // namespace WBICTrotPhase

template <typename T>
class WBICTrotTest : public Test<T> {
 public:
  WBICTrotTest(FloatingBaseModel<T>*, const RobotType&);
  virtual ~WBICTrotTest();

  DVec<T> _jpos_des_pre;

  DVec<T> _des_jpos;
  DVec<T> _des_jvel;
  DVec<T> _des_jacc;

  Vec3<T> _body_pos;
  Vec3<T> _body_vel;
  Vec3<T> _body_acc;

  Vec3<T> _body_ori_rpy;
  // This does not have any frame meaning. just derivation of rpy
  Vec3<T> _body_ang_vel;

  Vec3<T> _front_foot_loc;
  Vec3<T> _hind_foot_loc;

  std::string _folder_name;

  void handleVisionLCM(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                       const vision_data_t* msg);

  void visionLCMThread();

 protected:
  bool _b_remote_ctrl;
  Path<T>* _path;
  Vec3<T> _vision_loc;  // (x, y, yaw)
  lcm::LCM _visionLCM;
  std::thread _visionLCMThread;
  bool _b_loc_update;

  std::vector<filter<T>*> _filtered_input_vel;
  Vec3<T> _input_vel;
  T _target_body_height;
  virtual void _UpdateTestOneStep();
  virtual void _TestInitialization();
  virtual int _NextPhase(const int& phase);
  virtual void _UpdateExtraData(Cheetah_Extra_Data<T>* ext_data);

  void _SettingParameter();

  Controller<T>* body_up_ctrl_;
  Controller<T>* body_ctrl_;

  // Front Right and Hind Left leg swing
  Controller<T>* frhl_swing_start_trans_ctrl_;
  Controller<T>* frhl_swing_ctrl_;
  Controller<T>* frhl_swing_end_trans_ctrl_;

  // Front Left and Hind Right leg swing
  Controller<T>* flhr_swing_start_trans_ctrl_;
  Controller<T>* flhr_swing_ctrl_;
  Controller<T>* flhr_swing_end_trans_ctrl_;

  StateProvider<T>* _sp;

  lcm::LCM _wbcLCM;
  wbc_test_data_t _wbc_data_lcm;
};

#endif
