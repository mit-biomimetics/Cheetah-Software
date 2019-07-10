#ifndef KINEMATICS_BOUNDING_CTRL
#define KINEMATICS_BOUNDING_CTRL

#include <ParamHandler/ParamHandler.hpp>
#include <WBC_States/Controller.hpp>
#include <lcm-cpp.hpp>
#include <queue>
#include "ImpulseCurve.hpp"
#include "wbc_test_data_t.hpp"

template <typename T>
class ContactSpec;
template <typename T>
class KinWBC;
template <typename T>
class WBIC;
template <typename T>
class WBIC_ExtraData;
template <typename T>
class StateProvider;
template <typename T>
class BoundingTest;

template <typename T>
class KinBoundingCtrl : public Controller<T> {
 public:
  KinBoundingCtrl(BoundingTest<T>*, const FloatingBaseModel<T>*);
  virtual ~KinBoundingCtrl();

  virtual void OneStep(void* _cmd);
  virtual void FirstVisit();
  virtual void LastVisit();
  virtual bool EndOfPhase();

  /*!
   * Parameter setting from yaml file
   * @param category_name : name of this controller
   */
  virtual void CtrlInitialization(const std::string& category_name);
  virtual void SetTestParameter(const std::string& test_file);

 protected:
  std::vector<T> _Kp_joint, _Kd_joint;

  BoundingTest<T>* _bounding_test;

  bool _jump_signal;
  bool _b_jump_initiation;

  bool _b_first_stance;
  bool _b_front_swing;
  bool _b_hind_swing;

  bool _b_front_contact;
  bool _b_hind_contact;

  bool _b_front_contact_est;
  bool _b_hind_contact_est;

  std::queue<DVec<T> > _qdot_pre_queue;

  T _step_width;

  T _contact_vel_threshold;
  T _K_time;
  T _K_pitch;
  T _impact_amp;
  T _front_foot_offset;
  T _hind_foot_offset;
  T _swing_height;

  T _aerial_duration;
  T _front_swing_time;

  T _front_previous_stance;
  T _front_previous_swing;
  T _front_current_stance;

  T _hind_previous_stance;
  T _hind_previous_swing;
  T _hind_current_stance;

  T _nominal_gait_period;
  T _gait_period;

  T _end_time;
  T _swing_time;
  T _stance_time;
  T _nominal_stance_time;
  int _dim_contact;

  T _ctrl_start_time;
  T _front_start_time;
  T _hind_start_time;

  T _front_time;
  T _hind_time;

  T _step_length_lim;

  Task<T>* _local_head_pos_task;
  Task<T>* _local_tail_pos_task;

  Task<T>* _local_roll_task;
  Task<T>* _body_ryrz_task;

  Task<T>* _fr_foot_local_task;
  Task<T>* _fl_foot_local_task;
  Task<T>* _hr_foot_local_task;
  Task<T>* _hl_foot_local_task;

  Vec3<T> _fr_foot_pos;
  DVec<T> _fr_foot_vel;
  DVec<T> _fr_foot_acc;

  Vec3<T> _fl_foot_pos;
  DVec<T> _fl_foot_vel;
  DVec<T> _fl_foot_acc;

  Vec3<T> _hr_foot_pos;
  DVec<T> _hr_foot_vel;
  DVec<T> _hr_foot_acc;

  Vec3<T> _hl_foot_pos;
  DVec<T> _hl_foot_vel;
  DVec<T> _hl_foot_acc;

  ContactSpec<T>* _fr_contact;
  ContactSpec<T>* _fl_contact;
  ContactSpec<T>* _hr_contact;
  ContactSpec<T>* _hl_contact;

  WBIC<T>* _wbic;
  WBIC_ExtraData<T>* _wbic_data;

  T _target_leg_height;
  T _total_mass;

  ImpulseCurve<T> _front_z_impulse;
  ImpulseCurve<T> _hind_z_impulse;

  Vec3<T> _ini_front_body;
  Vec3<T> _ini_hind_body;

  Vec3<T> _ini_fr;
  Vec3<T> _ini_fl;
  Vec3<T> _ini_hr;
  Vec3<T> _ini_hl;

  Vec3<T> _fin_fr;
  Vec3<T> _fin_fl;
  Vec3<T> _fin_hr;
  Vec3<T> _fin_hl;

  void _body_task_setup();
  void _leg_task_setup();
  void _contact_update();
  void _compute_torque_wbic(DVec<T>& gamma);

  void _StatusCheck();
  void _setupTaskAndContactList();
  void _ContactUpdate();

  void _lcm_data_sending();
  void _save_file();

  DVec<T> _Fr_des;
  DVec<T> _Fr_result;

  DVec<T> _des_jpos;
  DVec<T> _des_jvel;
  DVec<T> _des_jacc;

  KinWBC<T>* _kin_wbc;
  std::vector<Task<T>*> _kin_task_list;
  std::vector<ContactSpec<T>*> _kin_contact_list;

  std::string _folder_name;
  ParamHandler* _param_handler;
  StateProvider<T>* _sp;
  lcm::LCM _wbcLCM;
  wbc_test_data_t _wbc_data_lcm;
};

#endif
